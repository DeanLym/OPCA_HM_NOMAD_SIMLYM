/*-----------------------------------------------------*/
/*  how to use the NOMAD library with a user function  */
/*-----------------------------------------------------*/
#include "nomad.hpp"
using namespace std;
using namespace NOMAD; //avoids putting  everywhere


#include "opca_bimodal.h"
#include "sim_model.hpp"
#include "generate_opca_model.hpp"
#include "util_funs.hpp"
// #include <armadillo>

// using namespace arma;

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public Evaluator {
public:
	My_Evaluator  ( const Parameters & p ) :
		Evaluator ( p ) {
		opca_bm_  = new OPCA_BIMODAL(1,1,1,1);
		dim_opca_ = 1;
		dim_kr_   = 1;
	}

	~My_Evaluator ( void ) {
		delete opca_bm_;
	}



	void update_iteration	(	NOMAD::success_type 	success,
			const NOMAD::Stats & 	stats,
			const NOMAD::Evaluator_Control & 	ev_control,
			const NOMAD::Barrier & 	true_barrier,
			const NOMAD::Barrier & 	sgte_barrier,
			const NOMAD::Pareto_Front & 	pareto_front,
			bool & 	stop
	)	{
		//		  num_iterations_ += 1;
		int num_iter  = stats.get_iterations();
		//		  cout << "Iteration # " << num_iter << "....";
		//		  if(success    == NOMAD::UNSUCCESSFUL ){
		//			  cout << "fails" << endl;
		//		  }else{
		//			  cout << "is successful." << endl;
		//		  }
		int num_eval  = stats.get_eval();
		//		  cout << "Number of function evaluations: " << num_eval << "..." << endl;
		int real_time = stats.get_real_time();
		//		  cout << "Wall clock time: " << real_time << "..." << endl;

		ofstream ofs;
		ofs.open("iter.txt",std::ofstream::out | std::ofstream::app);
		ofs << num_iter << "\t" << num_eval << "\t" << real_time << endl;
		ofs.close();

	}

	bool eval_x ( Eval_Point   & x          ,
			const Double & h_max      ,
			bool                & count_eval   ) const
	{
		vector<double> xi;
		double temp;
		for(int i=0;i<dim_opca_;i++){
			temp = x[i].value();
			xi.push_back(temp);
		}
		TranformUniform2Normal(dim_opca_,xi);

		double Nc = 3600;
		vector<double> m;
		m.resize(Nc);
		opca_bm_->GenerateOPCARealization(&(xi[0]) , &(m[0]));
		vector<double> perm;
		perm.resize(Nc);
		GeneratePerm(Nc, &(m[0]), &(perm[0]));

		// Transform kr parameters
		vector<double> kr;

		// Prior mean for kr parameters
		vector<double> kr_avg;
		kr_avg.push_back(0.16);
		kr_avg.push_back(0.16);
		kr_avg.push_back(2.6);
		kr_avg.push_back(2.4);
		kr_avg.push_back(0.38);
		kr_avg.push_back(0.7);

		vector<double> kr_std;
		kr_std.push_back(0.05);
		kr_std.push_back(0.05);
		kr_std.push_back(0.5);
		kr_std.push_back(0.5);
		kr_std.push_back(0.1);
		kr_std.push_back(0.1);

		for(int i=0; i<dim_kr_; i++){
			temp = x[i+dim_opca_].value();
			temp = TransformUniform2Normal_2(temp,kr_avg[i],kr_std[i]);
			kr.push_back(temp);
		}

#ifdef DEBUG
		SaveData("xi.debug",dim_opca_,&(xi[0]));
		SaveData("m.debug",Nc,&(m[0]));
		SaveData("perm.debug",Nc,&(perm[0]));
		SaveData("kr.debug",dim_kr_,&(kr[0]));
#endif

		SimCtrl* sim = GetSimulationModel(&perm[0],kr);

		sim->display_level_ = 0;
		sim->RunSim();
#ifdef DEBUG
		sim->OutputResult();
#endif
		sim->hm_ = new CHistoryMatching;
		ifstream in;
		in.open("HIST_FILE.DATA");
		string hist_file, temp_str;
		if(in.is_open())
			in >> temp_str >> hist_file;
		else{
			throw runtime_error("Can not open HIST_FILE.DATA");
			//			cout << "Can not open HIST_FILE.DATA" << endl;
		}
		in.close();
		//		cout << hist_file << endl;
		sim->hm_->SetHMTarget(hist_file.c_str());
		//		cout << "Set history matching target..." << endl;
		double Sd = sim->hm_->GetDataMismatch(sim->std_well_);

		// Read number of data points
		double Nd = 1.0;
		ifstream ifs;
		ifs.open(hist_file.c_str());
		ifs >> Nd;
		ifs.close();

		// Calculate model mismatch of xi
		double Sm_xi = 0.0;
		for(int i=0; i<dim_opca_;i++)
			Sm_xi += (xi[i]-xi_uc[i])*(xi[i]-xi_uc[i]);
		// Calculate model mismatch of kr parameters
		double Sm_kr = 0.0;
		for(int i=0; i<dim_kr_;  i++)
			Sm_kr += pow((kr[i]-kr_avg[i])/kr_std[i]-kr_uc[i],2);
		// Calculate normalized model mismatch + data mismatch
		double S = 0.5*(Sd+Sm_kr+Sm_xi)/Nd;
#ifdef DEBUG
		SaveData("Sd.debug",1,&(Sd));
		SaveData("Sm_xi.debug",1,&(Sm_xi));
		SaveData("Sm_kr.debug",1,&(Sm_kr));
#endif
		x.set_bb_output  ( 0 , S); // objective value
		count_eval = true; // count a black-box evaluation

		delete sim;
		return true;       // the evaluation succeeded
	}
public:
	OPCA_BIMODAL* opca_bm_;
	int 		  dim_opca_;
	int 	      dim_kr_;
	vector<double> xi_uc;
	vector<double> kr_uc;
};







/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

	// display:
	Display out ( std::cout );
	out.precision ( DISPLAY_PRECISION_STD );

	try {
		cout << "Begins" << endl;
		// NOMAD initializations:
		begin ( argc , argv );
		cout << "Begined" << endl;
		// parameters creation:
		Parameters p ( out );
		cout << "Initialize parameter" << endl;
		int dim_opca = 70;
		int dim_kr   = 6;
		int dim     = dim_opca + dim_kr;
		p.set_DIMENSION (dim);             // number of variables
		cout << "Set dimension" << endl;
		vector<bb_output_type> bbot (1); // definition of
		bbot[0] = OBJ;                   // output types
		p.set_BB_OUTPUT_TYPE ( bbot );

		//    p.set_DISPLAY_ALL_EVAL(true);   // displays all evaluations.
		p.set_DISPLAY_STATS ( "bbe obj" );
		cout << "Set X0" << endl;
#ifdef DEBUG
		cout << "Set X0" << endl;
#endif
		p.set_X0 ("ui_start.dat");  // starting point
#ifdef DEBUG
		cout << "Set X0 finished" << endl;
#endif
		p.set_LOWER_BOUND ( Point ( dim , 0.0001 ) ); // all var. >= -6
		p.set_UPPER_BOUND ( Point ( dim , 0.9999 ) ); // all var. >= -6

		p.set_INITIAL_MESH_SIZE(0.2);
		p.set_SNAP_TO_BOUNDS(0);
		p.set_DIRECTION_TYPE(NOMAD::ORTHO_2N);

		p.set_MODEL_SEARCH(0);
		p.set_ASYNCHRONOUS(0);
		p.set_MODEL_EVAL_SORT(0);
		p.set_MODEL_SEARCH_OPTIMISTIC(0);
		p.set_OPPORTUNISTIC_EVAL(0);
		p.set_SPECULATIVE_SEARCH(0);
#ifdef DEBUG
		p.set_MAX_BB_EVAL(1);
#endif
#ifndef DEBUG
		p.set_MAX_ITERATIONS (100);     // the algorithm terminates after
#endif
		// 100 black-box evaluations
		p.set_DISPLAY_DEGREE(2);
#ifndef DEBUG
		p.set_SOLUTION_FILE("solution.txt");
		p.set_STATS_FILE("stats.txt","eval bbe obj sol");
#endif
		p.set_ADD_SEED_TO_FILE_NAMES(0);

		// parameters validation:
		p.check();
#ifdef DEBUG
		cout << "Generate OPCA Model" << endl;
#endif
		// custom evaluator creation:
		My_Evaluator ev   ( p );
		ev.opca_bm_  = GenerateOPCAModel();
		ev.dim_opca_ = dim_opca;
		ev.dim_kr_   = dim_kr;

		// read unconditional realizations:
		ifstream ifs;
		ifs.open("UC_FILE.DATA");
		string xi_uc_file, kr_uc_file, temp_str;
		if(ifs.is_open()){
			ifs >> temp_str >> xi_uc_file;
			ifs >> temp_str >> kr_uc_file;
		}
#ifdef DEBUG
		cout << xi_uc_file << endl;
		cout << kr_uc_file << endl;
#endif
		ifs.close();
		double temp_data;
		// Read xi_uc
		ifs.open(xi_uc_file.c_str());
		for(int i = 0; i < ev.dim_opca_; i++){
			ifs >> temp_data;
			ev.xi_uc.push_back(temp_data);
		}
		ifs.close();
		// Read kr_uc
		ifs.open(kr_uc_file.c_str());
		for(int i = 0; i < ev.dim_kr_; i++){
			ifs >> temp_data;
			ev.kr_uc.push_back(temp_data);
		}
		ifs.close();

		// algorithm creation and execution:
		Mads mads ( p , &ev );
		mads.run();
	}
	catch ( exception & e ) {
		cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
	}

	Slave::stop_slaves ( out );
	end();

	return EXIT_SUCCESS;
}
