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
		Evaluator           ( p                                 ) ,
		_mesh_update_basis  ( p.get_mesh_update_basis().value() ) ,
		_initial_mesh_index ( p.get_initial_mesh_index()        ) ,
		_mesh_index         ( _initial_mesh_index               ) ,
		_initial_mesh_size  ( p.get_initial_mesh_size()         )
{
		opca_bm_  = new OPCA_BIMODAL(1,1,1,1);
		dim_opca_ = 1;
		dim_kr_   = 1;
}

	~My_Evaluator ( void ) {
		delete opca_bm_;
	}

	int get_mesh_index ( void ) const { return _mesh_index; }

	void get_mesh_size ( Point & mesh_size ) const
	{
		Mesh::get_delta_m ( mesh_size           ,
				_initial_mesh_size  ,
				_mesh_update_basis  ,
				_initial_mesh_index ,
				_mesh_index           );
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

		_mesh_index = Mesh::get_mesh_index();

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

		for(int i=0; i<dim_kr_; i++){
			temp = x[i+dim_opca_].value();
			temp = TransformUniform2Normal_2(temp,kr_avg[i],kr_std[i]);
			kr.push_back(temp);
		}
		// Boundary check
		kr[0] = kr[0]<0?0.01:kr[0];
		kr[1] = kr[1]<0?0.01:kr[1];
		kr[2] = kr[2]<1?1.0:kr[2];
		kr[3] = kr[3]<1?1.0:kr[3];
		kr[4] = kr[4]<0?0.01:kr[4];
		kr[4] = kr[4]>1?1:kr[4];
		if(kr[0]+kr[1]>1.0){
			x.set_bb_output  ( 0 , 1e20); // objective value
			count_eval = true; // count a black-box evaluation
			cout << "Sum of krw0 and kro0 bigger than 1." << endl;
			return true;       // the evaluation succeeded
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
#if defined(DEBUG) || defined(PRED)
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
		SaveData("Sd.debug"    , 1         ,&(Sd));
		SaveData("Sm_xi.debug" , 1         ,&(Sm_xi));
		SaveData("Sm_kr.debug" , 1         ,&(Sm_kr));
		SaveData("xi_uc.debug" , dim_opca_ ,&(xi_uc[0]));
		SaveData("kr_uc.debug" , dim_kr_   ,&(kr_uc[0]));
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
	vector<double> kr_avg;
	vector<double> kr_std;
	double _mesh_update_basis;
	int    _initial_mesh_index;
	int    _mesh_index;
	Point  _initial_mesh_size;
};







/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

	// display:
	Display out ( std::cout );
	out.precision ( DISPLAY_PRECISION_STD );
	int starting_level = 0;
	if(argc >=2)
		starting_level= char2num(argv[1]);

	try {
		//		cout << "Begins" << endl;
		// NOMAD initializations:
		begin ( argc , argv );
		//		cout << "Begined" << endl;
		// parameters creation:
		Parameters p ( out );
		//		cout << "Initialize parameter" << endl;
		int dim_opca = 70;
		int dim_kr   = 6;
		int dim     = dim_opca + dim_kr;
		p.set_DIMENSION (dim);             // number of variables
		//		cout << "Set dimension" << endl;
		vector<bb_output_type> bbot (1); // definition of
		bbot[0] = OBJ;                   // output types
		p.set_BB_OUTPUT_TYPE ( bbot );

		//    p.set_DISPLAY_ALL_EVAL(true);   // displays all evaluations.
		p.set_DISPLAY_STATS ( "bbe obj" );

		p.set_LOWER_BOUND ( Point ( dim , 0.0001 ) ); // all var. >= -6
		p.set_UPPER_BOUND ( Point ( dim , 0.9999 ) ); // all var. >= -6

		p.set_INITIAL_MESH_SIZE(0.2);
		p.set_SNAP_TO_BOUNDS(0);
		p.set_DIRECTION_TYPE(NOMAD::ORTHO_2N);
		p.set_MESH_COARSENING_EXPONENT(0);


		p.set_MODEL_SEARCH(0);
		p.set_ASYNCHRONOUS(0);
		p.set_MODEL_EVAL_SORT(0);
		p.set_MODEL_SEARCH_OPTIMISTIC(0);
		p.set_OPPORTUNISTIC_EVAL(0);
		p.set_SPECULATIVE_SEARCH(0);
#if defined(DEBUG) || defined(PRED)
		p.set_MAX_BB_EVAL(2);
#endif
#if !defined(DEBUG) && !defined(PRED)
		p.set_MAX_ITERATIONS (80);     // the algorithm terminates after
#endif
		// 100 black-box evaluations
		p.set_DISPLAY_DEGREE(2);
#if!defined(DEBUG) && !defined(PRED)
		if(starting_level>0){
			p.set_SOLUTION_FILE("solution1.txt");
			p.set_STATS_FILE("stats1.txt","eval bbe obj sol poll_size mesh_size");
		}

#endif
#ifdef DEBUG
		p.set_SOLUTION_FILE("dbg_solution1.txt");
		p.set_STATS_FILE("dbg_stats1.txt","eval bbe obj sol poll_size mesh_size");
#endif
		p.set_ADD_SEED_TO_FILE_NAMES(0);

		// parameters validation:
		double temp_data;
		vector<double> x0,xi_uc,kr_uc;
		// read unconditional realizations:
		ifstream ifs;
		ifs.open("UC_FILE.DATA");
		string xi_uc_file, kr_uc_file, temp_str;
		if(ifs.is_open()){
			ifs >> temp_str >> xi_uc_file;
			ifs >> temp_str >> kr_uc_file;
		}
		ifs.close();
#ifdef DEBUG
		cout << xi_uc_file << endl;
		cout << kr_uc_file << endl;
#endif
		ifs.open(xi_uc_file.c_str());
		for(int i = 0; i < dim_opca; i++){
			ifs >> temp_data;
			xi_uc.push_back(temp_data);
			x0.push_back(temp_data);
		}
		ifs.close();
		// Read kr_uc
		ifs.open(kr_uc_file.c_str());
		for(int i = 0; i < dim_kr; i++){
			ifs >> temp_data;
			kr_uc.push_back(temp_data);
			x0.push_back(temp_data);
		}
		ifs.close();
		// Read kr_avg and kr_std
		vector<double> kr_avg, kr_std;
		ifs.open("KR_PRIOR.DATA");
		for(int i = 0; i < dim_kr; i++){
			ifs >> temp_data;
			kr_avg.push_back(temp_data);
		}
		for(int i = 0; i < dim_kr; i++){
			ifs >> temp_data;
			kr_std.push_back(temp_data);
		}
		ifs.close();
#ifdef DEBUG
		SaveData("x0.debug",dim,&(x0[0]));
#endif
		TranformNormal2Uniform(dim,x0);
		CheckBound(0.9999,0.0001,dim,x0);
#ifdef DEBUG
		SaveData("u0.debug",dim,&(x0[0]));
#endif
		SaveData("ui_start.dat",dim,&x0[0]);
#ifndef PRED
		p.set_X0 ("ui_start.dat");  // starting point
#endif
#ifdef PRED
		p.set_X0 ("solution1.txt");
#endif
		p.set_FIXED_VARIABLE(70); // Fix swi
		p.set_FIXED_VARIABLE(71); // Fix sor
		p.set_FIXED_VARIABLE(75); // Fix kro_star
		p.check();
#ifdef DEBUG
		cout << "Generate OPCA Model" << endl;
#endif
		// custom evaluator creation:
		My_Evaluator ev   ( p );
		ev.opca_bm_  = GenerateOPCAModel();
		ev.dim_opca_ = dim_opca;
		ev.dim_kr_   = dim_kr;
		ev.xi_uc     = xi_uc;
		ev.kr_uc     = kr_uc;
		ev.kr_avg    = kr_avg;
		ev.kr_std    = kr_std;
		// algorithm creation and execution:
		Mads mads ( p , &ev );
	    // best solutions:
		// successive runs:
//		for ( int i = starting_level ; i < 2; ++i ) {

			// not for the first run:
			if ( starting_level > 0 )
			{
				// new starting points:
				p.reset_X0();
				/*
				string init_file = "solution";
				string sln_file  = "solution";
				string num;
				stringstream ss;
				ss << i;
				ss >> num;
				sln_file = sln_file + num + ".txt";
				*/
				p.set_X0 ( "solution1.txt" );
				p.set_SOLUTION_FILE("solution2.txt");
				p.set_STATS_FILE("stats2.txt","eval bbe obj sol poll_size mesh_size");
				// initial mesh:
				p.set_INITIAL_MESH_INDEX ( ev.get_mesh_index() );
				Point initial_mesh_size;
				ev.get_mesh_size ( initial_mesh_size );

				p.set_INITIAL_MESH_SIZE ( initial_mesh_size );
				for(int j=0;j<70;j++)
					p.set_FIXED_VARIABLE(j); // Fix O-PCA variable

				p.set_MAX_ITERATIONS (20);     // the algorithm terminates after
				// parameters validation:
				p.check();

				// reset the Mads object:
				mads.reset ( true , true );
			}
			// the run:
			mads.run();
			//cout << "run #" << i << endl;
//		}
	}
	catch ( exception & e ) {
		cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
	}

	Slave::stop_slaves ( out );
	end();

	return EXIT_SUCCESS;
}
