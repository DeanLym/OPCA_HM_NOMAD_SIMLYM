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
#include <mpi.h>
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
		dim_ = 1;
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
		for(int i=0;i<dim_;i++){
			temp = x[i].value();
			xi.push_back(temp);
		}
		TranformUniform2Normal(dim_,xi);

		double Nc = 3600;
		vector<double> m;
		m.resize(Nc);
		bool opca_flag;
		opca_flag = opca_bm_->GenerateOPCARealization(&(xi[0]) , &(m[0]));
		if(!opca_flag) // If O-PCA returns false, the solution is wrong
		{
			string log_file("opca_err_log.");
			log_file += num2str(_rank);
			LogOpcaError(log_file.c_str(), xi, dim_);
			return false;       // the evaluation failed
		}
		vector<double> perm;
		perm.resize(Nc);
		GeneratePerm(Nc, &(m[0]), &(perm[0]));

#ifdef DEBUG
		SaveData("xi.debug",dim_,&(xi[0]));
		SaveData("m.debug",Nc,&(m[0]));
		SaveData("perm.debug",Nc,&(perm[0]));
#endif

		SimCtrl* sim = GetSimulationModel(&perm[0]);

		sim->display_level_ = 0;
		//cout << "Begin to run simulation." <<endl;
		sim->RunSim();
		//cout << "Finished running simulation." << endl;
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
		//cout << "Finished setting history matching target." << endl;
		//		cout << "Set history matching target..." << endl;
		double Sd = sim->hm_->GetDataMismatch(sim->std_well_);
		//cout << "Finished calculating data mismatch" << endl;
		// Read number of data points
		double Nd = 1.0;
		ifstream ifs;
		ifs.open(hist_file.c_str());
		ifs >> Nd;
		ifs.close();

		// Calculate model mismatch of xi
		double Sm_xi = 0.0;
		for(int i=0; i<dim_;i++)
			Sm_xi += (xi[i]-xi_uc[i])*(xi[i]-xi_uc[i]);
		// Calculate normalized model mismatch + data mismatch
		double S = 0.5*(Sd+Sm_xi)/Nd;
#ifdef DEBUG
		SaveData("Sd.debug"    , 1         ,&(Sd));
		SaveData("Sm_xi.debug" , 1         ,&(Sm_xi));
		SaveData("xi_uc.debug" , dim_ ,&(xi_uc[0]));
#endif
		x.set_bb_output  ( 0 , S); // objective value
		count_eval = true; // count a black-box evaluation
		//cout << "Delete sim object." << endl;
		delete sim;
		//cout << "Finished evaluation!" << endl;
		return true;       // the evaluation succeeded
	}
public:
	OPCA_BIMODAL* opca_bm_;
	int 		  dim_;
	vector<double> xi_uc;
	double _mesh_update_basis;
	int    _initial_mesh_index;
	int    _mesh_index;
	Point  _initial_mesh_size;
	int    _rank;
};







/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

	// display:
	Display out ( std::cout );
	out.precision ( DISPLAY_PRECISION_STD );
	int dim = 80;
	int current_level = 0; //Default level is 0
	int force_init_zero = dim;
	if(argc > 1)
		current_level = char2num(argv[1]);
	if(argc > 2)
		force_init_zero = char2num(argv[2]);
	try {
		// NOMAD initializations:
		begin ( argc , argv );
		// parameters creation:
		int size, rank;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		//cout << "This is proc #" << world_rank << "in all " << world_size << " procs" << endl;

		Parameters p ( out );

		//===========================================================//
		// Read parameters from file
		double temp_data;
		vector<double> x0,xi_uc,kr_uc;
		// read unconditional realizations:
		int num_level = 4;
		vector<int> ml_dims;
		vector<int> ml_iters;
		ifstream ifs;
		// Read multilevel dimensions
		ifs.open("ML.DATA");
		ifs >> num_level;
		ml_dims.resize(num_level);
		ml_iters.resize(num_level);
		for(int i=0;i<num_level;i++){
			ifs>>ml_dims[i];
			ifs>>ml_iters[i];
		}
		ifs.close();

		// Read uc file path
		ifs.open("UC_FILE.DATA");
		string xi_uc_file, temp_str;
		if(ifs.is_open()){
			ifs >> temp_str >> xi_uc_file;
		}
		ifs.close();
#ifdef DEBUG
		//cout << xi_uc_file << endl;
#endif
		ifs.open(xi_uc_file.c_str());
		for(int i = 0; i < dim; i++){
			ifs >> temp_data;
			xi_uc.push_back(temp_data);
			if(i<force_init_zero)
			    x0.push_back(temp_data);
			else
				x0.push_back(0);
		}
		ifs.close();

		//===========================================================//

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
		p.set_MAX_BB_EVAL(1);
		//p.set_MAX_ITERATIONS (ml_iters[0]);     // the algorithm terminates after
#endif
#if !defined(DEBUG) && !defined(PRED)
		p.set_MAX_ITERATIONS (ml_iters[0]);     // the algorithm terminates after
#endif
		// 100 black-box evaluations
		p.set_DISPLAY_DEGREE(2);
		if(current_level==0){
#if!defined(DEBUG) && !defined(PRED)
		p.set_SOLUTION_FILE("solution1.txt");
		p.set_STATS_FILE("stats1.txt","eval bbe obj sol poll_size mesh_size");
#endif
		}
		p.set_ADD_SEED_TO_FILE_NAMES(0);


#ifdef DEBUG
		SaveData("x0.debug",dim,&(x0[0]));
#endif

		TranformNormal2Uniform(dim,x0);
		CheckBound(0.9999,0.0001,dim,x0);

#ifdef DEBUG
		SaveData("u0.debug",dim,&(x0[0]));
#endif

#ifndef DEBUG
		if( rank == 0)
			SaveData("ui_start1.dat",dim,&x0[0]);
#endif

		MPI_Barrier( MPI_COMM_WORLD );
#ifndef PRED
		p.set_X0 ("ui_start1.dat");  // starting point
#endif
#ifdef PRED
		string pred_init = "solution" + num2str(num_level) + ".txt";
		p.set_X0 (pred_init.c_str());
#endif
		for(int i=ml_dims[0];i<dim;i++)
			p.set_FIXED_VARIABLE(i); // Fix swi
		p.check();
#ifdef DEBUG
		//cout << "Generate OPCA Model" << endl;
#endif
		// custom evaluator creation:
		My_Evaluator ev   ( p );
		ev.opca_bm_  = GenerateOPCAModel();
		ev.dim_ = dim;
		ev.xi_uc     = xi_uc;
		ev._rank = rank;
		// algorithm creation and execution:
		Mads mads ( p , &ev );

	    // best solutions:
		// successive runs:
		int i = current_level;
		if ( i > 0 )
		{
			// new starting points:
			p.reset_X0();
			string init_file   = "solution" + num2str(i) + ".txt";
			string sln_file    = "solution" + num2str(i+1)   + ".txt";
			string stats_file  = "stats"    + num2str(i+1)   + ".txt";
			MPI_Barrier( MPI_COMM_WORLD ) ;
			p.set_X0 ( init_file.c_str() );

			p.set_SOLUTION_FILE(sln_file.c_str());
			p.set_STATS_FILE(stats_file.c_str(),"eval bbe obj sol poll_size mesh_size");
			// initial mesh:
			//p.set_INITIAL_MESH_INDEX ( ev.get_mesh_index() );
			//Point initial_mesh_size;
			//ev.get_mesh_size ( initial_mesh_size );
			//p.set_INITIAL_MESH_SIZE ( initial_mesh_size );
			for(int j=ml_dims[i-1];j<ml_dims[i];j++)
				p.set_FREE_VARIABLE(j);
			for(int j=0;j<ml_dims[i-1];j++)
				p.set_FIXED_VARIABLE(j); // Fix O-PCA variable
			for(int j=ml_dims[i];j<dim;j++)
				p.set_FIXED_VARIABLE(j);
			p.set_MAX_ITERATIONS (ml_iters[i]);     // the algorithm terminates after
			// parameters validation:
			p.check();

			// reset the Mads object:
			mads.reset ( true , true );
		}
		//cout << "Current level is :" << current_level << endl;
		//cout << "Begin to run MADS." << endl;
		// the run:
		mads.run();
		//cout << "run #" << i << endl;
	}
	catch ( exception & e ) {
		cerr << "\nNOMAD has been interrupted (" << e.what() << ")\n\n";
	}

	Slave::stop_slaves ( out );
	end();

	return EXIT_SUCCESS;
}
