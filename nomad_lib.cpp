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
#include <armadillo>


using namespace arma;

/*----------------------------------------*/
/*               The problem              */
/*----------------------------------------*/
class My_Evaluator : public Evaluator {
public:
	My_Evaluator  ( const Parameters & p ) :
		Evaluator ( p ) {
			opca_bm_ = new OPCA_BIMODAL(1,1,1,1);
			dim_ = 1;
	}

	~My_Evaluator ( void ) {
		delete opca_bm_;
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
		opca_bm_->GenerateOPCARealization(&(xi[0]) , &(m[0]));
		vector<double> perm;
		perm.resize(Nc);
		GeneratePerm(Nc, &(m[0]), &(perm[0]));

		SimCtrl* sim = GetSimulationModel(&perm[0]);

		sim->display_level_ = 0;
		sim->RunSim();
		sim->hm_ = new CHistoryMatching;
		sim->hm_->SetHMTarget("HIST.TXT");
		double Sd = sim->hm_->GetDataMismatch(sim->std_well_);
		double Nd = 9;
		double Sm = 0.0;
		for(int i=0;i<dim_;i++)
			Sm += xi[i]*xi[i];
		double S = 0.5*(Sd+Sm)/Nd;
		x.set_bb_output  ( 0 , S); // objective value
		count_eval = true; // count a black-box evaluation

		delete sim;
		return true;       // the evaluation succeeded
	}
public:
	OPCA_BIMODAL* opca_bm_;
	int dim_;
};




/*------------------------------------------*/
/*            NOMAD main function           */
/*------------------------------------------*/
int main ( int argc , char ** argv ) {

	// display:
	Display out ( std::cout );
	out.precision ( DISPLAY_PRECISION_STD );

	try {

		// NOMAD initializations:
		begin ( argc , argv );

		// parameters creation:
		Parameters p ( out );

		int dim = 70;
		p.set_DIMENSION (dim);             // number of variables

		vector<bb_output_type> bbot (1); // definition of
		bbot[0] = OBJ;                   // output types
		p.set_BB_OUTPUT_TYPE ( bbot );

		//    p.set_DISPLAY_ALL_EVAL(true);   // displays all evaluations.
		p.set_DISPLAY_STATS ( "bbe obj" );

		p.set_X0 ("ui_start.dat");  // starting point

		p.set_LOWER_BOUND ( Point ( 70 , 0.0001 ) ); // all var. >= -6
		p.set_UPPER_BOUND ( Point ( 70 , 0.9999 ) ); // all var. >= -6

		p.set_INITIAL_MESH_SIZE(0.2);
		p.set_SNAP_TO_BOUNDS(0);
		p.set_DIRECTION_TYPE(NOMAD::ORTHO_2N);
		p.set_ASYNCHRONOUS(0);
		p.set_MODEL_EVAL_SORT(0);
		p.set_MODEL_SEARCH_OPTIMISTIC(0);
		p.set_OPPORTUNISTIC_EVAL(0);
		p.set_SPECULATIVE_SEARCH(0);
		p.set_MAX_ITERATIONS (100);     // the algorithm terminates after
		// 100 black-box evaluations
		p.set_DISPLAY_DEGREE(2);
		p.set_SOLUTION_FILE("solution.txt");
		p.set_STATS_FILE("stats.txt","eval bbe obj sol");
		p.set_ADD_SEED_TO_FILE_NAMES(0);

		// parameters validation:
		p.check();

		// custom evaluator creation:
		My_Evaluator ev   ( p );
		ev.opca_bm_ = GenerateOPCAModel();
		ev.dim_ = dim;

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
