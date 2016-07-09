/*
 * sim_model.hpp
 *
 *  Created on: Mar 5, 2016
 *      Author: yiminliu
 */

#ifndef SIM_MODEL_HPP_
#define SIM_MODEL_HPP_

#include "CSimCtrl.h"

SimCtrl* GetSimulationModel(double *kx){
	// Model 1 Simulator file
	SimCtrl* sim;
	sim = new SimCtrl;
	sim->display_level_ = 1;
	// input and initialize grid
	sim->grid_ = new CCartGrid(60,60,1);
	sim->grid_->InputDx(164.042);
	sim->grid_->InputDy(164.042);
	sim->grid_->InputDz(32.81);
	sim->grid_->InputTops(25574.15);
	sim->grid_->InputKx(kx);
	sim->grid_->InputKy(kx);
	sim->grid_->InputKz(kx);
	sim->grid_->InputPoro(0.2);
	sim->InitializeGrid();

	// initialize PVT
	sim->pvt_ = new CPVT;
	sim->pvt_->SetDensity(49.1, 64.79, 0.06);
	sim->pvt_->set_pvdo_table(42,"PVDO.DAT");
	sim->pvt_->set_pvtw_table("PVTW.DAT");

	double swi  = 0.10 , sor  = 0.30 ;
	double aw   = 1.5  , ao   = 3.6  ;
	double krw0 = 0.7  , kro0 = 1.0  ;
	sim->sat_   = new CSAT_COREY(swi, sor , aw ,ao ,krw0 ,kro0 );

	// initialize SCH
	sim->sch_ = new CSchedule;

#if defined(DEBUG) || defined(PRED)
	sim->sch_->SetTEnd(2000.0);
#endif
#if !defined(DEBUG) && !defined(PRED)
	sim->sch_->SetTEnd(600.0);     // the algorithm terminates after
#endif

	sim->sch_->SetTCurrent(0.0);
	sim->sch_->SetTNext(sim->sch_->GetDt() + sim->sch_->GetTCurrent());
	sim->sch_->SetdTmax(100.0);
	vector<double> report_time;
	int num_report_time = 1;
#if !defined(DEBUG) && !defined(PRED)
	num_report_time = 6;
#endif
#if defined(DEBUG) || defined(PRED)
	num_report_time = 20;
#endif
	for(int i=0; i<num_report_time ;i++)
		report_time.push_back((i+1)*100);
	sim->sch_->SetReportTime(num_report_time, &report_time[0]);

	// initialize State
	sim->InitializeState();
	sim->SetInitPres(25590.55,4713.735);
	// Set initial saturation equal to swi
	sim->SetInitSat(swi, 0.0);

	// initialize Well
	double bhp_ = 2175.57;
	CStandardWell *p_1 = new CSTDProdWell("PROD-1", sim->grid_->GetIndex(12,6,0));
	p_1->set_r(0.5);
	p_1->set_ctrl_mode(CStandardWell::CBHP);
	p_1->set_TL_BHP(bhp_);
	p_1->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_1);

	CStandardWell *p_2 = new CSTDProdWell("PROD-2", sim->grid_->GetIndex(12,20,0));
	p_2->set_r(0.5);
	p_2->set_ctrl_mode(CStandardWell::CBHP);
	p_2->set_TL_BHP(bhp_);
	p_2->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_2);

	CStandardWell *p_3 = new CSTDProdWell("PROD-3", sim->grid_->GetIndex(12,32,0));
	p_3->set_r(0.5);
	p_3->set_ctrl_mode(CStandardWell::CBHP);
	p_3->set_TL_BHP(bhp_);
	p_3->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_3);

	CStandardWell *p_4 = new CSTDProdWell("PROD-4", sim->grid_->GetIndex(12,50,0));
	p_4->set_r(0.5);
	p_4->set_ctrl_mode(CStandardWell::CBHP);
	p_4->set_TL_BHP(bhp_);
	p_4->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_4);

	CStandardWell *p_5 = new CSTDProdWell("PROD-5", sim->grid_->GetIndex(24,6,0));
	p_5->set_r(0.5);
	p_5->set_ctrl_mode(CStandardWell::CBHP);
	p_5->set_TL_BHP(bhp_);
	p_5->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_5);

	CStandardWell *p_6 = new CSTDProdWell("PROD-6", sim->grid_->GetIndex(24,35,0));
	p_6->set_r(0.5);
	p_6->set_ctrl_mode(CStandardWell::CBHP);
	p_6->set_TL_BHP(bhp_);
	p_6->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_6);

	CStandardWell *p_7 = new CSTDProdWell("PROD-7", sim->grid_->GetIndex(40,6,0));
	p_7->set_r(0.5);
	p_7->set_ctrl_mode(CStandardWell::CBHP);
	p_7->set_TL_BHP(bhp_);
	p_7->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_7);

	CStandardWell *p_8 = new CSTDProdWell("PROD-8", sim->grid_->GetIndex(40,20,0));
	p_8->set_r(0.5);
	p_8->set_ctrl_mode(CStandardWell::CBHP);
	p_8->set_TL_BHP(bhp_);
	p_8->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_8);

	CStandardWell *p_9 = new CSTDProdWell("PROD-9", sim->grid_->GetIndex(40,35,0));
	p_9->set_r(0.5);
	p_9->set_ctrl_mode(CStandardWell::CBHP);
	p_9->set_TL_BHP(bhp_);
	p_9->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_9);

	CStandardWell *p_10 = new CSTDProdWell("PROD-10", sim->grid_->GetIndex(40,52,0));
	p_10->set_r(0.5);
	p_10->set_ctrl_mode(CStandardWell::CBHP);
	p_10->set_TL_BHP(bhp_);
	p_10->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_10);

	CStandardWell *p_11 = new CSTDProdWell("PROD-11", sim->grid_->GetIndex(52,6,0));
	p_11->set_r(0.5);
	p_11->set_ctrl_mode(CStandardWell::CBHP);
	p_11->set_TL_BHP(bhp_);
	p_11->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_11);

	CStandardWell *p_12 = new CSTDProdWell("PROD-12", sim->grid_->GetIndex(52,45,0));
	p_12->set_r(0.5);
	p_12->set_ctrl_mode(CStandardWell::CBHP);
	p_12->set_TL_BHP(bhp_);
	p_12->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_12);

	//=====================================
	double wrat_target_inj = -9447.0, bhp_limit_inj = 14503.8;
	CStandardWell *inj_1 = new CSTDWInjWell("INJ-1", sim->grid_->GetIndex(24, 20, 0));
	inj_1->set_r(0.5);
	inj_1->set_ctrl_mode(CStandardWell::CWRAT);
	inj_1->set_TL_WRAT(wrat_target_inj);
	inj_1->set_TL_BHP(bhp_limit_inj);
	inj_1->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(inj_1);

	CStandardWell *inj_2 = new CSTDWInjWell("INJ-2", sim->grid_->GetIndex(24, 52, 0));
	inj_2->set_r(0.5);
	inj_2->set_ctrl_mode(CStandardWell::CWRAT);
	inj_2->set_TL_WRAT(wrat_target_inj);
	inj_2->set_TL_BHP(bhp_limit_inj);
	inj_2->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(inj_2);

	CStandardWell *inj_3 = new CSTDWInjWell("INJ-3", sim->grid_->GetIndex(52, 25, 0));
	inj_3->set_r(0.5);
	inj_3->set_ctrl_mode(CStandardWell::CWRAT);
	inj_3->set_TL_WRAT(wrat_target_inj);
	inj_3->set_TL_BHP(bhp_limit_inj);
	inj_3->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(inj_3);

	CStandardWell *inj_4 = new CSTDWInjWell("INJ-4", sim->grid_->GetIndex(52, 56, 0));
	inj_4->set_r(0.5);
	inj_4->set_ctrl_mode(CStandardWell::CWRAT);
	inj_4->set_TL_WRAT(wrat_target_inj);
	inj_4->set_TL_BHP(bhp_limit_inj);
	inj_4->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(inj_4);
	// initialize Solver
	sim->InitializeSolver();
	return sim;

}


#endif /* SIM_MODEL_HPP_ */
