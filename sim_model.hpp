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
	// sim file master
	SimCtrl* sim;
	sim = new SimCtrl;
	sim->display_level_ = 0;
	// input and initialize grid
	sim->grid_ = new CCartGrid(60,60,1);
	sim->grid_->InputDx(82.02);
	sim->grid_->InputDy(82.02);
	sim->grid_->InputDz(32.81);
	sim->grid_->InputTops(11466.54);
	sim->grid_->InputKx(kx);
	sim->grid_->InputKy(kx);
	sim->grid_->InputKz(kx);
	sim->grid_->InputPoro(0.2);
	sim->InitializeGrid();
	//sim->grid_->OutputConnList();

	// initialize PVT
	sim->pvt_ = new CPVT;
	sim->pvt_->SetDensity(49.1, 64.79, 0.06);
	sim->pvt_->set_pvdo_table(42,"PVDO.DAT");
	sim->pvt_->set_pvtw_table("PVTW.DAT");

	// initialize SAT
//	sim->sat_ = new CSAT_TABLE;
//	sim->sat_->SetSWOF(9,"SCAL.DAT");

	sim->sat_ = new CSAT_COREY(0.1, 0.2 ,2.0 ,2.0 ,0.5 ,0.9 );

	// initialize SCH
	sim->sch_ = new CSchedule;
	sim->sch_->SetTEnd(540.0);
	sim->sch_->SetTCurrent(0.0);
	sim->sch_->SetTNext(sim->sch_->GetDt() + sim->sch_->GetTCurrent());
	sim->sch_->SetdTmax(30.0);
	vector<double> report_time;
	report_time.push_back(180.0); report_time.push_back(360.0);report_time.push_back(540.0);
	sim->sch_->SetReportTime(3, &report_time[0]);

	// initialize State
	sim->InitializeState();
	sim->SetInitPres(11482.94,5076.33);
	sim->SetInitSat(0.1, 0.0);


	// initialize Well
	CStandardWell *p_1 = new CSTDProdWell("PROD-1", sim->grid_->GetIndex(4,19,0));
	p_1->set_r(1);
	p_1->set_ctrl_mode(CStandardWell::CBHP);
	p_1->set_TL_BHP(2900.76);
	p_1->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_1);
	CStandardWell *p_2 = new CSTDProdWell("PROD-2", sim->grid_->GetIndex(9,9,0));
	p_2->set_r(1);
	p_2->set_ctrl_mode(CStandardWell::CBHP);
	p_2->set_TL_BHP(2900.76);
	p_2->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_2);
	CStandardWell *p_3 = new CSTDProdWell("PROD-3", sim->grid_->GetIndex(54,54,0));
	p_3->set_r(1);
	p_3->set_ctrl_mode(CStandardWell::CBHP);
	p_3->set_TL_BHP(2900.76);
	p_3->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(p_3);
	//=====================================
	CStandardWell *inj_1 = new CSTDWInjWell("INJ-1", sim->grid_->GetIndex(29,24, 0));
	inj_1->set_r(1);
	inj_1->set_ctrl_mode(CStandardWell::CBHP);
	inj_1->set_TL_BHP(7251.9);
	inj_1->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(inj_1);

	CStandardWell *inj_2 = new CSTDWInjWell("INJ-2", sim->grid_->GetIndex(29,44, 0));
	inj_2->set_r(1);
	inj_2->set_ctrl_mode(CStandardWell::CBHP);
	inj_2->set_TL_BHP(8702.28);
	inj_2->CalWellIndex(sim->grid_);
	sim->std_well_.push_back(inj_2);

	// initialize Solver
	sim->InitializeSolver();
	return sim;

}


#endif /* SIM_MODEL_HPP_ */
