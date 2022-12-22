/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/EMT/EMT_Ph3_ReducedOrderSynchronGeneratorVBR.h>
#include <dpsim-models/Solver/DAEInterface.h>

namespace CPS {
namespace EMT {
namespace Ph3 {
	/// @brief Voltage-Behind-Reactance (VBR) implementation
	/// of 4th order synchronous generator model
	class SynchronGenerator4OrderVBR :
		public ReducedOrderSynchronGeneratorVBR,
		public DAEInterface,
		public SharedFactory<SynchronGenerator4OrderVBR> {
	public:
		// ### Model specific elements ###
		/// transient voltage
		const Attribute<Matrix>::Ptr mEdq0_t;
	protected:
		/// history term of VBR
		Matrix mEhs_vbr;
		///
		Real mAd;
		///
		Real mAq;
		/// 
		Real mBd;
		///
		Real mBq;
		///
		Real mCq;

	public:
		///
		SynchronGenerator4OrderVBR(String uid, String name, Logger::Level logLevel = Logger::Level::off);
		///
		SynchronGenerator4OrderVBR(String name, Logger::Level logLevel = Logger::Level::off);
		///
		SimPowerComp<Real>::Ptr clone(String name);

		// #### General Functions ####
		///
		void specificInitialization();
		///
		void calculateAuxiliarConstants();
		///
		void stepInPerUnit();
		/// Setter 4th order parameters - extending base class setter by logging
		void setOperationalParametersPerUnit(Real nomPower, 
			Real nomVolt, Real nomFreq, Real H, Real Ld, Real Lq, Real L0,
			Real Ld_t, Real Lq_t, Real Td0_t, Real Tq0_t);

		// #### DAE Section ####
		/// Derivative of the current
		MatrixVar<Real> mIntfDerCurrent;
		///
		void daeInitialize(double time, double state[], double dstate_dt[], 
			double absoluteTolerances[], double stateVarTypes[], int& offset) override;
		/// Residual function for DAE Solver
		void daeResidual(double time, const double state[], const double dstate_dt[], 
			double resid[], std::vector<int>& off) override;
		/// Calculation of jacobian
		void daeJacobian(double current_time, const double state[], const double dstate_dt[], 
			SUNMatrix jacobian, double cj, std::vector<int>& off) override;
		///
		void daePostStep(double Nexttime, const double state[], 
			const double dstate_dt[], int& offset) override;
		///
		int getNumberOfStateVariables() override { return 6; }
	};
}
}
}
