/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <iostream>

#include <dpsim-models/SimPowerComp.h>
#include <dpsim-models/Solver/MNAInterface.h>
#include <dpsim-models/Solver/DAEInterface.h>
#include <dpsim-models/Base/Base_Ph1_Resistor.h>

namespace CPS {
namespace EMT {
namespace Ph1 {
	/// EMT Resistor
	class Resistor :
		public Base::Ph1::Resistor,
		public MNAInterface,
		public DAEInterface,
		public SimPowerComp<Real>,
		public SharedFactory<Resistor> {
	protected:
	public:
		/// Defines UID, name, component parameters and logging level
		Resistor(String uid, String name, Logger::Level logLevel = Logger::Level::off);
		/// Defines name, component parameters and logging level
		Resistor(String name, Logger::Level logLevel = Logger::Level::off)
			: Resistor(name, name, logLevel) { }

		SimPowerComp<Real>::Ptr clone(String name);

		// #### General ####
		/// Initializes component from power flow data
		void initializeFromNodesAndTerminals(Real frequency);

		// #### MNA section ####
		/// Initializes internal variables of the component
		void mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftSideVector);
		/// Stamps system matrix
		void mnaApplySystemMatrixStamp(Matrix& systemMatrix);
		/// Stamps right side (source) vector
		void mnaApplyRightSideVectorStamp(Matrix& rightVector) { }
		/// Update interface voltage from MNA system result
		void mnaUpdateVoltage(const Matrix& leftVector);
		/// Update interface current from MNA system result
		void mnaUpdateCurrent(const Matrix& leftVector);

		class MnaPostStep : public Task {
		public:
			MnaPostStep(Resistor& resistor, Attribute<Matrix>::Ptr leftSideVector) :
				Task(**resistor.mName + ".MnaPostStep"), mResistor(resistor), mLeftVector(leftSideVector)
			{
				mAttributeDependencies.push_back(mLeftVector);
				mModifiedAttributes.push_back(mResistor.attribute("v_intf"));
				mModifiedAttributes.push_back(mResistor.attribute("i_intf"));
			}

			void execute(Real time, Int timeStepCount);

		private:
			Resistor& mResistor;
			Attribute<Matrix>::Ptr mLeftVector;
		};

		// #### DAE Section ####
		///
		void daeInitialize(double time, double state[], double dstate_dt[], 
			double absoluteTolerances[], double stateVarTypes[], int& counter) override;
		///Residual Function for DAE Solver
		void daeResidual(double time, const double state[], const double dstate_dt[], 
			double resid[], std::vector<int>& off) override;
		/// Calculation of jacobian
		void daeJacobian(double current_time, const double state[], const double dstate_dt[], 
			SUNMatrix jacobian, double cj, std::vector<int>& off) override;
		///
		void daePostStep(double Nexttime, const double state[], 
			const double dstate_dt[], int& counter) override;
		///
		int getNumberOfStateVariables() override {return 0;}

	};
}
}
}
