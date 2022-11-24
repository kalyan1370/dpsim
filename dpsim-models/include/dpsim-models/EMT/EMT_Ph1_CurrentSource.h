/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/SimPowerComp.h>
#include <dpsim-models/Solver/MNAInterface.h>
#include <dpsim-models/Solver/DAEInterface.h>
#include <dpsim-models/Base/Base_Ph1_CurrentSource.h>

namespace CPS {
namespace EMT {
namespace Ph1 {
	/// \brief Ideal current source model
	///
	/// A positive current is flowing out of
	/// node1 and into node2.
	class CurrentSource :
		public MNAInterface,
		public DAEInterface,
		public SimPowerComp<Real>,
		public SharedFactory<CurrentSource> {
	public:
		const Attribute<Complex>::Ptr mCurrentRef;
		const Attribute<Real>::Ptr mSrcFreq;
		
		/// Defines UID, name and logging level
		CurrentSource(String uid, String name,
			Logger::Level logLevel = Logger::Level::off);
		///
		CurrentSource(String name, Logger::Level logLevel = Logger::Level::off)
			: CurrentSource(name, name, logLevel) { }

		SimPowerComp<Real>::Ptr clone(String name);

		void setParameters(Complex currentRef, Real srcFreq = -1);
		// #### General ####
		/// Initializes component from power flow data
		void initializeFromNodesAndTerminals(Real frequency) { }

		// #### MNA section ####
		/// Initializes internal variables of the component
		void mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector);
		/// Stamps system matrix
		void mnaApplySystemMatrixStamp(Matrix& systemMatrix) { }
		/// Stamps right side (source) vector
		void mnaApplyRightSideVectorStamp(Matrix& rightVector);
		///
		void mnaUpdateVoltage(const Matrix& leftVector);

		void updateState(Real time);

		class MnaPreStep : public Task {
		public:
			MnaPreStep(CurrentSource& currentSource) :
				Task(**currentSource.mName + ".MnaPreStep"), mCurrentSource(currentSource)
			{
				mAttributeDependencies.push_back(currentSource.attribute("I_ref"));
				mModifiedAttributes.push_back(currentSource.attribute("right_vector"));
				mModifiedAttributes.push_back(currentSource.attribute("i_intf"));
			}

			void execute(Real time, Int timeStepCount);

		private:
			CurrentSource& mCurrentSource;
		};

		class MnaPostStep : public Task {
		public:
			MnaPostStep(CurrentSource& currentSource, Attribute<Matrix>::Ptr leftSideVector) :
				Task(**currentSource.mName + ".MnaPostStep"), mCurrentSource(currentSource), mLeftVector(leftSideVector)
			{
				mAttributeDependencies.push_back(mLeftVector);
				mModifiedAttributes.push_back(mCurrentSource.attribute("v_intf"));
			}

			void execute(Real time, Int timeStepCount);
		private:
			CurrentSource& mCurrentSource;
			Attribute<Matrix>::Ptr mLeftVector;
		};

		///
		void daeInitialize(double time, double state[], double dstate_dt[],
			double absoluteTolerances[], double stateVarTypes[], int& offset) override;
		/// Residual function for DAE Solver
		void daeResidual(double time, const double state[], const double dstate_dt[], 
			double resid[], std::vector<int>& off) override;
		/// Calculation of jacobian
		void daeJacobian(double current_time, const double state[], const double dstate_dt[], 
			SUNMatrix jacobian, double cj, std::vector<int>& off) override {};
		///
		void daePostStep(double Nexttime, const double state[], 
			const double dstate_dt[], int& offset) override;
		///
		int getNumberOfStateVariables() override {return 0;}
	};
}
}
}
