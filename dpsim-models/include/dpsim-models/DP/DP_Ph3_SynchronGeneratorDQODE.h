/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#pragma once

#include <dpsim-models/DP/DP_Ph3_SynchronGeneratorDQ.h>

#ifdef WITH_SUNDIALS
	#include <sundials/sundials_types.h>		/* def. of type 'realtype' */
	#include <dpsim-models/Solver/ODEInterface.h>
	
	// Needed to print the computed to the console
	#if defined(SUNDIALS_EXTENDED_PRECISION)
		#define GSYM "Lg"
		#define ESYM "Le"
		#define FSYM "Lf"
	#else
		#define GSYM "g"
		#define ESYM "e"
		#define FSYM "f"
	#endif
#endif

namespace CPS {
namespace DP {
namespace Ph3 {
	class SynchronGeneratorDQODE :
		public SynchronGeneratorDQ,
		public ODEInterface,
		public SharedFactory<SynchronGeneratorDQODE> {
	public:
		SynchronGeneratorDQODE(String uid, String name, Logger::Level loglevel = Logger::Level::off);
		SynchronGeneratorDQODE(String name, Logger::Level loglevel = Logger::Level::off);

		// #### MNA Section ####
		void mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector);

		class MnaPreStep : public Task {
		public:
			MnaPreStep(SynchronGeneratorDQODE& synGen)
				: Task(**synGen.mName + ".MnaPreStep"), mSynGen(synGen) {
				mAttributeDependencies.push_back(synGen.attribute("ode_post_state"));
				mModifiedAttributes.push_back(synGen.attribute("right_vector"));
				mPrevStepDependencies.push_back(synGen.attribute("v_intf"));
			}

			void execute(Real time, Int timeStepCount);

		private:
			SynchronGeneratorDQODE& mSynGen;
		};

		class ODEPreStep : public Task {
		public:
			ODEPreStep(SynchronGeneratorDQODE& synGen)
			: Task(**synGen.mName + ".ODEPreStep"), mSynGen(synGen) {
				mModifiedAttributes.push_back(synGen.attribute("ode_pre_state"));
				mModifiedAttributes.push_back(synGen.attribute("i_intf"));
			}

			void execute(Real time, Int timeStepCount);

		private:
			SynchronGeneratorDQODE& mSynGen;
		};

	protected:
		// #### ODE Section ####

		/// Defines the ODE-System which shall be solved
		void odeStateSpace(double t, const double y[], double ydot[]); // ODE-Class approach

		/// Jacobian corresponding to the StateSpace system, needed for implicit solve
		void odeJacobian(double t, const double y[], double fy[], double J[],
		                 double tmp1[], double tmp2[], double tmp3[]);

		void odePreStep();
		void odePostStep();

		// ### Variables for ODE-Solver interaction ###
		/// Number of differential variables
		int mDim;
	};
}
}
}
