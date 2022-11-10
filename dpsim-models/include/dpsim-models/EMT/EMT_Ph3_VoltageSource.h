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
#include <dpsim-models/Signal/SignalGenerator.h>
#include <dpsim-models/Signal/SineWaveGenerator.h>
#include <dpsim-models/Signal/FrequencyRampGenerator.h>
#include <dpsim-models/Signal/CosineFMGenerator.h>

namespace CPS {
	namespace EMT {
		namespace Ph3 {
			/// \brief Ideal Voltage source model
			///
			/// This model uses modified nodal analysis to represent an ideal voltage source.
			/// For a voltage source between nodes j and k, a new variable (current across the voltage source)
			/// is added to the left side vector
			/// as unkown and it is taken into account for the equation of node j as positve and for the equation
			/// of node k as negative. Moreover
			/// a new equation ej - ek = V is added to the problem.
			class VoltageSource :
				public MNAInterface,
				public DAEInterface,
				public SimPowerComp<Real>,
				public SharedFactory<VoltageSource> {
			private:
				///
				CPS::Signal::SignalGenerator::Ptr mSrcSig;
			protected:
				// Updates voltage according to reference phasor and frequency
				void updateVoltage(Real time);
			public:
				const CPS::Attribute<MatrixComp>::Ptr mVoltageRef;
				const CPS::Attribute<Real>::Ptr mSrcFreq;

				/// Defines UID, name and logging level
				VoltageSource(String uid, String name, Logger::Level logLevel = Logger::Level::off);
				///
				VoltageSource(String name, Logger::Level logLevel = Logger::Level::off)
					: VoltageSource(name, name, logLevel) { }

				SimPowerComp<Real>::Ptr clone(String name);
				// #### General ####
				/// Initializes component from power flow data
				void initializeFromNodesAndTerminals(Real frequency);
				/// Setter for reference voltage and frequency with a sine wave generator
				void setParameters(MatrixComp voltageRef, Real srcFreq = 50.0);
				/// Setter for reference signal of type frequency ramp
				void setParameters(MatrixComp voltageRef, Real freqStart, Real rocof, Real timeStart, Real duration, bool smoothRamp = true);
				/// Setter for reference signal of type cosine frequency modulation
				void setParameters(MatrixComp voltageRef, Real modulationFrequency, Real modulationAmplitude, Real baseFrequency = 50.0, bool zigzag = false);

				// #### MNA section ####
				/// Initializes internal variables of the component
				void mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector);
				/// Stamps system matrix
				void mnaApplySystemMatrixStamp(Matrix& systemMatrix);
				/// Stamps right side (source) vector
				void mnaApplyRightSideVectorStamp(Matrix& rightVector);
				/// Returns current through the component
				void mnaUpdateCurrent(const Matrix& leftVector);
				/// MNA pre step operations
				void mnaPreStep(Real time, Int timeStepCount);
				/// MNA post step operations
				void mnaPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector);
				/// Add MNA pre step dependencies
				void mnaAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes);
				/// Add MNA post step dependencies
				void mnaAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector);

				class MnaPreStep : public Task {
				public:
					MnaPreStep(VoltageSource& voltageSource) :
						Task(**voltageSource.mName + ".MnaPreStep"), mVoltageSource(voltageSource) {
							mVoltageSource.mnaAddPreStepDependencies(mPrevStepDependencies, mAttributeDependencies, mModifiedAttributes);
						}
						void execute(Real time, Int timeStepCount) { mVoltageSource.mnaPreStep(time, timeStepCount); };
				private:
					VoltageSource& mVoltageSource;
				};

				class MnaPostStep : public Task {
				public:
					MnaPostStep(VoltageSource& voltageSource, Attribute<Matrix>::Ptr leftVector) :
						Task(**voltageSource.mName + ".MnaPostStep"),
						mVoltageSource(voltageSource), mLeftVector(leftVector) {
							mVoltageSource.mnaAddPostStepDependencies(mPrevStepDependencies, mAttributeDependencies, mModifiedAttributes, mLeftVector);
					}
					void execute(Real time, Int timeStepCount)  { mVoltageSource.mnaPostStep(time, timeStepCount, mLeftVector); };
				private:
					VoltageSource& mVoltageSource;
					Attribute<Matrix>::Ptr mLeftVector;
				};

				// #### DAE Section #### (Not yet implemented!!!)
				
				/// Derivative of the current
				MatrixVar<Real> mIntfDerCurrent;
				/// 
				Matrix intfDerCurrent() {return mIntfDerCurrent;}
				/// set init value of the current, calculate and set the 
				/// initial value of the derivative of the current
				void setInitialComplexIntfCurrent(Complex initCurrent);
				///
				void daePreStep(Real time);
				///
				void daeInitialize(double time, double state[], double dstate_dt[],
					double absoluteTolerances[], double stateVarTypes[], int& offset) {};
				/// Residual function for DAE Solver
				void daeResidual(double time, const double state[], const double dstate_dt[], double resid[], std::vector<int>& off) {};
				///
				void daePostStep(double Nexttime, const double state[], const double dstate_dt[], int& offset) {};
				///
				int getNumberOfStateVariables() {return 0;}

			};
		}
	}
}
