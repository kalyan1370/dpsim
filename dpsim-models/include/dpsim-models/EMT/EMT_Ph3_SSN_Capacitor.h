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
#include <dpsim-models/Base/Base_Ph3_Capacitor.h>

namespace CPS {
	namespace EMT {
		namespace Ph3 {
            namespace SSN{
			    class Capacitor :
				    public Base::Ph3::Capacitor,
				    public MNAInterface,
				    public SimPowerComp<Real>,
				    public SharedFactory<Capacitor> {
			    protected:
					Matrix Dufour_B_k_hat = Matrix::Zero(3, 3);
					Matrix Dufour_W_k_n = Matrix::Zero(3, 3);
                    Matrix mHistoricVoltage = Matrix::Zero(3, 1);
					
		        public:
				    /// Defines UID, name and logging level
				    Capacitor(String uid, String name, Logger::Level logLevel = Logger::Level::off);
				    /// Defines name and logging level
				    Capacitor(String name, Logger::Level logLevel = Logger::Level::off);

				    SimPowerComp<Real>::Ptr clone(String name);

				    // #### General ####
				    /// Initializes component from power flow data
				    void initializeFromNodesAndTerminals(Real frequency);

				    // #### MNA section ####
				    /// Initializes internal variables of the component
				    void mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector);
				    /// Stamps system matrix
				    void mnaApplySystemMatrixStamp(Matrix& systemMatrix);
				    /// Stamps right side (source) vector
				    void mnaApplyRightSideVectorStamp(Matrix& rightVector);
				    /// Update interface voltage from MNA system result
				    void mnaUpdateVoltage(const Matrix& leftVector);
				    /// Update interface current from MNA system result
				    void mnaUpdateCurrent(const Matrix& leftVector);
				    /// MNA pre step operations
				    void mnaPreStep(Real time, Int timeStepCount);
				    /// MNA post step operations
				    void mnaPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector);
				    /// Add MNA pre step dependencies
				    void mnaAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes);
				    /// Add MNA post step dependencies
				    void mnaAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector);

					bool isLinear() const
					{
						return true;
					}

				    class MnaPreStep : public Task {
				    public:
				    	MnaPreStep(Capacitor& capacitor)
						    : Task(**capacitor.mName + ".MnaPreStep"), mCapacitor(capacitor) {
						    	mCapacitor.mnaAddPreStepDependencies(mPrevStepDependencies, mAttributeDependencies, mModifiedAttributes);
					    }
					    void execute(Real time, Int timeStepCount) { mCapacitor.mnaPreStep(time, timeStepCount); };
				    private:
					    Capacitor& mCapacitor;
			    	};

				    class MnaPostStep : public Task {
				    public:
					    MnaPostStep(Capacitor& capacitor, Attribute<Matrix>::Ptr leftVector)
					    	: Task(**capacitor.mName + ".MnaPostStep"), mCapacitor(capacitor), mLeftVector(leftVector) {
						    	mCapacitor.mnaAddPostStepDependencies(mPrevStepDependencies, mAttributeDependencies, mModifiedAttributes, mLeftVector);
					    }
					    void execute(Real time, Int timeStepCount) { mCapacitor.mnaPostStep(time, timeStepCount, mLeftVector); };
				    private:
					    Capacitor& mCapacitor;
					    Attribute<Matrix>::Ptr mLeftVector;
				    };
			    };
            }
		}
	}
}
