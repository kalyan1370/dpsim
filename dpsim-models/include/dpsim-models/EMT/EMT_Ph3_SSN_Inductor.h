#pragma once

#include <dpsim-models/SimPowerComp.h>
#include <dpsim-models/Solver/MNAInterface.h>
#include <dpsim-models/Base/Base_Ph3_Inductor.h>

namespace CPS{
    namespace EMT{
        namespace Ph3{
            namespace SSN{
                class Inductor:
                    public Base::Ph3::Inductor,
				    public MNAInterface,
				    public SimPowerComp<Real>,
				    public SharedFactory<Inductor>
                    {
                protected:
                    //rightsideVector history term
                    Matrix historicCurrent =  Matrix::Zero(3, 1);
                    //dependency on latest Voltage, represented by Conductance in system matrix
					Matrix Dufour_B_k_hat = Matrix::Zero(3, 3);
					Matrix Dufour_W_k_n = Matrix::Zero(3, 3);
                public:
                    /// Defines UID, name, component parameters and logging level
				    Inductor(String uid, String name, Logger::Level logLevel = Logger::Level::off);
				    /// Defines name and logging level
				    Inductor(String name, Logger::Level logLevel = Logger::Level::off);

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
					    MnaPreStep(Inductor& inductor) :
						    Task(**inductor.mName + ".MnaPreStep"), mInductor(inductor) {
							    mInductor.mnaAddPreStepDependencies(mPrevStepDependencies, mAttributeDependencies, mModifiedAttributes);
					    }
					    void execute(Real time, Int timeStepCount) { mInductor.mnaPreStep(time, timeStepCount); };
				    private:
					    Inductor& mInductor;
				    };

				    class MnaPostStep : public Task {
				    public:
					    MnaPostStep(Inductor& inductor, Attribute<Matrix>::Ptr leftVector) :
					    	Task(**inductor.mName + ".MnaPostStep"),
						    mInductor(inductor), mLeftVector(leftVector) {
						    	mInductor.mnaAddPostStepDependencies(mPrevStepDependencies, mAttributeDependencies, mModifiedAttributes, mLeftVector);
					    }
					    void execute(Real time, Int timeStepCount) { mInductor.mnaPostStep(time, timeStepCount, mLeftVector); };
				    private:
				    	Inductor& mInductor;
				    	Attribute<Matrix>::Ptr mLeftVector;
                    };
                private:
                };    
            }
        }
    }
}