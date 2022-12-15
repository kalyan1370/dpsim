/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_RXLoad.h>

using namespace CPS;

EMT::Ph3::RXLoad::RXLoad(String uid, String name, Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel),
	mActivePower(Attribute<Matrix>::create("P", mAttributes)),
	mReactivePower(Attribute<Matrix>::create("Q", mAttributes)),
	mNomVoltage(Attribute<Real>::create("V_nom", mAttributes)) {
	mPhaseType = PhaseType::ABC;
	setTerminalNumber(1);

	mSLog->info("Create {} {}", this->type(), name);
	**mIntfVoltage = Matrix::Zero(3, 1);
	**mIntfCurrent = Matrix::Zero(3, 1);
	mSLog->flush();

	//DAE Solver
	mIntfDerCurrent = Matrix::Zero(3, 1);
}

EMT::Ph3::RXLoad::RXLoad(String name,
	Logger::Level logLevel)
	: RXLoad(name, name, logLevel) {
}

EMT::Ph3::RXLoad::RXLoad(String name,
	Matrix activePower, Matrix reactivePower, Real volt,
	Logger::Level logLevel)
	: RXLoad(name, logLevel) {
	**mActivePower = activePower;
	**mReactivePower = reactivePower;
	mPower = MatrixComp::Zero(3,3);
	mPower <<
		Complex((**mActivePower)(0, 0), (**mReactivePower)(0, 0)), Complex((**mActivePower)(0, 1), (**mReactivePower)(0, 1)), Complex((**mActivePower)(0, 2), (**mReactivePower)(0, 2)),
		Complex((**mActivePower)(1, 0), (**mReactivePower)(1, 0)), Complex((**mActivePower)(1, 1), (**mReactivePower)(1, 1)), Complex((**mActivePower)(1, 2), (**mReactivePower)(1, 2)),
		Complex((**mActivePower)(2, 0), (**mReactivePower)(2, 0)), Complex((**mActivePower)(2, 1), (**mReactivePower)(2, 1)), Complex((**mActivePower)(2, 2), (**mReactivePower)(2, 2));

	**mNomVoltage = volt;
	initPowerFromTerminal = false;
}

void EMT::Ph3::RXLoad::setParameters(Matrix activePower, Matrix reactivePower, Real volt) {
	**mActivePower = activePower;
	**mReactivePower = reactivePower;

	// complex power
	mPower = MatrixComp::Zero(3, 3);
	mPower(0, 0) = { (**mActivePower)(0, 0), (**mReactivePower)(0, 0) };
	mPower(1, 1) = { (**mActivePower)(1, 1), (**mReactivePower)(1, 1) };
	mPower(2, 2) = { (**mActivePower)(2, 2), (**mReactivePower)(2, 2) };

	**mNomVoltage = volt;

	mSLog->info("\nActive Power [W]: {}"
			"\nReactive Power [VAr]: {}",
			Logger::matrixToString(**mActivePower),
			Logger::matrixToString(**mReactivePower));
	mSLog->info("Nominal Voltage={} [V]", **mNomVoltage);

	initPowerFromTerminal = false;
}

SimPowerComp<Real>::Ptr EMT::Ph3::RXLoad::clone(String name) {
	// everything set by initializeFromNodesAndTerminals
	return RXLoad::make(name, mLogLevel);
}

void EMT::Ph3::RXLoad::initializeFromNodesAndTerminals(Real frequency) {

	if (initPowerFromTerminal) {
		**mActivePower = Matrix::Zero(3, 3);
		(**mActivePower)(0, 0) = mTerminals[0]->singleActivePower() / 3.;
		(**mActivePower)(1, 1) = mTerminals[0]->singleActivePower() / 3.;
		(**mActivePower)(2, 2) = mTerminals[0]->singleActivePower() / 3.;

		**mReactivePower = Matrix::Zero(3, 3);
		(**mReactivePower)(0, 0) = mTerminals[0]->singleReactivePower() / 3.;
		(**mReactivePower)(1, 1) = mTerminals[0]->singleReactivePower() / 3.;
		(**mReactivePower)(2, 2) = mTerminals[0]->singleReactivePower() / 3.;

		// complex power
		mPower = MatrixComp::Zero(3, 3);
		mPower(0, 0) = { (**mActivePower)(0, 0), (**mReactivePower)(0, 0) };
		mPower(1, 1) = { (**mActivePower)(1, 1), (**mReactivePower)(1, 1) };
		mPower(2, 2) = { (**mActivePower)(2, 2), (**mReactivePower)(2, 2) };

		**mNomVoltage = std::abs(mTerminals[0]->initialSingleVoltage());

		mSLog->info("\nActive Power [W]: {}"
					"\nReactive Power [VAr]: {}",
					Logger::matrixToString(**mActivePower),
					Logger::matrixToString(**mReactivePower));
		mSLog->info("Nominal Voltage={} [V]", **mNomVoltage);
	}

	if ((**mActivePower)(0,0) != 0) {
		mResistance = std::pow(**mNomVoltage / sqrt(3), 2) * (**mActivePower).inverse();
		mConductance = mResistance.inverse();
		mSubResistor = std::make_shared<EMT::Ph3::Resistor>(**mName + "_res", mLogLevel);
		mSubResistor->setParameters(mResistance);
		mSubResistor->connect({ SimNode::GND, mTerminals[0]->node() });
		mSubResistor->initialize(mFrequencies);
		mSubResistor->initializeFromNodesAndTerminals(frequency);
	}

	if ((**mReactivePower)(0, 0) != 0)
		mReactance = std::pow(**mNomVoltage / sqrt(3), 2) * (**mReactivePower).inverse();
	else
		mReactance = Matrix::Zero(1, 1);

	if (mReactance(0,0) > 0) {
		mInductance = mReactance / (2 * PI * frequency);

		mSubInductor = std::make_shared<EMT::Ph3::Inductor>(**mName + "_ind", mLogLevel);
		mSubInductor->setParameters(mInductance);
		mSubInductor->connect({ SimNode::GND, mTerminals[0]->node() });
		mSubInductor->initialize(mFrequencies);
		mSubInductor->initializeFromNodesAndTerminals(frequency);
	}
	else if (mReactance(0,0) < 0) {
		mCapacitance = -1 / (2 * PI * frequency) * mReactance.inverse();

		mSubCapacitor = std::make_shared<EMT::Ph3::Capacitor>(**mName + "_cap", mLogLevel);
		mSubCapacitor->setParameters(mCapacitance);
		mSubCapacitor->connect({ SimNode::GND, mTerminals[0]->node() });
		mSubCapacitor->initialize(mFrequencies);
		mSubCapacitor->initializeFromNodesAndTerminals(frequency);
	}

	MatrixComp vInitABC = MatrixComp::Zero(3, 1);
	vInitABC(0, 0) = RMS3PH_TO_PEAK1PH * mTerminals[0]->initialSingleVoltage();
	vInitABC(1, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_B;
	vInitABC(2, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_C;
	**mIntfVoltage = vInitABC.real();

	MatrixComp iInitABC = MatrixComp::Zero(3, 1);
	// v i^T* = S
	// v^T v i^T* = v^T S
	// i^T*= (|v|^2)^(-1) v^T S

	Complex v_ = vInitABC(0, 0)*vInitABC(0, 0) + vInitABC(1, 0)*vInitABC(1, 0) + vInitABC(2, 0)*vInitABC(2, 0);
	MatrixComp rhs_ = Complex(1, 0) / v_ * vInitABC.transpose() * mPower;
	iInitABC = rhs_.conjugate().transpose();
	**mIntfCurrent = iInitABC.real();

	mSLog->info(
		"\n--- Initialization from powerflow ---"
		"\nVoltage across: {:s}"
		"\nCurrent: {:s}"
		"\nTerminal 0 voltage: {:s}"
		"\nActive Power: {:s}"
		"\nReactive Power: {:s}"
		"\nResistance: {:s}"
		"\nReactance: {:s}"
		"\n--- Initialization from powerflow finished ---",
		Logger::matrixToString(**mIntfVoltage),
		Logger::matrixToString(**mIntfCurrent),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(0)),
		Logger::matrixToString(**mActivePower),
		Logger::matrixToString(**mReactivePower),
		Logger::matrixToString(mResistance),
		Logger::matrixToString(mReactance));
	mSLog->flush();

}

void EMT::Ph3::RXLoad::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();
	**mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
	if (mSubResistor) {
		mSubResistor->mnaInitialize(omega, timeStep, leftVector);
		for (auto task : mSubResistor->mnaTasks()) {
			mMnaTasks.push_back(task);
		}
	}
	if (mSubInductor) {
		mSubInductor->mnaInitialize(omega, timeStep, leftVector);
		for (auto task : mSubInductor->mnaTasks()) {
			mMnaTasks.push_back(task);
		}
	}
	if (mSubCapacitor) {
		mSubCapacitor->mnaInitialize(omega, timeStep, leftVector);
		for (auto task : mSubCapacitor->mnaTasks()) {
			mMnaTasks.push_back(task);
		}
	}
	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
}

void EMT::Ph3::RXLoad::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	if (mSubResistor)
		mSubResistor->mnaApplyRightSideVectorStamp(rightVector);
	if (mSubInductor)
		mSubInductor->mnaApplyRightSideVectorStamp(rightVector);
	if (mSubCapacitor)
		mSubCapacitor->mnaApplyRightSideVectorStamp(rightVector);
}

void EMT::Ph3::RXLoad::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	if (mSubResistor)
		mSubResistor->mnaApplySystemMatrixStamp(systemMatrix);
	if (mSubInductor)
		mSubInductor->mnaApplySystemMatrixStamp(systemMatrix);
	if (mSubCapacitor)
		mSubCapacitor->mnaApplySystemMatrixStamp(systemMatrix);
}

void EMT::Ph3::RXLoad::MnaPreStep::execute(Real time, Int timeStepCount) {
	mLoad.mnaApplyRightSideVectorStamp(**mLoad.mRightVector);
}

void EMT::Ph3::RXLoad::MnaPostStep::execute(Real time, Int timeStepCount) {
	mLoad.mnaUpdateVoltage(**mLeftVector);
	mLoad.mnaUpdateCurrent(**mLeftVector);
}

void EMT::Ph3::RXLoad::mnaUpdateVoltage(const Matrix& leftVector) {
	**mIntfVoltage = Matrix::Zero(3, 1);
	(**mIntfVoltage)(0, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
	(**mIntfVoltage)(1, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
	(**mIntfVoltage)(2, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
}

void EMT::Ph3::RXLoad::mnaUpdateCurrent(const Matrix& leftVector) {
	**mIntfCurrent = Matrix::Zero(3, 1);
	if (mSubResistor)
		**mIntfCurrent += mSubResistor->intfCurrent();
	if (mSubInductor)
		**mIntfCurrent += mSubInductor->intfCurrent();
	if (mSubCapacitor)
		**mIntfCurrent += mSubCapacitor->intfCurrent();
}

// #### DAE functions ####

void EMT::Ph3::RXLoad::daeInitialize(double time, double state[], double dstate_dt[], 
	double absoluteTolerances[], double stateVarTypes[], int& offset) {

	updateMatrixNodeIndices();

	// initial current
	**mIntfCurrent = Matrix::Zero(3,1);
	if ((**mActivePower)(0, 0) != 0)
		**mIntfCurrent += mSubResistor->intfCurrent();
	if ((**mReactivePower)(0, 0) == 0){
		// no state variable

		mSLog->info(
			"\n--- daeInitialize ---"
			"\nmReactivePower = 0  --> no state variables are needed"
			"\nInitial current through the resistor-ph1 of RXLoad '{:s}'={:f}A"
			"\nInitial current through the resistor-ph2 of RXLoad '{:s}'={:f}A"
			"\nInitial current through the resistor-ph3 of RXLoad '{:s}'={:f}A"
			"\n--- daeInitialize finished ---",
			this->name(), (**mIntfCurrent)(0,0),
			this->name(), (**mIntfCurrent)(1,0),
			this->name(), (**mIntfCurrent)(2,0)
		);
		mSLog->flush();
	
	} else if (mReactance(0,0) > 0) {
		// state variables are inductor currents
		**mIntfCurrent += mSubInductor->intfCurrent();

		// init current throw inductor: i_L = i-i_r
		state[offset] = (mSubInductor->intfCurrent())(0,0);
		dstate_dt[offset]   = (**mIntfVoltage)(0,0) / mInductance(0,0);
		state[offset+1] = (mSubInductor->intfCurrent())(1,0);
		dstate_dt[offset+1] = (**mIntfVoltage)(1,0) / mInductance(1,1);
		state[offset+2] = (mSubInductor->intfCurrent())(2,0);
		dstate_dt[offset+2] = (**mIntfVoltage)(2,0) / mInductance(2, 2);

		// init derivative current throw inductor:
		mIntfDerCurrent(0, 0) = dstate_dt[offset];
		mIntfDerCurrent(1, 0) = dstate_dt[offset+1];
		mIntfDerCurrent(2, 0) = dstate_dt[offset+2];

		//set state variable as differential variable
		stateVarTypes[offset]   = 0.0;
		stateVarTypes[offset+1] = 0.0;
		stateVarTypes[offset+2] = 0.0;

		//set absolute tolerance
		absoluteTolerances[offset] = mAbsTolerance;
		absoluteTolerances[offset+1] = mAbsTolerance;
		absoluteTolerances[offset+2] = mAbsTolerance;

		mSLog->info(
			"\n--- daeInitialize ---"
			"\nmReactance(0,0) > 0  --> state variable are inductor currents"
			"\nAdded current-phase1 through the inductor of RXLoad '{:s}' to state vector, initial value={:f}A"
			"\nAdded current-phase2 through the inductor of RXLoad '{:s}' to state vector, initial value={:f}A"
			"\nAdded current-phase3 through the inductor of RXLoad '{:s}' to state vector, initial value={:f}A"
			"\nAdded derivative of current-phase1 through the inductor of RXLoad '{:s}' to derivative state vector, initial value={:f}"
			"\nAdded derivative of current-phase2 through the inductor of RXLoad '{:s}' to derivative state vector, initial value={:f}"
			"\nAdded derivative of current-phase3 through the inductor of RXLoad '{:s}' to derivative state vector, initial value={:f}"
			"\nInitial current through the resistor-ph1 of RXLoad '{:s}'={:f}A"
			"\nInitial current through the resistor-ph2 of RXLoad '{:s}'={:f}A"
			"\nInitial current through the resistor-ph3 of RXLoadd '{:s}'={:f}A"
			"\nState variables set as differential"
			"\nAbsolute tolerances={:f}"
			"\n--- daeInitialize finished ---",
			this->name(), state[offset],
			this->name(), state[offset+1],
			this->name(), state[offset+2],
			this->name(), dstate_dt[offset],
			this->name(), dstate_dt[offset+1],
			this->name(), dstate_dt[offset+2],
			this->name(), (**mIntfVoltage)(0,0) * mConductance(0,0),
			this->name(), (**mIntfVoltage)(1,0) * mConductance(1,1),
			this->name(), (**mIntfVoltage)(2,0) * mConductance(2,2),
			absoluteTolerances[offset]
		);
		mSLog->flush();
		offset+=3;

	} else if (mReactance(0,0) < 0) {
		// no state variable are needed
		**mIntfCurrent += mSubCapacitor->intfCurrent();

		mSLog->info(
		"\n--- daeInitialize ---"
		"\nmReactance(0,0) < 0  --> no state variable are needed"
		"\nInitial current through RXLoad '{:s}'-ph1={:f}A"
		"\nInitial current through RXLoad '{:s}'-ph2={:f}A"
		"\nInitial current through RXLoad '{:s}'-ph3={:f}A"
		"\n--- daeInitialize finished ---",
		this->name(), (**mIntfCurrent)(0,0),
		this->name(), (**mIntfCurrent)(1,0),
		this->name(), (**mIntfCurrent)(2,0)
		);
		mSLog->flush();
	}
}

void EMT::Ph3::RXLoad::daeResidual(double sim_time, 
	const double state[], const double dstate_dt[], 
	double resid[], std::vector<int>& off) {

	// current offset for component	
	int c_offset = off[0] + off[1]; 

	if ((**mActivePower)(0, 0) != 0) {
		// add currents through resistor to nodal equations
		resid[matrixNodeIndex(0, 0)] += state[matrixNodeIndex(0, 0)] * mConductance(0,0);
		resid[matrixNodeIndex(0, 1)] += state[matrixNodeIndex(0, 1)] * mConductance(1,1);
		resid[matrixNodeIndex(0, 2)] += state[matrixNodeIndex(0, 2)] * mConductance(2,2);
	}
	if (mReactance(0,0) > 0) {
		// state variables are inductor currents
		resid[c_offset]   = state[matrixNodeIndex(0, 0)] - mInductance(0,0) * dstate_dt[c_offset];
		resid[c_offset+1] = state[matrixNodeIndex(0, 1)] - mInductance(1,1) * dstate_dt[c_offset+1];
		resid[c_offset+2] = state[matrixNodeIndex(0, 2)] - mInductance(2,2) * dstate_dt[c_offset+2];
		resid[matrixNodeIndex(0, 0)] += state[c_offset];
		resid[matrixNodeIndex(0, 1)] += state[c_offset+1];
		resid[matrixNodeIndex(0, 2)] += state[c_offset+2];
		off[1] += 3;
	}
	else if (mReactance(0,0) < 0) {
		resid[matrixNodeIndex(0, 0)] += mCapacitance(0,0)*dstate_dt[c_offset];
		resid[matrixNodeIndex(0, 1)] += mCapacitance(1,1)*dstate_dt[c_offset+1];
		resid[matrixNodeIndex(0, 2)] += mCapacitance(2,2)*dstate_dt[c_offset+2];
	}
}

void EMT::Ph3::RXLoad::daeJacobian(double current_time, const double state[], 
	const double dstate_dt[], SUNMatrix jacobian, double cj, std::vector<int>& off) {

	// current offset for component
	int c_offset = off[0] + off[1]; 

	if ((**mActivePower)(0, 0) != 0) {
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 0), matrixNodeIndex(0, 0)) += mConductance(0,0);
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 1), matrixNodeIndex(0, 1)) += mConductance(1,1);
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 2), matrixNodeIndex(0, 2)) += mConductance(2,2);
	}
	if (mReactance(0,0) > 0) {
		SM_ELEMENT_D(jacobian, c_offset,   c_offset)   += - cj * mInductance(0,0);
		SM_ELEMENT_D(jacobian, c_offset+1, c_offset+1) += - cj * mInductance(1,1);
		SM_ELEMENT_D(jacobian, c_offset+2, c_offset+2) += - cj * mInductance(2,2);

		SM_ELEMENT_D(jacobian, c_offset,   matrixNodeIndex(0, 0)) += 1.0;
		SM_ELEMENT_D(jacobian, c_offset+1, matrixNodeIndex(0, 1)) += 1.0;
		SM_ELEMENT_D(jacobian, c_offset+2, matrixNodeIndex(0, 2)) += 1.0;

		SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 0), c_offset)   += 1.0;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 1), c_offset+1) += 1.0;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 2), c_offset+2) += 1.0;

		off[1] += 3;
	}
	else if (mReactance(0,0) < 0) {
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 0), matrixNodeIndex(0, 0)) += cj * mConductance(0,0);
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 1), matrixNodeIndex(0, 1)) += cj * mConductance(1,1);
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 2), matrixNodeIndex(0, 2)) += cj * mConductance(2,2);
	}
}

void EMT::Ph3::RXLoad::daePostStep(double Nexttime, const double state[], 
	const double dstate_dt[], int& offset) {
	
	(**mIntfVoltage)(0, 0) = state[matrixNodeIndex(0, 0)];
	(**mIntfVoltage)(1, 0) = state[matrixNodeIndex(0, 1)];
	(**mIntfVoltage)(2, 0) = state[matrixNodeIndex(0, 2)];

	**mIntfCurrent = Matrix::Zero(3,1);
	if ((**mActivePower)(0, 0) != 0) 
		**mIntfCurrent += mConductance * **mIntfVoltage;
	if (mReactance(0,0) > 0) {
		(**mIntfCurrent)(0, 0) += state[offset];
		(**mIntfCurrent)(1, 0) += state[offset+1];
		(**mIntfCurrent)(2, 0) += state[offset+2];
		mIntfDerCurrent(0, 0) = dstate_dt[offset];
		mIntfDerCurrent(1, 0) = dstate_dt[offset+1];
		mIntfDerCurrent(2, 0) = dstate_dt[offset+2];

		offset+=3;
	}
	else if (mReactance(0,0) < 0) {
		(**mIntfCurrent)(0, 0) = mCapacitance(0,0) * dstate_dt[matrixNodeIndex(0, 0)];
		(**mIntfCurrent)(1, 0) = mCapacitance(1,1) * dstate_dt[matrixNodeIndex(0, 1)];
		(**mIntfCurrent)(2, 0) = mCapacitance(2,2) * dstate_dt[matrixNodeIndex(0, 2)];
	}
}


