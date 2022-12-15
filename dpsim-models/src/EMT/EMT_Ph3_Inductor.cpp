/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_Inductor.h>

using namespace CPS;

EMT::Ph3::Inductor::Inductor(String uid, String name, Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel),
	mIntfDerCurrent(Attribute<Matrix>::create("di_intf", mAttributes))  {

	mPhaseType = PhaseType::ABC;
	setTerminalNumber(2);
	mEquivCurrent = Matrix::Zero(3, 1);
	**mIntfVoltage = Matrix::Zero(3, 1);
	**mIntfCurrent = Matrix::Zero(3, 1);

	// initialize dae specific attributes
	**mIntfDerCurrent = Matrix::Zero(3,1);

	///FIXME: Initialization should happen in the base class declaring the attribute. However, this base class is currently not an AttributeList...
	mInductance = CPS::Attribute<Matrix>::create("L", mAttributes);
}

SimPowerComp<Real>::Ptr EMT::Ph3::Inductor::clone(String name) {
	auto copy = Inductor::make(name, mLogLevel);
	copy->setParameters(**mInductance);
	return copy;
}

void EMT::Ph3::Inductor::initializeFromNodesAndTerminals(Real frequency) {

	Real omega = 2 * PI * frequency;
	MatrixComp impedance = MatrixComp::Zero(3, 3);
	impedance <<
		Complex(0, omega * (**mInductance)(0, 0)), Complex(0, omega * (**mInductance)(0, 1)), Complex(0, omega * (**mInductance)(0, 2)),
		Complex(0, omega * (**mInductance)(1, 0)), Complex(0, omega * (**mInductance)(1, 1)), Complex(0, omega * (**mInductance)(1, 2)),
		Complex(0, omega * (**mInductance)(2, 0)), Complex(0, omega * (**mInductance)(2, 1)), Complex(0, omega * (**mInductance)(2, 2));

	MatrixComp vInitABC = Matrix::Zero(3, 1);
	vInitABC(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(1) - RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
	vInitABC(1, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_B;
	vInitABC(2, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_C;
	**mIntfVoltage = vInitABC.real();
	MatrixComp admittance = impedance.inverse();
	**mIntfCurrent = (admittance * vInitABC).real();

	mSLog->info("\nInductance [H]: {:s}"
				"\nImpedance [Ohm]: {:s}",
				Logger::matrixToString(**mInductance),
				Logger::matrixCompToString(impedance));
	mSLog->info(
		"\n--- Initialization from powerflow ---"
		"\nVoltage across: {:s}"
		"\nCurrent: {:s}"
		"\nTerminal 0 voltage: {:s}"
		"\nTerminal 1 voltage: {:s}"
		"\n--- Initialization from powerflow finished ---",
		Logger::matrixToString(**mIntfVoltage),
		Logger::matrixToString(**mIntfCurrent),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(0)),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(1)));
}

void EMT::Ph3::Inductor::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);

	updateMatrixNodeIndices();
	mEquivCond = timeStep / 2. * (**mInductance).inverse();
	// Update internal state
	mEquivCurrent = mEquivCond * **mIntfVoltage + **mIntfCurrent;

	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
	**mRightVector = Matrix::Zero(leftVector->get().rows(), 1);

	mSLog->info(
		"\n--- MNA initialization ---"
		"\nInitial voltage {:s}"
		"\nInitial current {:s}"
		"\nEquiv. current {:s}"
		"\n--- MNA initialization finished ---",
		Logger::matrixToString(**mIntfVoltage),
		Logger::matrixToString(**mIntfCurrent),
		Logger::matrixToString(mEquivCurrent));
	mSLog->flush();
}

void EMT::Ph3::Inductor::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	if (terminalNotGrounded(0)) {
		// set upper left block, 3x3 entries
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 0), mEquivCond(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 1), mEquivCond(0, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(0, 2), mEquivCond(0, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(0, 0), mEquivCond(1, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(0, 1), mEquivCond(1, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(0, 2), mEquivCond(1, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(0, 0), mEquivCond(2, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(0, 1), mEquivCond(2, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(0, 2), mEquivCond(2, 2));
	}
	if (terminalNotGrounded(1)) {
		// set buttom right block, 3x3 entries
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 0), mEquivCond(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 1), mEquivCond(0, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(1, 2), mEquivCond(0, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(1, 0), mEquivCond(1, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(1, 1), mEquivCond(1, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(1, 2), mEquivCond(1, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(1, 0), mEquivCond(2, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(1, 1), mEquivCond(2, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(1, 2), mEquivCond(2, 2));
	}
	// Set off diagonal blocks, 2x3x3 entries
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 0), -mEquivCond(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 1), -mEquivCond(0, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 0), matrixNodeIndex(1, 2), -mEquivCond(0, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(1, 0), -mEquivCond(1, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(1, 1), -mEquivCond(1, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 1), matrixNodeIndex(1, 2), -mEquivCond(1, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(1, 0), -mEquivCond(2, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(1, 1), -mEquivCond(2, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0, 2), matrixNodeIndex(1, 2), -mEquivCond(2, 2));


		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 0), -mEquivCond(0, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 1), -mEquivCond(0, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 0), matrixNodeIndex(0, 2), -mEquivCond(0, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(0, 0), -mEquivCond(1, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(0, 1), -mEquivCond(1, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 1), matrixNodeIndex(0, 2), -mEquivCond(1, 2));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(0, 0), -mEquivCond(2, 0));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(0, 1), -mEquivCond(2, 1));
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1, 2), matrixNodeIndex(0, 2), -mEquivCond(2, 2));
	}

	mSLog->info(
		"\nEquivalent Conductance: {:s}",
		Logger::matrixToString(mEquivCond));
}

void EMT::Ph3::Inductor::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	// Update internal state
	mEquivCurrent = mEquivCond * **mIntfVoltage + **mIntfCurrent;
	if (terminalNotGrounded(0)) {
		Math::setVectorElement(rightVector, matrixNodeIndex(0, 0), mEquivCurrent(0, 0));
		Math::setVectorElement(rightVector, matrixNodeIndex(0, 1), mEquivCurrent(1, 0));
		Math::setVectorElement(rightVector, matrixNodeIndex(0, 2), mEquivCurrent(2, 0));
	}
	if (terminalNotGrounded(1)) {
		Math::setVectorElement(rightVector, matrixNodeIndex(1, 0), -mEquivCurrent(0, 0));
		Math::setVectorElement(rightVector, matrixNodeIndex(1, 1), -mEquivCurrent(1, 0));
		Math::setVectorElement(rightVector, matrixNodeIndex(1, 2), -mEquivCurrent(2, 0));
	}
	mSLog->debug(
		"\nEquivalent Current (mnaApplyRightSideVectorStamp): {:s}",
		Logger::matrixToString(mEquivCurrent));
	mSLog->flush();
}

void EMT::Ph3::Inductor::mnaAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) {
	// actually depends on L, but then we'd have to modify the system matrix anyway
	prevStepDependencies.push_back(attribute("v_intf"));
	prevStepDependencies.push_back(attribute("i_intf"));
	modifiedAttributes.push_back(attribute("right_vector"));
}

void EMT::Ph3::Inductor::mnaPreStep(Real time, Int timeStepCount) {
	mnaApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::Inductor::mnaAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(attribute("v_intf"));
	modifiedAttributes.push_back(attribute("i_intf"));
}

void EMT::Ph3::Inductor::mnaPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	mnaUpdateVoltage(**leftVector);
	mnaUpdateCurrent(**leftVector);
}

void EMT::Ph3::Inductor::mnaUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	**mIntfVoltage = Matrix::Zero(3, 1);
	if (terminalNotGrounded(1)) {
		(**mIntfVoltage)(0, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 0));
		(**mIntfVoltage)(1, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 1));
		(**mIntfVoltage)(2, 0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1, 2));
	}
	if (terminalNotGrounded(0)) {
		(**mIntfVoltage)(0, 0) = (**mIntfVoltage)(0, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 0));
		(**mIntfVoltage)(1, 0) = (**mIntfVoltage)(1, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 1));
		(**mIntfVoltage)(2, 0) = (**mIntfVoltage)(2, 0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0, 2));
	}
	mSLog->debug(
		"\nUpdate Voltage: {:s}",
		Logger::matrixToString(**mIntfVoltage)
	);
}

void EMT::Ph3::Inductor::mnaUpdateCurrent(const Matrix& leftVector) {
	**mIntfCurrent = mEquivCond * **mIntfVoltage + mEquivCurrent;
	mSLog->debug(
		"\nUpdate Current: {:s}",
		Logger::matrixToString(**mIntfCurrent)
	);
	mSLog->flush();
}


// #### DAE Section ####

void EMT::Ph3::Inductor::daeInitialize(double time, double state[], double dstate_dt[], 
	double absoluteTolerances[], double stateVarTypes[], int &offset) {
	// state[offset] = current through inductor
	// dstate_dt[offset] = inductor current derivative
	updateMatrixNodeIndices();

	// set init values of states
	state[offset] 	= (**mIntfCurrent)(0,0);	
	state[offset+1] = (**mIntfCurrent)(1,0);	
	state[offset+2] = (**mIntfCurrent)(2,0);	
	dstate_dt[offset] 	= (**mIntfVoltage)(0,0) / (**mInductance)(0, 0);
	dstate_dt[offset+1] = (**mIntfVoltage)(1,0) / (**mInductance)(1, 1);	
	dstate_dt[offset+2] = (**mIntfVoltage)(2,0) / (**mInductance)(2, 2);		
	(**mIntfDerCurrent)(0,0) = dstate_dt[offset];
	(**mIntfDerCurrent)(1,0) = dstate_dt[offset+1];
	(**mIntfDerCurrent)(2,0) = dstate_dt[offset+2];
	
	// set state variable as differential variable
	stateVarTypes[offset]  	 = 1.0;
	stateVarTypes[offset+1]  = 1.0;
	stateVarTypes[offset+2]  = 1.0;

	// set absolute tolerance
	absoluteTolerances[offset] 	 = mAbsTolerance;
	absoluteTolerances[offset+1] = mAbsTolerance;
	absoluteTolerances[offset+2] = mAbsTolerance;

	mSLog->info(
		"\n--- daeInitialize ---"
		"\nAdded current (phase a) through inductor '{:s}' to state vector, initial value = {:f}A"
		"\nAdded current (phase b) through inductor '{:s}' to state vector, initial value = {:f}A"
		"\nAdded current (phase c) through inductor '{:s}' to state vector, initial value = {:f}A"
		"\nAdded derivative of current (phase a) inductor '{:s}' to derivative state vector, initial value = {:f}"
		"\nAdded derivative of current (phase b) inductor '{:s}' to derivative state vector, initial value = {:f}"
		"\nAdded derivative of current (phase c) inductor '{:s}' to derivative state vector, initial value = {:f}"
		"\nState variable set as differential"
		"\nAbsolute tolerance={:f}"
		"\n--- daeInitialize finished ---",
		this->name(), state[offset], 
		this->name(), state[offset+1], 
		this->name(), state[offset+2], 
		this->name(), dstate_dt[offset],
		this->name(), dstate_dt[offset+1],
		this->name(), dstate_dt[offset+2],
		absoluteTolerances[offset]
	);
	mSLog->flush();

	offset = offset + 3;
}

void EMT::Ph3::Inductor::daeResidual(double sim_time, 
	const double state[], const double dstate_dt[], 
	double resid[], std::vector<int>& off) {

	// current offset for component
	int c_offset = off[0] + off[1]; 

	resid[c_offset]   = - (**mInductance)(0, 0) * dstate_dt[c_offset];
	resid[c_offset+1] = - (**mInductance)(1, 1) * dstate_dt[c_offset+1];
	resid[c_offset+2] = - (**mInductance)(2, 2) * dstate_dt[c_offset+2];
	if (terminalNotGrounded(0)) {
		resid[c_offset]   -= state[matrixNodeIndex(0, 0)];
		resid[c_offset+1] -= state[matrixNodeIndex(0, 1)];
		resid[c_offset+2] -= state[matrixNodeIndex(0, 2)];
		resid[matrixNodeIndex(0, 0)] -= state[c_offset];
		resid[matrixNodeIndex(0, 1)] -= state[c_offset+1];
		resid[matrixNodeIndex(0, 2)] -= state[c_offset+2];
	}
	if (terminalNotGrounded(1)) {
		resid[c_offset]   += state[matrixNodeIndex(1, 0)];
		resid[c_offset+1] += state[matrixNodeIndex(1, 1)];
		resid[c_offset+2] += state[matrixNodeIndex(1, 2)];
		resid[matrixNodeIndex(1, 0)] += state[c_offset];
		resid[matrixNodeIndex(1, 1)] += state[c_offset+1];
		resid[matrixNodeIndex(1, 2)] += state[c_offset+2];
	}

	off[1] = off[1] + 3;
}

void EMT::Ph3::Inductor::daeJacobian(double current_time, const double state[], 
	const double dstate_dt[], SUNMatrix jacobian, double cj, std::vector<int>& off) {

	// current offset for component
	int c_offset = off[0] + off[1]; 

	SM_ELEMENT_D(jacobian, c_offset, c_offset)	   -= cj * (**mInductance)(0, 0);
	SM_ELEMENT_D(jacobian, c_offset+1, c_offset+1) -= cj * (**mInductance)(1, 1);
	SM_ELEMENT_D(jacobian, c_offset+2, c_offset+2) -= cj * (**mInductance)(2, 2);
	if (terminalNotGrounded(1)) {
		SM_ELEMENT_D(jacobian, c_offset, matrixNodeIndex(1, 0))   += 1.0;
		SM_ELEMENT_D(jacobian, c_offset+1, matrixNodeIndex(1, 1)) += 1.0;
		SM_ELEMENT_D(jacobian, c_offset+2, matrixNodeIndex(1, 2)) += 1.0;

		SM_ELEMENT_D(jacobian, matrixNodeIndex(1, 0), c_offset)   += 1.0;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(1, 1), c_offset+1) += 1.0;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(1, 2), c_offset+2) += 1.0;
	}

	if (terminalNotGrounded(0)) {
		SM_ELEMENT_D(jacobian, c_offset, matrixNodeIndex(0, 0))   -= 1.0;
		SM_ELEMENT_D(jacobian, c_offset+1, matrixNodeIndex(0, 1)) -= 1.0;
		SM_ELEMENT_D(jacobian, c_offset+2, matrixNodeIndex(0, 2)) -= 1.0;

		SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 0), c_offset)   -= 1.0;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 1), c_offset+1) -= 1.0;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0, 2), c_offset+2) -= 1.0;
	}

	off[1] = off[1] + 3;
}

void EMT::Ph3::Inductor::daePostStep(double Nexttime, const double state[], 
	const double dstate_dt[], int& offset) {
	
	**mIntfVoltage = Matrix::Zero(3,1);
	if (terminalNotGrounded(0)) {
		(**mIntfVoltage)(0,0) -= state[matrixNodeIndex(0,0)];
		(**mIntfVoltage)(1,0) -= state[matrixNodeIndex(0,1)];
		(**mIntfVoltage)(2,0) -= state[matrixNodeIndex(0,2)];
	}
	if (terminalNotGrounded(1)) {
		(**mIntfVoltage)(0,0) += state[matrixNodeIndex(1,0)];
		(**mIntfVoltage)(1,0) += state[matrixNodeIndex(1,1)];
		(**mIntfVoltage)(2,0) += state[matrixNodeIndex(1,2)];
	}

	(**mIntfCurrent)(0,0) = state[offset];
	(**mIntfCurrent)(1,0) = state[offset+1];
	(**mIntfCurrent)(2,0) = state[offset+2];

	(**mIntfDerCurrent)(0,0) = dstate_dt[offset];
	(**mIntfDerCurrent)(1,0) = dstate_dt[offset+1];
	(**mIntfDerCurrent)(2,0) = dstate_dt[offset+2];
	offset = offset + 3;
}


