/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph1_VoltageSource.h>

using namespace CPS;

EMT::Ph1::VoltageSource::VoltageSource(String uid, String name,	Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel),
	mVoltageRef(Attribute<Complex>::create("V_ref", mAttributes)),
	mSrcFreq(Attribute<Real>::create("f_src", mAttributes)) {
	setVirtualNodeNumber(1);
	setTerminalNumber(2);
	**mIntfVoltage = Matrix::Zero(1,1);
	**mIntfCurrent = Matrix::Zero(1,1);
}

void EMT::Ph1::VoltageSource::setParameters(Complex voltageRef, Real srcFreq) {
	**mVoltageRef = voltageRef;
	**mSrcFreq = srcFreq;

	mParametersSet = true;
}

SimPowerComp<Real>::Ptr EMT::Ph1::VoltageSource::clone(String name) {
	auto copy = VoltageSource::make(name, mLogLevel);
	copy->setParameters(**mVoltageRef, **mSrcFreq);
	return copy;
}

void EMT::Ph1::VoltageSource::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();
	(**mIntfVoltage)(0,0) = Math::abs(**mVoltageRef) * cos(Math::phase(**mVoltageRef));
	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
	**mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
}

void EMT::Ph1::VoltageSource::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	if (terminalNotGrounded(0)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), mVirtualNodes[0]->matrixNodeIndex(), -1);
		Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(0), -1);
	}
	if (terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), mVirtualNodes[0]->matrixNodeIndex(), 1);
		Math::addToMatrixElement(systemMatrix, mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(1), 1);
	}

	if (terminalNotGrounded(0)) {
		mSLog->info("Add {:f} to system at ({:d},{:d})", -1., matrixNodeIndex(0), mVirtualNodes[0]->matrixNodeIndex());
		mSLog->info("Add {:f} to system at ({:d},{:d})", -1., mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(0));
	}
	if (terminalNotGrounded(1)) {
		mSLog->info("Add {:f} to system at ({:d},{:d})", 1., matrixNodeIndex(1), mVirtualNodes[0]->matrixNodeIndex());
		mSLog->info("Add {:f} to system at ({:d},{:d})", 1., mVirtualNodes[0]->matrixNodeIndex(), matrixNodeIndex(1));
	}
}

void EMT::Ph1::VoltageSource::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	Math::setVectorElement(rightVector, mVirtualNodes[0]->matrixNodeIndex(), (**mIntfVoltage)(0,0));
}

void EMT::Ph1::VoltageSource::updateVoltage(Real time) {
	Complex voltageRef = mVoltageRef->get();
	Real srcFreq = mSrcFreq->get();
	if (srcFreq > 0)
		(**mIntfVoltage)(0,0) = Math::abs(voltageRef) * cos(time * 2.*PI*srcFreq + Math::phase(voltageRef));
	else
		(**mIntfVoltage)(0,0) = voltageRef.real();
}

void EMT::Ph1::VoltageSource::MnaPreStep::execute(Real time, Int timeStepCount) {
	mVoltageSource.updateVoltage(time);
	mVoltageSource.mnaApplyRightSideVectorStamp(**mVoltageSource.mRightVector);
}

void EMT::Ph1::VoltageSource::MnaPostStep::execute(Real time, Int timeStepCount) {
	mVoltageSource.mnaUpdateCurrent(**mLeftVector);
}

void EMT::Ph1::VoltageSource::mnaUpdateCurrent(const Matrix& leftVector) {
	(**mIntfCurrent)(0,0) = Math::realFromVectorElement(leftVector, mVirtualNodes[0]->matrixNodeIndex());
}


// #### DAE section ####

void EMT::Ph1::VoltageSource::setInitialComplexIntfCurrent(Complex initCurrent) {
	//set initial current 
    (**mIntfCurrent)(0, 0) = initCurrent.real();
	
	// Calculate initial derivative of current = -omega*Imag(Complex_voltage)
	Real omega = 2 * PI * attribute<Real>("f_src")->get();

	mIntfDerCurrent = Matrix::Zero(1, 1);
	mIntfDerCurrent(0,0) = -omega * initCurrent.imag();
}

void EMT::Ph1::VoltageSource::daeInitialize(double time, double state[], double dstate_dt[], 
	double absoluteTolerances[], double stateVarTypes[], int& offset) {
	// state[c_offset] = current through VoltageSource flowing into the voltage source. Current is positive when it flows into the positive voltage terminal of a voltage source
	// dstate_dt[offset] = derivative of current through voltage source

	updateMatrixNodeIndices();

	this->updateVoltage(time);

	state[offset] = (**mIntfCurrent)(0,0);
	dstate_dt[offset] = mIntfDerCurrent(0,0);

	//set state variable as algebraic variable
	stateVarTypes[offset]  = 0.0;

	//set absolute tolerance 
	absoluteTolerances[offset] = mAbsTolerance;

	mSLog->info(
		"\n--- daeInitialize ---"
		"\nAdded current through VoltageSource '{:s}' to state vector, initial value  = {:f}A"
		"\nAdded derivative of current through VoltageSource '{:s}' to derivative state vector, initial value= {:f}"
		"\nState variable set as algebraic"
		"\nAbsolute tolerance={:f}"
		"\n--- daeInitialize finished ---",
		this->name(), state[offset],
		this->name(), dstate_dt[offset],
		absoluteTolerances[offset]
	);
	mSLog->flush();
	offset++;
}

void EMT::Ph1::VoltageSource::daeResidual(double sim_time, 
	const double state[], const double dstate_dt[], 
	double resid[], std::vector<int>& off) {
	
	updateVoltage(sim_time);

	// current offset for component
	int c_offset = off[0] + off[1]; 

	resid[c_offset] = -(**mIntfVoltage)(0,0);
	if (terminalNotGrounded(0)) {
		resid[c_offset] -= state[matrixNodeIndex(0)];	
		resid[matrixNodeIndex(0)] -= state[c_offset];
	}
	if (terminalNotGrounded(1)) {
		resid[c_offset] += state[matrixNodeIndex(1)];
		resid[matrixNodeIndex(1)] += state[c_offset];
	}
	off[1] += 1;
}

void EMT::Ph1::VoltageSource::daeJacobian(double current_time, const double state[], 
			const double dstate_dt[], SUNMatrix jacobian, double cj, std::vector<int>& off) {

	// current offset for component
	int c_offset = off[0] + off[1]; 

	if (terminalNotGrounded(1)) {
		SM_ELEMENT_D(jacobian, c_offset, matrixNodeIndex(1)) += 1.0;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(1), c_offset) += 1.0;
	}

	if (terminalNotGrounded(0)) {
		SM_ELEMENT_D(jacobian, c_offset, matrixNodeIndex(0)) += -1.0;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0), c_offset) += -1.0;
	}

	off[1] += 1;
}

void EMT::Ph1::VoltageSource::daePostStep(double Nexttime, const double state[], 
	const double dstate_dt[], int& offset) {

	(**mIntfCurrent)(0,0) = state[offset];
	mIntfDerCurrent(0,0) = dstate_dt[offset];
	(**mIntfVoltage)(0,0) = 0.0;
	if (terminalNotGrounded(1)) {
		(**mIntfVoltage)(0,0) += state[matrixNodeIndex(1)];
	}
	if (terminalNotGrounded(0)) {
		(**mIntfVoltage)(0,0) -= state[matrixNodeIndex(0)];
	}
	offset++;
}

