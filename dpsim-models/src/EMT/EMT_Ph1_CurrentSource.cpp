/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph1_CurrentSource.h>

using namespace CPS;

EMT::Ph1::CurrentSource::CurrentSource(String uid, String name,	Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel),
	mCurrentRef(Attribute<Complex>::create("I_ref", mAttributes)),
	mSrcFreq(Attribute<Real>::create("f_src", mAttributes)) {
	setTerminalNumber(2);
	**mIntfVoltage = Matrix::Zero(1,1);
	**mIntfCurrent = Matrix::Zero(1,1);
}

SimPowerComp<Real>::Ptr EMT::Ph1::CurrentSource::clone(String name) {
	auto copy = CurrentSource::make(name, mLogLevel);
	copy->setParameters(**mCurrentRef, **mSrcFreq);
	return copy;
}

void EMT::Ph1::CurrentSource::setParameters(Complex currentRef, Real srcFreq) {
	**mCurrentRef = currentRef;
	**mSrcFreq = srcFreq;

	mParametersSet = true;
}

void EMT::Ph1::CurrentSource::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();
	(**mIntfCurrent)(0,0) = Math::abs(**mCurrentRef) * cos(Math::phase(**mCurrentRef));
	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
	**mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
}

void EMT::Ph1::CurrentSource::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	if (terminalNotGrounded(0))
		Math::setVectorElement(rightVector, matrixNodeIndex(0), -(**mIntfCurrent)(0,0));

	if (terminalNotGrounded(1))
		Math::setVectorElement(rightVector, matrixNodeIndex(1), (**mIntfCurrent)(0,0));
}

void EMT::Ph1::CurrentSource::updateState(Real time) {
	Complex currentRef = mCurrentRef->get();
	Real srcFreq = mSrcFreq->get();
	if (srcFreq > 0)
		(**mIntfCurrent)(0,0) = Math::abs(currentRef) * cos(time * 2.*PI*srcFreq + Math::phase(currentRef));
	else
		(**mIntfCurrent)(0,0) = currentRef.real();
}

void EMT::Ph1::CurrentSource::MnaPreStep::execute(Real time, Int timeStepCount) {
	mCurrentSource.updateState(time);
	mCurrentSource.mnaApplyRightSideVectorStamp(**mCurrentSource.mRightVector);
}

void EMT::Ph1::CurrentSource::MnaPostStep::execute(Real time, Int timeStepCount) {
	mCurrentSource.mnaUpdateVoltage(**mLeftVector);
}

void EMT::Ph1::CurrentSource::mnaUpdateVoltage(const Matrix& leftVector) {
	(**mIntfVoltage)(0,0) = 0;
	if (terminalNotGrounded(0))
		(**mIntfVoltage)(0,0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(0));
	if (terminalNotGrounded(1))
		(**mIntfVoltage)(0,0) = (**mIntfVoltage)(0,0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(1));
}

// #### DAE functions ####
void EMT::Ph1::CurrentSource::daeInitialize(double time, double state[], double dstate_dt[], 
	double absoluteTolerances[], double stateVarTypes[], int& offset) {

	updateMatrixNodeIndices();
	mSLog->info(
		"\n--- daeInitialize ---"
		"\nno variable was added by the ideal current source '{:s}' to the state vector"
		"\n--- daeInitialize finished ---",
		this->name());
	mSLog->flush();
}

void EMT::Ph1::CurrentSource::daeResidual(double sim_time, 
	const double state[], const double dstate_dt[], 
	double resid[], std::vector<int>& off) {
	// current source does not has states variables
	// Only the current flowing out of the component must be added to the node current equations

	int Pos1 = matrixNodeIndex(0);
    int Pos2 = matrixNodeIndex(1);

	if (terminalNotGrounded(0))
		resid[Pos1] += (**mIntfCurrent)(0,0);
	if (terminalNotGrounded(1))
		resid[Pos2] -= (**mIntfCurrent)(0,0);
}

void EMT::Ph1::CurrentSource::daePostStep(double Nexttime, const double state[], 
	const double dstate_dt[], int& offset) {

	(**mIntfVoltage)(0,0) = 0.0;
	if (terminalNotGrounded(1))
		(**mIntfVoltage)(0,0) += state[matrixNodeIndex(1)];
	if (terminalNotGrounded(0))
		(**mIntfVoltage)(0,0) -= state[matrixNodeIndex(0)];
}

