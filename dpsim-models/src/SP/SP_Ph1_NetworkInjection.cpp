/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/SP/SP_Ph1_NetworkInjection.h>

using namespace CPS;


SP::Ph1::NetworkInjection::NetworkInjection(String uid, String name,
    Logger::Level logLevel) : SimPowerComp<Complex>(uid, name, logLevel),
	mVoltageRef(Attribute<Complex>::createDynamic("V_ref", mAttributes)),
	mSrcFreq(Attribute<Real>::createDynamic("f_src", mAttributes)),
	mVoltageSetPoint(Attribute<Real>::create("V_set", mAttributes)),
	mVoltageSetPointPerUnit(Attribute<Real>::create("V_set_pu", mAttributes, 1.0)),
	mActivePowerInjection(Attribute<Real>::create("p_inj", mAttributes)),
	mReactivePowerInjection(Attribute<Real>::create("q_inj", mAttributes)) {

	mSLog->info("Create {} of type {}", **mName, this->type());
	mSLog->flush();
	**mIntfVoltage = MatrixComp::Zero(1, 1);
	**mIntfCurrent = MatrixComp::Zero(1, 1);
	setVirtualNodeNumber(0);
	setTerminalNumber(1);

	// Create electrical sub components
	mSubVoltageSource = std::make_shared<SP::Ph1::VoltageSource>(**mName + "_vs", mLogLevel);
	mSubComponents.push_back(mSubVoltageSource);
	mSLog->info("Electrical subcomponents: ");
	for (auto subcomp: mSubComponents)
		mSLog->info("- {}", subcomp->name());

	// MNA attributes
	mSubVoltageSource->mVoltageRef->setReference(mVoltageRef);
	mSubVoltageSource->mSrcFreq->setReference(mSrcFreq);
}

// #### Powerflow section ####

void SP::Ph1::NetworkInjection::setParameters(Real voltageSetPoint) {
	**mVoltageSetPoint = voltageSetPoint;

	mSLog->info("Voltage Set-Point ={} [V]", **mVoltageSetPoint);
	mSLog->flush();

	mParametersSet = true;
}

void SP::Ph1::NetworkInjection::setParameters(Complex initialPhasor, Real freqStart, Real rocof, Real timeStart, Real duration, bool smoothRamp) {
	mParametersSet = true;

	mSubVoltageSource->setParameters(initialPhasor, freqStart, rocof, timeStart, duration, smoothRamp);

	mSLog->info("\nVoltage Ref={:s} [V]"
				"\nFrequency={:s} [Hz]",
				Logger::phasorToString(initialPhasor),
				Logger::realToString(freqStart));
}

void SP::Ph1::NetworkInjection::setParameters(Complex initialPhasor, Real modulationFrequency, Real modulationAmplitude, Real baseFrequency /*= 0.0*/, bool zigzag /*= false*/) {
	mParametersSet = true;

	mSubVoltageSource->setParameters(initialPhasor, modulationFrequency, modulationAmplitude, baseFrequency, zigzag);

	mSLog->info("\nVoltage Ref={:s} [V]"
				"\nFrequency={:s} [Hz]",
				Logger::phasorToString(initialPhasor),
				Logger::realToString(baseFrequency));
}

void SP::Ph1::NetworkInjection::setBaseVoltage(Real baseVoltage) {
    mBaseVoltage = baseVoltage;
}

void SP::Ph1::NetworkInjection::calculatePerUnitParameters(Real baseApparentPower, Real baseOmega) {
    mSLog->info("#### Calculate Per Unit Parameters for {}", **mName);
	mSLog->info("Base Voltage={} [V]", mBaseVoltage);

    **mVoltageSetPointPerUnit = **mVoltageSetPoint / mBaseVoltage;

	mSLog->info("Voltage Set-Point ={} [pu]", **mVoltageSetPointPerUnit);
	mSLog->flush();
}

void SP::Ph1::NetworkInjection::modifyPowerFlowBusType(PowerflowBusType powerflowBusType) {
	mPowerflowBusType = powerflowBusType;
}

void SP::Ph1::NetworkInjection::updatePowerInjection(Complex powerInj) {
	**mActivePowerInjection = powerInj.real();
	**mReactivePowerInjection = powerInj.imag();
}

// #### MNA Section ####

void SP::Ph1::NetworkInjection::setParameters(Complex voltageRef, Real srcFreq) {
	mParametersSet = true;

	mSubVoltageSource->setParameters(voltageRef, srcFreq);

	mSLog->info("\nVoltage Ref={:s} [V]"
				"\nFrequency={:s} [Hz]",
				Logger::phasorToString(voltageRef),
				Logger::realToString(srcFreq));
}

SimPowerComp<Complex>::Ptr SP::Ph1::NetworkInjection::clone(String name) {
	auto copy = NetworkInjection::make(name, mLogLevel);
	copy->setParameters(**mVoltageRef);
	return copy;
}

void SP::Ph1::NetworkInjection::initializeFromNodesAndTerminals(Real frequency) {
	// Connect electrical subcomponents
	mSubVoltageSource->connect({ SimNode::GND, node(0) });

	// Initialize electrical subcomponents
	for (auto subcomp: mSubComponents) {
		subcomp->initialize(mFrequencies);
		subcomp->initializeFromNodesAndTerminals(frequency);
	}
}

// #### MNA functions ####

void SP::Ph1::NetworkInjection::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	// initialize electrical subcomponents
	for (auto subcomp: mSubComponents)
		if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
			mnasubcomp->mnaInitialize(omega, timeStep, leftVector);

	// collect right side vectors of subcomponents
	mRightVectorStamps.push_back(&**mSubVoltageSource->mRightVector);

	// collect tasks
	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));

	**mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
}

void SP::Ph1::NetworkInjection::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	for (auto subcomp: mSubComponents)
		if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
			mnasubcomp->mnaApplySystemMatrixStamp(systemMatrix);
}

void SP::Ph1::NetworkInjection::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	rightVector.setZero();
	for (auto stamp : mRightVectorStamps)
		rightVector += *stamp;

	mSLog->debug("Right Side Vector: {:s}",
				Logger::matrixToString(rightVector));
}

void SP::Ph1::NetworkInjection::mnaAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes) {
	// add pre-step dependencies of subcomponents
	for (auto subcomp: mSubComponents)
		if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
			mnasubcomp->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
	// add pre-step dependencies of component itself
	prevStepDependencies.push_back(mIntfCurrent);
	prevStepDependencies.push_back(mIntfVoltage);
	modifiedAttributes.push_back(mRightVector);
}

void SP::Ph1::NetworkInjection::mnaPreStep(Real time, Int timeStepCount) {
	// pre-step of subcomponents
	for (auto subcomp: mSubComponents)
		if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
			mnasubcomp->mnaPreStep(time, timeStepCount);
	// pre-step of component itself
	mnaApplyRightSideVectorStamp(**mRightVector);
}

void SP::Ph1::NetworkInjection::mnaAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	// add post-step dependencies of subcomponents
	for (auto subcomp: mSubComponents)
		if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
			mnasubcomp->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
	// add post-step dependencies of component itself
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(mIntfVoltage);
	modifiedAttributes.push_back(mIntfCurrent);
}

void SP::Ph1::NetworkInjection::mnaPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	// post-step of subcomponents
	for (auto subcomp: mSubComponents)
		if (auto mnasubcomp = std::dynamic_pointer_cast<MNAInterface>(subcomp))
			mnasubcomp->mnaPostStep(time, timeStepCount, leftVector);
	// post-step of component itself
	mnaUpdateCurrent(**leftVector);
	mnaUpdateVoltage(**leftVector);
}

void SP::Ph1::NetworkInjection::mnaUpdateVoltage(const Matrix& leftVector) {
	**mIntfVoltage = **mSubVoltageSource->mIntfVoltage;
}

void SP::Ph1::NetworkInjection::mnaUpdateCurrent(const Matrix& leftVector) {
	**mIntfCurrent = **mSubVoltageSource->mIntfCurrent;
}