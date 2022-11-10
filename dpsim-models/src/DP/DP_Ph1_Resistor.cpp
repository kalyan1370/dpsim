/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/DP/DP_Ph1_Resistor.h>

using namespace CPS;

DP::Ph1::Resistor::Resistor(String uid, String name, Logger::Level logLevel)
	: SimPowerComp<Complex>(uid, name, logLevel) {
	**mIntfVoltage = MatrixComp::Zero(1,1);
	**mIntfCurrent = MatrixComp::Zero(1,1);
	setTerminalNumber(2);

	///FIXME: Initialization should happen in the base class declaring the attribute. However, this base class is currently not an AttributeList...
	mResistance = CPS::Attribute<Real>::create("R", mAttributes);
}

SimPowerComp<Complex>::Ptr DP::Ph1::Resistor::clone(String name) {
	auto copy = Resistor::make(name, mLogLevel);
	copy->setParameters(**mResistance);
	return copy;
}

void DP::Ph1::Resistor::initializeFromNodesAndTerminals(Real frequency) {

	Complex impedance = { **mResistance, 0 };
	(**mIntfVoltage)(0,0) = initialSingleVoltage(1) - initialSingleVoltage(0);
	(**mIntfCurrent)(0,0) = (**mIntfVoltage)(0,0) / impedance;

	mSLog->info("\nResistance [Ohm]: {:s}"
				"\nImpedance [Ohm]: {:s}",
				Logger::realToString(**mResistance),
				Logger::complexToString(impedance));
	mSLog->info("\n--- Initialization from powerflow ---"
		"\nVoltage across: {:s}"
		"\nCurrent: {:s}"
		"\nTerminal 0 voltage: {:s}"
		"\nTerminal 1 voltage: {:s}"
		"\n--- Initialization from powerflow finished ---",
		Logger::phasorToString((**mIntfVoltage)(0,0)),
		Logger::phasorToString((**mIntfCurrent)(0,0)),
		Logger::phasorToString(initialSingleVoltage(0)),
		Logger::phasorToString(initialSingleVoltage(1)));
}

// #### MNA functions ####
void DP::Ph1::Resistor::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
	mSLog->info(
		"\n--- MNA initialization ---"
		"\nInitial voltage {:s}"
		"\nInitial current {:s}"
		"\n--- MNA initialization finished ---",
		Logger::phasorToString((**mIntfVoltage)(0, 0)),
		Logger::phasorToString((**mIntfCurrent)(0, 0)));
}

void DP::Ph1::Resistor::mnaInitializeHarm(Real omega, Real timeStep, std::vector<Attribute<Matrix>::Ptr> leftVectors) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mMnaTasks.push_back(std::make_shared<MnaPostStepHarm>(*this, leftVectors));
}

void DP::Ph1::Resistor::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	Complex conductance = Complex(1. / **mResistance, 0);

	for (UInt freq = 0; freq < mNumFreqs; freq++) {
		// Set diagonal entries
		if (terminalNotGrounded(0))
			Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(0), conductance, mNumFreqs, freq);
		if (terminalNotGrounded(1))
		// Set off diagonal entries
			Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(1), conductance, mNumFreqs, freq);
		if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
			Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(1), -conductance, mNumFreqs, freq);
			Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(0), -conductance, mNumFreqs, freq);
		}

		mSLog->info("-- Stamp frequency {:d} ---", freq);
		if (terminalNotGrounded(0))
			mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(conductance), matrixNodeIndex(0), matrixNodeIndex(0));
		if (terminalNotGrounded(1))
			mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(conductance), matrixNodeIndex(1), matrixNodeIndex(1));
		if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
			mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(-conductance), matrixNodeIndex(0), matrixNodeIndex(1));
			mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(-conductance), matrixNodeIndex(1), matrixNodeIndex(0));
		}
	}
}

void DP::Ph1::Resistor::mnaApplySystemMatrixStampHarm(Matrix& systemMatrix, Int freqIdx) {
	Complex conductance = Complex(1. / **mResistance, 0);
	// Set diagonal entries
	if (terminalNotGrounded(0))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(0), conductance);
	if (terminalNotGrounded(1))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(1), conductance);
	// Set off diagonal entries
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(1), -conductance);
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(0), -conductance);
	}

	mSLog->info("-- Stamp for frequency {:f} ---", mFrequencies(freqIdx,0));
	if (terminalNotGrounded(0))
		mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(conductance), matrixNodeIndex(0), matrixNodeIndex(0));
	if (terminalNotGrounded(1))
		mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(conductance), matrixNodeIndex(1), matrixNodeIndex(1));
	if (terminalNotGrounded(0)  &&  terminalNotGrounded(1)) {
		mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(-conductance), matrixNodeIndex(0), matrixNodeIndex(1));
		mSLog->info("Add {:s} to system at ({:d},{:d})", Logger::complexToString(-conductance), matrixNodeIndex(1), matrixNodeIndex(0));
	}
}

void DP::Ph1::Resistor::mnaAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(this->mIntfVoltage);
	modifiedAttributes.push_back(this->mIntfCurrent);
}

void DP::Ph1::Resistor::mnaPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	this->mnaUpdateVoltage(**leftVector);
	this->mnaUpdateCurrent(**leftVector);
}

void DP::Ph1::Resistor::MnaPostStepHarm::execute(Real time, Int timeStepCount) {
	for (UInt freq = 0; freq < mResistor.mNumFreqs; freq++)
		mResistor.mnaUpdateVoltageHarm(**mLeftVectors[freq], freq);
	mResistor.mnaUpdateCurrentHarm();
}

void DP::Ph1::Resistor::mnaUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	for (UInt freq = 0; freq < mNumFreqs; freq++) {
		(**mIntfVoltage)(0,freq) = 0;
		if (terminalNotGrounded(1))
			(**mIntfVoltage)(0,freq) = Math::complexFromVectorElement(leftVector, matrixNodeIndex(1), mNumFreqs, freq);
		if (terminalNotGrounded(0))
			(**mIntfVoltage)(0,freq) = (**mIntfVoltage)(0,freq) - Math::complexFromVectorElement(leftVector, matrixNodeIndex(0), mNumFreqs, freq);

		SPDLOG_LOGGER_DEBUG(mSLog, "Voltage {:s}", Logger::phasorToString((**mIntfVoltage)(0,freq)));
	}
}

void DP::Ph1::Resistor::mnaUpdateCurrent(const Matrix& leftVector) {
	for (UInt freq = 0; freq < mNumFreqs; freq++) {
		(**mIntfCurrent)(0,freq) = (**mIntfVoltage)(0,freq) / **mResistance;
		SPDLOG_LOGGER_DEBUG(mSLog, "Current {:s}", Logger::phasorToString((**mIntfCurrent)(0,freq)));
	}
}

void DP::Ph1::Resistor::mnaUpdateVoltageHarm(const Matrix& leftVector, Int freqIdx) {
	// v1 - v0
	(**mIntfVoltage)(0,freqIdx) = 0;
	if (terminalNotGrounded(1))
		(**mIntfVoltage)(0,freqIdx) = Math::complexFromVectorElement(leftVector, matrixNodeIndex(1));
	if (terminalNotGrounded(0))
		(**mIntfVoltage)(0,freqIdx) = (**mIntfVoltage)(0,freqIdx) - Math::complexFromVectorElement(leftVector, matrixNodeIndex(0));

	SPDLOG_LOGGER_DEBUG(mSLog, "Voltage {:s}", Logger::phasorToString((**mIntfVoltage)(0,freqIdx)));
}

void DP::Ph1::Resistor::mnaUpdateCurrentHarm() {
	for (UInt freq = 0; freq < mNumFreqs; freq++) {
		(**mIntfCurrent)(0,freq) = (**mIntfVoltage)(0,freq) / **mResistance;
		SPDLOG_LOGGER_DEBUG(mSLog, "Current {:s}", Logger::phasorToString((**mIntfCurrent)(0,freq)));
	}
}

void DP::Ph1::Resistor::mnaTearApplyMatrixStamp(Matrix& tearMatrix) {
	Math::addToMatrixElement(tearMatrix, mTearIdx, mTearIdx, Complex(**mResistance, 0));
}