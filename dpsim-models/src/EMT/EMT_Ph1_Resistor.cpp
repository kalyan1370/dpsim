/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph1_Resistor.h>

using namespace CPS;

EMT::Ph1::Resistor::Resistor(String uid, String name, Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel) {
	setTerminalNumber(2);
	**mIntfVoltage = Matrix::Zero(1,1);
	**mIntfCurrent = Matrix::Zero(1,1);

	///FIXME: Initialization should happen in the base class declaring the attribute. However, this base class is currently not an AttributeList...
	mResistance = CPS::Attribute<Real>::create("R", mAttributes);
}

SimPowerComp<Real>::Ptr EMT::Ph1::Resistor::clone(String name) {
	auto copy = Resistor::make(name, mLogLevel);
	copy->setParameters(**mResistance);
	return copy;
}

void EMT::Ph1::Resistor::initializeFromNodesAndTerminals(Real frequency) {

	(**mIntfVoltage)(0,0) = (initialSingleVoltage(1) - initialSingleVoltage(0)).real();
	(**mIntfCurrent)(0,0) = (**mIntfVoltage)(0,0) / **mResistance;

	mSLog->info(
		"\n--- Initialization from powerflow ---"
		"\nVoltage across: {:f}"
		"\nCurrent: {:f}"
		"\nTerminal 0 voltage: {:f}"
		"\nTerminal 1 voltage: {:f}"
		"\n--- Initialization from powerflow finished ---",
		(**mIntfVoltage)(0,0),
		(**mIntfCurrent)(0,0),
		initialSingleVoltage(0).real(),
		initialSingleVoltage(1).real());
}

void EMT::Ph1::Resistor::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
}

void EMT::Ph1::Resistor::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	Real conductance = 1. / **mResistance;
	// Set diagonal entries
	if (terminalNotGrounded(0))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(0), conductance);
	if (terminalNotGrounded(1))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(1), conductance);
	// Set off diagonal entries
	if (terminalNotGrounded(0)  &&  terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(1), -conductance);
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(0), -conductance);
	}

	if (terminalNotGrounded(0))
		mSLog->info("Add {:f} to system at ({:d},{:d})", conductance, matrixNodeIndex(0), matrixNodeIndex(0));
	if (terminalNotGrounded(1))
		mSLog->info("Add {:f} to system at ({:d},{:d})", conductance, matrixNodeIndex(1), matrixNodeIndex(1));
	if ( terminalNotGrounded(0)  &&  terminalNotGrounded(1) ) {
		mSLog->info("Add {:f} to system at ({:d},{:d})", -conductance, matrixNodeIndex(0), matrixNodeIndex(1));
		mSLog->info("Add {:f} to system at ({:d},{:d})", -conductance, matrixNodeIndex(1), matrixNodeIndex(0));
	}
}

void EMT::Ph1::Resistor::MnaPostStep::execute(Real time, Int timeStepCount) {
	mResistor.mnaUpdateVoltage(**mLeftVector);
	mResistor.mnaUpdateCurrent(**mLeftVector);
}

void EMT::Ph1::Resistor::mnaUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	(**mIntfVoltage)(0,0) = 0;
	if (terminalNotGrounded(1))
		(**mIntfVoltage)(0,0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1));
	if (terminalNotGrounded(0))
		(**mIntfVoltage)(0,0) = (**mIntfVoltage)(0,0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0));
}

void EMT::Ph1::Resistor::mnaUpdateCurrent(const Matrix& leftVector) {
	(**mIntfCurrent)(0,0) = (**mIntfVoltage)(0,0) / **mResistance;
}

// #### DAE Section ####

void EMT::Ph1::Resistor::daeInitialize(double time, double state[], double dstate_dt[],
	double absoluteTolerances[], double stateVarTypes[], int& offset) {
	updateMatrixNodeIndices();
	mSLog->info(
		"\n--- daeInitialize ---"
		"\nno variable was added by the resistor '{:s}' to the state vector"
		"\n--- daeInitialize finished ---",
		this->name());
	mSLog->flush();
}

void EMT::Ph1::Resistor::daeResidual(double sim_time,
	const double state[], const double dstate_dt[],
	double resid[], std::vector<int>& off) {
	// resid[Pos1] = nodal current equation of node matrixNodeIndex(1) --> add resistor current
	// resid[Pos2] = nodal current equation of node matrixNodeIndex(0) --> substract resistor current

	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		resid[matrixNodeIndex(0)] -= (state[matrixNodeIndex(1)] - state[matrixNodeIndex(0)]) / **mResistance;
		resid[matrixNodeIndex(1)] += (state[matrixNodeIndex(1)] - state[matrixNodeIndex(0)]) / **mResistance;
	} else if (terminalNotGrounded(0)) {
		resid[matrixNodeIndex(0)] += state[matrixNodeIndex(0)] / **mResistance;
	} else if (terminalNotGrounded(1)) {
		resid[matrixNodeIndex(1)] += state[matrixNodeIndex(1)] / **mResistance;
	}
}

void EMT::Ph1::Resistor::daeJacobian(double current_time, const double state[], 
			const double dstate_dt[], SUNMatrix jacobian, double cj, std::vector<int>& off) {

	if (terminalNotGrounded(1))
		SM_ELEMENT_D(jacobian, matrixNodeIndex(1), matrixNodeIndex(1)) += 1.0 / **mResistance;

	if (terminalNotGrounded(0))
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0), matrixNodeIndex(0)) += 1.0 / **mResistance;

	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		SM_ELEMENT_D(jacobian, matrixNodeIndex(1), matrixNodeIndex(0)) += -1.0 / **mResistance;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0), matrixNodeIndex(1)) += -1.0 / **mResistance;
	}
}

void EMT::Ph1::Resistor::daePostStep(double Nexttime, const double state[], const double dstate_dt[], 
	int& offset) {
	
	(**mIntfVoltage)(0,0) = 0.0;
	if (terminalNotGrounded(0))
		(**mIntfVoltage)(0,0) -= state[matrixNodeIndex(0)];
	if (terminalNotGrounded(1))
		(**mIntfVoltage)(0,0) += state[matrixNodeIndex(1)];
		
	(**mIntfCurrent)(0,0) = (**mIntfVoltage)(0,0) / **mResistance;
}