/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph1_Capacitor.h>

using namespace CPS;

EMT::Ph1::Capacitor::Capacitor(String uid, String name,	Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel),
	mIntfDerVoltage(Attribute<Matrix>::create("dv_intf", mAttributes))  {
	
	mEquivCurrent = 0;
	**mIntfVoltage = Matrix::Zero(1,1);
	**mIntfCurrent = Matrix::Zero(1,1);
	setTerminalNumber(2);

	// initialize dae specific attributes
	**mIntfDerVoltage = Matrix::Zero(1,1);

	///FIXME: Initialization should happen in the base class declaring the attribute. However, this base class is currently not an AttributeList...
	mCapacitance = CPS::Attribute<Real>::create("C", mAttributes);
}

SimPowerComp<Real>::Ptr EMT::Ph1::Capacitor::clone(String name) {
	auto copy = Capacitor::make(name, mLogLevel);
	copy->setParameters(**mCapacitance);
	return copy;
}

void EMT::Ph1::Capacitor::initializeFromNodesAndTerminals(Real frequency) {

	Real omega = 2 * PI * frequency;
	Complex impedance = { 0, - 1. / (omega * **mCapacitance) };
	(**mIntfVoltage)(0,0) = (initialSingleVoltage(1) - initialSingleVoltage(0)).real();
	(**mIntfCurrent)(0,0) = ((initialSingleVoltage(1) - initialSingleVoltage(0)) / impedance).real();

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

void EMT::Ph1::Capacitor::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();

	mEquivCond = (2.0 * **mCapacitance) / timeStep;
	// Update internal state
	mEquivCurrent = -(**mIntfCurrent)(0,0) + -mEquivCond * (**mIntfVoltage)(0,0);

	**mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
}

void EMT::Ph1::Capacitor::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	if (terminalNotGrounded(0))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(0), mEquivCond);
	if (terminalNotGrounded(1))
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(1), mEquivCond);
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(0), matrixNodeIndex(1), -mEquivCond);
		Math::addToMatrixElement(systemMatrix, matrixNodeIndex(1), matrixNodeIndex(0), -mEquivCond);
	}
}

void EMT::Ph1::Capacitor::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	mEquivCurrent = -(**mIntfCurrent)(0,0) + -mEquivCond * (**mIntfVoltage)(0,0);
	if (terminalNotGrounded(0))
		Math::setVectorElement(rightVector, matrixNodeIndex(0), mEquivCurrent);
	if (terminalNotGrounded(1))
		Math::setVectorElement(rightVector, matrixNodeIndex(1), -mEquivCurrent);
}

void EMT::Ph1::Capacitor::MnaPreStep::execute(Real time, Int timeStepCount) {
	mCapacitor.mnaApplyRightSideVectorStamp(**mCapacitor.mRightVector);
}

void EMT::Ph1::Capacitor::MnaPostStep::execute(Real time, Int timeStepCount) {
	mCapacitor.mnaUpdateVoltage(**mLeftVector);
	mCapacitor.mnaUpdateCurrent(**mLeftVector);
}

void EMT::Ph1::Capacitor::mnaUpdateVoltage(const Matrix& leftVector) {
	// v1 - v0
	(**mIntfVoltage)(0,0) = 0;
	if (terminalNotGrounded(1))
		(**mIntfVoltage)(0,0) = Math::realFromVectorElement(leftVector, matrixNodeIndex(1));
	if (terminalNotGrounded(0))
		(**mIntfVoltage)(0,0) = (**mIntfVoltage)(0,0) - Math::realFromVectorElement(leftVector, matrixNodeIndex(0));
}

void EMT::Ph1::Capacitor::mnaUpdateCurrent(const Matrix& leftVector) {
	(**mIntfCurrent)(0,0) = mEquivCond * (**mIntfVoltage)(0,0) + mEquivCurrent;
}

// #### DAE section ####

void EMT::Ph1::Capacitor::daeInitialize(double time, double state[], double dstate_dt[],
	double absoluteTolerances[], double stateVarTypes[], int &offset) {
	
	updateMatrixNodeIndices();

	mSLog->info(
		"\n--- daeInitialize ---"
		"\nNo state variables are needed"
		"\n--- daeInitialize finished ---"
	);
	mSLog->flush();
}

void EMT::Ph1::Capacitor::daeResidual(double sim_time, 
	const double state[], const double dstate_dt[], 
	double resid[], std::vector<int>& off) {

	int node_pos0 = matrixNodeIndex(0);
    int node_pos1 = matrixNodeIndex(1);
	int c_offset = off[0] + off[1]; //current offset for component

	// add currents to node equations
	if (terminalNotGrounded(0))
		resid[node_pos0] += **mCapacitance * dstate_dt[node_pos0];
	if (terminalNotGrounded(1))
		resid[node_pos1] += **mCapacitance * dstate_dt[node_pos1];
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		resid[node_pos1] -= **mCapacitance * dstate_dt[node_pos0];
		resid[node_pos0] -= **mCapacitance * dstate_dt[node_pos1];
	}
}

void EMT::Ph1::Capacitor::daeJacobian(double current_time, const double state[], const double dstate_dt[], 
			SUNMatrix jacobian, double cj, std::vector<int>& off) {

	int node_pos0 = matrixNodeIndex(0);
    int node_pos1 = matrixNodeIndex(1);

	if (terminalNotGrounded(1))
		SM_ELEMENT_D(jacobian, node_pos1, node_pos1) += cj * **mCapacitance;

	if (terminalNotGrounded(0))
		SM_ELEMENT_D(jacobian, node_pos0, node_pos0) += cj * **mCapacitance;
	
	if (terminalNotGrounded(0) && terminalNotGrounded(1)) {
		SM_ELEMENT_D(jacobian, node_pos1, node_pos0) += - cj * **mCapacitance;
		SM_ELEMENT_D(jacobian, node_pos0, node_pos1) += - cj * **mCapacitance;
	}
}

void EMT::Ph1::Capacitor::daePostStep(double Nexttime, const double state[], 
	const double dstate_dt[], int& offset) {

	(**mIntfCurrent)(0,0) = 0.0;
	(**mIntfVoltage)(0,0) = 0.0;
	(**mIntfDerVoltage)(0,0) = 0.0;

	if (terminalNotGrounded(1)) {
		(**mIntfVoltage)(0,0) += state[matrixNodeIndex(1)];
		(**mIntfDerVoltage)(0,0) += dstate_dt[matrixNodeIndex(1)]; 
	}
	if (terminalNotGrounded(0)) {
		(**mIntfVoltage)(0,0) -= state[matrixNodeIndex(0)];
		(**mIntfDerVoltage)(0,0) -= dstate_dt[matrixNodeIndex(0)]; 
	}

	**mIntfCurrent = **mCapacitance * **mIntfDerVoltage;
}


