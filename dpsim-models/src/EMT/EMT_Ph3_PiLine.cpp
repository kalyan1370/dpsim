/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_PiLine.h>

using namespace CPS;

EMT::Ph3::PiLine::PiLine(String uid, String name, Logger::Level logLevel)
	: SimPowerComp<Real>(uid, name, logLevel) {
	mPhaseType = PhaseType::ABC;
	setVirtualNodeNumber(1);
	setTerminalNumber(2);

	mSLog->info("Create {} {}", this->type(), name);
	**mIntfVoltage = Matrix::Zero(3, 1);
	**mIntfCurrent = Matrix::Zero(3, 1);

	///FIXME: Move initialization into base class
	mSeriesRes = Attribute<Matrix>::create("R_series", mAttributes);
	mSeriesInd = Attribute<Matrix>::create("L_series", mAttributes);
	mParallelCap = Attribute<Matrix>::create("C_parallel", mAttributes);
	mParallelCond = Attribute<Matrix>::create("G_parallel", mAttributes); 
	mSLog->flush();
}

/// DEPRECATED: Delete method
SimPowerComp<Real>::Ptr EMT::Ph3::PiLine::clone(String name) {
	auto copy = PiLine::make(name, mLogLevel);
	copy->setParameters(**mSeriesRes, **mSeriesInd, **mParallelCap, **mParallelCond);
	return copy;
}

void EMT::Ph3::PiLine::initializeFromNodesAndTerminals(Real frequency) {

	// By default there is always a small conductance to ground to
	// avoid problems with floating nodes.
	Matrix defaultParallelCond = Matrix::Zero(3, 3);
	defaultParallelCond <<
		1e-6, 0, 0,
		0, 1e-6, 0,
		0, 0, 1e-6;
	**mParallelCond = ((**mParallelCond)(0, 0) > 0) ? **mParallelCond : defaultParallelCond;

	// Static calculation
	Real omega = 2. * PI * frequency;
	MatrixComp impedance = MatrixComp::Zero(3, 3);
	impedance <<
		Complex((**mSeriesRes)(0, 0), omega * (**mSeriesInd)(0, 0)), Complex((**mSeriesRes)(0, 1), omega * (**mSeriesInd)(0, 1)), Complex((**mSeriesRes)(0, 2), omega * (**mSeriesInd)(0, 2)),
		Complex((**mSeriesRes)(1, 0), omega * (**mSeriesInd)(1, 0)), Complex((**mSeriesRes)(1, 1), omega * (**mSeriesInd)(1, 1)), Complex((**mSeriesRes)(1, 2), omega * (**mSeriesInd)(1, 2)),
		Complex((**mSeriesRes)(2, 0), omega * (**mSeriesInd)(2, 0)), Complex((**mSeriesRes)(2, 1), omega * (**mSeriesInd)(2, 1)), Complex((**mSeriesRes)(2, 2), omega * (**mSeriesInd)(2, 2));

	MatrixComp vInitABC = MatrixComp::Zero(3, 1);
	vInitABC(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(1) - RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
	vInitABC(1, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_B;
	vInitABC(2, 0) = vInitABC(0, 0) * SHIFT_TO_PHASE_C;
	MatrixComp iInit = impedance.inverse() * vInitABC;
	**mIntfCurrent = iInit.real();
	**mIntfVoltage = vInitABC.real();

	// Initialization of virtual node
	// Initial voltage of phase B,C is set after A
	MatrixComp vInitTerm0 = MatrixComp::Zero(3, 1);
	vInitTerm0(0, 0) = RMS3PH_TO_PEAK1PH * initialSingleVoltage(0);
	vInitTerm0(1, 0) = vInitTerm0(0, 0) * SHIFT_TO_PHASE_B;
	vInitTerm0(2, 0) = vInitTerm0(0, 0) * SHIFT_TO_PHASE_C;

	mVirtualNodes[0]->setInitialVoltage(PEAK1PH_TO_RMS3PH*(vInitTerm0 + **mSeriesRes * iInit));

	// Create series sub components
	mSubSeriesResistor = std::make_shared<EMT::Ph3::Resistor>(**mName + "_res", mLogLevel);
	mSubSeriesResistor->setParameters(**mSeriesRes);
	mSubSeriesResistor->connect({ mTerminals[0]->node(), mVirtualNodes[0] });
	mSubSeriesResistor->initialize(mFrequencies);
	mSubSeriesResistor->initializeFromNodesAndTerminals(frequency);

	mSubSeriesInductor = std::make_shared<EMT::Ph3::Inductor>(**mName + "_ind", mLogLevel);
	mSubSeriesInductor->setParameters(**mSeriesInd);
	mSubSeriesInductor->connect({ mVirtualNodes[0], mTerminals[1]->node() });
	mSubSeriesInductor->initialize(mFrequencies);
	mSubSeriesInductor->initializeFromNodesAndTerminals(frequency);

	// Create parallel sub components
	mSubParallelResistor0 = std::make_shared<EMT::Ph3::Resistor>(**mName + "_con0", mLogLevel);
	mSubParallelResistor0->setParameters(2. * (**mParallelCond).inverse());
	mSubParallelResistor0->connect(SimNode::List{ SimNode::GND, mTerminals[0]->node() });
	mSubParallelResistor0->initialize(mFrequencies);
	mSubParallelResistor0->initializeFromNodesAndTerminals(frequency);

	mSubParallelResistor1 = std::make_shared<EMT::Ph3::Resistor>(**mName + "_con1", mLogLevel);
	mSubParallelResistor1->setParameters(2. * (**mParallelCond).inverse());
	mSubParallelResistor1->connect(SimNode::List{ SimNode::GND, mTerminals[1]->node() });
	mSubParallelResistor1->initialize(mFrequencies);
	mSubParallelResistor1->initializeFromNodesAndTerminals(frequency);

	if ((**mParallelCap)(0,0) > 0) {
		mSubParallelCapacitor0 = std::make_shared<EMT::Ph3::Capacitor>(**mName + "_cap0", mLogLevel);
		mSubParallelCapacitor0->setParameters(**mParallelCap / 2.);
		mSubParallelCapacitor0->connect(SimNode::List{ SimNode::GND, mTerminals[0]->node() });
		mSubParallelCapacitor0->initialize(mFrequencies);
		mSubParallelCapacitor0->initializeFromNodesAndTerminals(frequency);

		mSubParallelCapacitor1 = std::make_shared<EMT::Ph3::Capacitor>(**mName + "_cap1", mLogLevel);
		mSubParallelCapacitor1->setParameters(**mParallelCap / 2.);
		mSubParallelCapacitor1->connect(SimNode::List{ SimNode::GND, mTerminals[1]->node() });
		mSubParallelCapacitor1->initialize(mFrequencies);
		mSubParallelCapacitor1->initializeFromNodesAndTerminals(frequency);
	}

	mSLog->debug(
		"\n--debug--"
		"\n seriesRes: {:s}"
		"\n seriesInd: {:s}"
		"\n Impedance: {:s}"
		"\n vInit: {:s}"
		"\n iInit: {:s}",
		Logger::matrixToString(**mSeriesRes),
		Logger::matrixToString(**mSeriesInd),
		Logger::matrixCompToString(impedance),
		Logger::matrixCompToString(vInitABC),
		Logger::matrixCompToString(iInit));

	mSLog->info(
		"\n--- Initialization from powerflow ---"
		"\nVoltage across: {:s}"
		"\nCurrent: {:s}"
		"\nTerminal 0 voltage: {:s}"
		"\nTerminal 1 voltage: {:s}"
		"\nVirtual Node 1 voltage: {:s}"
		"\n--- Initialization from powerflow finished ---",
		Logger::matrixToString(**mIntfVoltage),
		Logger::matrixToString(**mIntfCurrent),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(0)),
		Logger::phasorToString(RMS3PH_TO_PEAK1PH * initialSingleVoltage(1)),
		Logger::phasorToString(mVirtualNodes[0]->initialSingleVoltage()));
	mSLog->flush();
}

void EMT::Ph3::PiLine::mnaInitialize(Real omega, Real timeStep, Attribute<Matrix>::Ptr leftVector) {
	MNAInterface::mnaInitialize(omega, timeStep);
	updateMatrixNodeIndices();
	MNAInterface::List subComps({ mSubSeriesResistor, mSubSeriesInductor });

	mSubSeriesResistor->mnaInitialize(omega, timeStep, leftVector);
	mSubSeriesInductor->mnaInitialize(omega, timeStep, leftVector);
	mRightVectorStamps.push_back(&mSubSeriesInductor->attribute<Matrix>("right_vector")->get());

	mSubParallelResistor0->mnaInitialize(omega, timeStep, leftVector);
	mSubParallelResistor1->mnaInitialize(omega, timeStep, leftVector);
	subComps.push_back(mSubParallelResistor0);
	subComps.push_back(mSubParallelResistor1);

	if ((**mParallelCap)(0,0) > 0) {
		mSubParallelCapacitor0->mnaInitialize(omega, timeStep, leftVector);
		mSubParallelCapacitor1->mnaInitialize(omega, timeStep, leftVector);
		mRightVectorStamps.push_back(&mSubParallelCapacitor0->attribute<Matrix>("right_vector")->get());
		mRightVectorStamps.push_back(&mSubParallelCapacitor1->attribute<Matrix>("right_vector")->get());
		subComps.push_back(mSubParallelCapacitor0);
		subComps.push_back(mSubParallelCapacitor1);
	}
	mMnaTasks.push_back(std::make_shared<MnaPreStep>(*this));
	mMnaTasks.push_back(std::make_shared<MnaPostStep>(*this, leftVector));
	**mRightVector = Matrix::Zero(leftVector->get().rows(), 1);
}

void EMT::Ph3::PiLine::mnaApplySystemMatrixStamp(Matrix& systemMatrix) {
	mSubSeriesResistor->mnaApplySystemMatrixStamp(systemMatrix);
	mSubSeriesInductor->mnaApplySystemMatrixStamp(systemMatrix);

	mSubParallelResistor0->mnaApplySystemMatrixStamp(systemMatrix);
	mSubParallelResistor1->mnaApplySystemMatrixStamp(systemMatrix);

	if ((**mParallelCap)(0,0) > 0) {
		mSubParallelCapacitor0->mnaApplySystemMatrixStamp(systemMatrix);
		mSubParallelCapacitor1->mnaApplySystemMatrixStamp(systemMatrix);
	}
}

void EMT::Ph3::PiLine::mnaApplyRightSideVectorStamp(Matrix& rightVector) {
	rightVector.setZero();
	for (auto stamp : mRightVectorStamps)
		rightVector += *stamp;
}

void EMT::Ph3::PiLine::mnaAddPreStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes){
	// add pre-step dependencies of subcomponents
	mSubSeriesResistor->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
	mSubSeriesInductor->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
	mSubParallelResistor0->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
	mSubParallelResistor1->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);

	if ((**mParallelCap)(0, 0) > 0) {
		mSubParallelCapacitor0->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
		mSubParallelCapacitor1->mnaAddPreStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes);
	}
	// add pre-step dependencies of component itself
	prevStepDependencies.push_back(attribute("i_intf"));
	prevStepDependencies.push_back(attribute("v_intf"));
	modifiedAttributes.push_back(attribute("right_vector"));
}

void EMT::Ph3::PiLine::mnaPreStep(Real time, Int timeStepCount) {
	// pre-step of subcomponents
	mSubSeriesInductor->mnaPreStep(time, timeStepCount);
	if ((**mParallelCap)(0, 0) > 0) {
		mSubParallelCapacitor0->mnaPreStep(time, timeStepCount);
		mSubParallelCapacitor1->mnaPreStep(time, timeStepCount);
	}
	// pre-step of component itself
	mnaApplyRightSideVectorStamp(**mRightVector);
}

void EMT::Ph3::PiLine::mnaAddPostStepDependencies(AttributeBase::List &prevStepDependencies, AttributeBase::List &attributeDependencies, AttributeBase::List &modifiedAttributes, Attribute<Matrix>::Ptr &leftVector) {
	// add post-step dependencies of subcomponents
	mSubSeriesResistor->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
	mSubSeriesInductor->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
	if ((**mParallelCap)(0, 0) > 0) {
		mSubParallelCapacitor0->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
		mSubParallelCapacitor1->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
	}
	mSubParallelResistor0->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
	mSubParallelResistor1->mnaAddPostStepDependencies(prevStepDependencies, attributeDependencies, modifiedAttributes, leftVector);
	// add post-step dependencies of component itself
	attributeDependencies.push_back(leftVector);
	modifiedAttributes.push_back(attribute("v_intf"));
	modifiedAttributes.push_back(attribute("i_intf"));
}

void EMT::Ph3::PiLine::mnaPostStep(Real time, Int timeStepCount, Attribute<Matrix>::Ptr &leftVector) {
	// post-step of subcomponents
	mSubSeriesResistor->mnaPostStep(time, timeStepCount, leftVector);
	mSubSeriesInductor->mnaPostStep(time, timeStepCount, leftVector);
	if ((**mParallelCap)(0, 0) > 0) {
		mSubParallelCapacitor0->mnaPostStep(time, timeStepCount, leftVector);
		mSubParallelCapacitor1->mnaPostStep(time, timeStepCount, leftVector);
	}
	mSubParallelResistor0->mnaPostStep(time, timeStepCount, leftVector);
	mSubParallelResistor1->mnaPostStep(time, timeStepCount, leftVector);
	// post-step of component itself
	mnaUpdateVoltage(**leftVector);
	mnaUpdateCurrent(**leftVector);
}

void EMT::Ph3::PiLine::mnaUpdateVoltage(const Matrix& leftVector) {
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
}

void EMT::Ph3::PiLine::mnaUpdateCurrent(const Matrix& leftVector) {
	**mIntfCurrent = mSubSeriesInductor->intfCurrent();
}

// #### DAE functions ####

void EMT::Ph3::PiLine::daeInitialize(double time, double state[], double dstate_dt[], 
	double absoluteTolerances[], double stateVarTypes[], int& offset) {
	// state variables are: 3xinductor_current

	updateMatrixNodeIndices();
	mConductance = **mParallelCond * 1. / 2.;
	mCapacitance = **mParallelCap * 1. / 2.;

	// initialize inductor variables
	Matrix inductorVoltage = mSubSeriesInductor->intfVoltage();
	state[offset] = (**mIntfCurrent)(0,0);
	dstate_dt[offset]   = inductorVoltage(0,0) / (**mSeriesInd)(0,0);
	state[offset+1] = (**mIntfCurrent)(1,0);
	dstate_dt[offset+1] = inductorVoltage(1,0) / (**mSeriesInd)(1,1);
	state[offset+2] = (**mIntfCurrent)(2,0);
	dstate_dt[offset+2] = inductorVoltage(2,0) / (**mSeriesInd)(2, 2);

	// set inductor state variables as differential variable
	stateVarTypes[offset+0] = 0.0;
	stateVarTypes[offset+1] = 0.0;
	stateVarTypes[offset+2] = 0.0;

	// set absolute tolerance
	absoluteTolerances[offset]   = mAbsTolerance;
	absoluteTolerances[offset+1] = mAbsTolerance;
	absoluteTolerances[offset+2] = mAbsTolerance;

	mSLog->info(
		"\n--- daeInitialize ---"
		"\nAdded current-phase1 through the inductor of PiLine '{:s}' to state vector, initial value={:f}A"
		"\nAdded current-phase2 through the inductor of PiLine '{:s}' to state vector, initial value={:f}A"
		"\nAdded current-phase3 through the inductor of PiLine '{:s}' to state vector, initial value={:f}A"
		"\nAdded derivative of current-phase1 through the inductor of PiLine '{:s}' to derivative state vector, initial value={:f}"
		"\nAdded derivative of current-phase2 through the inductor of PiLine '{:s}' to derivative state vector, initial value={:f}"
		"\nAdded derivative of current-phase3 through the inductor of PiLine '{:s}' to derivative state vector, initial value={:f}"
		"\nState variables of inductors set as differential"
		"\nAbsolute tolerances={:f}"	
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
	offset+=3;
}

void EMT::Ph3::PiLine::daeResidual(double sim_time,
	const double state[], const double dstate_dt[],
	double resid[], std::vector<int>& off) {
	// offset+3, offset+4, offset+5 --> cap0 --> node1	(left node)
	// offset+3, offset+4, offset+5 --> cap1 --> node0  (right node)

	// current offset for component
	int c_offset = off[0]+off[1]; 

	// residual function of inductors: v1-(v0+vr) - L*di(t)/dt = 0
	resid[c_offset]    = -((**mSeriesInd)(0, 0) * dstate_dt[c_offset+0] + state[c_offset+0] * (**mSeriesRes)(0, 0));
	resid[c_offset+1]  = -((**mSeriesInd)(1, 1) * dstate_dt[c_offset+1] + state[c_offset+1] * (**mSeriesRes)(1, 1));
	resid[c_offset+2]  = -((**mSeriesInd)(2, 2) * dstate_dt[c_offset+2] + state[c_offset+2] * (**mSeriesRes)(2, 2));
	if (terminalNotGrounded(0)) {
		resid[c_offset]   -= state[matrixNodeIndex(0, 0)];
		resid[c_offset+1] -= state[matrixNodeIndex(0, 1)];
		resid[c_offset+2] -= state[matrixNodeIndex(0, 2)];

		// update residual equations of nodes
		resid[matrixNodeIndex(0, 0)] -= state[c_offset];
		resid[matrixNodeIndex(0, 1)] -= state[c_offset+1];
		resid[matrixNodeIndex(0, 2)] -= state[c_offset+2];
	}
	if (terminalNotGrounded(1)) {
		resid[c_offset]   += state[matrixNodeIndex(1, 0)];
		resid[c_offset+1] += state[matrixNodeIndex(1, 1)];
		resid[c_offset+2] += state[matrixNodeIndex(1, 2)];

		// update residual equations of nodes
		resid[matrixNodeIndex(1, 0)] += state[c_offset];
		resid[matrixNodeIndex(1, 1)] += state[c_offset+1];
		resid[matrixNodeIndex(1, 2)] += state[c_offset+2];
	}

	//add cap and cond currents to nodal equations
	
	resid[matrixNodeIndex(1, 0)] += dstate_dt[matrixNodeIndex(1, 0)] * mCapacitance(0,0) + state[matrixNodeIndex(1, 0)] * mConductance(0,0);
	resid[matrixNodeIndex(1, 1)] += dstate_dt[matrixNodeIndex(1, 1)] * mCapacitance(1,1) + state[matrixNodeIndex(1, 1)] * mConductance(1,1);
	resid[matrixNodeIndex(1, 2)] += dstate_dt[matrixNodeIndex(1, 2)] * mCapacitance(2,2) + state[matrixNodeIndex(1, 2)] * mConductance(2,2);
	resid[matrixNodeIndex(0, 0)] += dstate_dt[matrixNodeIndex(0, 0)] * mCapacitance(0,0) + state[matrixNodeIndex(0, 0)] * mConductance(0,0);
	resid[matrixNodeIndex(0, 1)] += dstate_dt[matrixNodeIndex(0, 1)] * mCapacitance(1,1) + state[matrixNodeIndex(0, 1)] * mConductance(1,1);
	resid[matrixNodeIndex(0, 2)] += dstate_dt[matrixNodeIndex(0, 2)] * mCapacitance(2,2) + state[matrixNodeIndex(0, 2)] * mConductance(2,2);
	
	off[1] += 3;
}

void EMT::Ph3::PiLine::daeJacobian(double current_time, const double state[], 
	const double dstate_dt[], SUNMatrix jacobian, double cj, std::vector<int>& off) {

	// current offset for component
	int c_offset = off[0] + off[1]; 

	// residual funcion (series part)
	SM_ELEMENT_D(jacobian, c_offset,   c_offset)   -= (**mSeriesRes)(0, 0) + cj * (**mSeriesInd)(0, 0);
	SM_ELEMENT_D(jacobian, c_offset+1, c_offset+1) -= (**mSeriesRes)(1, 1) + cj * (**mSeriesInd)(1, 1);
	SM_ELEMENT_D(jacobian, c_offset+2, c_offset+2) -= (**mSeriesRes)(2, 2) + cj * (**mSeriesInd)(2, 2);

	if (terminalNotGrounded(1)) {
		// residual funcion (series branch)
		SM_ELEMENT_D(jacobian, c_offset,   matrixNodeIndex(1,0)) += 1.0;
		SM_ELEMENT_D(jacobian, c_offset+1, matrixNodeIndex(1,1)) += 1.0;
		SM_ELEMENT_D(jacobian, c_offset+2, matrixNodeIndex(1,2)) += 1.0;

		// node 1 (series branch)
		SM_ELEMENT_D(jacobian, matrixNodeIndex(1,0), c_offset)   += 1.0;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(1,1), c_offset+1) += 1.0;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(1,2), c_offset+2) += 1.0;

		// node 1 (parallel branch)
		SM_ELEMENT_D(jacobian, matrixNodeIndex(1,0), matrixNodeIndex(1,0)) += mConductance(0,0) + cj * mCapacitance(0,0);
		SM_ELEMENT_D(jacobian, matrixNodeIndex(1,1), matrixNodeIndex(1,1)) += mConductance(1,1) + cj * mCapacitance(1,1);
		SM_ELEMENT_D(jacobian, matrixNodeIndex(1,2), matrixNodeIndex(1,2)) += mConductance(2,2) + cj * mCapacitance(2,2);
	}
	if (terminalNotGrounded(0)) {
		// residual funcion (series branch)
		SM_ELEMENT_D(jacobian, c_offset,   matrixNodeIndex(0,0)) += -1.0;
		SM_ELEMENT_D(jacobian, c_offset+1, matrixNodeIndex(0,1)) += -1.0;
		SM_ELEMENT_D(jacobian, c_offset+2, matrixNodeIndex(0,2)) += -1.0;

		// node 0 (series branch)
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0,0), c_offset)   += -1.0;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0,1), c_offset+1) += -1.0;
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0,2), c_offset+2) += -1.0;

		// node 0 (parallel branch)
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0,0), matrixNodeIndex(0,0)) += mConductance(0,0) + cj * mCapacitance(0,0);
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0,1), matrixNodeIndex(0,1)) += mConductance(1,1) + cj * mCapacitance(1,1);
		SM_ELEMENT_D(jacobian, matrixNodeIndex(0,2), matrixNodeIndex(0,2)) += mConductance(2,2) + cj * mCapacitance(2,2);
	}

	off[1] += 3;
}

void EMT::Ph3::PiLine::daePostStep(double Nexttime, const double state[], const double dstate_dt[], int& offset) {
	// update voltage
	(**mIntfVoltage)(0,0) = state[matrixNodeIndex(1, 0)] - state[matrixNodeIndex(0, 0)];
	(**mIntfVoltage)(1,0) = state[matrixNodeIndex(1, 1)] - state[matrixNodeIndex(0, 1)];
	(**mIntfVoltage)(2,0) = state[matrixNodeIndex(1, 2)] - state[matrixNodeIndex(0, 2)];

	// update current
	(**mIntfCurrent)(0, 0) = state[offset];
	(**mIntfCurrent)(1, 0) = state[offset+1];
	(**mIntfCurrent)(2, 0) = state[offset+2];

	offset+=3;
}


