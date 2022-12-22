/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/EMT/EMT_Ph3_SynchronGenerator3OrderVBR.h>

using namespace CPS;

EMT::Ph3::SynchronGenerator3OrderVBR::SynchronGenerator3OrderVBR
    (String uid, String name, Logger::Level logLevel)
	: ReducedOrderSynchronGeneratorVBR(uid, name, logLevel),
	mEdq0_t(Attribute<Matrix>::create("Edq0_t", mAttributes)) {

	// model specific variables
	**mEdq0_t = Matrix::Zero(3,1);
	mEhs_vbr = Matrix::Zero(3,1);
}

EMT::Ph3::SynchronGenerator3OrderVBR::SynchronGenerator3OrderVBR
	(String name, Logger::Level logLevel)
	: SynchronGenerator3OrderVBR(name, name, logLevel) {
}

SimPowerComp<Real>::Ptr EMT::Ph3::SynchronGenerator3OrderVBR::clone(String name) {
	
	auto copy = SynchronGenerator3OrderVBR::make(name, mLogLevel);
	return copy;
}

void EMT::Ph3::SynchronGenerator3OrderVBR::setOperationalParametersPerUnit(Real nomPower, 
			Real nomVolt, Real nomFreq, Real H, Real Ld, Real Lq, Real L0,
			Real Ld_t, Real Td0_t) {

	Base::ReducedOrderSynchronGenerator<Real>::setOperationalParametersPerUnit(nomPower, 
		nomVolt, nomFreq, H, Ld, Lq, L0,
		Ld_t, Td0_t);
	
	mSLog->info("Set base parameters: \n"
				"nomPower: {:e}\nnomVolt: {:e}\nnomFreq: {:e}\n",
				nomPower, nomVolt, nomFreq);

	mSLog->info("Set operational parameters in per unit: \n"
			"inertia: {:e}\n"
			"Ld: {:e}\nLq: {:e}\nL0: {:e}\n"
			"Ld_t: {:e}\nTd0_t: {:e}\n",
			H, Ld, Lq, L0, 
			Ld_t, Td0_t);
}

void EMT::Ph3::SynchronGenerator3OrderVBR::specificInitialization() {
	// initial voltage behind the transient reactance in the dq0 reference frame
	(**mEdq0_t)(1,0) = (**mVdq0)(1,0) + (**mIdq0)(0,0) * mLd_t;

	// calculate auxiliar VBR constants
	calculateAuxiliarConstants();

	// dq0 resistance matrix
	mResistanceMatrixDq0 = Matrix::Zero(3,3);
	mResistanceMatrixDq0 <<	0.0,			-mLq,	0.0,
							mLd_t - mAq,	0.0,	0.0,
					  		0.0,			0.0,	mL0;

	// initialize conductance matrix 
	mConductanceMatrix = Matrix::Zero(3,3);

	mSLog->info(
		"\n--- Model specific initialization  ---"
		"\nInitial Eq_t (per unit): {:f}"
		"\n--- Model specific initialization finished ---",
		(**mEdq0_t)(1,0)
	);
	mSLog->flush();
}

void EMT::Ph3::SynchronGenerator3OrderVBR::calculateAuxiliarConstants() {
	mAq = - mTimeStep * (mLd - mLd_t) / (2 * mTd0_t + mTimeStep);
	mBq = (2 * mTd0_t - mTimeStep) / (2 * mTd0_t + mTimeStep);
	mCq = 2 * mTimeStep * mEf / (2 * mTd0_t + mTimeStep);
}

void EMT::Ph3::SynchronGenerator3OrderVBR::stepInPerUnit() {

	if (mSimTime>0.0) {
		// calculate Eq_t at t=k
		(**mEdq0_t)(1,0) = (**mIdq0)(0,0) * mLd_t + (**mVdq0)(1,0);

		// calculate mechanical variables at t=k+1 with forward euler
		**mElecTorque = ((**mVdq0)(0,0) * (**mIdq0)(0,0) + (**mVdq0)(1,0) * (**mIdq0)(1,0));
		**mOmMech = **mOmMech + mTimeStep * (1. / (2. * mH) * (**mMechTorque - **mElecTorque));
		**mThetaMech = **mThetaMech + mTimeStep * (**mOmMech * mBase_OmMech);
		**mDelta = **mDelta + mTimeStep * (**mOmMech - 1.) * mBase_OmMech;
	}

	// get transformation matrix
	mAbcToDq0 = get_parkTransformMatrix();
	mDq0ToAbc = get_inverseParkTransformMatrix();

	// calculate resistance matrix at t=k+1
	calculateResistanceMatrix();

	// VBR history voltage
	mEhs_vbr(0,0) = 0.0;
	mEhs_vbr(1,0) = mAq * (**mIdq0)(0,0) + mBq * (**mEdq0_t)(1,0) + mCq;
	mEhs_vbr(2,0) = 0.0;

	// convert Edq_t into the abc reference frame
	**mEvbr = mDq0ToAbc * mEhs_vbr * mBase_V;
}


// #### DAE functions ####

void EMT::Ph3::SynchronGenerator3OrderVBR::daeInitialize(double time, double state[], double dstate_dt[], 
	double absoluteTolerances[], double stateVarTypes[], int& offset) {
	// state variables are: Eq_t, omega and mechanical theta
	updateMatrixNodeIndices();

	// initial voltage behind the transient reactance in the dq0 reference frame
	(**mEdq0_t)(1,0) = (**mVdq0)(1,0) + (**mIdq0)(0,0) * mLd_t;

	// init state variables
	state[offset] = (**mEdq0_t)(1,0);
	dstate_dt[offset] = 0.0;
	state[offset+1] = **mOmMech;
	dstate_dt[offset+1] = 0.0;
	state[offset+2] = **mThetaMech;
	dstate_dt[offset+2] = mBase_OmMech * state[offset+1];
	state[offset+3] = (**mIdq0)(0,0);
	dstate_dt[offset+3] = 0.0;
	state[offset+4] = (**mIdq0)(1,0);
	dstate_dt[offset+4] = 0.0;

	// set state variables as differential variable
	stateVarTypes[offset+0] = 0.0;
	stateVarTypes[offset+1] = 0.0;
	stateVarTypes[offset+2] = 0.0;
	stateVarTypes[offset+3] = 1.0;
	stateVarTypes[offset+4] = 1.0;

	// set absolute tolerance
	absoluteTolerances[offset]   = mAbsTolerance;
	absoluteTolerances[offset+1] = mAbsTolerance;
	absoluteTolerances[offset+2] = mAbsTolerance;
	absoluteTolerances[offset+3] = mAbsTolerance;
	absoluteTolerances[offset+4] = mAbsTolerance;

	mSLog->info(
		"\n--- daeInitialize ---"
		"\nInitial Eq_t = {:f} [p.u.]"
		"\nInitial mechanical omega = {:f} [p.u.]"
		"\nInitial mechanical theta = {:f} [p.u.]"
		"\nInitial derivative of mechanical theta = {:f} [p.u.]"
		"\n--- daeInitialize finished ---",
		state[offset],
		state[offset+1],
		state[offset+2],
		dstate_dt[offset+2]
	);
	mSLog->flush();
	offset+=5;
}

void EMT::Ph3::SynchronGenerator3OrderVBR::daeResidual(double sim_time, 
	const double state[], const double dstate_dt[], 
	double resid[], std::vector<int>& off) {

	// current offset for component	
	int c_offset = off[0] + off[1]; 

	// get transformation matrix
	**mThetaMech = state[c_offset+2];
	mAbcToDq0 = get_parkTransformMatrix();
	mDq0ToAbc = get_inverseParkTransformMatrix();

	// convert terminal voltage into dq0 reference frame
	(**mIntfVoltage)(0, 0) = state[matrixNodeIndex(0, 0)];
	(**mIntfVoltage)(1, 0) = state[matrixNodeIndex(0, 1)];
	(**mIntfVoltage)(2, 0) = state[matrixNodeIndex(0, 2)];
	**mVdq0 = mAbcToDq0 * **mIntfVoltage / mBase_V;

	// calculate dq0 and abc currents
	//(**mEdq0_t)(1,0) = state[c_offset];
	//(**mIdq0)(0,0) = ((**mEdq0_t)(1,0) - (**mVdq0)(1,0)) / mLd_t;
	//(**mIdq0)(1,0) = (**mVdq0)(0,0) / mLq;
	(**mIdq0)(0,0) = state[c_offset+3];
	(**mIdq0)(1,0) = state[c_offset+4];
	**mIntfCurrent = mBase_I * mDq0ToAbc * **mIdq0;
	
	// residual function for Eq_t
	resid[c_offset] = mTd0_t * dstate_dt[c_offset] + state[c_offset] - mEf + (mLd - mLd_t) * (**mIdq0)(0,0);

	// residual function for omega
	**mElecTorque = ((**mVdq0)(0,0) * (**mIdq0)(0,0) + (**mVdq0)(1,0) * (**mIdq0)(1,0));
	resid[c_offset+1] = 2 * mH * dstate_dt[c_offset+1] - **mMechTorque + **mElecTorque;

	// residual function for mechanical theta
	resid[c_offset+2] = dstate_dt[c_offset+2] - mBase_OmMech * state[c_offset+1];

	// residual function of currents
	resid[c_offset+3] = state[c_offset+3] * mLd_t - state[c_offset] + (**mVdq0)(1,0);
	resid[c_offset+4] = state[c_offset+4] * mLq - (**mVdq0)(0,0);

	// add currents to residual funcion of nodes
	resid[matrixNodeIndex(0, 0)] -= (**mIntfCurrent)(0, 0);
	resid[matrixNodeIndex(0, 1)] -= (**mIntfCurrent)(1, 0);
	resid[matrixNodeIndex(0, 2)] -= (**mIntfCurrent)(2, 0);
	off[1] += 5;
}

void EMT::Ph3::SynchronGenerator3OrderVBR::daeJacobian(double current_time, const double state[], 
	const double dstate_dt[], SUNMatrix jacobian, double cj, std::vector<int>& off) {

	// TODO

}

void EMT::Ph3::SynchronGenerator3OrderVBR::daePostStep(double Nexttime, const double state[], 
	const double dstate_dt[], int& offset) {
	
	(**mEdq0_t)(1,0) = state[offset];
	**mOmMech = state[offset+1];
	**mThetaMech = state[offset+2];

	// get transformation matrix
	mAbcToDq0 = get_parkTransformMatrix();
	mDq0ToAbc = get_inverseParkTransformMatrix();

	(**mIntfVoltage)(0, 0) = state[matrixNodeIndex(0, 0)];
	(**mIntfVoltage)(1, 0) = state[matrixNodeIndex(0, 1)];
	(**mIntfVoltage)(2, 0) = state[matrixNodeIndex(0, 2)];
	**mVdq0 = mAbcToDq0 * **mIntfVoltage / mBase_V;

	// calculate dq0 and abc currents
	//(**mIdq0)(0,0) = ((**mEdq0_t)(1,0) - (**mVdq0)(1,0)) / mLd_t;
	//(**mIdq0)(1,0) = (**mVdq0)(0,0) / mLq;
	(**mIdq0)(0,0) = state[offset+3];
	(**mIdq0)(1,0) = state[offset+4];
	**mIntfCurrent = mBase_I * mDq0ToAbc * **mIdq0;

	offset+=5;
}