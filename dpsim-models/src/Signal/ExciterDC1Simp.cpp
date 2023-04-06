/* Copyright 2017-2021 Institute for Automation of Complex Power Systems,
 *                     EONERC, RWTH Aachen University
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *********************************************************************************/

#include <dpsim-models/Definitions.h>
#include <dpsim-models/Signal/ExciterDC1Simp.h>
#include <dpsim-models/MathUtils.h>

using namespace CPS;
using namespace CPS::Signal;

ExciterDC1Simp::ExciterDC1Simp(const String & name, CPS::Logger::Level logLevel) 
	: SimSignalComp(name, name, logLevel) { 

    this->setExciterType(ExciterType::DC1Simp);
}

void ExciterDC1Simp::setParameters(Base::ExciterParameters parameters) {

	mTa = parameters.Ta;
	mTef = parameters.Tef;
	mTf = parameters.Tf;
	mTr = parameters.Tr;
	mKa = parameters.Ka;
	mKef = parameters.Kef;
	mKf = parameters.Kf;
	mAef = parameters.Aef;
	mBef = parameters.Bef;
	mMaxVa = parameters.MaxVa;
	mMinVa = parameters.MinVa;

	SPDLOG_LOGGER_INFO(mSLog, "Exciter parameters: \n"
				"Ta: {:e}"
				"\nKa: {:e}"
				"\nTef: {:e}"
				"\nKef: {:e}"
				"\nTf: {:e}"
				"\nKf: {:e}"
				"\nTr: {:e}"
				"\nAef: {:e}"
				"\nBef: {:e}"
				"\nMaximum amplifier Voltage: {:e}"
				"\nMinimum amplifier Voltage: {:e}\n",
				mTa, mKa, 
				mTef, mKef,
				mTf, mKf,
				mTr,
				mAef, mBef,
				mMaxVa, mMinVa);
	mSLog->flush();
}

void ExciterDC1Simp::initialize(Real Vh_init, Real Ef_init) {
	///
	mVh = Vh_init;
	mEf = Ef_init;
	SPDLOG_LOGGER_INFO(mSLog, "Initially set excitation system initial values: \n"
				"\ninit Vh: {:e}"
				"\ninit Ef: {:e}",
				mVh, mEf);

	/// init value of transducer output
	mVr = mVh;

	/// init value of stabilizing feedback output
	mVf = 0.0;

	/// ceiling function
	mVsat = mAef * exp(mBef * abs(mEf));

<<<<<<< HEAD
	//
	mVref = **mVr / mKa + **mVm;

	SPDLOG_LOGGER_INFO(mSLog, "Actually applied excitation system initial values:"
=======
	/// init value of amplifier output
	mVa = mKef * mEf + mVsat * mEf;
	if (mVa>mMaxVa)
		mVa = mMaxVa;
	if (mVa<mMinVa)
		mVa = mMinVa;

	/// init value of amplifier input 
	mVin = mVa / mKa;

	///
	mVref = mVr + mVin;

	/// check initial conditions
	if (mEf - mVa / (mVsat + mKef))
		mSLog->warn("Initial conditions are not consistent!!!");
	
	mSLog->info("Actually applied excitation system initial values:"
>>>>>>> 32920b70 (add exciter DC1)
				"\nVref : {:e}"
				"\ninit_Vr: {:e}"
				"\ninit_Ef: {:e}"
				"\ninit_Va: {:e}",
				mVref,
				mVr, 
				mEf,
				mVa);
	mSLog->flush();
}

Real ExciterDC1Simp::step(Real mVd, Real mVq, Real dt, Real Vpss) {
	// Voltage magnitude calculation
	mVh = sqrt(pow(mVd, 2.) + pow(mVq, 2.));

	// update state variables at time k-1
	mVr_prev = mVr;
	mVa_prev = mVa;
	mVf_prev = mVf;
	mEf_prev = mEf;

	// compute state variables at time k using euler forward

	// saturation function
	mVsat = mAef * exp(mBef * abs(mEf_prev));

	// Voltage Transducer equation
	mVr = mVr_prev + dt / mTr * (mVh - mVr_prev);

	// Voltage amplifier equation
	mVin = mVref + Vpss - mVr_prev - mVf_prev;
	mVa = mVa_prev + dt / mTa * (mVin * mKa - mVa_prev);
	if (mVa > mMaxVa)
		mVa = mMaxVa;
	else if (mVa < mMinVa)
		mVa = mMinVa;

	// Stabilizing feedback
	mVf = (1. - dt / mTf) * mVf_prev + dt * mKf / (mTf * mTef) * (mVa_prev - (mVsat + mKef) * mEf_prev);
	
	// Exciter output
	mEf = mEf_prev + dt / mTef * (mVa_prev - (mVsat + mKef) * mEf_prev); 

	return mEf;
}