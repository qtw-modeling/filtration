#ifndef __COMMON_PHASE_BALANCE_H__
#define __COMMON_PHASE_BALANCE_H__

#include "Props.h"


/** PhaseSolver termination reason */
enum Status 
{
	Ok = 0,
	StepBecameTooSmall = 1,
	MaximumNumberOfIterationsExceeded = 2,
	StepFactorIsZero = 3
};


#endif
