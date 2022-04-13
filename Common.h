#ifndef __COMMON_H__
#define __COMMON_H__


#include <vector>
#include <cmath>
#include "Array2D.hpp"

typedef double type;

enum Phases {Liquid = 1, Water = 2, Gas = 3};
enum Components {HO = 1, LO = 2, O2 = 3, CO2 = 4, H2O = 5};
enum Approximation {leftChoice, rightChoice, upwind, arithmeticAverage, geometricAverage};
enum Layer {current, middle, next};

// number of phases and components
const int NUM_COMPS  = 5;
const int NUM_PHASES = 3;

#endif
