#ifndef VARIABLES_HPP_
#define VARIABLES_HPP_

#include <stdlib.h>
#include <algorithm>
#include <cmath>
#include "Common.h"
#include "Filt.h"

using namespace std;

//constant parameters
const type ksi1  = 23;
const type ksi2  = 15.5;
const type heta1 = 15;
const type heta2 = 10;
const type nu1   = 16;
const type nu2   = 11;

const type alphaHO = 1;
const type betaHO  = 1;
const type alphaMO = 1;
const type betaMO  = 1;

const type R     = 8.31;
const type E_HO  = 40e3;
const type E_MO  = 40e3;     //change!
const type A_HO  = 1e2;
const type A_MO  = 1e2;      //change!
const type A_KER = 100;      //change!
const type E_KER = 5e4;

const type T_ACTIV_HO = 300;
const type T_ACTIV_MO = 300;

const type eps = 1e-7;


struct Densities
{
	vector<type> N;
	Densities() : N(NUM_COMPS + 1) { }

};


struct AbstractReaction {
	
	type predExponent, eActivation;

	AbstractReaction(type predExponent_, type eActivation_): predExponent(predExponent_), eActivation(eActivation_) {}

	virtual type SpeedOfReaction(type T) const = 0;
};


struct ReactionWithoutActivation : public AbstractReaction {

	ReactionWithoutActivation(type predExponent_, type eActivation_): AbstractReaction(predExponent_, eActivation_) {}

	type SpeedOfReaction(type T) const { return predExponent * exp(-eActivation / (R * T)); }
};



struct ReactionWithActivation : public AbstractReaction {

	const type Tactivation;

	ReactionWithActivation(type predExponent_, type eActivation_, type Tactivation_): Tactivation(Tactivation_), AbstractReaction(predExponent_, eActivation_) {}

	type SpeedOfReaction(type T) const {
		if (T >= Tactivation)
			return predExponent * exp(-eActivation / (R * T));
		else
			return 0;
	}
};


struct ChemicalReactions
{
	ReactionWithoutActivation HO;
	ReactionWithoutActivation LO;
	
	ChemicalReactions(): HO( /* pred exponent */ 1e-2, /* activation energy */ 40e3), LO(1e-2, 40e3) {}

	//Reaction Kerogene(TOSET, TOSET);	//!!!
};


struct ModelChemicalReactions
{
	ReactionWithoutActivation Fuel;

	ChemicalReactions(): Fuel (1e-2, 40e3) {}
};




struct AbstractOneStepChemicalSolver
{
	Densities oldDensities, newDensities;

	void SetInitialConditions(Densities& initialDensities);

	virtual Densities CalculateNewValuesInCell(type tau, type T, type P) = 0;
};


struct ExplicitImplicitChemicalSolver: public AbstractOneStepChemicalSolver
{
 	Densities CalculateNewValuesInCell(type tau, type T, type P);
};


struct PureImplicitChemicalSolver: public AbstractOneStepChemicalSolver
{
	Densities CalculateNewValuesInCell(type tau, type T, type P);
};


struct ExplicitImplicitChemicalSolverMod : public AbstractOneStepChemicalSolver
{
	Densities CalculateNewValuesInCell(type tau, type T, type P);
};

struct ExplicitImplicitChemicalSolverMod2 : public AbstractOneStepChemicalSolver
{
	Densities CalculateNewValuesInCell(type tau, type T, type P);
};


struct roots
{
     type x1, x2, x3;
     roots(type a, type b, type c)
     {
       if (a > b)
         swap(a, b);
       if (b > c)
         swap(b, c);
       if (a > b)
         swap(a, b);

           x1 = a;
           x2 = b;
           x3 = c;
     }

};


Densities CalculateTwoReactionsInCell(type tau, type nHO, type nMO, type nO2, type nCO2, type nH2O, type T, type P);
roots cubic_roots(type a, type b, type c);
type SpeedOfReactionKerogene(type T);
type alphaKerogene(type T);

// other functions
type Max2(const type& a, const type& b);
type Max3(const type& a, const type& b, const type& c);
#endif /* VARIABLES_HPP_ */
