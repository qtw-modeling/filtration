#ifndef FILT_H
#define FILT_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <algorithm>
#include <cstdlib>
#include "Util.h"
#include "assert.h"
#include "Common.h"
#include "PhaseSolver.h"
#include "ChemicalKinetics.h"
#include "Well.hpp"
#include "CompoundParameters.hpp"

using namespace std;

const type eps1 = 1e-7;
const type eps2 = 1e-3;

// parameters of the system
const type rhol = 1e3;
const type rhow = 1e3;
const type rhog = 1e3;
const type rhoS = 1e3;
const type f    = 1e-1;		// porosity
const type mul  = 1e-3;
const type muw  = 1e-3;
const type mug  = 1e-3;
//const type K    = 1e-13;	// TO CHANGE!!!
const type cl   = 1e3;
const type cw   = 1e3;
const type cg   = 1e3;		// TO CHANGE!!!
const type cPS  = 0;

// grid parameters
const type lengthOfAreaX = 1e2;
const type lengthOfAreaY = 1e2;


class Composition : public Array2D<type> {
public:
	Composition() : Array2D<type>(NUM_COMPS + 1, NUM_PHASES + 1) {}
};

struct Function1VariablePtr {
	Function1Variable *p;
	Function1VariablePtr() : p(0) { }
	Function1VariablePtr(Function1Variable *p) : p(p) { }
	Function1VariablePtr &operator=(const Function1VariablePtr &other) {
		p = other.p;
		return *this;
	}
	type operator()(type T) {
		return (*p)(T);
	}
};

//CLASS DEFINITION
class Filt2D
{
public:

		type tau;	// timestep
		const int numBlocksX, numBlocksY;
		const int numPointsX, numPointsY;
		const type hx;
		const type hy; 

//DEFINING THE STRUCTURES

		struct component
		{
			component(int numPointsX_, int numPointsY_, type vMolar_, type cPL_, type cPW_, type cPG_, type calorisity_) :
				omega(numPointsX_, numPointsY_),
				Nn(numPointsX_, numPointsY_), Nn1(numPointsX_, numPointsY_), Nmid(numPointsX_, numPointsY_), 
				Nwell(numPointsX_, numPointsY_), C(numPointsX_, numPointsY_),
				fluxXmid(numPointsX_ - 1, numPointsY_ - 1), fluxYmid(numPointsX_ - 1, numPointsY_ - 1)
			{
				vMolar     = vMolar_;
				calorisity = calorisity_;
				cP[Liquid] = cPL_;
				cP[Water]  = cPW_;
				cP[Gas]    = cPG_;
			}

			// default constructor
			component() {}

			Array2D <type> omega;
			type calorisity;
			type vMolar, cP[NUM_PHASES + 1];
			Array2D <type> Nn, Nn1, Nmid, Nwell, C;
			Array2D <type> fluxXmid, fluxYmid;
		};

		struct phase
		{
			phase(int numPointsX_, int numPointsY_): 
				rho(numPointsX_, numPointsY_), S(numPointsX_, numPointsY_),
				filtSpeedXmid(numPointsX_ - 1, numPointsY_ - 1), 
				filtSpeedYmid(numPointsX_ - 1, numPointsY_ - 1) 
				{}

			// default constructor
			phase() {}
			~phase() {
				delete mu.p;
			}

			void setViscosity(type A, type B, type C) { mu = new ConstantViscosity(A); }


			Array2D <type> rho, S;
			Function1VariablePtr mu;
			type phaseTransitionQ, cP, molarChemicalEnergyOfComponent[NUM_COMPS + 1];
			Array2D <type> filtSpeedXmid, filtSpeedYmid;
		};

//VARIABLES

		component Component[NUM_COMPS + 1];
		phase Phase[NUM_PHASES + 1];

		// matrix of the coeffs x(i)(alpha)
		Array2D <Composition> x; 
		// molar densities
		Array2D <type> lambda, omega, gamma;
		Array2D <type> fluxEnthalpyXmid, fluxEnthalpyYmid;
		Array2D <type> fluxHeatConductivityXmid, fluxHeatConductivityYmid;

		Array2D <type> P, Tn, Tn1, Hn, Hwell, Hn1, En, En1;
		Array2D <type> K;

		int counter; 	// number of timestep
		enum Approximation approximation;

//FUNCTIONS

		Filt2D(int, int);
		~Filt2D();
		void SetInitialConditions(type);
		//type etaMid(int, int, int, type*, type*);
		//type Tmiddle(int k, type* s, type* P);
		void CalculateSaturationsAndDensitiesOfPhasesInCell(int, int);
		//type CalculateNormOfArray(type* array, int l);
		//void CalculateDeltaN();
		//void CalculateMassTransferMod();
		//void CalculateDensitiesInCell(int l);
		void CalculateOmegasInCell(int, int);
		//void Progonka();
		
		void CalculatePressureByGaussSeidel(Array2D <type>& Pnew, Array2D <type>& Pprevious, type precision);
		void CalculatePhaseBalanceInCell(int l, int m, Array2D <Composition>& x, type eta, type P, type& T);
		void CalculateConcentrationOfComponentInCell(int l, int m, int i, type N, type NtotalInCell);
		int Iteration();
		void Calculate();
		int CalculatePhaseBalanceInCell(int l, int m, Array2D <Composition>& x, type eta, type P);
		type MolarEnthalpyOfComponentInPhaseInCell(int i, int alpha, type T);
		type CalculateTotalMolarEnthalpyInCell(int l, type T, type** x);
		type CalculateComponentFluxAtGridPoints(int i, int l);
		type CalculateCpOfPhase(int alpha, type** x);
		void CalculateTimeStepAuto();
		void CalculateMolarDensityOfComponentInCell(int, int, int, Array2D <type>&);
		void SetApproximation(enum Approximation);
		type MolarEnthalpyOfPhaseInCell(int alpha, type T, Composition& x);

// MASS TRANSFER FUNCTIONS
		void CalculateENOFluxesAtMidPoints();
		void CalculateENOFluxesAtMidPointsMod();
		type CalculateComponentPartFluxAtGridPoints(int i, int l);
		void CalculateRhoXAtMidPointsByENO();
		type rhoXPhi(int alpha, int i, int l);
		type Phi(int alpha, int l);
		type FiltSpeedAtMidPoints(int l);
		type RelativePermeability(type s) const;
		void CalculateFiltrationSpeeds();
		void CalculateMassFluxes();
		void CalculateMassTransfer();
		type UpwindApproximationAtMidPoint(type pL, type pR, type fL, type fR) const;

//APPROXIMATIONS FUNCTIONS

		type ApproximationAtMidPoint(type valueLeft, type valueRight, enum Approximation approx);

// TEMPERATURE FUNCTIONS
		void CalculateEnthalpyAtNextTimeLayer(Array2D <type>& Pcurrent, Array2D <type>& Pnext);
		type fFunction(int, type**);
		type gFunction(int, type**);
		//void IterationTemperature(type*, type*, type*, type*, type*, type*);
		//type CalcVelocityOfTheWave(type*, type*, const type, type);
		inline type CalculateCpOfPhase(int alpha, int l);
		type TemperatureUsingFullEnthalpy(int l, type** sNext, type H);
		void CalculateEnergyUsingEnthalpyInCell(type& E, const type& H, const type& P);
		void CalculateEnthalpyFluxes();
		type AverageMolarHeatCapacityInCell(int l);

// ARRAY FUNCTIONS
		void VisualizeVTK(const int step);

//DEBUG FUNCTIONS
		void DEBUG(type value);
		//OTHER FUNCTIONS PROTOTYPES
};


#endif
