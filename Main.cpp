#include "Filt.h"
#include "Well.hpp"

#define _GNU_SOURCE 1
#define CHEMICAL_REACTIONS
#define TEST

#include <fenv.h>
#include <cstdlib>



int main(int argc, char *argv[])
{
	//_clearfp();
	//_controlfp(0, _EM_OVERFLOW | _EM_ZERODIVIDE | _EM_INVALID); // _controlfp_s is the secure version of _controlfp

	feenableexcept(FE_OVERFLOW | FE_INVALID | FE_DIVBYZERO);

	// creating an object
	Filt2D filt(100, 100);
	ChemicalReactions reactions;
	ExplicitImplicitChemicalSolverMod2 chemSolver;
	Densities oldDensities, newDensities;

	InjectorWell iWell;
	RecoveryWell rWell;

	iWell.setPlace(1 / 5. * lengthOfAreaX, lengthOfAreaY / 2., filt);
	rWell.setPlace(4 / 5. * lengthOfAreaX, lengthOfAreaY / 2., filt);

	filt.SetInitialConditions(1e7);

	
	// setting initial conditions

	//filt.P.print();

	//filt.Phase[Gas].S.print();

	//filt.Component[O2].Nn.print();
	filt.SetApproximation(upwind);
	

	// calculating initial phase balance
	type totalMolarEnthalpyInCell, NtotalInCell;
	
	cycle(l, 0, filt.numPointsX - 1, m, 0, filt.numPointsY - 1) {
			// calculating NTotalInCell
			NtotalInCell = 0;
			for (int i = 1; i <= NUM_COMPS; i++)
			{
				NtotalInCell += filt.Component[i].Nn(l, m);
			}

			totalMolarEnthalpyInCell = filt.Hn(l, m) / NtotalInCell;
			
			filt.CalculatePhaseBalanceInCell(
												l, m, filt.x, totalMolarEnthalpyInCell, filt.P(l, m), filt.Tn(l, m)
											);
			
			filt.CalculateSaturationsAndDensitiesOfPhasesInCell(l, m);

			// calculating molar densities
			for (int i = 1; i <= NUM_COMPS; i++) 
			{
				filt.CalculateMolarDensityOfComponentInCell(l, m, i, filt.Component[i].Nn);
			}
		}
	}

	// tmp print
	//filt.Component[LO].C.print();


	//filt.Print(0);
	cout << "Initial phase balance has been calculated" << endl;
	//MALLOC ARRAYS FOR N AND N+1 LEVELS
	Array2D <type> Pn(filt.numPointsX, filt.numPointsY), Pn1(filt.numPointsX, filt.numPointsY);
	Array2D <type> NO2initial(filt.numPointsX, filt.numPointsY), Tinitial(filt.numPointsX, filt.numPointsY);


	// copying data to the arrays storing initial data
	filt.Component[O2].Nn.copy(NO2initial);
	filt.Tn.copy(Tinitial);


	const type tMax = 10e8;
	type t = 0;
	filt.tau = 2e5; 

	type a, deltaNO2, deltaK, cP;
	uint counter = 0;

	// main cycle ------------------------------------------------------------------------------------------------------------------------
	while (t <= tMax)
	{

		cycle (l, 0, filt.P.lastIndexX(),
			   m, 0, filt.P.lastIndexY()) {
			   		filt.CalculateEnergyUsingEnthalpyInCell(filt.En(l, m), filt.Hn(l, m), filt.P(l, m));
			}
		}

		// copying the data from n_level;
		filt.P.copy(Pn);
		//filt.CopyArray(filt.Tn1, filt.Tn);

		//filt.CalculateTimeStepAuto();


		// calculating chemical kinetics in the area
		for (int l = 0; l <= filt.P.lastIndexX(); l++)
			for (int m = 0; m <= filt.P.lastIndexY(); m++)
			{
				
				// setting oldDensities
				for (int i = 1; i <= NUM_COMPS; i++)
					oldDensities.N[i] = filt.Component[i].Nn(l, m);

				chemSolver.SetInitialConditions(oldDensities);

				for (int i = 1; i <= NUM_COMPS; i++)
					newDensities.N[i] = oldDensities.N[i];
				
				newDensities = chemSolver.CalculateNewValuesInCell(filt.tau, filt.Tn(l, m), filt.P(l, m));

				for (int i = 1; i <= NUM_COMPS; i++)
					filt.Component[i].Nmid(l, m) = newDensities.N[i];
			}

		iWell.injectMoles(filt, O2, 1e-1, t, 2e9);
		iWell.injectEnthalpy(filt, 350, t, 2e9);

		rWell.recoverMoles(filt, 1e-1);
		rWell.recoverEnthalpy(filt);
		



		filt.CalculatePressureByGaussSeidel(filt.P, Pn, 1e0);
		//filt.P.print();
		//cin >> a;
		filt.CalculateFiltrationSpeeds();
		filt.CalculateMassFluxes();

		//filt.Component[HO].fluxXmid.print();
		//cout << endl;
		// DEBUG Output
		//filt.Component[HO].fluxXmid.print();

		filt.CalculateMassTransfer();
		

		// copying the data from n+1_level
		filt.P.copy(Pn1);

		// calculating new enthalpy
		filt.CalculateEnthalpyFluxes();
		filt.CalculateEnthalpyAtNextTimeLayer(Pn, filt.P);

		// calculating phase balance & saturations
		for (int l = 0; l <= filt.P.lastIndexX(); l++)
			for (int m = 0; m <= filt.P.lastIndexY(); m++)
			{
				// calculating NTotalInCell
				NtotalInCell = 0;
				for (int i = 1; i <= NUM_COMPS; i++)
				{
					NtotalInCell += filt.Component[i].Nn1(l, m);
				}

				// calculating concentrations
				for (int i = 1; i <= NUM_COMPS; i++)
					filt.CalculateConcentrationOfComponentInCell(l, m, i, filt.Component[i].Nn1(l, m), NtotalInCell);

				totalMolarEnthalpyInCell = filt.Hn1(l, m) / NtotalInCell;
				filt.CalculatePhaseBalanceInCell(
													l, m, filt.x, totalMolarEnthalpyInCell, filt.P(l, m), filt.Tn1(l, m)
												);
				filt.CalculateSaturationsAndDensitiesOfPhasesInCell(l, m);
			}


		counter++;	

		type filtSpeedTotal = 0;
		for (int a = 1; a <= NUM_PHASES; a++)
			filtSpeedTotal += filt.Phase[a].filtSpeedXmid(iWell.getXindex(), iWell.getYindex());
		
		type courantNumber = filtSpeedTotal * filt.tau / filt.hx;
		// printing data in file
		if ((int)t % (int)1e7 == 0)
		{
			//cout << "courantNumber = " << courantNumber << "filtSpeedTotal = " << filtSpeedTotal << endl;
			filt.VisualizeVTK(counter);
			//cin >> a;
		}

		// changing n and n+1 layers
		swap(filt.Tn, filt.Tn1);
		swap(filt.Hn, filt.Hn1);
		swap(filt.En, filt.En1);
		for (int i = 1; i <= NUM_COMPS; i++)
			swap(filt.Component[i].Nn, filt.Component[i].Nn1);

		if ((int)t % (int)1e6 == 0)
			cout << "Timestep № " << counter << ", tau = " << filt.tau << ", t = " << t << endl;

		t += filt.tau;
		//counter++;
	}

	//filt.Print();


	return 0;
}
