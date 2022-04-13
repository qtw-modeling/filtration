/*
 * array_functions.cpp
 *
 *  Created on: 05.11.2012
 *      Author: Alexey Karpaev
 */

#include "Filt.h"


/*type** Filt1D::MemAlloc2D(int Y_SIZE, int X_SIZE)
{
	type **array = new type*[Y_SIZE];
	for (int m = 0; m < Y_SIZE; m++)
		array[m] = new type[X_SIZE];
	return array;
}

type* Filt1D::MemAlloc(int size)
{
	type* pArray = new type[size];
	return pArray;
}

type*** Filt1D::MemAlloc3D(int Z_SIZE, int Y_SIZE, int X_SIZE)
{
	type ***array = new type**[Z_SIZE];
	for (int k = 0; k < Z_SIZE; k++)
	{
  		array[k] = new type*[Y_SIZE];
		for (int m = 0; m < Y_SIZE; m++)
			array[k][m] = new type[X_SIZE];
	}
	return array;
}


void Filt1D::Print(int numIter)
{
		FILE *fp1;
		if ( (fp1 = fopen(util::format("../Output/OutputIERightMod_e%d.csv", numIter).c_str(), "w")) == NULL )
		{
		   		printf("Cannot open the file.\n");
		   		exit(1);
		}
		fprintf(fp1, "x,p,sLiquid,sWater,sGas,cHO,cLO,cH2O,cO2,cCO2,N_HO,N_LO,N_H2O,N_O2,N_CO2,T,H,E,WFull\n");
		int c, b = (int) (L / 10);
		type x;
		type filtSpeedFull;
		for (int l = 0; l <= L; l++)
		{
			c = b*l; x = l*h;
			filtSpeedFull = 0;
			for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
				filtSpeedFull += Phase[alpha].filtSpeedMid[l];
			fprintf(fp1, "%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e\n", x, P[l], Phase[Liquid].S[l], Phase[Water].S[l], Phase[Gas].S[l],
																  Component[HO].C[l], Component[LO].C[l], Component[H2O].C[l],
																  Component[O2].C[l],
																  Component[CO2].C[l], 
																  Component[HO].Nn1[l], Component[LO].Nn1[l], Component[H2O].Nn1[l],
																  Component[O2].Nn1[l],
																  Component[CO2].Nn1[l],
																  Tn[l], Hn[l], En[l], filtSpeedFull);
		}
		fclose(fp1);
}

void Filt1D::PrintArray(type* array)
{
	int c, b = L / 10;
	for (int l = 0; l <= 100; l++)
	{
		c = b * l;
		printf("%e ", array[l]);
	}
	printf("\n");
}

void Filt1D::fPrintArray(type* array)
{
	FILE *fp1;
	if ((fp1 = fopen("OUT_temperature.dat", "w")) == NULL)
	{
	   		printf("Cannot open the file.\n");
	    	//exit(1);
	}

	int c, b = L/10;
	type x;
	for (int l = 0; l <= L; l++)
	{
		c = b*l;
		x = l*h;
		fprintf(fp1, "%e %e\n ", x, array[l]);
	}
	printf("\n");
}


void Filt1D::SummArray(type* init, type* delta)
{
	for (int l = 0; l <= L; l++)
		init[l] += delta[l];
}

void Filt1D::CopyArray(type* dest, type* source)
{
	for (int l = 0; l <= L; l++)
		dest[l] = source[l];
}


//DEBUG: PRINT VALUE;
void Filt1D::DEBUG(type value)
{
	printf("n = %d, Value = %e\n", counter, value);
}
*/
