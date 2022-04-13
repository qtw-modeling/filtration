#include "Filt.h"
//#define heatConductivity

static const type W = 5e-7;
static const type coeffHeatConductivity = 7e-7;

void Filt2D::CalculateEnthalpyAtNextTimeLayer(Array2D <type>& Pcurrent, Array2D <type>& Pnext)
{
	// the main calculation
	for (int l = 1; l <= P.lastIndexX() - 1; l++)
		for (int m = 1; m <= P.lastIndexY() - 1; m++)
		{
			

			Hn1(l, m) = Hwell(l, m) - tau  * ( 
											(fluxEnthalpyXmid(l, m) - fluxEnthalpyXmid(l - 1, m)) / hx +
											(fluxEnthalpyYmid(l, m) - fluxEnthalpyYmid(l, m - 1)) / hy
										  )
										+ (Pnext(l, m) - Pcurrent(l, m));
			if (Hn1(l, m) < 0)
			{
				Hn1(l, m) = 0;
			}
		}
	
	// setting boundary conditions
	for (int l = 1; l <= P.lastIndexX() - 1; l++)
	{
		Hn1(l, 0) = Hn1(l, 1);
		Hn1(l, P.lastIndexY()) = Hn1(l, P.lastIndexY() - 1);
	}

	for (int m = 0; m <= P.lastIndexY(); m++)
	{
		Hn1(0, m) = Hn1(1, m);
		Hn1(P.lastIndexX(), m) = Hn1(P.lastIndexX() - 1, m);
	}

}


void Filt2D::CalculateEnthalpyFluxes()
{
	type Tmid, hAlphaMid, rhohAlphaXmid, rhohAlphaYmid;

	// calculating fluxes at all other cells
	for (int l = 0; l <= P.lastIndexX() - 1; l++)
	{
		for (int m = 1; m <= P.lastIndexY() - 1; m++)
		{
			fluxEnthalpyXmid(l, m) = 0;

			for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
			{
				rhohAlphaXmid = UpwindApproximationAtMidPoint(
																P(l, m), 
																P(l + 1, m), 
																Phase[alpha].rho(l, m) 	* MolarEnthalpyOfPhaseInCell(alpha, Tn(l, m), x(l, m)), 
																Phase[alpha].rho(l + 1, m) * MolarEnthalpyOfPhaseInCell(alpha, Tn(l + 1, m), x(l + 1, m))
															);
					
				fluxEnthalpyXmid(l, m)  += rhohAlphaXmid * Phase[alpha].filtSpeedXmid(l, m);
			}
		}	

			
	}

	for (int l = 1; l <= P.lastIndexX() - 1; l++)
	{
		for (int m = 0; m <= P.lastIndexY() - 1; m++)
		{
			fluxEnthalpyYmid(l, m) = 0;

			for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
			{
				rhohAlphaYmid = UpwindApproximationAtMidPoint(
																P(l, m), 
																P(l, m + 1), 
																Phase[alpha].rho(l, m) * MolarEnthalpyOfPhaseInCell(alpha, Tn(l, m), x(l, m)), 
																Phase[alpha].rho(l, m + 1) * MolarEnthalpyOfPhaseInCell(alpha, Tn(l, m + 1), x(l, m + 1))
															);
				fluxEnthalpyYmid(l, m) += rhohAlphaYmid * Phase[alpha].filtSpeedYmid(l, m);
			}	
		}
	}

}


/*type Filt1D::AverageMolarHeatCapacityInCell(int l)
{
	type result = 0, Ntotal = 0;
	for (int i = 1; i <= NUM_COMPS; i++)
	{
		result += Component[i].cP[Liquid] * Component[i].Nn1[l];
		Ntotal += Component[i].Nn1[l];
	}

	result = result / Ntotal;
	return result;
}
*/


void Filt2D::CalculateEnergyUsingEnthalpyInCell(type& E, const type& H, const type& P)
{
	E = H - P;
}

/*
type Filt1D::TemperatureUsingFullEnthalpy(int l, type** sNext, type H)
{
	type sum1 = 0, sum2 = 0;
	for (int i = 1; i <= NUM_COMPS; i++)
		for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
			{
				sum1 += Phase[alpha].rho[l] * sNext[alpha][l] * x[l][i][alpha] * Component[i].calorisity;
				sum2 += Phase[alpha].rho[l] * sNext[alpha][l] * x[l][i][alpha] * Component[i].cP[alpha];
			}
			sum1 *= f; sum2 *= f;

			type T = (H - sum1) / sum2;
			return T;
}
*/