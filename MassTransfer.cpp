#include "Filt.h"

/*void Filt1D::CalculateENOFluxesAtMidPoints()
{

	for (int i = 1; i <= NUM_COMPS; i++)
		for (int l = 0; l <= L; l++)
		{
			Component[i].flux[l] = CalculateComponentFluxAtGridPoints(i, l);
		}
	//--------------------------------------------------------------------------------------------------
	type fluxMidTmp;
	for (int i = 1; i <= NUM_COMPS; i++)
	{
		// calculate flux at point (1 + 1/2)
		Component[i].fluxMid[1] = Component[i].flux[1] + (Component[i].flux[1] - Component[i].flux[0]) / 2 + (Component[i].flux[2] - 2*Component[i].flux[1] + Component[i].flux[0]) * 3.0 / 8;

		// calculate flux at point (l + 1/2) by ENO-procedure
		for (int l = 2; l <= L-1; l++)
		{
			if ( fabs(Component[i].flux[l] - 2 * Component[i].flux[l-1] + Component[i].flux[l-2])
				 < fabs(Component[i].flux[l+1] - 2 * Component[i].flux[l] + Component[i].flux[l-1]) )
				fluxMidTmp  = Component[i].flux[l] + (Component[i].flux[l] - Component[i].flux[l-1]) / 2 + (Component[i].flux[l] - 2*Component[i].flux[l-1] + Component[i].flux[l-2]) * 3.0 / 8;
			else
				fluxMidTmp  = Component[i].flux[l] + (Component[i].flux[l] - Component[i].flux[l-1]) / 2 + (Component[i].flux[l+1] - 2*Component[i].flux[l] + Component[i].flux[l-1]) * 3.0 / 8;

			Component[i].fluxMid[l] = fluxMidTmp;
		}
	}
}


type Filt1D::CalculateComponentFluxAtGridPoints(int i, int l)
{
	// launch this function only after the calculating the filtSpeeds!!!
	type flux = 0;
	for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
		flux += Phase[alpha].rho[l] * x[l][i][alpha] * Phase[alpha].filtSpeed[l];
	return flux;
}


void Filt1D::CalculateRhoXAtMidPointsByENO()
{
	for (int i = 1; i <= NUM_COMPS; i++)
		for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
		{
			// calculate rhoXPhi at (1+1/2) point
			A[i][alpha][1] = rhoXPhi(alpha, i, 1) + (rhoXPhi(alpha, i, 1)  - rhoXPhi(alpha, i, 0) ) / 2
							     	 	 	 	 	 + (rhoXPhi(alpha, i, 2) - 2*rhoXPhi(alpha, i, 1)  + rhoXPhi(alpha, i, 0)) * 3.0 / 8;

			for (int l = 2; l <= L-1; l++)
			{
				// comparing 2nd order finite differences
				if ( fabs(rhoXPhi(alpha, i, l)  - 2*rhoXPhi(alpha, i, l-1)  + rhoXPhi(alpha, i, l-2))
					 < fabs(rhoXPhi(alpha, i, l+1)  - 2*rhoXPhi(alpha, i, l)  + rhoXPhi(alpha, i, l-1)) )

					A[i][alpha][l] = rhoXPhi(alpha, i, l) + (rhoXPhi(alpha, i, l)  - rhoXPhi(alpha, i, l-1) ) / 2
									     + (rhoXPhi(alpha, i, l) - 2*rhoXPhi(alpha, i, l-1)  + rhoXPhi(alpha, i, l-2)) * 3.0 / 8;
				else
					A[i][alpha][l] = rhoXPhi(alpha, i, l) + (rhoXPhi(alpha, i, l)  - rhoXPhi(alpha, i, l-1) ) / 2
				     	 	 	 	 	 + (rhoXPhi(alpha, i, l+1) - 2*rhoXPhi(alpha, i, l)  + rhoXPhi(alpha, i, l-1)) * 3.0 / 8;
			}
		}
}


type Filt1D::rhoXPhi(int alpha, int i, int l)
{
	type value = Phase[alpha].rho[l] * x[l][i][alpha] * Phi(alpha, l);
	return value;
}


type Filt1D::Phi(int alpha, int l)
{
	type totalMobility = 0;
	for (int a = 1; a <= NUM_PHASES; a++)
		totalMobility += RelativePermeability(Phase[a].S[l]) / Phase[a].mu;

	type result = (RelativePermeability(Phase[alpha].S[l]) / Phase[alpha].mu) / totalMobility;
	return result;
}

type Filt1D::FiltSpeedAtMidPoints(int l)
{
	type totalMobilityL = 0, totalMobilityR = 0;
	type value = 0;
	for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
	{
		totalMobilityL += RelativePermeability(Phase[alpha].S[l]) / Phase[alpha].mu;
		totalMobilityR += RelativePermeability(Phase[alpha].S[l+1]) / Phase[alpha].mu;
	}

	type B = ApproximationAtMidPoints(totalMobilityL, totalMobilityR, geometricAverage);
	value = - B * (P[l+1] - P[l]) / h;
	return value;
}
*/


type Filt2D::RelativePermeability(type s) const
{
	type result = s * s;
	return result;
}

type Filt2D::UpwindApproximationAtMidPoint(type pL, type pR,  type fL, type fR) const
{
	type result;
	if (pL > pR)
		result = fL;
	else
		result = fR;
	return result;
}


void Filt2D::CalculateFiltrationSpeeds()
{
	type fPhaseXMid, fPhaseYMid;
	for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
	{
		for (int l = 0; l <= P.lastIndexX() - 1; l++)
			for (int m = 1; m <= P.lastIndexY() - 1; m++)
			{
				fPhaseXMid = UpwindApproximationAtMidPoint(
															P(l, m), 
															P(l + 1, m), 
															RelativePermeability(Phase[alpha].S(l, m)), 
															RelativePermeability(Phase[alpha].S(l + 1, m))
														  );
				// filtration speed between l and l+1 gridpoints with 2nd order approximation
				type KXmid = ApproximationAtMidPoint(K(l, m), K(l + 1, m), geometricAverage);
				type muXmid = UpwindApproximationAtMidPoint(P(l, m), 
															P(l + 1, m), 
															Phase[alpha].mu(Tn(l, m)), 
															Phase[alpha].mu(Tn(l + 1, m)));
				
				Phase[alpha].filtSpeedXmid(l, m) = -KXmid * fPhaseXMid / muXmid * (P(l + 1, m) - P(l, m)) / (hx);
			}
		
		for (int l = 1; l <= P.lastIndexX() - 1; l++)
			for (int m = 0; m <= P.lastIndexY() - 1; m++)
			{
				fPhaseYMid = UpwindApproximationAtMidPoint(
															P(l, m), 
															P(l, m + 1), 
															RelativePermeability(Phase[alpha].S(l, m)),
															RelativePermeability(Phase[alpha].S(l, m + 1))
														  );
				type KYmid = ApproximationAtMidPoint(K(l, m), K(l, m + 1), geometricAverage);
				type muYmid = UpwindApproximationAtMidPoint(P(l, m), 
															P(l, m + 1), 
															Phase[alpha].mu(Tn(l, m)), 
															Phase[alpha].mu(Tn(l, m + 1)));
				
				Phase[alpha].filtSpeedYmid(l, m) = -KYmid * fPhaseYMid / muYmid * (P(l, m + 1) - P(l, m)) / (hy);
			}
	}

}


void Filt2D::CalculateMassFluxes()
{
	type rhoXalphaXmid, rhoXalphaYmid;
	for (int i = 1; i <= NUM_COMPS; i++)
	{
		for (int l = 0; l <= P.lastIndexX() - 1; l++)
			for (int m = 1; m <= P.lastIndexY() - 1; m++)
			{
				Component[i].fluxXmid(l, m) = 0;

				for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
				{
					rhoXalphaXmid = UpwindApproximationAtMidPoint(
																	P(l, m), 
																	P(l + 1, m), 
																	Phase[alpha].rho(l, m) * x(l, m)(i, alpha), 
																	Phase[alpha].rho(l + 1, m) * x(l + 1, m)(i, alpha)
																);
					Component[i].fluxXmid(l, m) += rhoXalphaXmid * Phase[alpha].filtSpeedXmid(l, m);
				}		
			}
		
		for (int l = 1; l <= P.lastIndexX() - 1; l++)
			for (int m = 0; m <= P.lastIndexY() - 1; m++)
			{
				Component[i].fluxYmid(l, m) = 0;

				for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
				{
					rhoXalphaYmid = UpwindApproximationAtMidPoint(
																	P(l, m), 
																	P(l, m + 1), 
																	Phase[alpha].rho(l, m) * x(l, m)(i, alpha), 
																	Phase[alpha].rho(l, m + 1) * x(l, m + 1)(i, alpha)
																);
					Component[i].fluxYmid(l, m) += rhoXalphaYmid * Phase[alpha].filtSpeedYmid(l, m);
				}	
			}
	}
}


void Filt2D::CalculateMassTransfer()
{
	type tmp;
	// calculating mass transfer
	for (int i = 1; i <= NUM_COMPS; i++)
	{
		for (int l = 1; l <= P.lastIndexX() - 1; l++)
		{
			for (int m = 1; m <= P.lastIndexY() - 1; m++)
			{
					// the main calculation -------------------------------------------------------
					Component[i].Nn1(l, m) = Component[i].Nwell(l, m) - tau * ( 
																(Component[i].fluxXmid(l, m) - Component[i].fluxXmid(l - 1, m)) / hx + 
									 							(Component[i].fluxYmid(l, m) - Component[i].fluxYmid(l, m - 1)) / hy
									 	  									);
					// ----------------------------------------------------------------------------
					//cout << additive << " ";
					if(Component[i].Nn1(l, m) < 0)
					{
						//cerr << "N[" << i << "](" << l << ", " << m << ") = " << Component[i].Nn1(l, m) << endl;
						Component[i].Nn1(l, m) = 0;
					}
			}
			//cout << endl;
		}

		// setting values at the borders
		for (int l = 1; l <= P.lastIndexX() - 1; l++)
		{
			Component[i].Nn1(l, 0) = Component[i].Nn1(l, 1);
			Component[i].Nn1(l, P.lastIndexY()) = Component[i].Nn1(l, P.lastIndexY() - 1);
		}

		for (int m = 0; m <= P.lastIndexY(); m++)
		{
			Component[i].Nn1(P.lastIndexX(), m) = Component[i].Nn1(P.lastIndexX() - 1, m);
			Component[i].Nn1(0, m) = Component[i].Nn1(1, m);
		}

	}
}


type Filt2D::ApproximationAtMidPoint(type valueLeft, type valueRight, enum Approximation approx)
{
	type value;

	switch(approx)
	{
		case arithmeticAverage:
			value = 0.5 * (valueLeft + valueRight);
			break;
		case leftChoice:
			value = valueLeft;
			break;
		case rightChoice:
			value = valueRight;
			break;
		case geometricAverage:
			value = 2 * valueLeft * valueRight / (valueLeft + valueRight + 1e-13);
			break;
		default: // return algebraicAverage
			value = 0.5 * (valueLeft + valueRight);
			break;
	}

	return value;
}

