#include "Filt.h"
#include "Well.hpp"


//constructor
Filt2D::Filt2D(int numBlocksX_, int numBlocksY_) : 

numBlocksX(numBlocksX_), numBlocksY(numBlocksY_), 
numPointsX(numBlocksX_ + 1), numPointsY(numBlocksY_ + 1),
hx(lengthOfAreaX / numBlocksX_),
hy(lengthOfAreaY / numBlocksY_),
P(numPointsX, numPointsY),
lambda(numPointsX, numPointsY), 
omega(numPointsX, numPointsY),
gamma(numPointsX, numPointsY),
K(numPointsX, numPointsY),
Tn(numPointsX, numPointsY), Tn1(numPointsX, numPointsY), 
Hn(numPointsX, numPointsY), Hn1(numPointsX, numPointsY), Hwell(numPointsX, numPointsY),
En(numPointsX, numPointsY), En1(numPointsX, numPointsY),
fluxEnthalpyXmid(numPointsX - 1, numPointsY - 1), fluxEnthalpyYmid(numPointsX - 1, numPointsY - 1),
fluxHeatConductivityXmid(numPointsX - 1, numPointsY - 1), fluxHeatConductivityYmid(numPointsX - 1, numPointsY - 1),
x(numPointsX, numPointsY)

{
		cycle(l, 0, numPointsX - 1, 
				m, 0, numPointsY - 1)
			{
				x(l, m)(HO, Liquid)  = 0.8;
				x(l, m)(LO, Liquid)  = 0.2;
				x(l, m)(H2O, Liquid) = 0.0;
				x(l, m)(CO2, Liquid) = 0.0;
				x(l, m)(O2,Liquid)   = 0.0;

				x(l, m)(HO, Water)   = 0.0;
				x(l, m)(LO, Water)   = 0.0;
				x(l, m)(H2O, Water)  = 1.0;
				x(l, m)(CO2, Water)  = 0.0;
				x(l, m)(O2, Water)   = 0.0;

				x(l, m)(HO, Gas) 	 = 0.0;
				x(l, m)(LO, Gas)     = 0.0;
				x(l, m)(H2O, Gas)    = 0.0;
				x(l, m)(CO2, Gas)    = 0.0;
				x(l, m)(O2, Gas)     = 1.0;
			}
		}

		const type tmp = 1e3 + 0e0;
		const type calTmp = 2e6;
		const type vMolarTmp = 1e-5;
		const type cpTmp = 4.2e3;	/// setting component data
		
		Component[HO]  = component(numPointsX, numPointsY, /* volumeMolar */ vMolarTmp, cpTmp, cpTmp, cpTmp, /* calorisity */ 0);
		Component[LO]  = component(numPointsX, numPointsY, vMolarTmp, cpTmp, cpTmp, cpTmp, 0);
		Component[H2O] = component(numPointsX, numPointsY, vMolarTmp, cpTmp, cpTmp, cpTmp, 0);
		Component[CO2] = component(numPointsX, numPointsY, vMolarTmp, cpTmp, cpTmp, cpTmp, 0);
		Component[O2]  = component(numPointsX, numPointsY, vMolarTmp, cpTmp, cpTmp, cpTmp, 2e6);


		// setting phase data
		Phase[Liquid] = phase(numPointsX, numPointsY);
		Phase[Water]  = phase(numPointsX, numPointsY);
		Phase[Gas] 	  = phase(numPointsX, numPointsY);

		Phase[Liquid].setViscosity(1e-3, 0, 0);
		Phase[Water].setViscosity(1e-3, 0, 0);
		Phase[Gas].setViscosity(1e-3, 0, 0);
}


// destructor
Filt2D::~Filt2D()
{
	std::cout << "~Filt2D" << std::endl;
}


void Filt2D::SetApproximation(enum Approximation approximation_)
{
	approximation = approximation_;
}


void Filt2D::SetInitialConditions(type Pinside)
{
	// setting initial  pressure
	type Ksmall = 1e-13;
	type Kbig = 1e-13;
	int x1, y1;
	for (int l = 0; l <= P.lastIndexX(); l++)
		for (int m = 0; m <= P.lastIndexY(); m++)
		{
			x1 = l * hx;
			y1 = m * hy;
			
			if (y1 <= lengthOfAreaY / 3. || y1 >= lengthOfAreaY * 2. / 3. )
			{
				K(l, m) = Ksmall;
			}
			else
			{
				K(l, m) = Kbig;
			}
			
			P(l, m) = Pinside;

			lambda(l, m) = 1.0;
			omega(l, m) = 0.0;
			gamma(l, m) = 0.0;

			for (int i = 1; i <= NUM_COMPS; i++) 
			{
				Component[i].C(l, m) = omega(l, m) * x(l, m)(i, Water) + lambda(l, m) * x(l, m)(i, Liquid) + gamma(l, m) * x(l, m)(i, Gas);
			}

			CalculateSaturationsAndDensitiesOfPhasesInCell(l, m);

			// calculating densities
			for (int i = 1; i <= NUM_COMPS; i++) {
				CalculateMolarDensityOfComponentInCell(l, m, i, Component[i].Nn);
			}

			Tn(l, m) = 300;
			Hn(l, m) = 0;
			for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
			{
				Hn(l, m) += f * Phase[alpha].rho(l, m) * Phase[alpha].S(l, m) * MolarEnthalpyOfPhaseInCell(alpha, Tn(l, m), x(l, m));
			}
			Hn(l, m) += (1 - f) * rhoS * cPS * Tn(l, m);
			CalculateEnergyUsingEnthalpyInCell(En(l, m), Hn(l, m), P(l, m));

			CalculateOmegasInCell(l, m);
		}

		// setting boundary permeability
		/*for (int l = 0; l <= P.lastIndexX(); l++)
		{
			K(l, 0) = K(l, P.lastIndexY()) = 0;
		}

		for (int m = 0; m <= P.lastIndexY(); m++)
		{
			K(0, m) = K(P.lastIndexX(), m) = 0;
		}

		K(0, 2) = Kleft;*/

	// BORDERS----------------------------------------------------------------------
}


void Filt2D::CalculateOmegasInCell(int l, int m)
{
	Component[HO].omega(l, m)  = Component[HO].vMolar;
	Component[LO].omega(l, m)  = Component[LO].vMolar;
	Component[H2O].omega(l, m) = Component[H2O].vMolar;
	Component[CO2].omega(l, m) = Component[CO2].vMolar;
	Component[O2].omega(l, m)  = Component[O2].vMolar;
}

void Filt2D::CalculateSaturationsAndDensitiesOfPhasesInCell(int l, int m)
{
	type vMolarL = 0, vMolarW = 0, vMolarG = 0;
	type vL, vW, vG, vTotal;
	for (int i = 1; i <= NUM_COMPS; i++)
	{
		vMolarL += x(l, m)(i, Liquid) * Component[i].vMolar;
		vMolarW += x(l, m)(i, Water) * Component[i].vMolar;
		vMolarG += x(l, m)(i, Gas) * Component[i].vMolar;
	}

	vL = lambda(l, m) * vMolarL;
	vW = omega(l, m) * vMolarW;
	vG = gamma(l, m) * vMolarG;

	vTotal = vL + vW + vG;

	assert(vMolarL >= 0);
	assert(vMolarW >= 0);
	assert(vMolarG >= 0);
	assert(vTotal > 0);

	// calculating densities
	Phase[Liquid].rho(l, m) = 1.0 / vMolarL;
	Phase[Water].rho(l, m)  = 1.0 / vMolarW;
	Phase[Gas].rho(l, m)    = 1.0 / vMolarG;
	//cout << vMolarG << endl;

	// calculating saturations
	Phase[Liquid].S(l, m) = vL / vTotal;
	Phase[Water].S(l, m)  = vW / vTotal;
	Phase[Gas].S(l, m)    = vG / vTotal;

	for (int a = 1; a <= NUM_PHASES; a++)
	{
		if (Phase[a].S(l, m) < 0)
			Phase[a].S(l, m) = eps2; // NOT FAIR!
	}

}

void Filt2D::CalculateMolarDensityOfComponentInCell(int l, int m, int i, Array2D <type>& N)
{
			N(l, m) = 0;
			for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
				N(l, m) += f * Phase[alpha].rho(l, m) * Phase[alpha].S(l, m) * x(l, m)(i, alpha);
}


/*type Filt2D::etaMid(int l, int i, int alpha, type* S, type* P)
{
	type fPhase, eta, xRhoMid;
	switch(approximation)
	{
		case upwind:
			xRhoMid = UpwindApproximationAtMidPoint(P[l], P[l + 1], x[l][i][alpha] * Phase[alpha].rho[l], x[l + 1][i][alpha] * Phase[alpha].rho[l + 1]);
			fPhase = UpwindApproximationAtMidPoint(P[l], P[l + 1], RelativePermeability(Phase[alpha].S[l]), RelativePermeability(Phase[alpha].S[l + 1]));
			eta =  xRhoMid * K * fPhase / Phase[alpha].mu;
			break;
		case arithmeticAverage:
			xRhoMid = ApproximationAtMidPoint(x[l][i][alpha] * Phase[alpha].rho[l], x[l + 1][i][alpha] * Phase[alpha].rho[l + 1], approximation);
			fPhase = ApproximationAtMidPoint(RelativePermeability(Phase[alpha].S[l]), RelativePermeability(Phase[alpha].S[l + 1]), approximation);
			eta =  xRhoMid * K * fPhase / Phase[alpha].mu;
			break;
	}

	return eta;
}*/


/*type Filt2D::Tmiddle(int k, type* T, type* P)
{
	if (P[k] > P[k+1])
		return T[k];
	else
		return T[k+1];
}*/

/*type Filt2D::CalculateNormOfArray(type* array, int length)
{
		type delta, max = 0;
		for (int k = 0; k <= length; k++)
		{
			delta = fabs(array[k]);
				max = delta;
		}
		return max;
}*/


void Filt2D::CalculateConcentrationOfComponentInCell(int l, int m, int i, type N, type NtotalInCell)
{
		assert(NtotalInCell > 0);
		// the main calculation
		Component[i].C(l, m) = N / NtotalInCell;

		assert(Component[i].C(l, m) >= 0);
}


void Filt2D::CalculatePressureByGaussSeidel(Array2D <type>& Pnew, Array2D <type>& Pprevious, type precision)
{


	// initializing arrays
	type B[6];

	Array2D<vector<type> > mobilityOfPhaseXmid(Pnew.getXsize() - 1, Pnew.getYsize() - 1);
	Array2D<vector<type> > mobilityOfPhaseYmid(Pnew.getXsize() - 1, Pnew.getYsize() - 1);

	Array2D<type> totalMobilityXmid(Pnew.getXsize() - 1, Pnew.getYsize() - 1);
	Array2D<type> totalMobilityYmid(Pnew.getXsize() - 1, Pnew.getYsize() - 1);

	for (int l = 0; l <= Pnew.lastIndexX() - 1; l++)
		for (int m = 0; m <= Pnew.lastIndexY() - 1; m++)
		{
			mobilityOfPhaseXmid(l, m) = vector<type>(NUM_PHASES + 1);
			mobilityOfPhaseYmid(l, m) = vector<type>(NUM_PHASES + 1);
		}


	// setting intial values for the iterations to start with
	type x, y;
	for (int l = 0; l <= Pnew.lastIndexX(); l++)
		for (int m = 0; m <= Pnew.lastIndexY(); m++)
			Pnew(l, m) = Pprevious(l, m);


	type max, rhs, Ptmp, PGaussSeidel, delta, omega = 1; // MAX norma raznosti; temprorary variable; norma raznosti
	uint counter = 0;

	for (int l = 0; l <= Pnew.lastIndexX() - 1; l++)
			for (int m = 1; m <= Pnew.lastIndexY() - 1; m++)
			{
				for (int alpha = 1; alpha <= NUM_PHASES; alpha++) 
				{
					mobilityOfPhaseXmid(l, m)[alpha] = UpwindApproximationAtMidPoint(P(l, m),
																					 P(l + 1, m),
																					 RelativePermeability(Phase[alpha].S(l, m)) / Phase[alpha].mu(Tn(l, m)),
																					 RelativePermeability(Phase[alpha].S(l + 1, m)) / Phase[alpha].mu(Tn(l + 1, m)));
				}

				totalMobilityXmid(l, m) = 0;
				for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
					totalMobilityXmid(l, m) += mobilityOfPhaseXmid(l, m)[alpha];
			}
		

		for (int l = 1; l <= Pnew.lastIndexX() - 1; l++)
			for (int m = 0; m <= Pnew.lastIndexY() - 1; m++)
			{
				for (int alpha = 1; alpha <= NUM_PHASES; alpha++) 
				{
					mobilityOfPhaseYmid(l, m)[alpha] = UpwindApproximationAtMidPoint(P(l, m),
																					 P(l, m + 1),
																					 RelativePermeability(Phase[alpha].S(l, m)) / Phase[alpha].mu(Tn(l, m)),
																					 RelativePermeability(Phase[alpha].S(l, m + 1)) / Phase[alpha].mu(Tn(l, m + 1)));
				}

				totalMobilityYmid(l, m) = 0;
				for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
					totalMobilityYmid(l, m) += mobilityOfPhaseYmid(l, m)[alpha];
			}


	do
	{
		max = 0;
		// the main cycle
		for (int l = 1; l <= Pnew.lastIndexX() - 1; l++)
			for (int m = 1; m <= Pnew.lastIndexY() - 1; m++)
			{
				// setting up the coefficients
				B[1] =  ApproximationAtMidPoint(K(l, m), K(l + 1, m), geometricAverage) * totalMobilityXmid(l, m) / (hx * hx);
				B[2] =  ApproximationAtMidPoint(K(l, m), K(l, m + 1), geometricAverage) * totalMobilityYmid(l, m) / (hy * hy);	//w1*A1*k1MidPlusY + w2*A2*k2MidPlusY +  w3*A3*k3MidPlusY )/(hy*hy);
				B[3] =  ApproximationAtMidPoint(K(l, m), K(l - 1, m), geometricAverage) * totalMobilityXmid(l - 1, m) / (hx * hx);	//( w1*A1*k1MidMinusX +  w2*A2*k2MidMinusX +  w3*A3*k3MidMinusX )/(hx*hx);
				B[4] =  ApproximationAtMidPoint(K(l, m), K(l, m - 1), geometricAverage) * totalMobilityYmid(l, m - 1) / (hy * hy);	//( w1*A1*k1MidMinusY +  w2*A2*k2MidMinusY +  w3*A3*k3MidMinusY )/(hy*hy);
				B[5] = -(B[1] + B[2] + B[3] + B[4]);
				rhs = 0;
				for (int i = 1; i <= NUM_COMPS; i++)
					rhs += Component[i].omega(l, m) * (Component[i].Nwell(l, m) - Component[i].Nn(l, m)) / tau;

				/*if ((m == 1 && l == 1) || (m == 1 && l == 2))
				{
					printf("l = %d, m = %d\n", l, m);
					for (int j = 1; j <= 5; j++)
						printf("B[%d] = %e\n", j, B[j]);
				}*/



				Ptmp = Pnew(l, m);
				
				// calculating new values
				PGaussSeidel = (-(B[1] * Pnew(l + 1, m) + B[2] * Pnew(l, m + 1) + B[3] * Pnew(l - 1, m) + B[4] * Pnew(l, m - 1)) - rhs) / B[5];

				Pnew(l, m) = Ptmp + omega * (PGaussSeidel - Ptmp);
				
				delta = fabs(Pnew(l, m) - Ptmp);
				if (delta > max)
					max = delta;
			}
		
			counter++;
			/*if (counter > (int)2e4)
			{
				cerr << "Pressure iterations haven't converged!" << endl;
				//break;
			}*/

			if ((counter % (int)1e3) == 0)
				cout << "Residual = " << max << endl;

	} while(max > precision);
	
	//printf("Iteration = %d , deltaMax = %e\n", counter, max);
}


/*type Filt2D::CalculateTotalMolarEnthalpyInCell(int l, type T, type** x)
{
	type hL = MolarEnthalpyOfPhaseInCell(Liquid, T, x);
	type hW = MolarEnthalpyOfPhaseInCell(Water, T, x);
	type hG = MolarEnthalpyOfPhaseInCell(Gas, T, x);

	type hTotal = lambda[l] * hL + omega[l] * hW + gamma[l] * hG;

	return hTotal;
}*/


inline type Filt2D::MolarEnthalpyOfComponentInPhaseInCell(int i, int alpha, type T)
{
	type value = Component[i].cP[alpha] * T + Component[i].calorisity;
	return value;
}


/*type Filt2D::CalculateCpOfPhase(int alpha, type** x)
{
	type value = 0;
	for (int i = 1; i <= NUM_COMPS; i++)
	{
		value += Component[i].cP[alpha] * x[i][alpha];
	}
	return value;
}*/

type Filt2D::MolarEnthalpyOfPhaseInCell(int alpha, type T, Composition& x)
{
	type hAlpha = 0;	// molar enthalpy of phase alpha
	type hia[NUM_COMPS + 1][NUM_PHASES + 1];
	for (int i = 1; i <= NUM_COMPS; i++)
		for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
			hia[i][alpha] = MolarEnthalpyOfComponentInPhaseInCell(i, alpha, T);

	for (int i = 1; i <= NUM_COMPS; i++)
		hAlpha += hia[i][alpha] * x(i, alpha);

	return hAlpha;
}


/*void Filt2D::CalculateTimeStepAuto()
{
	type wLMax = CalculateNormOfArray(Phase[Liquid].filtSpeedMid, L-1);
	type wWMax = CalculateNormOfArray(Phase[Water].filtSpeedMid, L-1);
	type wGMax = CalculateNormOfArray(Phase[Gas].filtSpeedMid, L-1);

	type wFullMax = wLMax + wWMax + wGMax;	// max speed of filtration of the mixture

	type uMax = wFullMax / f;	// max speed of the mixture
	if (uMax <= 1e-13)			// if uMax == 0, but computer arithmetics sets it almost to 0
		tau = 5e3;
	else
		tau = 0.1 * h / uMax;
}*/


void Filt2D::CalculatePhaseBalanceInCell(int l, int m, Array2D <Composition>& x, type eta, type P, type& Tnew)
{
	la::vector <double> c(NUM_COMPS);

	c[0] = Component[H2O].C(l, m);
	c[1] = Component[LO].C(l, m);
	c[2] = Component[HO].C(l, m);
	c[3] = Component[CO2].C(l, m);
	c[4] = Component[O2].C(l, m);

	PhaseSolver ps(1, 1, 2, 0, 1e-6);	// LO , HO , G , S , eps 
	ps.setEquilibriumConstantLog(0, new Kappa(11.8572e9, 3816.44, 46.13)); // H2O
	ps.setEquilibriumConstantLog(1, new Kappa(1.0026e9, 2477.07, 39.94)); // C5H12

	//type tmp = 1 + 1e-1;
	//ps.setEquilibriumConstantLog(0, new KappaConst(tmp)); // H2O 
	//ps.setEquilibriumConstantLog(1, new KappaConst(tmp));

	ps.setStandartEnthalpy(0, new GasEnthalpy(Component[H2O].cP[Gas], 0, Component[H2O].calorisity)); // H2O
	ps.setStandartEnthalpy(1, new GasEnthalpy(Component[LO].cP[Gas], 0, Component[LO].calorisity)); // C5H12
	ps.setStandartEnthalpy(2, new GasEnthalpy(Component[HO].cP[Gas], 0, Component[HO].calorisity));
	ps.setStandartEnthalpy(3, new GasEnthalpy(Component[CO2].cP[Gas], 0, Component[CO2].calorisity));
	ps.setStandartEnthalpy(4, new GasEnthalpy(Component[O2].cP[Gas], 0, Component[O2].calorisity));


	ps.setVaporizationEnthalpy(0, new VaporizationEnthalpy1(Component[H2O].cP[Gas], Component[H2O].cP[Water]));
	ps.setVaporizationEnthalpy(1, new VaporizationEnthalpy1(Component[LO].cP[Gas], Component[LO].cP[Liquid]));

	ps.setTempRange(100, 2000); // Range for which correlations are valid

	type Tapprox = 300;	// initial value for iteration
	ps.setParams(c, P, Tapprox, eta);	// eta = H / N (molar enthalpy)

	la::vector<double> solution = ps.solve();

	la::vector<double> cNew = ps.getAdjustedC();
	double res;

	Status s = ps.getStatus(res);

	if (s != Ok)
		std::cerr << s << " " << res << std::endl;

	Kappa tmp1(11.8572e9, 3816.44, 46.13);
	Kappa tmp2(1.0026e9, 2477.07, 39.94);

	// setting phase balance
	omega(l, m) = solution[0];
	type lambda1 = solution[1];
	Tnew = solution[2];

	lambda(l, m) = lambda1 + cNew[2];
	gamma(l, m) = (cNew[0] - omega(l, m)) + (cNew[1] - lambda1) + cNew[3] + cNew[4];

	//if (l == L/2)
	//	cout << "c = " << c[0] << " cNew = " << cNew[0] << endl;


	assert(gamma(l, m) > 0);
	assert(lambda(l, m) > 0);
	assert(omega(l, m) > 0);

	x(l, m)(HO, Liquid) = cNew[2] / lambda(l, m);
	x(l, m)(LO, Liquid) = lambda1 / lambda(l, m);
	x(l, m)(H2O, Liquid) = 0;
	x(l, m)(O2, Liquid) = 0;
	x(l, m)(CO2, Liquid) = 0;

	x(l, m)(HO, Water) = 0;
	x(l, m)(LO, Water) = 0;
	x(l, m)(H2O, Water) = 1;
	x(l, m)(O2, Water) = 0;
	x(l, m)(CO2, Water) = 0;

	x(l, m)(HO, Gas) = 0;
	x(l, m)(LO, Gas) = (cNew[1] - lambda1) / gamma(l, m);
	x(l, m)(H2O, Gas) = (cNew[0] - omega(l, m)) / gamma(l, m);
	x(l, m)(O2, Gas) = cNew[4] / gamma(l, m);
	x(l, m)(CO2, Gas) = cNew[3] / gamma(l, m);

	//cout << "eta was = " << eta << ", new eta = " << etaNew << endl;

	for (int i = 1; i <= NUM_COMPS; i++)
		for (int alpha = 1; alpha <= NUM_PHASES; alpha++)
		{
			if (x(l, m)(i, alpha) < 0 || x(l, m)(i, alpha) > 1)
				x(l, m)(i, alpha) = eps2;
		}
}


//CALCULATE dP
/*
void Filt_1D::Progonka()
{
	type a[L+1], b[L+1], c[L+1], d[L+1], alpha[L], beta[L];
	type etaMidPlus[NUM_PHASES+1], etaMidMinus[NUM_PHASES+1];

	for(int l = 1; l <= L-1; l++)
	{
		a[l] = b[l] = c[l] = d[l] = 0;

		for(int i = 1; i <= NUM_COMPS; i++)
		{
			//CALCULATE THE MASSIVES OF hiMiddles;
			for(int alpha1 = 1; alpha1 <= NUM_PHASES; alpha1++)
			{
				etaMidPlus[alpha1] = etaMiddle(l, i, alpha1,  Phase[alpha1].S, P);
				etaMidMinus[alpha1] = etaMiddle(l-1, i, alpha1, Phase[alpha1].S, P);
			}


			a[l] += -tau/(h*h) * Comp[i].omega * (etaMidPlus[1] + etaMidPlus[2] + etaMidPlus[3]);

			b[l] += -tau/(h*h) * Comp[i].omega * (
										   -(etaMidMinus[1] + etaMidMinus[2] + etaMidMinus[3])
										   - (etaMidPlus[1] + etaMidPlus[2] + etaMidPlus[3])
										   );

			c[l] += -tau/(h*h) * Comp[i].omega * (etaMidMinus[1] + etaMidMinus[2] + etaMidMinus[3]);

			d[l] +=	Comp[i].omega * (Comp[i].Nn[l] - Comp[i].N_s[l])
					+ tau/(h*h) *  Comp[i].omega * (
											(etaMidPlus[1] + etaMidPlus[2] + etaMidPlus[3]) * (P[l+1] - P[l])
											- (etaMidMinus[1] + etaMidMinus[2] + etaMidMinus[3]) * (P[l] - P[l-1])
										  );
		}пушной

	}
		a[0] = 0; b[0] = 1; c[0] = 0; d[0] = 0;
		a[L] = 0; b[L] = 1; c[L] = 0; d[L] = 0;

		alpha[0] = -a[0] / b[0];
		beta[0] = d[0] / b[0];

		for(int l = 1; l <= L-1; l++)
		{
			alpha[l] = -a[l] / (b[l] + c[l] * alpha[l-1]);
			beta[l] = (d[l] - c[l] * beta[l-1]) / (b[l] + c[l] * alpha[l-1]);
		}

		dP[L] = (d[L] - c[L] * beta[L-1]) / (b[L] + c[L] * alpha[L-1]);

		for(int l = L - 1; l >= 0; l--)
		{
			dP[l] = alpha[l] * dP[l+1] + beta[l];
		 << "dp = " << dP[l] << " ";
		}
}
*/

//one iteration's actions
/*
int Filt_1D::Iteration()
{

	for (int i = 1; i <= NUM_COMPS; i++)
		CopyArray(Comp[i].N_s, Comp[i].Nn); //Copy_array(N2_s, N2);  // copy the data from the prev level

	int count = 0;
	type delta;

	do
	{
		Progonka();
		delta = CalculateNormOfArray(dP);
		//DEBUG(delta);
		SummArray(P, dP);
		CalculateDeltaN();
		count++;
	} while(delta >= eps);
			//PrintArray(P);

	for(int i = 1; i <= NUM_COMPS; i++)
		CopyArray(Comp[i].Nn, Comp[i].N_s);

	return 0;
}
*/

/*type Filt_1D::A(int i){ //calc A(i)
			if( i == 1 )
				return po1*K/mu1;
			else
				return po2*K/mu2;
		}

void Filt_1D::Summ_array(type* init, type* delta){
		for(int k = 0; k <= L; k++)
			init[k] += delta[k];
}

void Filt_1D::Copy_array(type* dest, type* source){
		for(int k = 0; k <= L; k++)
			dest[k] = source[k];
}


void Filt_1D::Print(){ // print array i need
		FILE *fp1;
		if((fp1 = fopen("OUT_temperature.dat", "w"))== NULL){
		   		printf("Cannot open the file.\n");
		    	//exit(1);
		}
		type* f_skach = new type[L+1];
		int c, b = L/10; type x;
		for(int k = 0; k <= L; k++){
			c = b*k; x = k*h;
			f_skach[k] = S2[k]*S2[k]/mu2/(S1[k]*S1[k]/mu1 + S2[k]*S2[k]/mu2);
			//printf(" %e  %e  %e  %e %e %e %e %e\n ", x, array1[k], array2[k], array3[k], array4[k], array5[k], array6[k], array7[k], array8[k]);
			fprintf(fp1, "%e %e %e %e %e %e %e %e %e %e\n", x, P1[k], C1[k], C2[k], S1[k], S2[k], N1[k], N2[k], Tn[k], f_skach[k]);
		}
		//printf("\n\n");
		//fprintf(fp, "\n\n");
		//
}

/*void Filt_1D::Print_array(const type* array){
	int c, b = L/10;
	for(int k = 0; k <= 100; k++){
		c = b*k;
		printf("%e ", array[k]);
	}
	printf("\n");
}

void Filt_1D::fPrint_array(const type* array){
	FILE *fp1;
	if((fp1 = fopen("OUT_temperature.dat", "w"))== NULL){
	   		printf("Cannot open the file.\n");
	    	//exit(1);
	}
	int c, b = L/10;
	type x;
	for(int k = 0; k <= L; k++){
		c = b*k;
		x = k*h;
		fprintf(fp1, " %e %e\n ", x, array[k]);
	}
	printf("\n");
}





type Filt_1D::Calc_filt_speed(){ // calculates the speed of filtration
	type speed, tmp;
	speed = (-K*(S1[1]*S1[1]/mu1 + S2[1]*S2[1]/mu2))*(P1[2] - P1[1])/h;
	for(int k = 2; k <= L-1; k++){
		tmp = (-K*(S1[k]*S1[k]/mu1 + S2[k]*S2[k]/mu2))*(P1[k+1] - P1[k])/h;
		/*if(tmp != speed){
			printf(" Filtration's speed is not a constant! ");

		}
		speed = tmp;
			//continue;
		//printf("%f ", speed);
	}
	//printf(" %e ", speed);
	//printf("\n");

	return speed;
}
*/
