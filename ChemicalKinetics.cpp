#include "ChemicalKinetics.h"


type Max2(const type& a, const type& b)
{
	if (a > b)
		return a;
	else
		return b;
}


type Max3(const type& a, const type& b, const type& c)
{
	return Max2(a, Max2(b, c));
}

void AbstractOneStepChemicalSolver::SetInitialConditions(Densities& initialDensities)
{
	oldDensities.N = initialDensities.N;
}


Densities ExplicitImplicitChemicalSolverMod::CalculateNewValuesInCell(type tau, type T, type P)
{
	ChemicalReactions iReactions;

	// initializing first integrals
	vector<type> c(4);
	c[0] = 0; // not needed at all
	c[1] = 2 * oldDensities.N[O2] + 2 * oldDensities.N[CO2] + oldDensities.N[H2O];
	c[2] = 15 * oldDensities.N[HO] + 10 * oldDensities.N[LO] + oldDensities.N[CO2];
	c[3] = 32 * oldDensities.N[HO] + 22 * oldDensities.N[LO] + 2 * oldDensities.N[H2O];

	// computing new values
	newDensities.N[HO] = oldDensities.N[HO] / (1 + tau * iReactions.HO.SpeedOfReaction(T) * oldDensities.N[O2]);

	newDensities.N[LO] = oldDensities.N[LO] / (1 + tau * iReactions.LO.SpeedOfReaction(T) * oldDensities.N[O2]);

	newDensities.N[CO2] = c[2] - 15 * newDensities.N[HO] - 10 * newDensities.N[LO];

	newDensities.N[H2O] = (c[3] - 32 * newDensities.N[HO] - 22 * newDensities.N[LO]) / 2.;

	newDensities.N[O2] = oldDensities.N[O2] / (
													1 + tau * (23 * iReactions.HO.SpeedOfReaction(T) * newDensities.N[HO] +
													15.5 * iReactions.LO.SpeedOfReaction(T) * newDensities.N[LO])
												);
	return newDensities;
}


Densities ExplicitImplicitChemicalSolverMod2::CalculateNewValuesInCell(type tau, type T, type P)
{
	ChemicalReactions iReactions;

	// initializing first integrals
	vector<type> c(4);
	c[0] = 0; // not needed at all
	c[1] = 2 * oldDensities.N[O2] + 2 * oldDensities.N[CO2] + oldDensities.N[H2O];
	c[2] = 15 * oldDensities.N[HO] + 10 * oldDensities.N[LO] + oldDensities.N[CO2];
	c[3] = 32 * oldDensities.N[HO] + 22 * oldDensities.N[LO] + 2 * oldDensities.N[H2O];

	// computing new values
	newDensities.N[HO] = oldDensities.N[HO] / (1 + tau * iReactions.HO.SpeedOfReaction(T) * oldDensities.N[O2]);

	newDensities.N[LO] = oldDensities.N[LO] / (1 + tau * iReactions.LO.SpeedOfReaction(T) * oldDensities.N[O2]);

	newDensities.N[CO2] = c[2] - 15 * newDensities.N[HO] - 10 * newDensities.N[LO];

	newDensities.N[H2O] = (c[3] - 32 * newDensities.N[HO] - 22 * newDensities.N[LO]) / 2.;

	newDensities.N[O2] = (c[1] - 2 * newDensities.N[CO2] - newDensities.N[H2O]) / 2.;

	// NO2 can possibly be < 0; checking this case
	if (newDensities.N[O2] < 0)
	{
		newDensities.N[O2] = 0;

		type s, s1, s2;		// the speeds of reactions
		s = -oldDensities.N[O2] / tau;

		type a1, a2;

		a1 = iReactions.HO.SpeedOfReaction(T) * oldDensities.N[HO];
		a2 = iReactions.LO.SpeedOfReaction(T) * oldDensities.N[LO];

		/* s1 = gamma a1, s2 = gamma a2 */
		type gamma = s / (23 * a1 + 15.5 * a2);

		s1 = gamma * a1;
		s2 = gamma * a2;

		newDensities.N[HO] = oldDensities.N[HO] + tau * s1;
		
		newDensities.N[LO] = oldDensities.N[LO] + tau * s2;

		newDensities.N[CO2] = c[2] - 15 * newDensities.N[HO] - 10 * newDensities.N[LO];

		newDensities.N[H2O] = (c[3] - 32 * newDensities.N[HO] - 22 * newDensities.N[LO]) / 2.;
	}

	return newDensities;

}


Densities ExplicitImplicitChemicalSolver::CalculateNewValuesInCell(type tau, type T, type P)
{
	ChemicalReactions iReactions;

	newDensities.N[HO] = oldDensities.N[HO] / (1 + tau * iReactions.HO.SpeedOfReaction(T) * oldDensities.N[O2]);

	newDensities.N[LO] = oldDensities.N[LO] / (1 + tau * iReactions.LO.SpeedOfReaction(T) * oldDensities.N[O2]);

	newDensities.N[O2] = oldDensities.N[O2] / (
													1 + tau * (23 * iReactions.HO.SpeedOfReaction(T) * newDensities.N[HO] +
													15.5 * iReactions.LO.SpeedOfReaction(T) * newDensities.N[LO])
												);

	newDensities.N[CO2] = oldDensities.N[CO2] + tau * (15 * iReactions.HO.SpeedOfReaction(T) * newDensities.N[HO] * newDensities.N[O2] +
																		10 * iReactions.LO.SpeedOfReaction(T) * newDensities.N[LO] * newDensities.N[O2]);

	newDensities.N[H2O]  = oldDensities.N[H2O]  + tau * (16 * iReactions.HO.SpeedOfReaction(T) * newDensities.N[HO] * newDensities.N[O2] +
																		11 * iReactions.LO.SpeedOfReaction(T) * newDensities.N[LO] * newDensities.N[O2]);

	return newDensities;
}


Densities PureImplicitChemicalSolver::CalculateNewValuesInCell(type tau, type T, type P)
{
  //cout << "Launching chemical kinetics..." << endl;

	ChemicalReactions iReactions;

	type A = iReactions.HO.SpeedOfReaction(T);
	type B = iReactions.LO.SpeedOfReaction(T);
	//cout << " A = " << A << " B = " << B << endl;

	type nHONew, nMONew, nO2New, nCO2New, nH2ONew;

  // then solve cubic equation
  if (A != 0 && B != 0)
  {
      type a = (tau * (A + B) - tau * tau * A * B * oldDensities.N[O2]) / (tau * tau * A * B);
      type b = (tau * (ksi1 * A * oldDensities.N[HO] + ksi2 * B * oldDensities.N[LO]) - tau * A * oldDensities.N[O2] - tau * B * oldDensities.N[O2] + 1) / (tau * tau * A * B);
      type c =  - oldDensities.N[O2] / (tau * tau * A * B);

      // solve cubic equation with coeffs {a,b,c}
      roots root = cubic_roots(a, b, c);

      //cout << "rootA = " << root.x1 << " rootB = " << root.x2 << " rootC = " << root.x3 << endl;

      // block that chooses the right root
      /*if (root.x1 < 0 && root.x2 < 0 && root.x3 >= 0)
      {
      	cout << "nO2New = " << root.x3 << "\n";
         nO2New = root.x3;
      }
      else if (root.x1 >= 0 && root.x2 < 0 && root.x3 < 0)
      {
         cout << "nO2New = " << root.x1 << "\n";
         nO2New = root.x1;
      }
      else if (root.x1 < 0 && root.x2 >= 0 && root.x3 < 0)
      {
      	cout << "nO2New = " << root.x2 << '\n';
      	nO2New = root.x2;
      }
      else
      {
      	cerr << "Two positive  roots!" << '\n';
      	nO2New = Max3(root.x1, root.x2, root.x3);
      	//exit(1);
      }*/
      if (root.x1 >= 0)
      	nO2New = root.x1;
      else if (root.x2 >= 0)
      	nO2New = root.x2;
      else
		nO2New = root.x3;
		//cout << root.x1 << " " << root.x2 << " " << root.x3 << '\n';
      //assert(nO2New >= 0);


  }
  // then solve square equation
  else if ((A + B) != 0)
  {
    type a = tau * (A + B);
    type b = tau * (ksi1 * A * oldDensities.N[HO] + ksi2 * B * oldDensities.N[LO]) - tau * A * oldDensities.N[O2] - tau * B * oldDensities.N[O2] + 1;
	type c =  -oldDensities.N[O2];

   	type Discriminant = b * b - 4 * a * c;

   	assert(a != 0);
    type x1 = (-b + sqrt(Discriminant)) / (2*a);
    nO2New = x1;
    assert(nO2New >= 0);
  }
  // then reaction doesn't start
  else
      nO2New = oldDensities.N[O2];

  // if timestep is too large
  if (nO2New < 0)
	  nO2New = 0;

  // set new values
	assert ((1 + tau * A * nO2New) != 0);
  nHONew = oldDensities.N[HO] / (1 + tau * A * nO2New);
	assert ((1 + tau * B * nO2New) != 0);
  nMONew = oldDensities.N[LO] / (1 + tau * B * nO2New);
  nCO2New = oldDensities.N[CO2] + tau * (heta1 * A * nHONew * nO2New + heta2 * B * nMONew * nO2New);
  nH2ONew = oldDensities.N[H2O] + tau * (nu1 * A * nHONew * nO2New + nu2 * B * nMONew * nO2New);

  const type epsChemical = 0;

  if(nHONew < 0)
    nHONew = epsChemical;
  if(nMONew < 0)
    nMONew = epsChemical;
  if(nH2ONew < 0)
    nH2ONew = epsChemical;
  if(nO2New < 0)
    nO2New = epsChemical;
  if(nCO2New < 0)
    nCO2New = epsChemical;

  Densities returnValues;

  returnValues.N[HO] = nHONew;
  returnValues.N[LO] = nMONew;
  returnValues.N[O2] = nO2New;
  returnValues.N[H2O] = nH2ONew;
  returnValues.N[CO2] = nCO2New;

  //cout << "Chemical kinetics was calculated." << endl;

  return returnValues;
}


roots cubic_roots(type a, type b, type c)
{
                type Q = (a * a - 3 * b) / 9;
                type R = (2 * a * a * a - 9 * a * b + 27 * c) / 54;
                type Q3 = Q * Q * Q;
                type R2 = R * R;

                if (R2 < Q3) {
                        type sQ = sqrt(Q);
								type v = R / (sQ * Q);
								if (v < -1) v = -1;
								if (v > 1)  v = 1;
                        type theta = acos(v);
                        type phi_min = theta / 3;
                        type phi_max = phi_min + 6.2831853071795864769 / 3;
                        type phi_med = phi_min - 6.2831853071795864769 / 3;

                        return roots(
                                -a / 3 - 2 * sQ * cos(phi_min),
                                -a / 3 - 2 * sQ * cos(phi_med),
                                -a / 3 - 2 * sQ * cos(phi_max));
                } else {
                        type sR = R > 0 ? 1 : -1;
                        type R2mQ3 = R2 - Q3;
                        if (R2mQ3 < 0)
                        	R2mQ3 = 0;
                        type A = -sR * pow(fabs(R) + sqrt(R2mQ3), 1. / 3);
                        type B = (A != 0) ? Q / A : 0;
                        type x3 = -a / 3 + (A + B);
                        return roots(x3, x3, x3);
                }
}



/*
type ChemicalReactions::RhsHO(int l, type T, const Filt1D& filt, enum Layer layer) const
{
	type speed  = ChemicalReactions::HO.SpeedOfReaction(T);
	type result;
	switch(layer)
	{
		case current:
			result = -speed * filt.Component[HO].Nn[l] * filt.Component[O2].Nn[l];
			break;
		case middle:
			result = -speed * filt.Component[HO].Nmid[l] * filt.Component[O2].Nmid[l];
			break;
		case next:
			result = -speed * filt.Component[HO].Nn1[l] * filt.Component[O2].Nn1[l];
			break;
	}
	return result;
}


type ChemicalReactions::RhsMO(int l, type T, const Filt1D& filt, enum Layer layer) const
{
	type speed  = ChemicalReactions::LO.SpeedOfReaction(T);
	type result;
	switch(layer)
	{
		case current:
			result = -speed * filt.Component[LO].Nn[l] * filt.Component[O2].Nn[l];
			break;
		case middle:
			result = -speed * filt.Component[LO].Nmid[l] * filt.Component[O2].Nmid[l];
			break;
		case next:
			result = -speed * filt.Component[LO].Nn1[l] * filt.Component[O2].Nn1[l];
			break;
	}
	return result;
}


type ChemicalReactions::RhsO2(int l, type T, const Filt1D& filt, enum Layer layer) const
{
	type speedHO = ChemicalReactions::HO.SpeedOfReaction(T);
	type speedMO = ChemicalReactions::LO.SpeedOfReaction(T);
	type result;

	switch(layer)
	{
		case current:
			result = -23 * speedHO * filt.Component[HO].Nn[l] * filt.Component[O2].Nn[l] - 15.5 * speedMO * filt.Component[LO].Nn[l] * filt.Component[O2].Nn[l];
			break;
		case middle:
			result = -23 * speedHO * filt.Component[HO].Nmid[l] * filt.Component[O2].Nmid[l] - 15.5 * speedMO * filt.Component[LO].Nmid[l] * filt.Component[O2].Nmid[l];
			break;
		case next:
			result = -23 * speedHO * filt.Component[HO].Nn1[l] * filt.Component[O2].Nn1[l] - 15.5 * speedMO * filt.Component[LO].Nn1[l] * filt.Component[O2].Nn1[l];
			break; // TOFIX!!
	}
	return result;
}


type ChemicalReactions::RhsCO2(int l, type T, const Filt1D& filt, enum Layer layer) const
{
	type speedHO = ChemicalReactions::HO.SpeedOfReaction(T);
	type speedMO = ChemicalReactions::LO.SpeedOfReaction(T);
	type result;

	switch(layer)
	{
		case current:
			result = 15 * speedHO * filt.Component[HO].Nn[l] * filt.Component[O2].Nn[l] + 10 * speedMO * filt.Component[LO].Nn[l] * filt.Component[O2].Nn[l];
			break;
		case middle:
			result = 15 * speedHO * filt.Component[HO].Nmid[l] * filt.Component[O2].Nmid[l] + 10 * speedMO * filt.Component[LO].Nmid[l] * filt.Component[O2].Nmid[l];
			break;
		case next:
			result = 15 * speedHO * filt.Component[HO].Nn1[l] * filt.Component[O2].Nn1[l] + 10 * speedMO * filt.Component[LO].Nn1[l] * filt.Component[O2].Nn1[l];
			break; // TOFIX!!
	}
	return result;
}


type ChemicalReactions::RhsH2O(int l, type T, const Filt1D& filt, enum Layer layer) const
{
	type speedHO = ChemicalReactions::HO.SpeedOfReaction(T);
	type speedMO = ChemicalReactions::LO.SpeedOfReaction(T);
	type result;

	switch(layer)
	{
		case current:
			result = 16 * speedHO * filt.Component[HO].Nn[l] * filt.Component[O2].Nn[l] + 11 * speedMO * filt.Component[LO].Nn[l] * filt.Component[O2].Nn[l];
			break;
		case middle:
			result = 16 * speedHO * filt.Component[HO].Nmid[l] * filt.Component[O2].Nmid[l] + 11 * speedMO * filt.Component[LO].Nmid[l] * filt.Component[O2].Nmid[l];
			break;
		case next:
			result = 16 * speedHO * filt.Component[HO].Nn1[l] * filt.Component[O2].Nn1[l] + 11 * speedMO * filt.Component[LO].Nn1[l] * filt.Component[O2].Nn1[l];
			break;// TOFIX!!
	}
	return result;
}
*/
