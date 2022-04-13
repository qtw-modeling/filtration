#include "Well.hpp"
#include "Filt.h"

void Well::setPlace(type x, type y, Filt2D& filt) 
{ 
	xCoordinate = x; 
	xIndex      = (int)(x / filt.hx);
		
	yCoordinate = y;
	yIndex		= (int) (y / filt.hy);
}


void InjectorWell::injectMoles(Filt2D& filt, enum Components num, type speed, type t, type tInjection)
{	
	for (int i = 1; i <= NUM_COMPS; i++)
		filt.Component[i].Nmid.copy(filt.Component[i].Nwell);

	if (t <= tInjection)
	{
		filt.Component[num].Nwell(getXindex(), getYindex()) += filt.tau * speed;
	}
}


void InjectorWell::injectEnthalpy(Filt2D& filt, type T, type t, type tInjection)
{
	if (t <= tInjection)
	{
		filt.Hn.copy(filt.Hwell);

		vector<type> NinjectedInCell(NUM_COMPS + 1);
		
		for (int i = 1; i <= NUM_COMPS; i++)
		{
			
			NinjectedInCell[i] = filt.Component[i].Nwell(getXindex(), getYindex()) - filt.Component[i].Nmid(getXindex(), getYindex());
			filt.Hwell(getXindex(), getYindex()) += NinjectedInCell[i] * filt.MolarEnthalpyOfComponentInPhaseInCell(i, 1, T);
		}
	}
}


void RecoveryWell::recoverMoles(Filt2D& filt, type speed) 
{
	type newMolarDensity;
	for (int i = 1; i <= NUM_COMPS; i++)
	{
		if ((newMolarDensity = filt.Component[i].Nmid(getXindex(), getYindex()) - 
					filt.tau * speed * filt.Component[i].C(getXindex(), getYindex())) < 0)
			filt.Component[i].Nwell(getXindex(), getYindex()) = 1e1;
		else
			filt.Component[i].Nwell(getXindex(), getYindex()) = newMolarDensity;
	}
}


void RecoveryWell::recoverEnthalpy(Filt2D& filt)
{

	type NtotalInCellWell = 0;
	type NtotalInCellOld = 0;
	for (int i = 1; i <= NUM_COMPS; i++)
	{
		NtotalInCellWell += filt.Component[i].Nwell(getXindex(), getYindex());
		NtotalInCellOld  += filt.Component[i].Nmid(getXindex(), getYindex());
	}

	type alpha = (NtotalInCellOld - NtotalInCellWell) / NtotalInCellOld;			 
	filt.Hwell(getXindex(), getYindex()) *= (1 - alpha);
}