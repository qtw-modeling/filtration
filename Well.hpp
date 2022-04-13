#ifndef _WELL_H_
#define _WELL_H_

#include "Array2D.hpp"
#include "Common.h"

class Filt2D;

class Well {
protected:
	type xCoordinate, yCoordinate;
	int xIndex, yIndex;

public:
	void setPlace(type x, type y, Filt2D& filt);

	type getCoordinateX() { return xCoordinate; }
	type getCoordinateY() { return yCoordinate; }

	int getXindex() { return xIndex; }
	int getYindex() { return yIndex; }
};


class InjectorWell : public Well {
public:
	void injectMoles(Filt2D& filt, enum Components comp, type speed, type t, type tInjection);
	void injectEnthalpy(Filt2D& filt, type speed, type t, type tInjection);
};


class RecoveryWell : public Well {
public:
	void recoverMoles(Filt2D& filt, type speed);
	void recoverEnthalpy(Filt2D& filt);
};


#endif
