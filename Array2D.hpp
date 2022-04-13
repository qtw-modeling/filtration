#ifndef ARRAY_2D
#define ARRAY_2D

#define cycle(i, iBegin, iEnd, j, jBegin, jEnd) for (int i = iBegin; i <= iEnd; i++) { \
												for (int j = jBegin; j <= jEnd; j++)

#include <vector>
#include <string>
#include <iostream>
#include <cassert>



template <class T>
class Array2D {
private:
	std::vector <T> array;
	int xSize, ySize;

public:
	Array2D(int xSize_, int ySize_): 
	xSize(xSize_), ySize(ySize_), array(xSize_ * ySize_) {}

	Array2D() : xSize(0), ySize(0) {}

	T& operator()(int i, int j) 
	{
		checkIJ(i, j);
		return array[i + j * xSize]; 
	}

	const T& operator()(int i, int j) const 
	{ 
		checkIJ(i, j);
		return array[i + j * xSize]; 
	}

	void checkIJ(int i, int j) const
	{
		assert(i >= 0 && i <= lastIndexX() && j >= 0 && j <= lastIndexY()); 
	}

	int lastIndexX() const { return (int)(xSize - 1); }
	int lastIndexY() const { return (int)(ySize - 1); }

	const int getXsize() const { return xSize; }
	const int getYsize() const { return ySize; }


	void print() const
	{
		cycle(j, 0, lastIndexY(), i, 0, lastIndexX()) {
				std::cout << (*this)(i, j) << " ";
			}
			std::cout << std::endl;
		}
	}

	void copy (Array2D <T>& destination) {
		cycle(j, 0, lastIndexY(), i, 0, lastIndexX()) {
				destination(i, j) = (*this)(i, j);
			}
		}
	}
};


#endif