#include "Filt.h"
#include <fstream>

void Filt2D::VisualizeVTK(const int step)
{
	char fn[256];
	sprintf(fn, "../Output/resWith.%d.vtk", step);

	std::fstream f(fn, std::ios::out);
	f << "# vtk DataFile Version 3.0" << std::endl;
	f << "Solution" << std::endl;
	f << "ASCII" << std::endl;
	f << "DATASET RECTILINEAR_GRID" << std::endl;
	f << "DIMENSIONS " << P.getXsize() - 1 << " " << P.getYsize() - 1 << " 1" << std::endl;
	
	f << "X_COORDINATES " << P.getXsize() - 1 << " double" << std::endl;
	
	for (int i = 1; i <= P.getXsize() - 1; i++)
		f << i * hx << " ";
	f << std::endl;
	
	f << "Y_COORDINATES " << P.getYsize() - 1 << " double" << std::endl;
	
	for (int i = 1; i <= P.getYsize() - 1; i++)
		f << i * hy << " ";
	f << std::endl;
	f << "Z_COORDINATES 1 double\n0" << std::endl;
	f << "CELL_DATA " << (P.getXsize() - 2) * (P.getYsize() - 2) << std::endl;
	
	
	// printing all the data in the file
	f << "SCALARS Pressure double\nLOOKUP_TABLE default" << std::endl;
	for (int j = 1; j <= P.lastIndexY() - 1; j++) 
	{
		for (int i = 1; i <= P.lastIndexX() - 1; i++)
			f << P(i, j) << " ";
		f << std::endl;
	}
	

	f << "SCALARS Sat_Liquid double\nLOOKUP_TABLE default" << std::endl;
	for (int j = 1; j <= P.lastIndexY() - 1; j++) 
	{
		for (int i = 1; i <= P.lastIndexX() - 1; i++)
			f << Phase[Liquid].S(i, j) << " ";
		f << std::endl;
	}


	f << "SCALARS Sat_Water double\nLOOKUP_TABLE default" << std::endl;
	for (int j = 1; j <= P.lastIndexY() - 1; j++) 
	{
		for (int i = 1; i <= P.lastIndexX() - 1; i++)
			f << Phase[Water].S(i, j) << " ";
		f << std::endl;
	}


	f << "SCALARS Sat_Gas double\nLOOKUP_TABLE default" << std::endl;
	for (int j = 1; j <= P.lastIndexY() - 1; j++) 
	{
		for (int i = 1; i <= P.lastIndexX() - 1; i++)
			f << Phase[Gas].S(i, j) << " ";
		f << std::endl;
	}


	f << "SCALARS C_HO double\nLOOKUP_TABLE default" << std::endl;
	for (int j = 1; j <= P.lastIndexY() - 1; j++) 
	{
		for (int i = 1; i <= P.lastIndexX() - 1; i++)
			f << Component[HO].C(i, j) << " ";
		f << std::endl;
	}
	

	f << "SCALARS C_LO double\nLOOKUP_TABLE default" << std::endl;
	for (int j = 1; j <= P.lastIndexY() - 1; j++) 
	{
		for (int i = 1; i <= P.lastIndexX() - 1; i++)
			f << Component[LO].C(i, j) << " ";
		f << std::endl;
	}


	f << "SCALARS C_H2O double\nLOOKUP_TABLE default" << std::endl;
	for (int j = 1; j <= P.lastIndexY() - 1; j++) 
	{
		for (int i = 1; i <= P.lastIndexX() - 1; i++)
			f << Component[H2O].C(i, j) << " ";
		f << std::endl;
	}


	f << "SCALARS C_CO2 double\nLOOKUP_TABLE default" << std::endl;
	for (int j = 1; j <= P.lastIndexY() - 1; j++) 
	{
		for (int i = 1; i <= P.lastIndexX() - 1; i++)
			f << Component[CO2].C(i, j) << " ";
		f << std::endl;
	}


	f << "SCALARS C_O2 double\nLOOKUP_TABLE default" << std::endl;
	for (int j = 1; j <= P.lastIndexY() - 1; j++) 
	{
		for (int i = 1; i <= P.lastIndexX() - 1; i++)
			f << Component[O2].C(i, j) << " ";
		f << std::endl;
	}


	f << "SCALARS Temperature double\nLOOKUP_TABLE default" << std::endl;
	for (int j = 1; j <= P.lastIndexY() - 1; j++) 
	{
		for (int i = 1; i <= P.lastIndexX() - 1; i++)
			f << Tn1(i, j) << " ";
		f << std::endl;
	}


	f << "SCALARS N_H2O double\nLOOKUP_TABLE default" << std::endl;
	for (int j = 1; j <= P.lastIndexY() - 1; j++) 
	{
		for (int i = 1; i <= P.lastIndexX() - 1; i++)
			f << Component[H2O].Nn1(i, j) << " ";
		f << std::endl;
	}





	f.close();
}