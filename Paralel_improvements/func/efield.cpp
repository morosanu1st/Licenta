// This function computes the ELECTRIC FIELD at a given time for all grid cells in the simulation domain

#include "../utils/nr3.h"

void efield(VecDoub &ex, VecDoub &ey, VecDoub &ez, VecDoub by, VecDoub bz, VecDoub jxe, VecDoub jye, VecDoub jze,
			VecDoub jxi, VecDoub jyi, VecDoub jzi, int m, double c)
{
	int i;
	int numThreads=m/800;
	if (numThreads==0){
		numThreads=1;
	}
	// compute electric field from Ampere's law
	// components of the Efield are computed at a given time by knowing their values at a previous time-step
	// see Notebook 5, page 15 for more details
	VecDoub *ex1, *ey1, *ez1, *by1, *bz1, *jxe1, *jye1, *jze1, *jxi1, *jyi1, *jzi1;
#pragma omp parallel private(ex1, ey1, ez1, jxe1, jye1, jze1, jxi1, jyi1, jzi1, bz1, by1) num_threads(numThreads)
	{
		ex1 = &ex;
		ey1 = &ey;
		ez1 = &ez;
		by1 = &by;
		bz1 = &bz;
		jxe1 = &jxe;
		jye1 = &jye;
		jze1 = &jze;
		jxi1 = &jxi;
		jyi1 = &jyi;
		jzi1 = &jzi;

#pragma omp for
		for (i = 2; i < m - 3; i++)
		{
			ex1[0][i] = ex1[0][i] - (jxe1[0][i] + jxi1[0][i]);
			ey1[0][i] = ey1[0][i] - (jye1[0][i] + jyi1[0][i]) - c * (bz1[0][i + 1] - bz1[0][i]);
			ez1[0][i] = ez1[0][i] - (jze1[0][i] + jzi1[0][i]) + c * (by1[0][i + 1] - by1[0][i]);
		}
	}
	// apply periodic boundary conditions for the electric field
	// see Notebook 5, page 29 for more details
	ex[m - 3] = ex[2];
	ey[m - 3] = ey[2];
	ez[m - 3] = ez[2];

	// apply periodicity to the guard cells too
	ex[1] = ex[m - 4];
	ey[1] = ey[m - 4];
	ez[1] = ez[m - 4];

	ex[0] = ex[m - 5];
	ey[0] = ey[m - 5];
	ez[0] = ez[m - 5];

	ex[m - 2] = ex[3];
	ey[m - 2] = ey[3];
	ez[m - 2] = ez[3];

	ex[m - 1] = ex[4];
	ey[m - 1] = ey[4];
	ez[m - 1] = ez[4];
}
