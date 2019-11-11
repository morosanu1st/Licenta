// This function computes the MAGNETIC FIELD at a given time for all grid cells in the simulation domain

#include "../utils/nr3.h"

void bfield(VecDoub &by, VecDoub &bz, VecDoub ey, VecDoub ez, int m, double c)
{
	int i;
	VecDoub *by1, *bz1, *ey1, *ez1;
	// compute magnetic field from Faraday's law
	// components of the Bfield are computed at a given time by knowing their values at a previous time-step
	// see Notebook 5, page 15 for more details
	// advance the magnetic field only for a half time-step: delta_t=0.5

#pragma omp parallel private(by1, bz1, ey1, ez1) num_threads(4)
	{
		by1 = &by;
		bz1 = &bz;
		ey1 = &ey;
		ez1 = &ez;
#pragma omp for
		for (i = 3; i <= m - 3; i++)
		{
			by1[0][i] = by1[0][i] + 0.5 * c * (ez1[0][i] - ez1[0][i - 1]);
			bz1[0][i] = bz1[0][i] - 0.5 * c * (ey1[0][i] - ey1[0][i - 1]);
		}
	}

	// apply periodic boundary conditions for the magnetic field
	// see Notebook 5, page 30 for more details
	by[2] = by[m - 3];
	bz[2] = bz[m - 3];

	// apply periodicity to the guard cells too
	by[1] = by[m - 4];
	bz[1] = bz[m - 4];

	by[0] = by[m - 5];
	bz[0] = bz[m - 5];

	by[m - 2] = by[3];
	bz[m - 2] = bz[3];

	by[m - 1] = by[4];
	bz[m - 1] = bz[4];
}
