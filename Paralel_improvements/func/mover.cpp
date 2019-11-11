// This function MOVE all the particles over one time-step

#include <cmath>
#include "../utils/nr3.h"

typedef struct t
{
	VecDoub *x;
	VecDoub *vx;
	VecDoub *vy;
	VecDoub *vz;
	VecDoub *ex;
	VecDoub *ey;
	VecDoub *ez;
	VecDoub *by;
	VecDoub *bz;
	double ex0;
	double ey0;
	double ez0;
	double bx0;
	double by0;
	double bz0;
	double qme;
	double qmi;
	double c;
	int start;
	int end;
	int np;
} param_struct;

void *mover_thread(void *param)
{
	param_struct *p = (param_struct *)param;

	VecDoub x = *(p->x);
	VecDoub vx = *(p->vx);
	VecDoub vy = *(p->vy);
	VecDoub vz = *(p->vz);
	VecDoub ex = *(p->ex);
	VecDoub ey = *(p->ey);
	VecDoub ez = *(p->ez);
	VecDoub by = *(p->by);
	VecDoub bz = *(p->bz);
	double ex0 = p->ex0;
	double ey0 = p->ey0;
	double ez0 = p->ez0;
	double bx0 = p->bx0;
	double by0 = p->by0;
	double bz0 = p->bz0;
	double qme = p->qme;
	double qmi = p->qmi;
	double c = p->c;
	int start = p->start;
	int end = p->end;
	int np = p->np;

	int k, i, j;
	double ex_p, ey_p, ez_p, by_p, bz_p, qm, gamma_p, eps_x, eps_y, eps_z, beta_x, beta_y, beta_z, ux1, uy1, uz1, wx, wy, wz,
		ux2, uy2, uz2, ux, uy, uz;
	double c_square = c * c;
	double expr_1;
	for (k = start; k < end; k++)
	{
		// particle 'k' is located between i and i+1 and between j+1/2 and j+3/2
		// see diagram on page 22 of Notebook 5
		i = x[k];
		j = x[k] - 0.5;

		// compute the electric and magnetic fields in the actual position of each particle
		// see Notebook 5, page 45 for details
		// relocation from half-integer to full-integer grid points is needed for Ex: interpolarion from full-integer grid-points
		ex_p = (i + 1 - x[k]) * 0.5 * (ex[i - 1] + ex[i]) + (x[k] - i) * 0.5 * (ex[i] + ex[i + 1]);

		// no relocation is needed for Ey and Ez: interpolation from half-integer grid-points
		ey_p = (j + 1.5 - x[k]) * ey[j] + (x[k] - j - 0.5) * ey[j + 1];
		ez_p = (j + 1.5 - x[k]) * ez[j] + (x[k] - j - 0.5) * ez[j + 1];

		// relocation from full- to half-integer grid points is needed for By and Bz: interpolation from half-integer grid-points
		by_p = (j + 1.5 - x[k]) * 0.5 * (by[j] + by[j + 1]) + (x[k] - j - 0.5) * 0.5 * (by[j + 1] + by[j + 2]);
		bz_p = (j + 1.5 - x[k]) * 0.5 * (bz[j] + bz[j + 1]) + (x[k] - j - 0.5) * 0.5 * (bz[j + 1] + bz[j + 2]);

		// check if particle 'k' is electron or ion
		if (k < np)
		{
			// electron
			qm = qme;
		}
		else
		{
			// ion
			qm = qmi;
		}

		// factor proportional to E-field
		// see Notebook 5, pg. 26 for details
		// note that the external electric field is added to the internal one at this step
		eps_x = qm * 0.5 * (ex_p + ex0);
		eps_y = qm * 0.5 * (ey_p + ey0);
		eps_z = qm * 0.5 * (ez_p + ez0);

		// compute the relativistic factor gamma
		// see page 26 on Notebook 5 for more details
		gamma_p = 1.0 / sqrt(1.0 - (vx[k] * vx[k] + vy[k] * vy[k] + vz[k] * vz[k]) / c_square);

		// compute u1 as shown in Notebook 5, page 26
		ux1 = gamma_p * vx[k] + eps_x;
		uy1 = gamma_p * vy[k] + eps_y;
		uz1 = gamma_p * vz[k] + eps_z;

		// compute the new relativistic factor
		gamma_p = sqrt(1.0 + (ux1 * ux1 + uy1 * uy1 + uz1 * uz1) / c_square);

		// factor proportional to B-field
		// see Notebook 5, pg. 26 for details
		// note that the external magnetic field is added to the internal one at this step
		beta_x = (qm * 0.5 / gamma_p) * (bx0 / c);
		beta_y = (qm * 0.5 / gamma_p) * ((by_p + by0) / c);
		beta_z = (qm * 0.5 / gamma_p) * ((bz_p + bz0) / c);

		// intermediate quantity: Notebook 5, page 25
		wx = ux1 + uy1 * beta_z - uz1 * beta_y;
		wy = uy1 + uz1 * beta_x - ux1 * beta_z;
		wz = uz1 + ux1 * beta_y - uy1 * beta_x;

		// compute u2 as shown in Notebook 5, page 26

		expr_1 = 2 / (1 + beta_x * beta_x + beta_y * beta_y + beta_z * beta_z);
		ux2 = ux1 + expr_1 * (wy * beta_z - wz * beta_y);
		uy2 = uy1 + expr_1 * (wz * beta_x - wx * beta_z);
		uz2 = uz1 + expr_1 * (wx * beta_y - wy * beta_x);

		// compute u as shown in Notebook 5, page 26
		ux = ux2 + eps_x;
		uy = uy2 + eps_y;
		uz = uz2 + eps_z;

		// compute the new relativistic factor
		gamma_p = sqrt(1.0 + (ux * ux + uy * uy + uz * uz) / c_square);

		// compute the new velocity after one time-step: page 26, Notebook 5
		vx[k] = ux / gamma_p;
		vy[k] = uy / gamma_p;
		vz[k] = uz / gamma_p;

		// move the particle over one time-step: page 25 on Notebook 5
		x[k] = x[k] + vx[k];
	}

	return NULL;
}

void mover(VecDoub &x, VecDoub &vx, VecDoub &vy, VecDoub &vz, VecDoub ex, VecDoub ey, VecDoub ez, VecDoub by, VecDoub bz,
		   double ex0, double ey0, double ez0, double bx0, double by0, double bz0, double qme, double qmi, double c, int np, int m)
{
	int numThreads = 4;
	int chunk = np * 2 / numThreads;
	pthread_t *threads = (pthread_t *)malloc(sizeof(pthread_t) * numThreads);
	param_struct **params = (param_struct **)malloc(sizeof(param_struct *)*numThreads);
	for (int i = 0; i < numThreads; i++)
	{
		params[i] = (param_struct *)malloc(sizeof(param_struct));
		params[i]->x = &x;
		params[i]->vx = &vx;
		params[i]->vy = &vy;
		params[i]->vz = &vz;
		params[i]->ex = &ex;
		params[i]->ey = &ey;
		params[i]->ez = &ez;
		params[i]->by = &by;
		params[i]->bz = &bz;
		params[i]->ex0 = ex0;
		params[i]->ey0 = ey0;
		params[i]->ez0 = ez0;
		params[i]->bx0 = bx0;
		params[i]->by0 = by0;
		params[i]->bz0 = bz0;
		params[i]->qme = qme;
		params[i]->qmi = qmi;
		params[i]->c = c;
		params[i]->np = np;
		params[i]->start = i * chunk;
		if (i == numThreads - 1)
		{
			params[i]->end = np * 2;
		}
		else
		{
			params[i]->end = (i + 1) * chunk;
		}

		pthread_create(&threads[i], NULL, mover_thread, (void *)params[i]);
	}
	for (int i = 0; i < numThreads; i++)
	{
		pthread_join(threads[i], NULL);
	}
	for (int i = 0; i < numThreads; i++)
	{
		free(params[i]);
	}
	free(params);
	free(threads);
}
