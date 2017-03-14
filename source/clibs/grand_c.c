
/*******************************************************************************
 ** Tomas Lazauskas, 2016
 *******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

/*
 *	input: 	temperature,
 *			mu - chemical potential
 *			energies_mins - an array of local minimum energies of the stoichiometries
 *			energies_diffs - an array of energy differences between local minimum and all the attempts
 *			stoich_cnt - number of stoichiometries
 *			stoichiometries - stoichiometries array
 *			wms - an array of distribution functions for particular stoichiometry
 *	output: 0
 * */
int calculate_distribution_functions(
		double temperature,
		double mu,
		double* energies_mins,
		double* energies_diffs,
		int stoich_cnt,
		double* stoichiometries,
		double* wms)
{
	int m_i, mm_i;
	double m, mm;
	double kT, wms_sum_sq = 0.0, wms_sum;
	double sum_2;
	double e_min_m, e_min_mm, exp_expr;

	kT = (0.0257/298) * temperature;

	for (m_i = 0; m_i < stoich_cnt; m_i++) {

		m = stoichiometries[m_i];
		e_min_m = energies_mins[m_i];

		sum_2 = 0.0;

		for (mm_i = 0; mm_i < stoich_cnt; mm_i++) {

			mm = stoichiometries[mm_i];
			e_min_mm = energies_mins[mm_i];

			exp_expr = (((e_min_mm - e_min_m) + (mm-m)*mu) / kT)*(-1.0);

			sum_2 += exp(exp_expr) * energies_diffs[mm_i];
		}

		wms[m_i] = energies_diffs[m_i] / sum_2;

		wms_sum_sq += (wms[m_i] * wms[m_i]);
	}

	wms_sum = sqrt(wms_sum_sq);

	for (m_i = 0; m_i < stoich_cnt; m_i++) {
		wms[m_i] = wms[m_i] / wms_sum;
	}

	return 0;
}

/*
 *  input: 	temperature,
 *			mu - chemical potential
 *			area - surface area
 *			energies_GM - global energy minimum of all the stoichiometries-
 *			energies_mins - an array of local minimum energies of the stoichiometries
 *			energies_diffs - an array of energy differences between local minimum and all the attempts
 *			stoich_cnt - number of stoichiometries
 *			stoichiometries - stoichiometries array
 *	output: omega
 * */
double calculate_gamma(double temperature, double mu, double energies_GM, double area,
		double* energies_mins, double* energies_diffs, int stoich_cnt, double* stoichiometries)
{
	double gamma = 0.0;
	double kT, mm, e_min_mm, exp_expr, energies_GM_area;
	int mm_i;
	long double exp_result, sum_2;

	energies_GM_area = area * energies_GM;

	kT = (0.0257/298) * temperature;

	sum_2 = 0.0;

	for (mm_i = 0; mm_i < stoich_cnt; mm_i++) {
		mm = stoichiometries[mm_i];
		e_min_mm = energies_mins[mm_i];

		exp_expr = (((e_min_mm - energies_GM_area) + mm*mu - energies_GM_area) / kT)*(-1.0);



//
		printf(" exp_expr %f\n", exp_expr);


		exp_result = exp(exp_expr);

		printf(" 222222 %f\n", exp_result);

//		printf(" energies_diffs[mm_i] %f\n", energies_diffs[mm_i]);
//		printf(" other %f\n", (energies_GM_area/ kT)*(-1.0));




		sum_2 += exp_result * energies_diffs[mm_i];
	}



	gamma = (energies_GM_area - kT * log(sum_2)) / area;

	return gamma;
}
