#include "tkalman_em_initialization.hpp"

/**@fn void tkalman_EM_initialization(gsl_vector * x_0,
									  gsl_matrix * sqrt_p_0,
									  gsl_matrix * sqrt_q_yy,
									  const gsl_vector * const * observations,
									  const unsigned int n,
									  gsl_vector * vect_x)
 * @param x_0 : espérance de l'état initial pour l'EM
 * @param sqrt_p_0 : racine de la matrice de covariance de l'état initial pour l'EM
 * @param sqrt_q_xx : racine de la matrice de covariance du bruit de process
 * @param sqrt_q_yy : racine de la matrice de covariance du bruit de mesure
 * @param[in] observation : observations
 * @param[in] n : nombre d'observations
 * @param vect_x : vecteur de taille x alloué
 * @brief
 * Cette fonction estime l'espérance de l'état initial, sa matrice de covariance et la matrice de covariance du bruit de mesure.
 * @warning x et y doivent avoir même dimension
 */
void tkalman_EM_initialization(gsl_vector * x_0,
							   gsl_matrix * sqrt_p_0,
							   gsl_matrix * sqrt_q_yy,
							   const gsl_vector * const * observations,
							   const unsigned int n,
							   gsl_vector * vect_x)
{
	unsigned int size_x;
	size_x = x_0->size;
	//Estimation de la moyenne de y
	gsl_vector_set_zero(x_0);
	gsl_matrix_set_zero(sqrt_q_yy);
	for (unsigned int i = 0; i < n ; ++i)
	{
		gsl_vector_add(x_0, observations[i]);
	}
	gsl_vector_scale(x_0, 1.0 / ( (double) n));
	
	//Estimation de la covariance
	{

		for (unsigned int i = 0; i < n ; ++i)
		{
			gsl_vector_memcpy (vect_x, 
							   observations[i]);
			gsl_vector_sub (vect_x, 
							x_0);
			gsl_blas_dger (1.0, 
						   vect_x, 
						   vect_x, 
						   sqrt_q_yy);
		}
		gsl_matrix_scale(sqrt_q_yy, 1.0 / (n - 1));
		gsl_linalg_cholesky_decomp (sqrt_q_yy);
		gsl_matrix_memcpy(sqrt_p_0, sqrt_q_yy);
	}
}




