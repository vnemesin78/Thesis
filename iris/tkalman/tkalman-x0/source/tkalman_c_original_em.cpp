#include "tkalman_c_original_em.hpp"
/**@fn virtual void tkalman_c_orignal_em  :: filter(const gsl_vector * const * observations,
                                                  unsigned int nb_observations) = 0
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette méthode effectue le filtrage des données par le filtre de Kalman Triple non supervisé.
 */
void tkalman_c_original_em  :: filter(const gsl_vector * const * observations,
                                     unsigned int nb_observations)
{
	if (nb_observations != _n)
	{
		tkalman_base :: free_moments();
		tkalman_base :: initialize_moments();
		_n = nb_observations;
		tkalman_base :: alloc_moments();
	}
	//Algorithme EM
	tkalman_c_original_em :: do_em_algorithm(observations);
	//Filtrage
	tkalman_original_filter :: filter_without_equivalents(observations);
	//Restauration des données
    tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
    tkalman_original_filter :: compute_equivalents_p_p_and_p_f();

}

/**@fn virtual void tkalman_c_orignal_em :: smooth(const gsl_vector * const * observations,
                                                  unsigned int nb_observations)
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette méthode effectue le lissage des données par le filtre de Kalman Triple non supervisé.
 */
void tkalman_c_original_em  :: smooth(const gsl_vector * const * observations,
                                     unsigned int nb_observations)
{
	if (nb_observations != _n)
	{
		tkalman_base :: free_moments();
		tkalman_base :: initialize_moments();
		_n = nb_observations;
		tkalman_base :: alloc_moments();
	}
	//Algorithme EM
	tkalman_original_em :: do_em_algorithm(observations);
	//Lissage
	tkalman_original_filter :: smooth_without_equivalents(observations);
	//Restauration des données
    tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
    tkalman_original_filter :: compute_equivalents_p_p_and_p_f();
    tkalman_base :: compute_equivalents_x_s(observations);
    tkalman_original_filter :: compute_equivalents_p_s();
}

/**@fn bool tkalman_c_orignal_em :: operator!()
 * @return
 - 0 si l'objet est valide
 - 1 sinon
 * @brief
 Cette méthode teste la validité de chaque attribut.
**/
bool tkalman_c_original_em  :: operator!() const
{
    return (tkalman_original_em  :: operator !() || _size_x == _size_y);
}


/**@fn void tkalman_c_orignal_em  :: do_em_algorithm(const gsl_vector * const * observations);
 * @param[in] observations : observations
 * @brief
 * Cette méthode estime les paramètres du filtre de Kalman couple à partir d'un jeu d'observations.
 */
void tkalman_c_original_em  :: do_em_algorithm(const gsl_vector * const * observations)
{
	for (unsigned int i = 0; i < nb_iter; ++ i)
	{
		tkalman_original_em :: estimate_parameters(observations,
												   i);
	}
    //Paramètrage de la matrice p
	int signum;
	//Recopie de p
	gsl_matrix_memcpy(&p_xx, &f_yx);
	gsl_matrix_memcpy(&p_xy, &f_yy);

	//Inversion de p_xx
	gsl_matrix_memcpy(mat_xx_1, &p_xx);
	gsl_linalg_LU_decomp(mat_xx_1, perm_x_1, &signum);
	gsl_linalg_LU_invert(mat_xx_1, perm_x_1, &p_inv_xx);

	//Calcul de p_inv_xy = -p_inv_xx¯¹.p_xy
	gsl_blas_dgemm(CblasNoTrans,
				   CblasNoTrans,
				   1.0,
				   &p_inv_xx,
				   &p_xy,
				   0,
				   &p_inv_xy);
	gsl_matrix_scale(&p_inv_xy,
					 -1);

	//Correction des erreurs sur p (p_yx = 0, p_yy = Id)
	gsl_matrix_set_identity(&p_yy);
	gsl_matrix_set_zero(&p_yx);

	gsl_matrix_set_identity(&p_inv_yy);
	gsl_matrix_set_zero(&p_inv_yx);
}
/**@fn tkalman_c_original_em :: tkalman_x_original_em(const gsl_vector * x0,
                                                  const gsl_matrix * p0,
                                                  const gsl_matrix * f,
                                                  const gsl_matrix * q,
                                                  unsigned int n = 0,
                                                  unsigned int p = 0,
                                                  bool data = false);
 * @param[in] x0 : Espérance de l'état initial
 * @param[in] p0 : Matrice de covariance de l'état initial (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
 * @param[in] f : Matrice d'évolution
 * @param[in] q : Matrice de covariance (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
 * @param[in] n : Nombre d'observations (0 par défaut)
 * @param[in] p : Nombre d'itérations de l'EM (0 par défaut et dans ce cas, cela équivaut à un filtrage simple)
 * @param[in] data : Booléen (True = suivi de l'EM et stockage des paramètres et de la vraisemblance à chaque itération, False = Pas de suivi de l'EM)
 * @brief
 * Constructeur
*/
tkalman_c_original_em ::tkalman_c_original_em(const gsl_vector * x0,
                                              const gsl_matrix * p0,
                                              const gsl_matrix * f,
                                              const gsl_matrix * q,
                                              unsigned int n,
                                              unsigned int p,
                                              bool data)
: tkalman_original_em :: tkalman_original_em(x0, p0, f, q, n, p, data)
{}

