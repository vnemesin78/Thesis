#include "kalman_c_robust_em.hpp"
/**@fn virtual void tkalman_c_orignal_em  :: filter(const gsl_vector * const * observations,
                                                  unsigned int nb_observations) = 0
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette m�thode effectue le filtrage des donn�es par le filtre de Kalman Triple non supervis�.
 */
void kalman_c_robust_em  :: filter(const gsl_vector * const * observations,
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
	kalman_c_robust_em :: do_em_algorithm(observations);
	//Filtrage
	tkalman_robust_filter :: filter_without_equivalents(observations);
	//Restauration des donn�es
    tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
    tkalman_robust_filter :: compute_equivalents_sqrt_p_p_and_sqrt_p_f();

}

/**@fn virtual void kalman_c_robust_em :: smooth(const gsl_vector * const * observations,
                                                  unsigned int nb_observations)
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette m�thode effectue le lissage des donn�es par le filtre de Kalman Triple non supervis�.
 */
void kalman_c_robust_em :: smooth(const gsl_vector * const * observations,
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
	kalman_c_robust_em :: do_em_algorithm(observations);
	//Lissage
	tkalman_robust_filter :: smooth_without_equivalents(observations);
	//Restauration des donn�es
    tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
    tkalman_robust_filter :: compute_equivalents_sqrt_p_p_and_sqrt_p_f();
    tkalman_base :: compute_equivalents_x_s(observations);
    tkalman_robust_filter :: compute_equivalents_sqrt_p_s();
}

/**@fn bool kalman_c_robust_em  :: operator!()
 * @return
 - 0 si l'objet est valide
 - 1 sinon
 * @brief
 Cette m�thode teste la validit� de chaque attribut.
**/
bool kalman_c_robust_em :: operator!() const
{
    return (kalman_robust_em  :: operator !() || _size_x == _size_y);
}


/**@fn void kalman_c_robust_em  do_em_algorithm(const gsl_vector * const * observations);
 * @param[in] observations : observations
 * @brief
 * Cette m�thode estime les param�tres du filtre de Kalman couple � partir d'un jeu d'observations.
 */
void kalman_c_robust_em   :: do_em_algorithm(const gsl_vector * const * observations)
{
	for (unsigned int i = 0; i < nb_iter; ++ i)
	{
		kalman_robust_em :: estimate_parameters(observations,
												   i);
	}
    //Param�trage de la matrice p
	int signum;
	gsl_matrix_set_identity(_p);
	//Construction de p
	gsl_matrix_memcpy(&p_ux, &f_yx);
	gsl_matrix_memcpy(&p_uy, &f_yy);
	//Inversion de p_xx
	gsl_matrix_memcpy(mat_xx_1, &p_xx);
	gsl_linalg_LU_decomp(mat_xx_1, perm_x_1, &signum);
	gsl_linalg_LU_invert(mat_xx_1, perm_x_1, &p_inv_xx);

	//Calcul de p_inv_xy = -p_inv_xx��.p_xy
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
/**@fn kalman_c_robust_em  :: kalman_c_robust_em (const gsl_vector * x0,
                                                  const gsl_matrix * p0,
                                                  const gsl_matrix * f,
                                                  const gsl_matrix * q,
                                                  unsigned int n = 0,
                                                  unsigned int p = 0,
                                                  bool data = false);
 * @param[in] x0 : Esp�rance de l'�tat initial
 * @param[in] p0 : Matrice de covariance de l'�tat initial (remplac�e par sa d�composition de Cholesky dans certaines des classes filles)
 * @param[in] f : Matrice d'�volution
 * @param[in] q : Matrice de covariance (remplac�e par sa d�composition de Cholesky dans certaines des classes filles)
 * @param[in] n : Nombre d'observations (0 par d�faut)
 * @param[in] p : Nombre d'it�rations de l'EM (0 par d�faut et dans ce cas, cela �quivaut � un filtrage simple)
 * @param[in] data : Bool�en (True = suivi de l'EM et stockage des param�tres et de la vraisemblance � chaque it�ration, False = Pas de suivi de l'EM)
 * @brief
 * Constructeur
*/
kalman_c_robust_em  ::kalman_c_robust_em (const gsl_vector * x0,
                                              const gsl_matrix * p0,
                                              const gsl_matrix * f,
                                              const gsl_matrix * q,
                                              unsigned int n,
                                              unsigned int p,
                                              bool data)
: kalman_robust_em :: kalman_robust_em(x0, p0, f, q, n, p, data)
{};

