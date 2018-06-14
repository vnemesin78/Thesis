#include "tkalman_argmax.hpp"
/**@fn void tkalman_compute_argmax(gsl_vector * t_0,
								   gsl_matrix * sqrt_q_0,
								   gsl_matrix * f,
								   gsl_matrix * sqrt_q,
								   const gsl_matrix * sqrt_c_00,
								   const gsl_matrix * sqrt_c_01,
								   const gsl_matrix * sqrt_c_11,
								   unsigned int n,
								   const gsl_vector * t_s_0,
								   const gsl_matrix * sqrt_q_s_0,
								   gsl_matrix * mat_tt,
								   gsl_permutation * perm_t)
 * @param t_0 : Espérance de t_0 (NULL pour ne pas l'estimer)
 * @param sqrt_q_0 : Matrice de covariance de t_0 (NULL pour ne pas l'estimer)
 * @param f : Matrice d'évolution
 * @param q : Matrice de covariance du bruit
 * @param sqrt_c_00: Somme
 * @param sqrt_c_10 : Somme
 * @param sqrt_c_11 : Somme
 * @param n : nombre d'observations
 * @param t_s_0 : t_{0|N} (NULL pour ne pas l'estimer)
 * @param sqrt_q_s_0 : racine de Q_{0|N}
 * @param mat_tt : matrice de taille (n_t, n_t) allouée
 * @param perm_t : permutation de taille (n_t) préallouée
 * @brief
 * Cette fonction calcule les paramètres optimaux selon la fonction auxiliaire pour le filtre de Kalman triplet.
 */
void tkalman_compute_argmax(gsl_vector * t_0,
						    gsl_matrix * sqrt_q_0,
						    gsl_matrix * f,
						    gsl_matrix * sqrt_q,
						    const gsl_matrix * sqrt_c_00,
						    const gsl_matrix * sqrt_c_01,
						    const gsl_matrix * sqrt_c_11,
						    unsigned int n,
						    const gsl_vector * t_s_0,
						    const gsl_matrix * sqrt_q_s_0,
						    gsl_matrix * mat_tt,
						    gsl_permutation * perm_t)
{
	//Réestimation falcutative
		//Calcul de x0, p0
			if (t_0 != NULL && t_s_0 != NULL)
				gsl_vector_memcpy(t_0,
								  t_s_0);


			if (sqrt_q_0 != NULL && sqrt_q_s_0 != NULL)
				gsl_matrix_memcpy(sqrt_q_0,
								  sqrt_q_s_0);

		//Calcul de F
			//Inversion de sqrt_c_00
				gsl_permutation_init(perm_t);
				gsl_linalg_LU_invert(sqrt_c_00,
									 perm_t,
									 sqrt_q);
			//Calcul de F = sqrt_c01^T [C00^{1/2}]^{-1}
				gsl_blas_dgemm(CblasTrans,
							   CblasTrans,
							   1.0,
							   sqrt_c_01,
							   sqrt_q,
							   0.0,
							   f);

		//Calcul de Q
		gsl_matrix_memcpy(sqrt_q,
						  sqrt_c_11);

		gsl_matrix_scale(sqrt_q,
						 sqrt(1.0/ ((double) n)) );
}



//Objet de prédiction
/**@fn tkalman_argmax :: tkalman_argmax(const gsl_matrix * sqrt_c_00,
									    const gsl_matrix * sqrt_c_01,
									    const gsl_matrix * sqrt_c_11);
* @param f2_x : [F2xx, F2xy] 
* @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
* @param q2_xy : Qxy.Qyy
* @brief
* Ce constructeur alloue les variables temp. de la prédiction du filtre de Kalman couple.
*/
tkalman_argmax :: tkalman_argmax(const gsl_matrix * sqrt_c_00,
						         const gsl_matrix * sqrt_c_01,
						         const gsl_matrix * sqrt_c_11)
{
	tkalman_argmax :: initialize();
	tkalman_argmax :: setup (sqrt_c_00,
						     sqrt_c_01,
						     sqrt_c_11);
}
/**@fn int tkalman_argmax :: setup(const gsl_matrix * sqrt_c_00,
						           const gsl_matrix * sqrt_c_01,
						           const gsl_matrix * sqrt_c_11);
 * @param f2_x : [F2xx, F2xy]
 * @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
 * @param q2_xy : Qxy.Qyy
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction permet de modifier les paramètres de la prédiction (size_x et size_y)
**/
int tkalman_argmax :: setup (const gsl_matrix * sqrt_c_00,
						     const gsl_matrix * sqrt_c_01,
						     const gsl_matrix * sqrt_c_11)
	
{
	unsigned int size_t = sqrt_c_00->size1;
	
	if (size_t != _size_t)
	{
		tkalman_argmax :: free();
		tkalman_argmax :: initialize();;
		if (tkalman_argmax :: set_params(sqrt_c_00,
								 sqrt_c_01,
								 sqrt_c_11))
		{
			tkalman_argmax :: free();
			tkalman_argmax :: initialize();
			return 1;
		}
		_size_t = size_t;
		if ( tkalman_argmax :: alloc() )
		{
			tkalman_argmax :: free();
			tkalman_argmax :: initialize();
			return 1;
		}
		tkalman_argmax :: create_views();
	}	
	return 0;
}

/**@fn tkalman_argmax :: ~ tkalman_argmax()
 * @brief
 * Destructeur
 */
tkalman_argmax :: ~tkalman_argmax()
{
	tkalman_argmax :: free();
	tkalman_argmax :: initialize();
}

/**@fn bool tkalman_argmax :: operator!() const;
 * @return
 * - 0 si l'objet est correctement alloué
 * - 1 sinon
 * @brief
 * Check de l'objet.
 */
bool tkalman_argmax :: operator!() const
{
	return (! (mat_tt && perm_t && _sqrt_c_00 && _sqrt_c_01 && _sqrt_c_11) );
}

/**@fn void tkalman_argmax :: free();
 * @brief
 * Cette fonction désalloue la mémoire utilisée par les variables tmp.
 */
void tkalman_argmax :: free()
{
	if (mat_tt)
		gsl_matrix_free(mat_tt);
	if (perm_t)
		gsl_permutation_free(perm_t);
}

/**@fn int tkalman_argmax :: alloc();
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction alloue la mémoire utilisée par les variables tmp.
 */
int tkalman_argmax :: alloc()
{
	if (!mat_tt)
	{
		mat_tt = gsl_matrix_alloc(_size_t, 
										 _size_t);
	}
	if (!perm_t)
	{
		perm_t = gsl_permutation_alloc(_size_t);
	}
	return (!(mat_tt && perm_t) );
}		
/**@fn void tkalman_argmax :: initialize();
 * @brief
 * Cette fonction met les variables tmp à NULL.
 */
void tkalman_argmax :: initialize()
{
	_size_t = 0;
	mat_tt= NULL;
	perm_t = NULL;
	_sqrt_c_00 = NULL;
	_sqrt_c_01 = NULL;
	_sqrt_c_11 = NULL;
}

/**@fn virtual int tkalman_argmax :: set_params(const gsl_matrix * sqrt_c_00,
												const gsl_matrix * sqrt_c_01,
												const gsl_matrix * sqrt_c_11);
 * @param f2_x : [F2xx, F2xy]
 * @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
 * @param q2_xy : Qxy.Qyy
 * @brief
 * Cette méthode change les paramètres de la prédiction.
 */
int tkalman_argmax :: set_params(const gsl_matrix * sqrt_c_00,
						         const gsl_matrix * sqrt_c_01,
						         const gsl_matrix * sqrt_c_11)
{

	if (!sqrt_c_00 || !sqrt_c_01 || !sqrt_c_11)
		return 1;
	_sqrt_c_00 = sqrt_c_00;
	_sqrt_c_01 = sqrt_c_01;
	_sqrt_c_11 = sqrt_c_11;
	return 0;
}
