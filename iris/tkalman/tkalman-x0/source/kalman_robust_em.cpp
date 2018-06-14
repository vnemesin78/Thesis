#include "kalman_robust_em.hpp"
/**@fn bool kalman_robust_em :: check_tmp()
 * @return
 - 0 si les tmp ont été bien alloués
 - 1 sinon
 * @brief
 Cette méthode teste la bonne allocation des variables temporaires utilisées par l'objet.
**/
bool kalman_robust_em :: check_tmp() const
{
    return (! (mat_4x4x_1 && mat_4xp12x_1 && mat_2tp1t_1 && vect_4x_1 && perm_y_1) );
}

/**@fn bool kalman_robust_em :: check_data()
 * @return
 - 0 si les données de suivi de l'EM ont été bien alloués
 - 1 sinon
 * @brief
 Cette méthode teste la bonne allocation des données de suivi de l'EM.
**/
bool kalman_robust_em :: check_data() const
{
	return (! (_log_likelihood && _x0_est && _p0_est && _f_est && _q_est) );
}



/**@fn void kalman_robust_em :: free();
* @brief
 Cette méthode libère la mémoire utilisée par les attributs de l'objet.
**/
void kalman_robust_em :: free()
{
    kalman_robust_em :: free_tmp();
    kalman_robust_em :: free_data();
}

/**@fn void kalman_robust_em :: free_tmp();
 * @brief
 Cette méthode libère la mémoire utilisée par les temporaires de l'objet.
**/
void kalman_robust_em :: free_tmp()
{
    if (mat_4x4x_1)
        gsl_matrix_free(mat_4x4x_1);

    if (mat_4xp12x_1)
        gsl_matrix_free(mat_4xp12x_1);

    if (mat_2tp1t_1)
        gsl_matrix_free(mat_2tp1t_1);

    if (vect_4x_1)
        gsl_vector_free(vect_4x_1);
        
    if (perm_y_1)
        gsl_permutation_free(perm_y_1);
}


/**@fn void kalman_robust_em :: free_data();
 *  @brief
 Cette méthode libère la mémoire utilisée par les données de suivi de l'EM.
**/
void kalman_robust_em :: free_data()
{
	tkalman_esperance_unref(_x0_est, nb_iter);
	tkalman_covariance_unref(_p0_est, nb_iter);
	tkalman_covariance_unref(_f_est, nb_iter);
	tkalman_covariance_unref(_q_est, nb_iter);
	if (_log_likelihood)
		delete[] _log_likelihood;
}


/**@fn int kalman_robust_em :: alloc(unsigned int p, bool data);
 * @param p : nombre d'itérations de l'EM
 * @param data : booléen
 * @return
  - 0 si l'allocation des attributs de l'objet s'est bien déroulée
  - 1 sinon
  *@brief
  Cette méthode alloue les attributs de l'objet.
 */
int kalman_robust_em :: alloc(unsigned int p, bool data)
{
	if (kalman_robust_em :: alloc_tmp())
		return 1;
	nb_iter = p;
	if (data)
	{
		if (kalman_robust_em :: alloc_data())
			return 1;
	}
	return 0;
}

/**@fn int kalman_robust_em :: alloc_tmp();
 * @return
  - 0 si l'allocation des temporaires de l'objet s'est bien déroulée
  - 1 sinon
  *@brief
  Cette méthode alloue les temporaires de l'objet.
 */
int kalman_robust_em :: alloc_tmp()
{
    gsl_matrix_view view_1;
    gsl_vector_view view_2;
    if (!mat_4x4x_1)
        mat_4x4x_1 = gsl_matrix_alloc(4 * _size_x, 4 * _size_x);

    if (!mat_4xp12x_1)
        mat_4xp12x_1 = gsl_matrix_alloc(4 * _size_x + 1, 2 * _size_x);


    if (!mat_2tp1t_1)
        mat_2tp1t_1 = gsl_matrix_alloc(2 * _size_t + 1, _size_t);

    if (!vect_4x_1)
        vect_4x_1 = gsl_vector_alloc(4 * _size_x);

    if (!perm_y_1)
        perm_y_1 = gsl_permutation_alloc(_size_y);

    if (kalman_robust_em :: check_tmp())
        return 1;

    //Création des vues
        //Mat_4x4x
        view_1 = gsl_matrix_submatrix(mat_4x4x_1, 0, 0, _size_x, _size_x);
        mat_4x4x_1_view_00 = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4x4x_1, _size_x, _size_x, _size_x, _size_x);
        mat_4x4x_1_view_11 = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4x4x_1, 0, 2 * _size_x, _size_x, _size_x);
        mat_4x4x_1_view_02 = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4x4x_1, 2 * _size_x, 0, _size_x, _size_x);
        mat_4x4x_1_view_20 = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4x4x_1, 3 * _size_x, 2 * _size_x, _size_x, _size_x);
        mat_4x4x_1_view_32 = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4x4x_1, 3 * _size_x, 3 * _size_x, _size_x, _size_x);
        mat_4x4x_1_view_33 = view_1.matrix;

		view_1 = gsl_matrix_submatrix(mat_4x4x_1, 2 * _size_x, 2 * _size_x, 2 * _size_x, 2 * _size_x);
		mat_4x4x_1_view_corr = view_1.matrix;

		//mat_4xp12x_1

		view_1 = gsl_matrix_submatrix(mat_4xp12x_1, 0, 0, _size_x, _size_x);
		mat_4xp12x_1_view_00 = view_1.matrix;
		
		view_1 = gsl_matrix_submatrix(mat_4xp12x_1, 0, _size_x, _size_x, _size_x);
		mat_4xp12x_1_view_01 = view_1.matrix;
		
		view_1 = gsl_matrix_submatrix(mat_4xp12x_1, _size_x, _size_x, _size_x, _size_x);
		mat_4xp12x_1_view_11 = view_1.matrix;
		
		view_1 = gsl_matrix_submatrix(mat_4xp12x_1, 2 * _size_x, 0, 2 * _size_x + 1, 2 * _size_x);
		mat_4xp12x_1_view_2xp12x = view_1.matrix;
		
		view_1 = gsl_matrix_submatrix(mat_4xp12x_1, 2 * _size_x, 0, 2 * _size_x, 2 * _size_x);
		mat_2xp12x_view_0 = view_1.matrix;
		
        view_2 = gsl_matrix_subrow(mat_4xp12x_1, 4 * _size_x, 0, _size_x);
        mat_2xp12x_view_10 = view_2.vector;
		
        view_2 = gsl_matrix_subrow(mat_4xp12x_1, 4 * _size_x, _size_x, _size_x);
        mat_2xp12x_view_11 = view_2.vector;	
					
		//mat_2tp1t_1
 
		view_1 = gsl_matrix_submatrix(mat_2tp1t_1, 0, 0, _size_x, _size_x);
		mat_2tp1t_1_view_00 = view_1.matrix;
		
		view_1 = gsl_matrix_submatrix(mat_2tp1t_1, 0, _size_x, _size_x, _size_y);
		mat_2tp1t_1_view_01 = view_1.matrix;

		view_1 = gsl_matrix_submatrix(mat_2tp1t_1, _size_x, _size_x, _size_y, _size_y);
		mat_2tp1t_1_view_11 = view_1.matrix;

		view_1 = gsl_matrix_submatrix(mat_2tp1t_1, _size_t, 0, _size_t + 1, _size_t);
		mat_2tp1t_1_view_tp1t = view_1.matrix;
		
		view_1 = gsl_matrix_submatrix(mat_2tp1t_1, _size_t, 0, _size_x, _size_x);
		mat_tp1t_view_xx = view_1.matrix;

        view_2 = gsl_matrix_subrow(mat_2tp1t_1, 2 * _size_t, 0, _size_x);
        mat_tp1t_view_10 = view_2.vector;
		
        view_2 = gsl_matrix_subrow(mat_2tp1t_1, 2 * _size_t, _size_x, _size_y);
        mat_tp1t_view_11 = view_2.vector;	
					
		
		
    return 0;
}

/**@fn int kalman_robust_em :: alloc_data()
 * @brief
 * Cette méthode alloue les données
 *
 */
int kalman_robust_em :: alloc_data()
{
	if (!_log_likelihood)
		_log_likelihood = new double[nb_iter];
	tkalman_esperance_ref(_x0_est, nb_iter, _size_x);
	tkalman_covariance_ref(_p0_est, nb_iter, _size_x, _size_x);
	tkalman_covariance_ref(_f_est, nb_iter, _size_t, _size_t);
	tkalman_covariance_ref(_q_est, nb_iter, _size_t, _size_t);
	if (kalman_robust_em :: check_data())
		return 1;
	return 0;
}

/**@fn void kalman_robust_em :: initialize();
 * @brief
 Cette méthode met tous les attributs à zéro.
**/
void kalman_robust_em :: initialize()
{
	kalman_robust_em  :: initialize_tmp();
	nb_iter  = 0;
	_log_likelihood = NULL;
	_x0_est = NULL;
	_p0_est = NULL;
	_f_est = NULL;
	_q_est = NULL;

}

/**@fn void kalman_robust_em :: initialize_tmp();
 * @brief
 Cette méthode met tous les temporaires à zéros
**/
void kalman_robust_em :: initialize_tmp()
{
	perm_y_1 = NULL;
	mat_4x4x_1 = NULL,
	mat_4xp12x_1 = NULL;
	mat_2tp1t_1 = NULL;
	vect_4x_1 = NULL;
}



/**@fn void kalman_robust_em :: follow(const gsl_vector * const * observations, unsigned int i);
 * @param i : itération de l'EM
 * @brief
 * Cette méthode enregistre les paramètres de l'itération i de l'algorithme EM.
 **/
void kalman_robust_em :: follow(unsigned int i)
{
	_log_likelihood[i] = log_likelihood(&mat_tt_1_yy, &mat_tp1t_view_11, perm_y_1);
	gsl_vector_memcpy(_x0_est[i], _x0);
	gsl_matrix_memcpy(_p0_est[i], _p0);
	gsl_matrix_memcpy(_f_est[i], _f);
	gsl_matrix_memcpy(_q_est[i], _q);
}

/**@fn void kalman_robust_em :: estimate_parameters(const gsl_vector * const * observations,
                                                       unsigned int i);
 * @param[in] observations : observations
 * @param[in] i : itération de l'EM
 * @brief
 * Cette méthode ré-estime les paramètres du filtre de Kalman couple à partir d'un jeu d'observations.
 */
void kalman_robust_em :: estimate_parameters(const gsl_vector * const * observations,
                                              unsigned int i)
{
	tkalman_robust_filter :: smooth_without_equivalents(observations);
	//Calcul des sommes
	kalman_get_sqrt_sum_corr_xx(_p_f,
							    _x_s,
							    _p_s,
								_c_s,
							    _n,
							    f2_xx,
							    q2_xx,
								mat_4x4x_1,
								&mat_4x4x_1_view_00,
								&mat_4x4x_1_view_11,
								&mat_4x4x_1_view_02,
								&mat_4x4x_1_view_20,
								&mat_4x4x_1_view_32,
								&mat_4x4x_1_view_33,
								&mat_4x4x_1_view_corr,
								mat_4xp12x_1,
								&mat_4xp12x_1_view_2xp12x,
								&mat_2xp12x_view_0,
								&mat_2xp12x_view_10,
								&mat_2xp12x_view_11,
								vect_2x_1,
								vect_4x_1);
	
	kalman_get_sqrt_sum_corr_xy(_x_s,
							    _p_s,
								observations,
								_n,
								mat_2tp1t_1,
								&mat_2tp1t_1_view_tp1t,
								&mat_tp1t_view_xx,
								&mat_tp1t_view_10,
								&mat_tp1t_view_11,
								vect_t_1);
	
	//Calcul des paramètres optimaux
	kalman_get_arg_max( NULL,
					    NULL,
					    &f_xx,
					    &f_yx,						 
					    &q_xx,
					    &q_yy,
					    &mat_4xp12x_1_view_00,
					    &mat_4xp12x_1_view_01,
					    &mat_4xp12x_1_view_11,
					    &mat_2tp1t_1_view_00,
					    &mat_2tp1t_1_view_01,
					    &mat_2tp1t_1_view_11,
					    _n,
					    NULL,
					    NULL,
					    perm_x_1
				        );
	
	//Calcul des constantes
	tkalman_robust_filter :: compute_constants();

	
	if (_log_likelihood != NULL)
		kalman_robust_em :: follow(i);

}


/**@fn void kalman_robust_em :: do_em_algorithm(const gsl_vector * const * observations);
 * @param[in] observations : observations
 * @brief
 * Cette méthode estime les paramètres du filtre de Kalman couple à partir d'un jeu d'observations.
 */
void kalman_robust_em :: do_em_algorithm(const gsl_vector * const * observations)
{
	for (unsigned int i = 0; i < nb_iter; ++ i)
	{
		kalman_robust_em :: estimate_parameters(observations,
												   i);
	}
}

/**@fn bool kalman_robust_em :: operator!() const
 * @return
 - 0 si l'objet est valide
 - 1 sinon
 * @brief
 Cette méthode teste la validité de chaque attribut.
**/
bool kalman_robust_em ::operator!() const
{
    if (tkalman_robust_filter :: operator!() || kalman_robust_em :: check_tmp() )
        return true;
    if (_log_likelihood)
        return kalman_robust_em :: check_data();
    return true;
}


/**@fn kalman_robust_em  :: ~kalman_robust_em()
 * @brief
 * Destructeur
 *
 */
kalman_robust_em :: ~kalman_robust_em(void)
{
    kalman_robust_em :: free();
    kalman_robust_em :: initialize();
}


/**@fn virtual void kalman_robust_em  :: filter(const gsl_vector * const * observations,
                                                  unsigned int nb_observations) = 0
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette méthode effectue le filtrage des données par le filtre de Kalman Triple non supervisé.
 */
void kalman_robust_em :: filter(const gsl_vector * const * observations,
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
	kalman_robust_em :: do_em_algorithm(observations);
	//Filtrage
	tkalman_robust_filter :: filter_without_equivalents(observations);
	//Restauration des données
    tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
    tkalman_robust_filter :: compute_equivalents_sqrt_p_p_and_sqrt_p_f();

}

/**@fn virtual void kalman_robust_em :: smooth(const gsl_vector * const * observations,
                                                  unsigned int nb_observations)
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette méthode effectue le lissage des données par le filtre de Kalman Triple non supervisé.
 */
void kalman_robust_em :: smooth(const gsl_vector * const * observations,
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
	kalman_robust_em :: do_em_algorithm(observations);
	//Lissage
	tkalman_robust_filter :: smooth_without_equivalents(observations);
	//Restauration des données
    tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
    tkalman_robust_filter :: compute_equivalents_sqrt_p_p_and_sqrt_p_f();
    tkalman_base :: compute_equivalents_x_s(observations);
    tkalman_robust_filter :: compute_equivalents_sqrt_p_s();
}

/**@fn kalman_robust_em :: kalman_robust_em (const gsl_vector * x0,
                                                  const gsl_matrix * sqrt_p0,
                                                  const gsl_matrix * f,
                                                  const gsl_matrix * sqrt_q,
                                                  unsigned int n = 0,
                                                  unsigned int p = 0,
                                                  bool data = false);
 * @param[in] x0 : Espérance de l'état initial
 * @param[in] sqrt_p0 : Matrice de covariance de l'état initial (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
 * @param[in] f : Matrice d'évolution
 * @param[in] sqrt_q : Matrice de covariance (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
 * @param[in] n : Nombre d'observations (0 par défaut)
 * @param[in] p : Nombre d'itérations de l'EM (0 par défaut et dans ce cas, cela équivaut à un filtrage simple)
 * @param[in] data : Booléen (True = suivi de l'EM et stockage des paramètres et de la vraisemblance à chaque itération, False = Pas de suivi de l'EM)
 * @brief
 * Constructeur
*/
kalman_robust_em ::  kalman_robust_em(const gsl_vector * x0,
                                        const gsl_matrix * p0,
                                        const gsl_matrix * f,
                                        const gsl_matrix * q,
                                        unsigned int n,
                                        unsigned int p,
                                        bool data)
: tkalman_robust_filter :: tkalman_robust_filter(x0,
                                                   p0,
                                                   f,
                                                   q,
                                                   n)
{
	kalman_robust_em:: initialize();
	kalman_robust_em :: alloc(p, data);

}

/**@fn virtual int kalman_robust_em  :: setup(const gsl_vector * x0,
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
 * Cette méthode réinitialise l'objet.
 */
int kalman_robust_em :: setup(const gsl_vector * x0,
                               const gsl_matrix * p0,
                               const gsl_matrix * f,
                               const gsl_matrix * q,
                               unsigned int n,
                               unsigned int p,
                               bool data)
{
	if (! tkalman_robust_filter :: setup(x0, p0, f, q, n) )
		return 1;
	else
	{
		kalman_robust_em :: free();
		kalman_robust_em :: initialize();
		if (! kalman_robust_em :: alloc(p, data))
			return 1;
	}


	return 0;
}

