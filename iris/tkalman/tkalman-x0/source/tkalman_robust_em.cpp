#include "tkalman_robust_em.hpp"
#include <iostream>
using namespace std;
/**@fn bool tkalman_robust_em :: check_tmp()
 * @return
 - 0 si les tmp ont été bien alloués
 - 1 sinon
 * @brief
 Cette méthode teste la bonne allocation des variables temporaires utilisées par l'objet.
**/
bool tkalman_robust_em :: check_tmp() const
{
    return (! (mat_4x4x_1 && mat_4tp12t_1 && vect_2t_1 && vect_4x_1 && perm_t_1 && perm_y_1) );
}

/**@fn bool tkalman_robust_em :: check_data()
 * @return
 - 0 si les données de suivi de l'EM ont été bien alloués
 - 1 sinon
 * @brief
 Cette méthode teste la bonne allocation des données de suivi de l'EM.
**/
bool tkalman_robust_em :: check_data() const
{
	return (! (_log_likelihood && _x0_est && _p0_est && _f_est && _q_est) );
}



/**@fn void tkalman_robust_em :: free();
* @brief
 Cette méthode libère la mémoire utilisée par les attributs de l'objet.
**/
void tkalman_robust_em :: free()
{
    tkalman_robust_em :: free_tmp();
    tkalman_robust_em :: free_data();
}

/**@fn void tkalman_robust_em :: free_tmp();
 * @brief
 Cette méthode libère la mémoire utilisée par les temporaires de l'objet.
**/
void tkalman_robust_em :: free_tmp()
{
    if (mat_4x4x_1)
        gsl_matrix_free(mat_4x4x_1);

    if (mat_4tp12t_1)
        gsl_matrix_free(mat_4tp12t_1);

    if (vect_2t_1)
        gsl_vector_free(vect_2t_1);

    if (vect_4x_1)
        gsl_vector_free(vect_4x_1);

    if (perm_t_1)
        gsl_permutation_free(perm_t_1);

    if (perm_y_1)
        gsl_permutation_free(perm_y_1);
}


/**@fn void tkalman_robust_em :: free_data();
 *  @brief
 Cette méthode libère la mémoire utilisée par les données de suivi de l'EM.
**/
void tkalman_robust_em :: free_data()
{
	tkalman_esperance_unref(_x0_est, nb_iter);
	tkalman_covariance_unref(_p0_est, nb_iter);
	tkalman_covariance_unref(_f_est, nb_iter);
	tkalman_covariance_unref(_q_est, nb_iter);
	if (_log_likelihood)
		delete[] _log_likelihood;
}


/**@fn int tkalman_robust_em :: alloc(unsigned int p, bool data);
 * @param p : nombre d'itérations de l'EM
 * @param data : booléen
 * @return
  - 0 si l'allocation des attributs de l'objet s'est bien déroulée
  - 1 sinon
  *@brief
  Cette méthode alloue les attributs de l'objet.
 */
int tkalman_robust_em :: alloc(unsigned int p, bool data)
{
	if (tkalman_robust_em :: alloc_tmp())
		return 1;
	nb_iter = p;
	if (data)
	{
		if (tkalman_robust_em :: alloc_data())
			return 1;
	}
	return 0;
}

/**@fn int tkalman_robust_em :: alloc_tmp();
 * @return
  - 0 si l'allocation des temporaires de l'objet s'est bien déroulée
  - 1 sinon
  *@brief
  Cette méthode alloue les temporaires de l'objet.
 */
int tkalman_robust_em :: alloc_tmp()
{
    gsl_matrix_view view_1;
    gsl_vector_view view_2;
    if (!mat_4x4x_1)
        mat_4x4x_1 = gsl_matrix_alloc(4 * _size_x, 4 * _size_x);

    if (!mat_4tp12t_1)
        mat_4tp12t_1 = gsl_matrix_alloc(4 * _size_t + 1, 2 * _size_t);

    if (!vect_2t_1)
        vect_2t_1 = gsl_vector_alloc(2 * _size_t);

    if (!vect_4x_1)
        vect_4x_1 = gsl_vector_alloc(4 * _size_x);

    if (!perm_t_1)
        perm_t_1 = gsl_permutation_alloc(_size_t);

    if (!perm_y_1)
        perm_y_1 = gsl_permutation_alloc(_size_y);

    if (tkalman_robust_em :: check_tmp())
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

        view_1 = gsl_matrix_submatrix(mat_4x4x_1, 2 * _size_x, 2 * _size_x, _size_x, _size_x);
        mat_4x4x_1_view_22 = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4x4x_1, 2 * _size_x, 3 * _size_x, _size_x, _size_x);
        mat_4x4x_1_view_23 = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4x4x_1, 3 * _size_x, 2 * _size_x, _size_x, _size_x);
        mat_4x4x_1_view_32 = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4x4x_1, 3 * _size_x, 3 * _size_x, _size_x, _size_x);
        mat_4x4x_1_view_33 = view_1.matrix;

        //mat_4tp12t_1
        view_1 = gsl_matrix_submatrix(mat_4tp12t_1, 0, 0, 4 * _size_t, 2 * _size_t);
        mat_4tp12t_1_view_mat_4t2t = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4tp12t_1,
                                      2 * _size_t,
                                      0,
                                      2 * _size_t + 1,
                                      2 * _size_t);
        mat_4tp12t_1_view_10b = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4tp12t_1,
                                      2 * _size_t,
                                      0,
                                      _size_x,
                                      _size_x);
        mat_4tp12t_1_view_10b_view_00 = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4tp12t_1,
                                      2 * _size_t,
                                      _size_t,
                                      _size_x,
                                      _size_x);
        mat_4tp12t_1_view_10b_view_02 = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4tp12t_1,
                                      3 * _size_t,
                                      _size_t,
                                      _size_x,
                                      _size_x);
        mat_4tp12t_1_view_10b_view_22 = view_1.matrix;


        view_2 = gsl_matrix_subrow(mat_4tp12t_1, 4 * _size_t, 0, _size_x);
        mat_4tp12t_1_view_10b_view_30 = view_2.vector;

        view_2 = gsl_matrix_subrow(mat_4tp12t_1, 4 * _size_t, _size_x, _size_y);
        mat_4tp12t_1_view_10b_view_31 = view_2.vector;

        view_2 = gsl_matrix_subrow(mat_4tp12t_1, 4 * _size_t, _size_t, _size_x);
        mat_4tp12t_1_view_10b_view_32 = view_2.vector;

        view_2 = gsl_matrix_subrow(mat_4tp12t_1, 4 * _size_t, _size_t + _size_x, _size_y);
        mat_4tp12t_1_view_10b_view_33 = view_2.vector;

        view_1 = gsl_matrix_submatrix(mat_4tp12t_1,
                                      0,
                                      0,
                                      _size_t,
                                      _size_t);
        mat_4tp12t_1_view_sqrt_sum_view_00 = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4tp12t_1,
                                      0,
                                      _size_t,
                                      _size_t,
                                      _size_t);
        mat_4tp12t_1_view_sqrt_sum_view_01 = view_1.matrix;

        view_1 = gsl_matrix_submatrix(mat_4tp12t_1,
                                      _size_t,
                                      _size_t,
                                      _size_t,
                                      _size_t);
        mat_4tp12t_1_view_sqrt_sum_view_11 = view_1.matrix;

    return 0;
}

/**@fn int tkalman_robust_em :: alloc_data()
 * @brief
 * Cette méthode alloue les données
 *
 */
int tkalman_robust_em :: alloc_data()
{
	if (!_log_likelihood)
		_log_likelihood = new double[nb_iter];
	tkalman_esperance_ref(_x0_est, nb_iter, _size_x);
	tkalman_covariance_ref(_p0_est, nb_iter, _size_x, _size_x);
	tkalman_covariance_ref(_f_est, nb_iter, _size_t, _size_t);
	tkalman_covariance_ref(_q_est, nb_iter, _size_t, _size_t);
	if (tkalman_robust_em :: check_data())
		return 1;
	return 0;
}

/**@fn void tkalman_robust_em :: initialize();
 * @brief
 Cette méthode met tous les attributs à zéro.
**/
void tkalman_robust_em :: initialize()
{
	tkalman_robust_em  :: initialize_tmp();
	nb_iter  = 0;
	_log_likelihood = NULL;
	_x0_est = NULL;
	_p0_est = NULL;
	_f_est = NULL;
	_q_est = NULL;

}

/**@fn void tkalman_robust_em :: initialize_tmp();
 * @brief
 Cette méthode met tous les temporaires à zéros
**/
void tkalman_robust_em :: initialize_tmp()
{
    mat_4x4x_1 = NULL;
    mat_4tp12t_1 = NULL;
    vect_2t_1 = NULL;
    vect_4x_1 = NULL;
    perm_t_1 = NULL;
    perm_y_1 = NULL;
}



/**@fn void tkalman_robust_em :: follow(const gsl_vector * const * observations, unsigned int i);
 * @param i : itération de l'EM
 * @brief
 * Cette méthode enregistre les paramètres de l'itération i de l'algorithme EM.
 **/
void tkalman_robust_em :: follow(unsigned int i)
{
	_log_likelihood[i] = log_likelihood(&mat_tt_1_yy, &mat_4tp12t_1_view_10b_view_33, perm_y_1);
	gsl_vector_memcpy(_x0_est[i], _x0);
	gsl_matrix_memcpy(_p0_est[i], _p0);
	gsl_matrix_memcpy(_f_est[i], _f);
	gsl_matrix_memcpy(_q_est[i], _q);
}

/**@fn void tkalman_robust_em :: estimate_parameters(const gsl_vector * const * observations,
                                                       unsigned int i);
 * @param[in] observations : observations
 * @param[in] i : itération de l'EM
 * @brief
 * Cette méthode ré-estime les paramètres du filtre de Kalman couple à partir d'un jeu d'observations.
 */
void tkalman_robust_em :: estimate_parameters(const gsl_vector * const * observations,
                                              unsigned int i)
{
	tkalman_robust_filter :: smooth_without_equivalents(observations);
	//Calcul des sommes
    tkalman_robust_get_sqrt_sum (_p_f,
                                 _x_s,
                                 _p_s,
                                 _c_s,
                                 observations,
                                 _n,
                                 f2_xx,
                                 q2_xx,
                                 vect_y_zero,
                                 mat_4x4x_1,
                                 &mat_4x4x_1_view_00,
                                 &mat_4x4x_1_view_11,
                                 &mat_4x4x_1_view_02,
                                 &mat_4x4x_1_view_20,
                                 &mat_4x4x_1_view_22,
                                 &mat_4x4x_1_view_23,
                                 &mat_4x4x_1_view_32,
                                 &mat_4x4x_1_view_33,
                                 mat_4tp12t_1,
                                 &mat_4tp12t_1_view_mat_4t2t,
                                 &mat_4tp12t_1_view_10b,
                                 &mat_4tp12t_1_view_10b_view_00,
                                 &mat_4tp12t_1_view_10b_view_02,
                                 &mat_4tp12t_1_view_10b_view_22,
                                 &mat_4tp12t_1_view_10b_view_30,
                                 &mat_4tp12t_1_view_10b_view_31,
                                 &mat_4tp12t_1_view_10b_view_32,
                                 &mat_4tp12t_1_view_10b_view_33,
                                 vect_2t_1,
                                 vect_4x_1 );
     /*
    gsl_matrix * mat = gsl_matrix_alloc( _size_t, _size_t );    
    
    gsl_blas_dgemm(	CblasTrans,
					CblasNoTrans,
					1.0,
					&mat_4tp12t_1_view_sqrt_sum_view_00,
					&mat_4tp12t_1_view_sqrt_sum_view_00,
					0.0,
					mat );
					
	cout << "C00" << endl;
	gsl_matrix_fprintf_(stdout, API_DEFAULT_MATRIX_FORMAT, mat );
	cout << endl;
	
    gsl_blas_dgemm(	CblasTrans,
					CblasNoTrans,
					1.0,
					&mat_4tp12t_1_view_sqrt_sum_view_01,
					&mat_4tp12t_1_view_sqrt_sum_view_00,
					0.0,
					mat );
	
	
	
	
	
	cout << "C10" << endl;
	gsl_matrix_fprintf_(stdout, API_DEFAULT_MATRIX_FORMAT, mat );
	cout << endl;
	
    gsl_blas_dgemm(	CblasTrans,
					CblasNoTrans,
					1.0,
					&mat_4tp12t_1_view_sqrt_sum_view_01,
					&mat_4tp12t_1_view_sqrt_sum_view_01,
					0.0,
					mat );
					
    gsl_blas_dgemm(	CblasTrans,
					CblasNoTrans,
					1.0,
					&mat_4tp12t_1_view_sqrt_sum_view_11,
					&mat_4tp12t_1_view_sqrt_sum_view_11,
					1.0,
					mat );
					
					
	
	cout << "C11" << endl;
	gsl_matrix_fprintf_(stdout, API_DEFAULT_MATRIX_FORMAT, mat );
	cout << endl;
                                 
                                 
    gsl_matrix_free( mat );
                        */
	//Calcul des paramètres optimaux
    tkalman_robust_argmax(NULL,
						  NULL,
						  _f,
						  _q,
						  &mat_4tp12t_1_view_sqrt_sum_view_00,
                          &mat_4tp12t_1_view_sqrt_sum_view_01,
                          &mat_4tp12t_1_view_sqrt_sum_view_11,
                          _n,
						  NULL,
                          NULL,
						  perm_t_1);
	
	//Calcul des constantes
	tkalman_robust_filter :: compute_constants();

	
	if (_log_likelihood != NULL)
		tkalman_robust_em :: follow(i);

}


/**@fn void tkalman_robust_em :: do_em_algorithm(const gsl_vector * const * observations);
 * @param[in] observations : observations
 * @brief
 * Cette méthode estime les paramètres du filtre de Kalman couple à partir d'un jeu d'observations.
 */
void tkalman_robust_em :: do_em_algorithm(const gsl_vector * const * observations)
{
	for (unsigned int i = 0; i < nb_iter; ++ i)
	{
		tkalman_robust_em :: estimate_parameters(observations,
												   i);
	}
}

/**@fn bool tkalman_robust_em :: operator!() const
 * @return
 - 0 si l'objet est valide
 - 1 sinon
 * @brief
 Cette méthode teste la validité de chaque attribut.
**/
bool tkalman_robust_em ::operator!() const
{
    if (tkalman_robust_filter :: operator!() || tkalman_robust_em :: check_tmp() )
        return true;
    if (_log_likelihood)
        return tkalman_robust_em :: check_data();
    return true;
}


/**@fn tkalman_robust_em  :: ~tkalman_robust_em()
 * @brief
 * Destructeur
 *
 */
tkalman_robust_em :: ~tkalman_robust_em(void)
{
    tkalman_robust_em :: free();
    tkalman_robust_em :: initialize();
}


/**@fn virtual void tkalman_robust_em  :: filter(const gsl_vector * const * observations,
                                                  unsigned int nb_observations) = 0
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette méthode effectue le filtrage des données par le filtre de Kalman Triple non supervisé.
 */
void tkalman_robust_em :: filter(const gsl_vector * const * observations,
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
	tkalman_robust_em :: do_em_algorithm(observations);
	//Filtrage
	tkalman_robust_filter :: filter_without_equivalents(observations);
	//Restauration des données
    tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
    tkalman_robust_filter :: compute_equivalents_sqrt_p_p_and_sqrt_p_f();

}

/**@fn virtual void tkalman_robust_em :: smooth(const gsl_vector * const * observations,
                                                  unsigned int nb_observations)
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette méthode effectue le lissage des données par le filtre de Kalman Triple non supervisé.
 */
void tkalman_robust_em :: smooth(const gsl_vector * const * observations,
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
	tkalman_robust_em :: do_em_algorithm(observations);
	//Lissage
	tkalman_robust_filter :: smooth_without_equivalents(observations);
	//Restauration des données
    tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
    tkalman_robust_filter :: compute_equivalents_sqrt_p_p_and_sqrt_p_f();
    tkalman_base :: compute_equivalents_x_s(observations);
    tkalman_robust_filter :: compute_equivalents_sqrt_p_s();
}

/**@fn tkalman_robust_em :: tkalman_robust_em (const gsl_vector * x0,
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
tkalman_robust_em ::  tkalman_robust_em(const gsl_vector * x0,
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
	tkalman_robust_em:: initialize();
	tkalman_robust_em :: alloc(p, data);

}

/**@fn virtual int tkalman_robust_em  :: setup(const gsl_vector * x0,
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
int tkalman_robust_em :: setup(const gsl_vector * x0,
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
		tkalman_robust_em :: free();
		tkalman_robust_em :: initialize();
		if (! tkalman_robust_em :: alloc(p, data))
			return 1;
	}


	return 0;
}

