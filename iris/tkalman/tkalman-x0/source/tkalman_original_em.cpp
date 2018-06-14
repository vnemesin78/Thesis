#include "tkalman_original_em.hpp"
#include <iostream>
using namespace std;
/**@fn tkalman_original_em :: tkalman_original_em(const gsl_vector * x0,
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
tkalman_original_em :: tkalman_original_em(const gsl_vector * x0,
										   const gsl_matrix * p0,
										   const gsl_matrix * f,
										   const gsl_matrix * q,
										   unsigned int n,
										   unsigned int p,
										   bool data)
: tkalman_original_filter :: tkalman_original_filter(x0,
													 p0,
													 f,
													 q,
													 n)
{
	tkalman_original_em:: initialize();
	tkalman_original_em :: alloc(p, data);
}

/**@fn virtual int tkalman_original_em :: setup(const gsl_vector * x0,
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
int tkalman_original_em :: setup(const gsl_vector * x0,
								 const gsl_matrix * p0,
								 const gsl_matrix * f,
								 const gsl_matrix * q,
							     unsigned int n,
							     unsigned int p,
								 bool data)
{
	if (! tkalman_original_filter :: setup(x0, p0, f, q, n) )
		return 1;
	else
	{
		tkalman_original_em :: free();
		tkalman_original_em :: initialize();
		if (! tkalman_original_em :: alloc(p, data))
			return 1;
	}


	return 0;
}

/**@fn virtual void tkalman_original_em :: filter(const gsl_vector * const * observations,
												  unsigned int nb_observations) = 0
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette méthode effectue le filtrage des données par le filtre de Kalman Triple non supervisé.
 */
void tkalman_original_em :: filter(const gsl_vector * const * observations,
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
	//Filtrage
	if (tkalman_original_filter :: check_positivity())
	{
		return;
	}
	tkalman_original_filter :: filter_without_equivalents(observations);
	//Restauration des données
    tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
    tkalman_original_filter :: compute_equivalents_p_p_and_p_f();
}

/**@fn virtual void tkalman_original_em :: smooth(const gsl_vector * const * observations,
												  unsigned int nb_observations)
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette méthode effectue le lissage des données par le filtre de Kalman Triple non supervisé.
 */
void tkalman_original_em :: smooth(const gsl_vector * const * observations,
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
	if (tkalman_original_filter :: check_positivity())
	{
		return;
	}
	tkalman_original_filter :: smooth_without_equivalents(observations);
	//Restauration des données
	tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
	tkalman_original_filter :: compute_equivalents_p_p_and_p_f();
	tkalman_base :: compute_equivalents_x_s(observations);
	tkalman_original_filter :: compute_equivalents_p_s();
}

/**@fn tkalman_original_em :: ~tkalman_original_em()
 * @brief
 * Destructeur
 *
 */
tkalman_original_em :: ~tkalman_original_em(void)
{
	tkalman_original_em :: free();
	tkalman_original_em :: initialize();
}

/**@fn void tkalman_original_em :: do_em_algorithm(const gsl_vector * const * observations);
 * @param[in] observations : observations
 * @brief
 * Cette méthode estime les paramètres du filtre de Kalman couple à partir d'un jeu d'observations.
 */
void tkalman_original_em :: do_em_algorithm(const gsl_vector * const * observations)
{
	for (unsigned int i = 0; i < nb_iter; ++ i)
	{	
		tkalman_original_em :: estimate_parameters(observations,
												   i);
		//if (tkalman_original_filter :: check_positivity())
		//{
			//fprintf(stderr, "Warning : Parameter matrices are not inversible - step %d!\n", i);
		//	break;
		//}
	}
}

/**@fn void tkalman_original_em :: estimate_parameters(const gsl_vector * const * observations,
													   unsigned int nb_observations,
													   unsigned int i);
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @param[in] i : itération de l'EM
 * @brief
 * Cette méthode ré-estime les paramètres du filtre de Kalman couple à partir d'un jeu d'observations.
 */
void tkalman_original_em :: estimate_parameters(const gsl_vector * const * observations,
												unsigned int i)
{
	tkalman_original_filter :: smooth_without_equivalents(observations);
	
	//Calcul des sommes
	tkalman_original_get_sums(c_00,
							  &c_00_xx,
							  c_10,
							  &c_10_xx,
							  c_11,
							  &c_11_xx,
							  _x_s,
							  observations,
							  _p_s,
							  _c_s,
							  _n,
							  vect_t_1,
							  &vect_t_1_x,
							  &vect_t_1_y,
							  vect_t_2,
							  &vect_t_2_x,
							  &vect_t_2_y);
							  
	
	/*cout << "C00" << endl;
	/gsl_matrix_fprintf_(stdout, API_DEFAULT_MATRIX_FORMAT, c_00 );
	cout << endl;
	**/
	/**
	cout << "C10" << endl;
	gsl_matrix_fprintf_(stdout, API_DEFAULT_MATRIX_FORMAT, c_10 );
	cout << endl;
	cout << "C11" << endl;
	gsl_matrix_fprintf_(stdout, API_DEFAULT_MATRIX_FORMAT, c_11 );
	cout << endl;
	**/
	//Calcul des paramètres optimaux
	tkalman_original_argmax(NULL,
							NULL,
							_f,
							_q,
							c_00,
							c_10,
							c_11,
							_n,
							_x_s[0],
							_p_s[0]);

	//Calcul des constantes
	tkalman_original_filter :: compute_constants();
	if (_log_likelihood != NULL)
		tkalman_original_em :: follow(observations, i);
}

/**@fn void tkalman_original_em :: follow(const gsl_vector * const * observations, unsigned int i);
 * @param i : itération de l'EM
 * @brief
 * Cette méthode enregistre les paramètres à l'itération i de l'algorithme EM.
 **/
void tkalman_original_em :: follow(const gsl_vector * const * observations,
								   unsigned int i)
{
	_log_likelihood[i] = log_likelihood(observations, &mat_tt_1_yy, &vect_t_1_y, perm_y_1, mat_yy_1);
	gsl_vector_memcpy(_x0_est[i], _x0);
	gsl_matrix_memcpy(_p0_est[i], _p0);
	gsl_matrix_memcpy(_f_est[i], _f);
	gsl_matrix_memcpy(_q_est[i], _q);
}


/**@fn void tkalman_original_em :: initialize();
 * @brief
 Cette méthode met tous les attributs à zéro.
**/
void tkalman_original_em :: initialize()
{
	tkalman_original_em :: initialize_tmp();
	nb_iter  = 0;
	_log_likelihood = NULL;
	_x0_est = NULL;
	_p0_est = NULL;
	_f_est = NULL;
	_q_est = NULL;
}

/**@fn void tkalman_original_em :: initialize_tmp();
 * @brief
 Cette méthode met tous les temporaires à zéros
**/
void tkalman_original_em :: initialize_tmp()
{
	c_00 = NULL;
	c_10 = NULL;
	c_11 = NULL;
	vect_t_1 = NULL;
	vect_t_2 = NULL;
	perm_y_1 = NULL;
	mat_yy_1 = NULL;
}

/**@fn int tkalman_original_em :: alloc(unsigned int p, bool data);
 * @param p : nombre d'itérations de l'EM
 * @param data : booléen
 * @return
  - 0 si l'allocation des attributs de l'objet s'est bien déroulée
  - 1 sinon
  *@brief
  Cette méthode alloue les attributs de l'objet.
 */
int  tkalman_original_em :: alloc(unsigned int p, bool data)
{
	if (tkalman_original_em :: alloc_tmp())
		return 1;
	nb_iter = p;
	if (data)
	{
		if (tkalman_original_em :: alloc_data())
			return 1;
	}
	return 0;

}
/**@fn int tkalman_original_em :: alloc_data()
 * @brief
 * Cette méthode alloue les données
 *
 */
int tkalman_original_em :: alloc_data()
{
	if (!_log_likelihood)
		_log_likelihood = new double[nb_iter];
	tkalman_esperance_ref(_x0_est, nb_iter, _size_x);
	tkalman_covariance_ref(_p0_est, nb_iter, _size_x, _size_x);
	tkalman_covariance_ref(_f_est, nb_iter, _size_t, _size_t);
	tkalman_covariance_ref(_q_est, nb_iter, _size_t, _size_t);
	if (tkalman_original_em :: check_data())
		return 1;
	return 0;
}
/**@fn int tkalman_original_em :: alloc_tmp();
 * @return
  - 0 si l'allocation des temporaires de l'objet s'est bien déroulée
  - 1 sinon
  *@brief
  Cette méthode alloue les temporaires de l'objet.
 */
int tkalman_original_em :: alloc_tmp()
{
	gsl_matrix_view view;
	gsl_vector_view view2;
	if (!c_00)
		c_00 = gsl_matrix_alloc(_size_t, _size_t);

	if (!c_10)
		c_10 = gsl_matrix_alloc(_size_t, _size_t);

	if (!c_11)
		c_11 = gsl_matrix_alloc(_size_t, _size_t);

	if (!vect_t_1)
		vect_t_1 = gsl_vector_alloc(_size_t);

	if (!vect_t_2)
		vect_t_2 = gsl_vector_alloc(_size_t);

	if (!mat_yy_1)
		mat_yy_1 = gsl_matrix_alloc(_size_y, _size_y);

    if (!perm_y_1)
        perm_y_1 = gsl_permutation_alloc(_size_y);




	if (tkalman_original_em :: check_tmp())
	{
		return 1;
	}

	view = gsl_matrix_submatrix(c_00,
								0,
								0,
								_size_x,
								_size_x);
	c_00_xx = view.matrix;


	view = gsl_matrix_submatrix(c_10,
								0,
								0,
								_size_x,
								_size_x);
	c_10_xx = view.matrix;



	view = gsl_matrix_submatrix(c_11,
								0,
								0,
								_size_x,
								_size_x);
	c_11_xx = view.matrix;

	view2 = gsl_vector_subvector (vect_t_1, 0, _size_x);
	vect_t_1_x = view2.vector;

	view2 = gsl_vector_subvector (vect_t_1, _size_x, _size_y);
	vect_t_1_y = view2.vector;

	view2 = gsl_vector_subvector (vect_t_2, 0, _size_x);
	vect_t_2_x = view2.vector;

	view2 = gsl_vector_subvector (vect_t_2, _size_x, _size_y);
	vect_t_2_y = view2.vector;

	return 0;
}

/**@fn void tkalman_original_em :: free();
 * @brief
 Cette méthode libère la mémoire utilisée par les attributs de l'objet.
**/
void tkalman_original_em :: free()
{
	tkalman_original_em :: free_data();
	tkalman_original_em :: free_tmp();
}

/**@fn void tkalman_original_em :: free_tmp();
 * @brief
 Cette méthode libère la mémoire utilisée par les temporaires de l'objet.
**/
void tkalman_original_em :: free_tmp()
{
	if (c_00)
		gsl_matrix_free(c_00);

	if (c_10)
		gsl_matrix_free(c_10);

	if (c_11)
		gsl_matrix_free(c_11);

	if (vect_t_1)
		gsl_vector_free(vect_t_1);

	if (vect_t_2)
		gsl_vector_free(vect_t_2);

    if (mat_yy_1)
		gsl_matrix_free(mat_yy_1);

    if (perm_y_1)
        gsl_permutation_free(perm_y_1);

}

/**@fn void tkalman_original_em :: free_data();
 *  @brief
 Cette méthode libère la mémoire utilisée par les données de suivi de l'EM.
**/
void tkalman_original_em :: free_data()
{
	tkalman_esperance_unref(_x0_est, nb_iter);
	tkalman_covariance_unref(_p0_est, nb_iter);
	tkalman_covariance_unref(_f_est, nb_iter);
	tkalman_covariance_unref(_q_est, nb_iter);
	if (_log_likelihood)
		delete[] _log_likelihood;
}

/**@fn bool tkalman_original_em :: check_tmp()
 * @return
 - 0 si les tmp ont été bien alloués
 - 1 sinon
 * @brief
 Cette méthode teste la bonne allocation des variables temporaires utilisées par l'objet.
**/
bool tkalman_original_em :: check_tmp() const
{
	return (! (c_00 && c_10 && c_11 && vect_t_1 && vect_t_2 && perm_y_1 && mat_yy_1) );

}

/**@fn bool tkalman_original_em :: check_data()
 * @return
 - 0 si les données de suivi de l'EM ont été bien alloués
 - 1 sinon
 * @brief
 Cette méthode teste la bonne allocation des données de suivi de l'EM.
**/
bool tkalman_original_em :: check_data() const
{
	return (! (_log_likelihood && _x0_est && _p0_est && _f_est && _q_est) );
}

/**@fn bool tkalman_original_em :: operator!()
 * @return
 - 0 si l'objet est valide
 - 1 sinon
 * @brief
 Cette méthode teste la validité de chaque attribut.
**/
bool tkalman_original_em :: operator!() const
{
	if (tkalman_original_filter :: operator!())
		return true;
	else
	{
		if ( tkalman_original_em :: check_tmp() )
			return true;
		else
		{
			if (_log_likelihood != NULL)
			{
				if (! tkalman_original_em :: check_data() )
					return true;
				else
					return false;
			}

		}
	}
	return false;
}














