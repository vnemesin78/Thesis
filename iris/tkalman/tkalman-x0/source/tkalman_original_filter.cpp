#include "tkalman_original_filter.hpp"
/**@fn  tkalman_original_filter :: tkalman_original_filter(const gsl_vector * x0,
														 const gsl_matrix * p0,
														 const gsl_matrix * f,
														 const gsl_matrix * q,
														 unsigned int n = 0);
 * @param[in] x0 : Espérance de l'état initial
 * @param[in] p0 : Matrice de covariance de l'état initial (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
 * @param[in] f : Matrice d'évolution
 * @param[in] q : Matrice de covariance (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
 * @param[in] n : Nombre d'observations (0 par défaut)
 * @brief
 constructeur de l'objet.
 */
tkalman_original_filter :: tkalman_original_filter(const gsl_vector * x0,
												   const gsl_matrix * p0,
												   const gsl_matrix * f,
												   const gsl_matrix * q,
												   unsigned int n)
:
tkalman_base(x0,
			 p0,
			 f,
		     q,
			 n)
{
	tkalman_original_filter :: initialize();
	tkalman_original_filter :: alloc();
	//Calcul des constantes
	tkalman_original_filter :: compute_constants();

}

/**@fn virtual int tkalman_original_filter :: setup(const gsl_vector * x0,
												   const gsl_matrix * p0,
												   const gsl_matrix * f,
												   const gsl_matrix * q,
												   unsigned int n = 0);
 * @param[in] x0 : Espérance de l'état initial
 * @param[in] p0 : Matrice de covariance de l'état initial (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
 * @param[in] f : Matrice d'évolution
 * @param[in] q : Matrice de covariance (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
 * @param[in] n : Nombre d'observations (0 par défaut)
 * @return
 * - 0 si l'objet est valide
 * - 1 en cas de problème
 * @brief
 * Cette méthode initialise l'objet.
 */
int tkalman_original_filter :: setup(const gsl_vector * x0,
									const gsl_matrix * p0,
									const gsl_matrix * f,
									const gsl_matrix * q,
									unsigned int n)
{
	unsigned int size_x,
				 size_y;

	//Test des args. d'entrée
	if (!x0 || !p0 || !f || !q)
		return 1;

	size_x = x0->size;
	size_y = f->size1 - size_x;

	//Test si les dim. sont différentes
	if (size_x != _size_x || size_y != _size_y)
	{
		tkalman_base :: free();
		tkalman_base :: initialize();
		tkalman_original_filter :: free();
		tkalman_original_filter :: initialize();
		_size_x = size_x;
		_size_y = size_y;
		_size_t = size_x + size_y;
		_n = n;
		if (tkalman_base :: alloc())
			return 1;
		if (tkalman_original_filter :: alloc())
			return 1;
	}
	else
	{
		if (n != _n)
		{
			_n = n;
			if (tkalman_base :: alloc_moments())
				return 1;
		}
	}

	gsl_vector_memcpy(_x0, x0);
	gsl_matrix_memcpy(_p0, p0);
	gsl_matrix_memcpy(_f, f);
	gsl_matrix_memcpy(_q, q);
	tkalman_original_filter :: compute_constants();

	return 0;
}

/**@fn virtual void tkalman_original_filter :: filter(const gsl_vector * const * observations,
                                                      unsigned int nb_observations)
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette méthode effectue le filtrage des données par le filtre de Kalman Triple.
 */
void tkalman_original_filter :: filter(const gsl_vector * const * observations,
									  unsigned int nb_observations)
{
	if (nb_observations != _n)
	{
		tkalman_base :: free_moments();
		tkalman_base :: initialize_moments();
		_n = nb_observations;
		tkalman_base :: alloc_moments();
	}
	//Passe-avant
        filter_without_equivalents(observations);
	//Restauration des données
		tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
		tkalman_original_filter :: compute_equivalents_p_p_and_p_f();
}

/**@fn virtual void tkalman_original_filter :: smooth(const gsl_vector * const * observations,
										   unsigned int nb_observations)
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette méthode effectue le lissage des données par le filtre de Kalman Triple.
 */
void tkalman_original_filter :: smooth(const gsl_vector * const * observations,
									  unsigned int nb_observations)
{
	if (nb_observations != _n)
	{
		tkalman_base :: free_moments();
		tkalman_base :: initialize_moments();
		_n = nb_observations;
		tkalman_base :: alloc_moments();
	}
	//Lissage
        smooth_without_equivalents(observations);
	//Restauration des données
		tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
		tkalman_original_filter :: compute_equivalents_p_p_and_p_f();

		tkalman_base :: compute_equivalents_x_s(observations);
		tkalman_original_filter :: compute_equivalents_p_s();

}

/**@fn void tkalman_original_filter :: filter_without_equivalents(const gsl_vector * const * observations);
 * @param observations : observations
 * @brief
 * Cette méthode effectue le filtrage avec les paramètres réduits.
 */
void tkalman_original_filter :: filter_without_equivalents(const gsl_vector * const * observations)
{
	//Première prédiction
	gsl_vector_memcpy(_x_p[0],
					  _x0);

	gsl_matrix_memcpy(_p_p[0],
					  _p0);

	//Premier filtrage
	tkalman_original_filtering(_x_f[0],
							   _p_f[0],
							   _innovation[0],
							   _s[0],
							   _x_p[0],
							   _p_p[0],
							   vect_y_zero,
							   observations[0],
							   &f_yx,
							   &f_yy,
							   &q_yy,
							   &mat_tt_1_xy, // Gain
							   &mat_tt_1_yx,
							   &mat_tt_1_yy);

	//Prédiction 1
	tkalman_original_prediction(_x_p[1],
							    _p_p[1],
								_x_f[0],
								_p_f[0],
								vect_y_zero,
								observations[0],
								f2_xx, //F2_xx
								f2_xy, //F2_xy
								q2_xx, //Q2_xx
								q2_xy, //Q2_xy
								mat_xx_1);


	//Boucle

	for (unsigned int i = 1; i < _n; ++ i)
	{

		//Filtrage i
		tkalman_original_filtering(_x_f[i],
								   _p_f[i],
								   _innovation[i],
								   _s[i],
								   _x_p[i],
								   _p_p[i],
								   observations[i - 1],
								   observations[i],
								   &f_yx,
								   &f_yy,
								   &q_yy,
								   &mat_tt_1_xy, // Gain
								   &mat_tt_1_yx,
								   &mat_tt_1_yy);

		//Prédiction i + 1
		tkalman_original_prediction(_x_p[i + 1],
									_p_p[i + 1],
									_x_f[i],
									_p_f[i],
									observations[i - 1],
									observations[i],
									f2_xx, //F2_xx
									f2_xy, //F2_xy
									q2_xx, //Q2_xx
									q2_xy, //Q2_xy
									mat_xx_1);
	}
}


/**@fn void tkalman_original_filter :: smooth_without_equivalents(const gsl_vector * const * observations);
 * @param observations : observations
 * @brief
 * Cette méthode effectue le lissage avec les paramètres réduits.
 */
void tkalman_original_filter :: smooth_without_equivalents(const gsl_vector * const * observations)
{
    //Passe avant
    filter_without_equivalents(observations);
	//Passe arrière
	//Dernier lissage
	gsl_vector_memcpy(_x_s[_n],
					  _x_p[_n]);

	gsl_matrix_memcpy(_p_s[_n],
					  _p_p[_n]);

	for (unsigned int i = _n - 1; i != ((unsigned int) -1) ; -- i)
	{
		tkalman_original_smoothing(_x_s[i],
								   _p_s[i],
								   _c_s[i],
								   _x_f[i],
								   _p_f[i],
								   _x_p[i + 1],
								   _p_p[i + 1],
								   _x_s[i + 1],
								   _p_s[i + 1],
								   f2_xx,
								   mat_xx_1,
								   &mat_tt_1_xx);
	}

}


/**@fn tkalman_original_filter :: ~tkalman_original_filter();
 * @brief
 destructeur de l'objet.
 */
tkalman_original_filter :: ~tkalman_original_filter()
{
	tkalman_original_filter :: free();
	tkalman_original_filter :: initialize();
}

/**@fn void tkalman_original_filter :: compute_equivalents_p_p_and_p_f();
 * @brief
 * Cette méthode calcule les p_p et p_f équivalents.
 */
void tkalman_original_filter :: compute_equivalents_p_p_and_p_f()
{
	for ( unsigned int i = 0; i < _n; ++ i )
	{
		gsl_blas_dgemm(CblasNoTrans,
					   CblasTrans,
					   1.0,
					   _p_p[i],
					   &p_xx,
					   0.0,
					   mat_xx_1);
		gsl_blas_dgemm(CblasNoTrans,
					   CblasNoTrans,
					   1.0,
					   &p_xx,
					   mat_xx_1,
					   0.0,
					   _p_p[i]);

		gsl_blas_dgemm(CblasNoTrans,
					   CblasTrans,
					   1.0,
					   _p_f[i],
					   &p_xx,
					   0.0,
					   mat_xx_1);
		gsl_blas_dgemm(CblasNoTrans,
					   CblasNoTrans,
					   1.0,
					   &p_xx,
					   mat_xx_1,
					   0.0,
					   _p_f[i]);

	}
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   _p_p[_n],
				   &p_xx,
				   0.0,
				   mat_xx_1);
	gsl_blas_dgemm(CblasNoTrans,
				   CblasNoTrans,
				   1.0,
				   &p_xx,
				   mat_xx_1,
				   0.0,
				   _p_p[_n]);
}

/**@fn void tkalman_original_filter :: compute_equivalents_p_s()
 * @brief
 * Cette méthode calcule les p_s équivalents.
 */
void tkalman_original_filter :: compute_equivalents_p_s()
{

	for ( unsigned int i = 0; i <= _n; ++ i )
	{
		gsl_blas_dgemm(CblasNoTrans,
					   CblasTrans,
					   1.0,
					   _p_s[i],
					   &p_xx,
					   0.0,
					   mat_xx_1);
		gsl_blas_dgemm(CblasNoTrans,
					   CblasNoTrans,
					   1.0,
					   &p_xx,
					   mat_xx_1,
					   0.0,
					   _p_s[i]);

	}
}


/**@fn void tkalman_original_filter :: compute_constants()
 * @brief
 * Cette méthode calcule les constantes associées aux paramètres du filtre de Kalman couple.
 *
 */
void tkalman_original_filter :: compute_constants()
{
	tkalman_original_get_constants(f2_xx,
                                   f2_xy,
                                   q2_xx,
                                   q2_xy,
                                   &f_xx,
                                   &f_xy,
                                   &f_yx,
                                   &f_yy,
                                   &q_xx,
                                   &q_xy,
                                   &q_yy,
                                   &mat_tt_1_yy);
}

/**@fn void tkalman_original_filter :: initialize();
 * @brief
 Cette méthode met tous les attributs à zéro.
**/
void tkalman_original_filter :: initialize()
{
	tkalman_original_filter :: initialize_tmp();
}

/**@fn void tkalman_original_filter :: initialize_tmp();
 * @brief
 Cette méthode met tous les temporaires à zéros
**/
void tkalman_original_filter :: initialize_tmp()
{
	mat_tt_1 = NULL;
}

/**@fn int tkalman_original_filter :: alloc();
 * @return
  - 0 si l'allocation des attributs de l'objet s'est bien déroulée
  - 1 sinon
  *@brief
  Cette méthode alloue les attributs de l'objet.
 */
int tkalman_original_filter :: alloc()
{
	return tkalman_original_filter :: alloc_tmp();
}


/**@fn int tkalman_original_filter :: alloc_tmp();
 * @return
  - 0 si l'allocation des temporaires de l'objet s'est bien déroulée
  - 1 sinon
  *@brief
  Cette méthode alloue les temporaires de l'objet.
 */
int tkalman_original_filter :: alloc_tmp()
{
	gsl_matrix_view view;

	if (!mat_tt_1)
	{
		mat_tt_1 = gsl_matrix_alloc(_size_t, _size_t);
	}
	if (!mat_tt_1)
		return 1;


	view = gsl_matrix_submatrix(mat_tt_1,
								0,
								0,
								_size_x,
								_size_x);
	mat_tt_1_xx = view.matrix;

	view = gsl_matrix_submatrix(mat_tt_1,
								_size_x,
								0,
								_size_y,
								_size_x);
	mat_tt_1_yx = view.matrix;

	view = gsl_matrix_submatrix(mat_tt_1,
								_size_x,
								_size_x,
								_size_y,
								_size_y);
	mat_tt_1_yy = view.matrix;

	view = gsl_matrix_submatrix(mat_tt_1,
								0,
								_size_x,
								_size_x,
								_size_y);
	mat_tt_1_xy = view.matrix;

	return 0;

}


/**@fn void tkalman_original_filter :: free();
 * @brief
 Cette méthode libère la mémoire utilisée par les attributs de l'objet.
**/
void tkalman_original_filter :: free()
{
	tkalman_original_filter :: free_tmp();

}

/**@fn void tkalman_original_filter :: free_tmp();
 * @brief
 Cette méthode libère la mémoire utilisée par les temporaires de l'objet.
**/
void tkalman_original_filter :: free_tmp()
{
	if (mat_tt_1)
		gsl_matrix_free(mat_tt_1);
}

/**@fn bool tkalman_original_filter :: check_tmp() const;
 * @return
 * - 0 si les tmps sont bien alloués
 * - 1 sinon
 * @brief
 * Cette fonction controle l'allocation mémoire des tmp.
 */
bool tkalman_original_filter :: check_tmp() const
{
	return (!mat_tt_1);
}

/**@fn virtual bool tkalman_original_filter :: operator!() const;
 * @return
 * - 0 si l'objet est bien initialisé.
 * - 1 sinon
 * @brief
 * Cette fonction controle les allocations mémoires des différents attributs de l'objet.
*/
bool tkalman_original_filter :: operator!() const
{
	return (tkalman_base :: operator !() || tkalman_original_filter :: check_tmp());
}


/**@fn int tkalman_original_filter :: check_positivity() const;
 * @return
 * - 0 si Q et P0 sont positives
 * - 1 sinon.
 */
int tkalman_original_filter :: check_positivity() const
{
	return (!(gsl_matrix_ispos (_q) && gsl_matrix_ispos (_p0)));
}






