#include "tkalman_robust_filter.hpp"
/**@fn void tkalman_robust_filter :: get_p0(gsl_vector * p0) const;
 * @param p0 : Matrice de covariance de l'�tat initial (Pr�allou�e)
 * @param mat_xx : matrice temporaire de taille (x.x) (Pr�allou�e)
 * @brief
 Cette fonction calcule la matrice de covariance de l'�tat initial. (P^{x,x}\;P_0\;[P^{x,x}]^{-1})
 */
void tkalman_robust_filter :: get_p0(gsl_matrix * p0,
                                             gsl_matrix * mat_xx) const
{
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   _p0,
				   &p_xx,
				   0,
				   mat_xx);

	gsl_blas_dgemm(CblasTrans,
				   CblasNoTrans,
				   1.0,
				   mat_xx,
				   mat_xx,
				   0,
				   p0);
}


/**@fn void  tkalman_robust_filter :: get_q(gsl_matrix * f,
                                  gsl_matrix * mat_tt) const;
 * @param q : Matrice de covariance du bruit (Pr�allou�e)
 * @param mat_tt : matrice de taille (t.t) (Pr�allou�e)
 * @brief
 Cette fonction calcule la matrice de covariance. ((P\;Q\;P^{T})
 */
void tkalman_robust_filter :: get_q(gsl_matrix * q,
                                    gsl_matrix * mat_tt) const
{
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   _q,
				   _p,
				   0,
				   mat_tt);

	gsl_blas_dgemm(CblasTrans,
				   CblasNoTrans,
				   1.0,
				   mat_tt,
				   mat_tt,
				   0,
				   q);

}


/**@fn  tkalman_robust_filter :: tkalman_robust_filter(const gsl_vector * x0,
                                                         const gsl_matrix * sqrt_p0,
                                                         const gsl_matrix * f,
                                                         const gsl_matrix * sqrt_q,
                                                         unsigned int n = 0);
 * @param[in] x0 : Esp�rance de l'�tat initial
 * @param[in] sqrt_p0 : racine de la matrice de covariance de l'�tat initial
 * @param[in] f : Matrice d'�volution
 * @param[in] sqrt_q : racine de la matrice de covariance
 * @param[in] n : Nombre d'observations (0 par d�faut)
 * @brief
 constructeur de l'objet.
 */
tkalman_robust_filter :: tkalman_robust_filter(const gsl_vector * x0,
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
	tkalman_robust_filter :: initialize();
	tkalman_robust_filter :: alloc();
	//Calcul des constantes
	tkalman_robust_filter :: compute_constants();

}


/**@fn virtual int tkalman_robust_filter :: setup(const gsl_vector * x0,
                                                   const gsl_matrix * p0,
                                                   const gsl_matrix * f,
                                                   const gsl_matrix * q,
                                                   unsigned int n = 0);
 * @param[in] x0 : Esp�rance de l'�tat initial
 * @param[in] sqrt_p0 : racine de la matrice de covariance de l'�tat initial
 * @param[in] f : Matrice d'�volution
 * @param[in] sqrt_q : racine de la matrice de covariance
 * @param[in] n : Nombre d'observations (0 par d�faut)
 * @return
 * - 0 si l'objet est valide
 * - 1 en cas de probl�me
 * @brief
 * Cette m�thode initialise l'objet.
 */
int tkalman_robust_filter :: setup(const gsl_vector * x0,
                                           const gsl_matrix * p0,
                                           const gsl_matrix * f,
                                           const gsl_matrix * q,
                                           unsigned int n)
{
	unsigned int size_x,
				 size_y;

	//Test des args. d'entr�e
	if (!x0 || !p0 || !f || !q)
		return 1;

	size_x = x0->size;
	size_y = f->size1 - size_x;

	//Test si les dim. sont diff�rentes
	if (size_x != _size_x || size_y != _size_y)
	{
		tkalman_base :: free();
		tkalman_base :: initialize();
		tkalman_robust_filter :: free();
		tkalman_robust_filter :: initialize();
		_size_x = size_x;
		_size_y = size_y;
		_size_t = size_x + size_y;
		_n = n;
		if (tkalman_base :: alloc())
			return 1;
		if (tkalman_robust_filter :: alloc())
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
	tkalman_robust_filter :: compute_constants();
	return 0;

}

/**@fn virtual void tkalman_robust_filter :: filter(const gsl_vector * const * observations,
										   unsigned int nb_observations)
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette m�thode effectue le filtrage des donn�es par le filtre de Kalman Triple.
 */
void tkalman_robust_filter :: filter(const gsl_vector * const * observations,
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
        tkalman_robust_filter :: filter_without_equivalents(observations);
	//Restauration des donn�es
		tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
		tkalman_robust_filter :: compute_equivalents_sqrt_p_p_and_sqrt_p_f();
}

/**@fn virtual void tkalman_robust_filter :: smooth(const gsl_vector * const * observations,
										   unsigned int nb_observations) = 0
 * @param[in] observations : observations
 * @param[in] nb_observations : nombre d'observations
 * @brief
 Cette m�thode effectue le lissage des donn�es par le filtre de Kalman Triple.
 */
void tkalman_robust_filter :: smooth(const gsl_vector * const * observations,
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
        tkalman_robust_filter :: smooth_without_equivalents(observations);
	//Restauration des donn�es
		tkalman_base :: compute_equivalents_x_f_and_x_p(observations);
		tkalman_robust_filter :: compute_equivalents_sqrt_p_p_and_sqrt_p_f();

		tkalman_base :: compute_equivalents_x_s(observations);
		tkalman_robust_filter :: compute_equivalents_sqrt_p_s();

}


/**@fn virtual bool tkalman_original_filter :: operator!() const;
 * @return
 * - 0 si l'objet est bien initialis�.
 * - 1 sinon
 * @brief
 * Cette fonction controle les allocations m�moires des diff�rents attributs de l'objet.
 */
bool tkalman_robust_filter ::  operator!() const
{
    return ( tkalman_base :: operator !() || tkalman_robust_filter :: check_params() || tkalman_robust_filter :: check_tmp() );
}


/**@fn tkalman_robust_filter :: ~tkalman_robust_filter();
 * @brief
 destructeur de l'objet.
 */
tkalman_robust_filter :: ~tkalman_robust_filter()
{
    tkalman_robust_filter :: free();
    tkalman_robust_filter :: initialize();
}





/**@fn void tkalman_robust_filter :: filter_without_equivalents(const gsl_vector * const * observations);
 * @param observations : observations
 * @brief
 * Cette m�thode effectue le filtrage avec les param�tres r�duits.
 */
void tkalman_robust_filter :: filter_without_equivalents(const gsl_vector * const * observations)
{
   //Premi�re pr�diction
        gsl_vector_memcpy(_x_p[0], _x0);
        gsl_matrix_memcpy(_p_p[0], _p0);
        //Premier filtrage
        tkalman_robust_filtering( _x_f[0],
                                  _p_f[0],
                                  _innovation[0],
                                  _s[0],
                                  _x_p[0],
                                  _p_p[0],
                                  observations[0],
                                  vect_y_zero,
                                  &f_yx,
                                  &f_yy,
                                  sqrt_q_yy,
                                  mat_tt_1,
                                  &mat_tt_1_yy,
                                  &mat_tt_1_yx,
                                  &mat_tt_1_xy,
                                  &mat_tt_1_xx,
                                  mat_xy_1,
                                  perm_y_1,
                                  vect_t_1);

        //Pr�diction 2
        tkalman_robust_prediction(_x_p[1],
                                  _p_p[1],
                                  _x_f[0],
                                  _p_f[0],
                                  vect_y_zero,
                                  observations[0],
                                  f2_xx,
                                  f2_xy,
                                  q2_xx,
                                  q2_xy,
                                  &mat_3x2x_1_view_mat_2xx,
                                  &mat_3x2x_1_view_00,
                                  &mat_3x2x_1_view_10,
                                  vect_x_1);

        //Boucle
        for (unsigned i = 1; i < _n; ++ i)
        {
            tkalman_robust_filtering( _x_f[i],
                                      _p_f[i],
                                      _innovation[i],
                                      _s[i],
                                      _x_p[i],
                                      _p_p[i],
                                      observations[i],
                                      observations[i - 1],
                                      &f_yx,
                                      &f_yy,
                                      sqrt_q_yy,
                                      mat_tt_1,
                                      &mat_tt_1_yy,
                                      &mat_tt_1_yx,
                                      &mat_tt_1_xy,
                                      &mat_tt_1_xx,
                                      mat_xy_1,
                                      perm_y_1,
                                      vect_t_1);

            tkalman_robust_prediction(_x_p[i + 1],
                                      _p_p[i + 1],
                                      _x_f[i],
                                      _p_f[i],
                                      observations[i - 1],
                                      observations[i],
                                      f2_xx,
                                      f2_xy,
                                      q2_xx,
                                      q2_xy,
                                      &mat_3x2x_1_view_mat_2xx,
                                      &mat_3x2x_1_view_00,
                                      &mat_3x2x_1_view_10,
                                      vect_x_1);
        }







}


/**@fn void tkalman_robust_filter :: smooth_without_equivalents(const gsl_vector * const * observations);
 * @param observations : observations
 * @brief
 * Cette m�thode effectue le lissage avec les param�tres r�duits.
 */
void tkalman_robust_filter :: smooth_without_equivalents(const gsl_vector * const * observations)
{
    //Passe avant
        tkalman_robust_filter :: filter_without_equivalents(observations);
    //Passe arri�re
        //Dernier lissage
        gsl_vector_memcpy(_x_s[_n], _x_p[_n]);
        gsl_matrix_memcpy(_p_s[_n], _p_p[_n]);

    for (unsigned int i = _n - 1; i != ((unsigned int) -1); --i)
    {
        tkalman_robust_smoothing(_x_s[i],
                                 _p_s[i],
                                 _c_s[i],
                                 _x_f[i],
                                 _p_f[i],
                                 _x_p[i + 1],
                                 _p_p[i + 1],
                                 _x_s[i + 1],
                                 _p_s[i + 1],
                                 f2_xx,
                                 q2_xx,
                                 &mat_tt_1_xx,
                                 mat_3x2x_1,
                                 &mat_3x2x_1_view_00,
                                 &mat_3x2x_1_view_01,
                                 &mat_3x2x_1_view_10,
                                 &mat_3x2x_1_view_11,
                                 &mat_3x2x_1_view_20,
                                 &mat_3x2x_1_view_21,
                                 perm_x_1,
                                 vect_2x_1);
    }
}

/**@fn void tkalman_robust_filter :: compute_equivalents_sqrt_p_p_and_sqrt_p_f();
 * @brief
 * Cette m�thode calcule les p_p et p_f �quivalents.
 */
void tkalman_robust_filter :: compute_equivalents_sqrt_p_p_and_sqrt_p_f()
{
	for ( unsigned int i = 0; i < _n; ++ i )
	{
	    gsl_matrix_memcpy(mat_xx_1, _p_p[i]);
        gsl_blas_dgemm(CblasNoTrans,
                       CblasTrans,
                       1.0,
                       mat_xx_1,
                       &p_xx,
                       0.0,
                       _p_p[i]);
        gsl_linalg_QR_decomp(_p_p[i], vect_x_1);
        gsl_triangle_matrix(_p_p[i]);

	    gsl_matrix_memcpy(mat_xx_1, _p_f[i]);
        gsl_blas_dgemm(CblasNoTrans,
                       CblasTrans,
                       1.0,
                       mat_xx_1,
                       &p_xx,
                       0.0,
                       _p_f[i]);
        gsl_linalg_QR_decomp(_p_f[i], vect_x_1);
        gsl_triangle_matrix(_p_f[i]);
	}

    gsl_matrix_memcpy(mat_xx_1, _p_p[_n]);
    gsl_blas_dgemm(CblasNoTrans,
                   CblasTrans,
                   1.0,
                   mat_xx_1,
                   &p_xx,
                   0.0,
                   _p_p[_n]);
    gsl_linalg_QR_decomp(_p_p[_n], vect_x_1);
    gsl_triangle_matrix(_p_p[_n]);
}

/**@fn void tkalman_robust_filter :: compute_equivalents_p_s()
 * @brief
 * Cette m�thode calcule les p_s �quivalents.
 */
void tkalman_robust_filter :: compute_equivalents_sqrt_p_s()
{

	for ( unsigned int i = 0; i <= _n; ++ i )
	{
	    gsl_matrix_memcpy(mat_xx_1, _p_s[i]);
        gsl_blas_dgemm(CblasNoTrans,
                       CblasTrans,
                       1.0,
                       mat_xx_1,
                       &p_xx,
                       0.0,
                       _p_s[i]);
        gsl_linalg_QR_decomp(_p_s[i], vect_x_1);
        gsl_triangle_matrix(_p_s[i]);

	}
}

/**@fn void tkalman_robust_filter :: compute_constants()
 * @brief
 * Cette m�thode calcule les constantes associ�es aux param�tres du filtre de Kalman couple.
 *
 */
void tkalman_robust_filter :: compute_constants()
{
	tkalman_robust_get_constants ( f2_xx,
                                   f2_xy,
                                   q2_xx,
                                   q2_xy,
                                   sqrt_q_yy,
                                   &f_xx,
                                   &f_xy,
                                   &f_yx,
                                   &f_yy,
                                   &q_xx,
                                   &q_xy,
                                   &q_yy,
                                   mat_tt_1,
                                   &mat_tt_1_yy,
                                   &mat_tt_1_yx,
                                   &mat_tt_1_xy,
                                   &mat_tt_1_xx,
                                   vect_t_1,
                                   perm_y_1);                         
}

/**@fn void tkalman_robust_filter :: initialize();
 * @brief
 Cette m�thode met tous les attributs � z�ro.
**/
void tkalman_robust_filter :: initialize()
{
    tkalman_robust_filter :: initialize_params();
    tkalman_robust_filter :: initialize_tmp();
}

/**@fn void tkalman_robust_filter :: initialize_tmp();
 * @brief
 Cette m�thode met tous les temporaires � z�ros
**/
void tkalman_robust_filter :: initialize_tmp()
{
    mat_tt_1 = NULL;
    mat_xy_1 = NULL;
    mat_3x2x_1 = NULL;
    perm_y_1 = NULL;
    vect_t_1 = NULL;
    vect_2x_1 = NULL;

}

/**@fn void tkalman_robust_filter :: initialize_params()
 * @brief
 Cette m�thode met les param�tres � z�ro.
**/
void tkalman_robust_filter :: initialize_params()
{
    sqrt_q_yy = NULL;
}

/**@fn int tkalman_robust_filter :: alloc();
 * @return
  - 0 si l'allocation des attributs de l'objet s'est bien d�roul�e
  - 1 sinon
  *@brief
  Cette m�thode alloue les attributs de l'objet.
 */
int tkalman_robust_filter :: alloc()
{
    return (tkalman_robust_filter :: alloc_params() || tkalman_robust_filter :: alloc_tmp());
}

/**@fn int tkalman_robust_filter :: alloc_tmp();
 * @return
  - 0 si l'allocation des temporaires de l'objet s'est bien d�roul�e
  - 1 sinon
  *@brief
  Cette m�thode alloue les temporaires de l'objet.
 */
int tkalman_robust_filter :: alloc_tmp()
{
    gsl_matrix_view view;
    //Alloc
    if (! mat_tt_1)
        mat_tt_1 = gsl_matrix_alloc(_size_t, _size_t);
    if (!mat_xy_1)
        mat_xy_1 = gsl_matrix_alloc(_size_x, _size_y);
    if (!mat_3x2x_1)
        mat_3x2x_1 = gsl_matrix_alloc(_size_x * 3, 2 * _size_x);
    if (!perm_y_1)
        perm_y_1 = gsl_permutation_alloc(_size_y);
    if (!vect_t_1)
        vect_t_1 = gsl_vector_alloc(_size_t);
    if (!vect_2x_1)
        vect_2x_1 = gsl_vector_alloc(2 * _size_x);
    //Test de l'allocation
    if (tkalman_robust_filter :: check_tmp())
        return 1;

    //Cr�ation des vues
        view = gsl_matrix_submatrix ( mat_tt_1,
                                      0,
                                      0,
                                      _size_y,
                                      _size_y);
        mat_tt_1_yy = view.matrix;

        view = gsl_matrix_submatrix ( mat_tt_1,
                                      0,
                                      _size_y,
                                      _size_y,
                                      _size_x);
        mat_tt_1_yx = view.matrix;

        view = gsl_matrix_submatrix ( mat_tt_1,
                                      _size_y,
                                      0,
                                      _size_x,
                                      _size_y);
        mat_tt_1_xy = view.matrix;

        view = gsl_matrix_submatrix ( mat_tt_1,
                                      _size_y,
                                      _size_y,
                                      _size_x,
                                      _size_x);
        mat_tt_1_xx = view.matrix;

        view = gsl_matrix_submatrix ( mat_3x2x_1,
                                      0,
                                      0,
                                      _size_x * 2,
                                      _size_x);
        mat_3x2x_1_view_mat_2xx = view.matrix;

        //00
        view = gsl_matrix_submatrix ( mat_3x2x_1,
                                      0,
                                      0,
                                      _size_x,
                                      _size_x);
        mat_3x2x_1_view_00 = view.matrix;

        //01
        view = gsl_matrix_submatrix ( mat_3x2x_1,
                                      0,
                                      _size_x,
                                      _size_x,
                                      _size_x);
        mat_3x2x_1_view_01 = view.matrix;



        //10
        view = gsl_matrix_submatrix ( mat_3x2x_1,
                                      _size_x,
                                      0,
                                      _size_x,
                                      _size_x);
        mat_3x2x_1_view_10 = view.matrix;



        //11
        view = gsl_matrix_submatrix ( mat_3x2x_1,
                                      _size_x,
                                      _size_x,
                                      _size_x,
                                      _size_x);
        mat_3x2x_1_view_11 = view.matrix;


        //20
        view = gsl_matrix_submatrix ( mat_3x2x_1,
                                      2 * _size_x,
                                      0,
                                      _size_x,
                                      _size_x);
        mat_3x2x_1_view_20 = view.matrix;



        //21
        view = gsl_matrix_submatrix ( mat_3x2x_1,
                                      2 * _size_x,
                                      _size_x,
                                      _size_x,
                                      _size_x);
        mat_3x2x_1_view_21 = view.matrix;

    return 0;
}

/**@fn int tkalman_robust_filter :: alloc_params()
 * @brief
 Cette m�thode alloue les attributs stockant les param�tres du filtre de Kalman triplet.
**/
int tkalman_robust_filter :: alloc_params()
{
    if (! sqrt_q_yy)
        sqrt_q_yy = gsl_matrix_calloc(_size_y, _size_y);

    if (tkalman_robust_filter :: check_params())
        return 1;
    return 0;
}

/**@fn void tkalman_robust_filter :: free();
 * @brief
 Cette m�thode lib�re la m�moire utilis�e par les attributs de l'objet.
**/
void  tkalman_robust_filter :: free()
{
    tkalman_robust_filter :: free_tmp();
    tkalman_robust_filter :: free_params();
}

/**@fn void tkalman_robust_filter :: free_tmp();
 * @brief
 Cette m�thode lib�re la m�moire utilis�e par les temporaires de l'objet.
**/
void  tkalman_robust_filter :: free_tmp()
{
    if (mat_tt_1)
        gsl_matrix_free(mat_tt_1);
    if (mat_xy_1)
        gsl_matrix_free(mat_xy_1);
    if (mat_3x2x_1)
        gsl_matrix_free(mat_3x2x_1);
    if (perm_y_1)
        gsl_permutation_free(perm_y_1);
    if (vect_t_1)
        gsl_vector_free(vect_t_1);
    if (vect_2x_1)
        gsl_vector_free(vect_2x_1);
}

/**@fn void tkalman_robust_filter :: free_params()
 * @brief
 Cette m�thode d�salloue les param�tres
**/
void  tkalman_robust_filter :: free_params()
{
    if (sqrt_q_yy)
        gsl_matrix_free(sqrt_q_yy);
}

/**@fn bool tkalman_robust_filter :: check_tmp() const;
 * @return
 * - 0 si les tmps sont bien allou�s
 * - 1 sinon
 * @brief
 * Cette fonction controle l'allocation m�moire des tmp.
 */
bool tkalman_robust_filter :: check_tmp() const
{
    return ( !(mat_tt_1 && mat_xy_1 && mat_3x2x_1 && perm_y_1 && vect_t_1 && vect_2x_1) );
}


/**@fn bool tkalman_robust_filter :: check_params() const;
 * @return
 * - 0 si les tmps sont bien allou�s
 * - 1 sinon
 * @brief
 * Cette fonction controle l'allocation m�moire des param�tres.
 */
bool tkalman_robust_filter :: check_params() const
{
    return ( !(sqrt_q_yy));
}


/**@fn virtual double tkalman_robust_filter :: log_likelihood(gsl_matrix * mat_yy_1,
															  gsl_vector * vect_y,
															  gsl_permutation * perm_y) const
 * @param mat_yy_1 : matrice temporaire de taille (y.y) (pr�allou�e)
 * @param vect_y : vecteur temporaire de taille (y) (pr�allou�)
 * @param perm_y : permutation de taille (y) pr�allou�e
 * @param mat_yy_2 : matrice temporaire de taille (y.y) (pr�allou�e)
 * @return
 Valeur du log-vraisemblance.
 */
double tkalman_robust_filter :: log_likelihood(gsl_matrix * mat_yy_1,
                                               gsl_vector * vect_y_1,
                                               gsl_permutation * perm_y) const
{
	double l = 0;
	for (unsigned int i = 0; i < _n; ++ i)
	{

		//Calcul du d�terminant
		l -= gsl_linalg_LU_lndet (_s[i]);

		//Inversion
		gsl_permutation_init(perm_y);
        gsl_linalg_LU_invert(_s[i], perm_y, mat_yy_1);


						  
		//y - inn
		gsl_blas_dgemv (CblasTrans, 
						1.0, 
						mat_yy_1, 
						_innovation[i], 
						0.0, 
						vect_y_1);
		l -= gsl_blas_dnrm2(vect_y_1);
	}
	l -= (log(2.0 * M_PI) * (_size_y * _n)) / 2.0;
	return l;
}



