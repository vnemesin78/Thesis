#include "tkalman_base.hpp"

/**@fn  tkalman_base :: tkalman_base(const gsl_vector * x0,
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
tkalman_base :: tkalman_base(const gsl_vector * x0,
					 		 const gsl_matrix * p0,
					 		 const gsl_matrix * f,
					 		 const gsl_matrix * q,
							 unsigned int n)
{
	tkalman_base :: initialize();
	tkalman_base :: setup(x0,
				 	      p0,
				  	      f,
				 	      q,
				  		  n);
}

/**@fn int tkalman_base :: setup(const gsl_vector * x0,
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
int tkalman_base :: setup(const gsl_vector * x0,
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
		_size_x = size_x;
		_size_y = size_y;
		_size_t = size_x + size_y;
		_n = n;
		if (tkalman_base :: alloc())
		{
			return 1;
		}
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
   
	return 0;
}

/**@fn void tkalman_base :: get_equivalent(const gsl_matrix * p);
 * @param[in] p : matrice de passage entre les système
 * @brief
 Cette fonction modifie la matrice de passage du filtre de Kalman triplet. Le filtre utilisé pour débruiter sera le filtre original. Cependant la méthode restaura les moments du filtre souhaité à partir du filtre original.
 */
void tkalman_base :: get_equivalent(const gsl_matrix * p)
{
	int signum;
	//Recopie de p
	gsl_matrix_memcpy(_p, p);

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

/**@fn virtual void tkalman_base :: get_equivalent(const gsl_matrix * p,
												   const gsl_matrix * p_inv);
 * @param[in] p : matrice de passage entre les système
 * @param[in] p_inv : inverse de la matrice de passage
 * @brief
 Cette fonction modifie la matrice de passage du filtre de Kalman triplet. Le filtre utilisé pour débruiter sera le filtre original. Cependant la méthode restaura les moments du filtre souhaité à partir du filtre original.
 */
void tkalman_base :: get_equivalent(const gsl_matrix * p,
									const gsl_matrix * p_inv)
{
	//Recopie de p
	gsl_matrix_memcpy(_p, p);
	//Recopie de p
	gsl_matrix_memcpy(_p_inv, p_inv);

	//Correction des erreurs sur p (p_yx = 0, p_yy = Id)
	gsl_matrix_set_identity(&p_yy);
	gsl_matrix_set_zero(&p_yx);

	gsl_matrix_set_identity(&p_inv_yy);
	gsl_matrix_set_zero(&p_inv_yx);


}

/**@fn void tkalman_base :: get_x0(gsl_vector * x0) const;
 * @param x0 : Espérance de l'état initial(Préalloué)
 * @brief
 Cette fonction calcule l'espérance de l'état initial. (P^{x,x}\;\hat{x}_0)
 */
void tkalman_base :: get_x0(gsl_vector * x0) const
{
	//trop dur!
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   &p_xx,
				   _x0,
				   0,
				   x0);
}

/**@fn void tkalman_base :: get_p0(gsl_vector * p0) const;
 * @param p0 : Matrice de covariance de l'état initial (Préallouée)
 * @param mat_xx : matrice temporaire de taille (x.x) (Préallouée)
 * @brief
 Cette fonction calcule la matrice de covariance de l'état initial. (P^{x,x}\;P_0\;[P^{x,x}]^{T})
 */
void tkalman_base :: get_p0(gsl_matrix * p0,
							gsl_matrix * mat_xx) const
{

	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   _p0,
				   &p_xx,
				   0,
				   mat_xx);

	gsl_blas_dgemm(CblasNoTrans,
				   CblasNoTrans,
				   1.0,
				   &p_xx,
				   mat_xx,
				   0,
				   p0);
}

/**@fn void tkalman_base :: get_f(gsl_matrix * f,
								  gsl_matrix * mat_tt) const;
 * @param f : Matrice d'évolution (Préallouée)
 * @param mat_tt : matrice de taille (t.t) (Préallouée)
 * @brief
 Cette fonction calcule la matrice d'évolution. ((P\;F\;P^{-1})
 */
void tkalman_base :: get_f(gsl_matrix * f,
					   	   gsl_matrix * mat_tt) const
{
	gsl_blas_dgemm(CblasNoTrans,
				   CblasNoTrans,
				   1.0,
				   _f,
				   _p_inv,
				   0,
				   mat_tt);

	gsl_blas_dgemm(CblasNoTrans,
				   CblasNoTrans,
				   1.0,
				   _p,
				   mat_tt,
				   0,
				   f);
}

/**@fn void tkalman_base :: get_f(gsl_matrix * f,
								  gsl_matrix * mat_tt) const;
 * @param q : Matrice de covariance du bruit (Préallouée)
 * @param mat_tt : matrice de taille (t.t) (Préallouée)
 * @brief
 Cette fonction calcule la matrice de covariance. ((P\;Q\;P^{T})
 */
void tkalman_base :: get_q(gsl_matrix * q,
		   				   gsl_matrix * mat_tt) const
{
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   _q,
				   _p,
				   0,
				   mat_tt);

	gsl_blas_dgemm(CblasNoTrans,
				   CblasNoTrans,
				   1.0,
				   _p,
				   mat_tt,
				   0,
				   q);
}

/**@fn bool tkalman_base :: operator ! () const;
 * @return
 * - 1 si l'objet est invalide
 * - 0 sinon
 * @brief
 Cette méthode contrôle la validité de l'objet.
 */
bool tkalman_base :: operator ! () const
{
	return !(tkalman_base :: check_params() || tkalman_base :: check_moments() || tkalman_base :: check_tmp() || !vect_y_zero);
}

/**@fn tkalman_base :: ~tkalman_base();
 * @brief
 destructeur de l'objet.
 */
tkalman_base :: ~tkalman_base()
{
	tkalman_base :: free();
	tkalman_base :: initialize();

}

/**@fn void tkalman_base :: initialize();
 * @brief
 Cette méthode met tous les attributs à zéro.
**/
void tkalman_base :: initialize()
{
	tkalman_base :: initialize_params();
	tkalman_base :: initialize_moments();
	tkalman_base :: initialize_tmp();
	vect_y_zero = NULL;

}

/**@fn void tkalman_base :: initialize_params();
 * @brief
 Cette méthode met tous les paramètres à zéros
**/
void tkalman_base :: initialize_params()
{
	_p = NULL;
	_p_inv = NULL;
	_x0 = NULL;
	_p0 = NULL;
	_f = NULL;
	_q = NULL;
	q2_xx = NULL;
	q2_xy = NULL;
	f2_xx = NULL;
	f2_xy = NULL;
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;

}

/**@fn void tkalman_base :: initialize_moments();
 * @brief
 Cette méthode met tous les moments à zéros
**/
void tkalman_base :: initialize_moments()
{
	_n = 0;
	_x_p = NULL;
	_x_f = NULL;
	_x_s = NULL;
	_innovation = NULL;
	_p_p = NULL;
	_p_f = NULL;
	_p_s = NULL;
	_c_s = NULL;
	_s = NULL;



}

/**@fn void tkalman_base :: initialize_tmp();
 * @brief
 Cette méthode met tous les temporaires à zéros
**/
void tkalman_base :: initialize_tmp()
{
	mat_xx_1 = NULL;
	perm_x_1 = NULL;
	vect_x_1 = NULL;
}

/**@fn void tkalman_base :: free();
 * @brief
 Cette méthode libère la mémoire utilisée par les attributs de l'objet.
**/
void tkalman_base :: free()
{
	if (vect_y_zero)
		gsl_vector_free(vect_y_zero);
	 tkalman_base :: free_moments();
	 tkalman_base :: free_params();
	 tkalman_base :: free_tmp();
}

/**@fn void tkalman_base :: free_params();
 * @brief
 Cette  méthode libère la mémoire utilisée par les paramètres de l'objet.
**/
void tkalman_base :: free_params()
{
	if (_p)
		gsl_matrix_free(_p);
	if (_p_inv)
		gsl_matrix_free(_p_inv);
	if (_x0)
		gsl_vector_free(_x0);
	if (_p0)
		gsl_matrix_free(_p0);
	if (_f)
		gsl_matrix_free(_f);
	if (_q)
		gsl_matrix_free(_q);
	if (q2_xx)
		gsl_matrix_free(q2_xx);
	if (q2_xy)
		gsl_matrix_free(q2_xy);
	if (f2_xx)
		gsl_matrix_free(f2_xx);
	if (f2_xy)
		gsl_matrix_free(f2_xy);
}

/**@fn void tkalman_base :: free_moments();
 * @brief
 Cette méthode libère la mémoire utilisée par les moments de l'objet.
**/
void tkalman_base :: free_moments()
{
	if ( _n > 0 )
	{
		unsigned int p = _n + 1;
		tkalman_esperance_unref(_x_p, p)
		tkalman_esperance_unref(_x_f, _n)
		tkalman_esperance_unref(_x_s, p)
		tkalman_esperance_unref(_innovation, _n)
		tkalman_covariance_unref(_c_s, _n)
		tkalman_covariance_unref(_p_p, p)
		tkalman_covariance_unref(_p_f, _n)
		tkalman_covariance_unref(_p_s, p)
		tkalman_covariance_unref(_s, _n)
	}
}

/**@fn void tkalman_base :: free_tmp();
 * @brief
 Cette méthode libère la mémoire utilisée par les temporaires de l'objet.
**/
void tkalman_base :: free_tmp()
{
	if (mat_xx_1)
		gsl_matrix_free(mat_xx_1);
	if (perm_x_1)
		gsl_permutation_free(perm_x_1);
	if (vect_x_1)
		gsl_vector_free(vect_x_1);
}

/**@fn int tkalman_base :: alloc();
 * @return
  - 0 si l'allocation des attributs de l'objet s'est bien déroulée
  - 1 sinon
  *@brief
  Cette méthode alloue les attributs de l'objet.
 */
int tkalman_base :: alloc()
{
	if (!vect_y_zero)
		vect_y_zero = gsl_vector_calloc(_size_y);
	if (!vect_y_zero)
		return 1;

	if (tkalman_base :: alloc_moments() || tkalman_base :: alloc_params() || tkalman_base :: alloc_tmp())
	{
		return 1;
	}
	return 0;
}

/**@fn int tkalman_base :: alloc_params();
 * @return
  - 0 si l'allocation des paramètres de l'objet s'est bien déroulée
  - 1 sinon
  *@brief
  Cette méthode alloue les paramètres de l'objet.
 */
int tkalman_base :: alloc_params()
{
	if (!_p)
		_p = gsl_matrix_calloc(_size_t, _size_t);
	if (!_p_inv)
		_p_inv = gsl_matrix_calloc(_size_t, _size_t);
	if (!_x0)
		_x0 = gsl_vector_calloc(_size_x);
	if (!_p0)
		_p0 = gsl_matrix_calloc(_size_x, _size_x);
	if (!_f)
		_f = gsl_matrix_calloc(_size_t, _size_t);
	if (!_q)
		_q = gsl_matrix_calloc(_size_t, _size_t);


	//Création des vues sur f et q
		gsl_matrix_view view;
		view = gsl_matrix_submatrix(_f, 0, 0, _size_x, _size_x);
		f_xx = view.matrix;
		view = gsl_matrix_submatrix(_f, _size_x, 0, _size_y, _size_x);
		f_yx = view.matrix;
		view = gsl_matrix_submatrix(_f, 0, _size_x, _size_x, _size_y);
		f_xy = view.matrix;
		view = gsl_matrix_submatrix(_f, _size_x, _size_x, _size_y, _size_y);
		f_yy = view.matrix;

		view = gsl_matrix_submatrix(_p, 0, 0, _size_x, _size_x);
		p_xx = view.matrix;
		view = gsl_matrix_submatrix(_p, _size_x, 0, _size_y, _size_x);
		p_yx = view.matrix;
		view = gsl_matrix_submatrix(_p, 0, _size_x, _size_x, _size_y);
		p_xy = view.matrix;
		view = gsl_matrix_submatrix(_p, _size_x, _size_x, _size_y, _size_y);
		p_yy = view.matrix;

		view = gsl_matrix_submatrix(_p, 0, 0, _size_y, _size_x);
		p_ux = view.matrix;
		view = gsl_matrix_submatrix(_p, 0, _size_x, _size_y, _size_y);
		p_uy = view.matrix;

		view = gsl_matrix_submatrix(_p_inv, 0, 0, _size_x, _size_x);
		p_inv_xx = view.matrix;
		view = gsl_matrix_submatrix(_p_inv, _size_x, 0, _size_y, _size_x);
		p_inv_yx = view.matrix;
		view = gsl_matrix_submatrix(_p_inv, 0, _size_x, _size_x, _size_y);
		p_inv_xy = view.matrix;
		view = gsl_matrix_submatrix(_p_inv, _size_x, _size_x, _size_y, _size_y);
		p_inv_yy = view.matrix;

		gsl_matrix_set_identity(_p);
		gsl_matrix_set_identity(_p_inv);










		view = gsl_matrix_submatrix(_q, 0, 0, _size_x, _size_x);
		q_xx = view.matrix;
		view = gsl_matrix_submatrix(_q, _size_x, 0, _size_y, _size_x);
		q_yx = view.matrix;
		view = gsl_matrix_submatrix(_q, 0, _size_x, _size_x, _size_y);
		q_xy = view.matrix;
		view = gsl_matrix_submatrix(_q, _size_x, _size_x, _size_y, _size_y);
		q_yy = view.matrix;

	if (!q2_xx)
		q2_xx = gsl_matrix_calloc(_size_x, _size_x);
	if (!q2_xy)
		q2_xy = gsl_matrix_calloc(_size_x, _size_y);
	if (!f2_xx)
		f2_xx = gsl_matrix_calloc(_size_x, _size_x);
	if (!f2_xy)
		f2_xy = gsl_matrix_calloc(_size_x, _size_y);

	return (tkalman_base :: check_params());
}

/**@fn int tkalman_base :: alloc_moments();
 * @return
  - 0 si l'allocation des moments de l'objet s'est bien déroulée
  - 1 sinon
  *@brief
  Cette méthode alloue les moments de l'objet.
 */
int tkalman_base :: alloc_moments()
{
	if (_n > 0)
	{
		unsigned int p = _n + 1;
		tkalman_esperance_ref(_x_p, p, _size_x)
		tkalman_esperance_ref(_x_f, _n, _size_x)
		tkalman_esperance_ref(_x_s, p, _size_x)
		tkalman_esperance_ref(_innovation, _n, _size_y)

		tkalman_covariance_ref(_p_p, p, _size_x, _size_x)
		tkalman_covariance_ref(_p_f, _n, _size_x, _size_x)
		tkalman_covariance_ref(_p_s, p, _size_x, _size_x)
		tkalman_covariance_ref(_s, _n, _size_y, _size_y)
        tkalman_covariance_ref(_c_s, _n, _size_x, _size_x)
		return (tkalman_base :: check_moments());
	}
	return 0;
}

/**@fn int tkalman_base :: alloc_tmp();
 * @return
  - 0 si l'allocation des temporaires de l'objet s'est bien déroulée
  - 1 sinon
  *@brief
  Cette méthode alloue les temporaires de l'objet.
 */
int tkalman_base :: alloc_tmp()
{
	if (!mat_xx_1)
		mat_xx_1 = gsl_matrix_alloc(_size_x, _size_x);

	if (!perm_x_1)
		perm_x_1 = gsl_permutation_alloc(_size_x);

	if (!vect_x_1)
		vect_x_1 = gsl_vector_alloc(_size_x);

	return (tkalman_base :: check_tmp());
}
/**@fn virtual double tkalman_base :: log_likelihood(gsl_matrix * mat_yy_1,
                                                     gsl_vector * vect_y,
                                                     gsl_permutation * perm_y,
                                                     gsl_matrix * mat_yy_2) const
 * @param mat_yy_1 : matrice temporaire de taille (y.y) (préallouée)
 * @param vect_y : vecteur temporaire de taille (y) (préalloué)
 * @param perm_y : permutation de taille (y) préallouée
 * @param mat_yy_2 : matrice temporaire de taille (y.y) (préallouée)
 * @return
 Valeur du log-vraisemblance.
 */
double tkalman_base :: log_likelihood(const gsl_vector * const * observations,
              						  gsl_matrix * mat_yy_1,
									  gsl_vector * vect_y,
									  gsl_permutation * perm_y,
									  gsl_matrix * mat_yy_2) const
{
	double l = 0;
	for (unsigned int i = 0; i < _n; ++ i)
	{
		//Copie
		gsl_matrix_memcpy (mat_yy_2, _s[i]);

		//Cholesky
		gsl_linalg_cholesky_decomp (mat_yy_2);

		//Calcul du déterminant
		gsl_triangle_matrix(mat_yy_2);
		l -= gsl_linalg_LU_lndet (mat_yy_2);

		//Inversion

		gsl_permutation_init(perm_y);
        gsl_linalg_LU_invert(mat_yy_2, perm_y, mat_yy_1);

		//Recopie de vect_y
		gsl_vector_memcpy(vect_y,
						  _innovation[i]);
		//y - inn
		gsl_vector_sub(vect_y, observations[i]);

		gsl_blas_dtrmv(CblasLower,
					   CblasNoTrans,
					   CblasNonUnit,
					   mat_yy_1,
					   vect_y);
		l -= gsl_blas_dnrm2(vect_y) / 2;
	}
	l -= log(2 * M_PI) * _size_x * _n / 2;
	return l;
}

/**@fn bool tkalman_base :: check_params() const;
 * @return
 * - 0 si les paramètres sont valides
 * - 1 sinon
 */
bool tkalman_base :: check_params() const
{
	return !(_p && _p_inv && _p0 && _f && _q && q2_xx && q2_xy &&  f2_xx && f2_xy);
}

/**@fn bool tkalman_base :: check_moments() const;
 * @return
 * - 0 si les moments sont valides.
 * - 1 sinon
 */
bool tkalman_base :: check_moments() const
{
	if (_x_p && _x_f && _x_s && _innovation && _p_p && _p_f && _p_s && _s && _c_s)
	{
		for (unsigned int i = 0; i < _n ; ++ i)
		{
			if ( !(_x_p[i] && _x_f[i] && _x_s[i] && _innovation[i] && _p_p[i] && _p_f[i] && _p_s[i] && _s[i] && _c_s[i]) )
				return true;
		}
		if ( !(_x_p[_n] && _x_s[_n] && _p_p[_n] && _p_s[_n]) )
			return true;
	}
	else
	{
		return true;
	}
	return false;
}

/**@fn bool tkalman_base :: check_tmp() const;
 * @return
 * - 0 si les tmp sont valides.
 * - 1 sinon
 */
bool tkalman_base :: check_tmp() const
{
	return !(mat_xx_1 && vect_x_1 && perm_x_1);
}

/**@fn void tkalman_base :: compute_equivalents_x_f_and_x_p(const gsl_vector * const * observations);
 * @brief
 * Cette méthode calcule les équivalents avant le lissage.
**/
void tkalman_base :: compute_equivalents_x_f_and_x_p(const gsl_vector * const * observations)
{
	gsl_vector_memcpy(vect_x_1, _x_f[0]);
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   &p_xx,
				   vect_x_1,
				   0.0,
				   _x_f[0]);
	gsl_vector_memcpy(vect_x_1, _x_p[0]);
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   &p_xx,
				   vect_x_1,
				   0.0,
				   _x_p[0]);

	for (unsigned int i = 1; i < _n ; ++ i)
	{
		gsl_vector_memcpy(vect_x_1, _x_f[i]);
		gsl_blas_dgemv(CblasNoTrans,
					   1.0,
					   &p_xx,
					   vect_x_1,
					   0.0,
					   _x_f[i]);
		gsl_blas_dgemv(CblasNoTrans,
					   1.0,
					   &p_xy,
					   observations[i - 1],
					   1.0,
					   _x_f[i]);


		gsl_vector_memcpy(vect_x_1, _x_p[i]);
		gsl_blas_dgemv(CblasNoTrans,
					   1.0,
					   &p_xx,
					   vect_x_1,
					   0.0,
					   _x_p[i]);
		gsl_blas_dgemv(CblasNoTrans,
					   1.0,
					   &p_xy,
					   observations[i - 1],
					   1.0,
					   _x_p[i]);
	}

	gsl_vector_memcpy(vect_x_1, _x_p[_n]);
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   &p_xx,
				   vect_x_1,
				   0.0,
				   _x_p[_n]);
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   &p_xy,
				   observations[_n - 1],
				   1.0,
				   _x_p[_n]);

}

/**@fn void tkalman_base :: compute_equivalents_x_s(const gsl_vector * const * observations);
 * @brief
 * Cette méthode calcule les équivalents après le lissage.
**/
void tkalman_base :: compute_equivalents_x_s(const gsl_vector * const * observations)
{
	gsl_vector_memcpy(vect_x_1, _x_s[0]);
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   &p_xx,
				   vect_x_1,
				   0.0,
				   _x_s[0]);
	for (unsigned int i = 1; i <= _n ; ++ i)
	{
		gsl_vector_memcpy(vect_x_1, _x_s[i]);
		gsl_blas_dgemv(CblasNoTrans,
					   1.0,
					   &p_xx,
					   vect_x_1,
					   0.0,
					   _x_s[i]);
		gsl_blas_dgemv(CblasNoTrans,
					   1.0,
					   &p_xy,
					   observations[i - 1],
					   1.0,
					   _x_s[i]);
	}
}
