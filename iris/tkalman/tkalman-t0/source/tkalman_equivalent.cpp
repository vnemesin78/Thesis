#include "tkalman_equivalent.hpp"


/**@fn void tkalman_compute_equivalent_f(gsl_matrix * f_eq,
										 const gsl_matrix * f,
										 const gsl_matrix * m,
										 const gsl_matrix * m_inv,
										 gsl_matrix * mat_tt)
 * @param f_eq : f équivalent (MFM⁻¹)
 * @param f : f original
 * @param m : matrice de passage
 * @param m_inv : inverse de la matrice de passage
 * @param mat_tt : matrice de taille (n_t, n_t) préallouée
 * @brief
 * Cette fonction calcule le F équivalent par la matrice de passage M. \\
  (MFM⁻¹)
 * 
 */
void tkalman_compute_equivalent_f(gsl_matrix * f_eq,
								  const gsl_matrix * f,
								  const gsl_matrix * m,
								  const gsl_matrix * m_inv,
								  gsl_matrix * mat_tt)
{
	
	gsl_blas_dgemm(CblasNoTrans,
				   CblasNoTrans,
				   1.0,
				   f,
				   m_inv,
				   0.0,
				   mat_tt);
				   
	
	gsl_blas_dgemm(CblasNoTrans,
				   CblasNoTrans,
				   1.0,
				   m,
				   mat_tt,
				   0.0,
				   f_eq);	   
	
}

/**@fn void tkalman_compute_equivalent_sqrt_q(gsl_matrix * sqrt_q_eq,
											  const gsl_matrix * sqrt_q,
											  const gsl_matrix * m,
											  gsl_vector * vect_t)
 * @param sqrt_q_eq : racine de q équivalent R(QM')
 * @param sqrt_q : racine de la matrice de covariance originale
 * @param m : matrice de passage
 * @param vect_t : vecteur de taille (n_t) préalloué
 * @brief
 * Cette fonction calcule la  racine de q équivalent par la matrice de passage M. \\
  R(QM')
 * 
 */
void tkalman_compute_equivalent_sqrt_q(gsl_matrix * sqrt_q_eq,
									   const gsl_matrix * sqrt_q,
									   const gsl_matrix * m,
									   gsl_vector * vect_t)
{
	
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   sqrt_q,
				   m,
				   0.0,
				   sqrt_q_eq);
				   
	gsl_linalg_QR_decomp(sqrt_q_eq,
						 vect_t);

	gsl_triangle_matrix(sqrt_q_eq);
}


/**@fn void tkalman_compute_equivalent_t_0(gsl_vector * t_0_eq,
										   const gsl_vector * t_0,
										   const gsl_matrix * m)
 * @param t_0_eq : espérance du vecteur t initial équivalent (M t0)
 * @param t_0 : espérance du vecteur t initial
 * @param m : matrice de passe M
 * @brief
 * Cette fonction calcule l'espérance du vecteur t équivalent.
*/
void tkalman_compute_equivalent_t_0(gsl_vector * t_0_eq,
									const gsl_vector * t_0,
									const gsl_matrix * m)
{
	
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   m,
				   t_0,
				   0.0,
				   t_0_eq);
}


/**@fn void tkalman_compute_equivalent_expectation(gsl_vector * x,
												   const gsl_vector * _y,
												   const gsl_matrix * m_xx,
												   const gsl_matrix * m_xy,
												   gsl_vector * vect_x)
 * @param x : in <- espérance de l'état original out <- espérance de l'état équivalent
 * @param _y : espérance de l'observation précédente (remplacée par sa valeur si elle est connue)
 * @param m_xx : vue sur la matrice M, (0,0) à (n_x - 1, n_x - 1)
 * @param m_xy : vue sur la matrice M, (0, n_x) à (n_x, n_t)
 * @param vect_x : vecteur de taille (n_x) préalloué.
 * Cette fonction calcule l'espérance de l'état caché équivalent. X_n(eq) = Mxx.X_n + Mxy.Y_{n-1}
 */
void tkalman_compute_equivalent_expectation(gsl_vector * x,
											const gsl_vector * _y,
											const gsl_matrix * m_xx,
											const gsl_matrix * m_xy,
											gsl_vector * vect_x)
{
	gsl_vector_memcpy(vect_x, x);
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   m_xx,
				   vect_x,
				   0.0,
				   x);
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   m_xy,
				   _y,
				   1.0,
				   x);
}

/**@fn void tkalman_compute_equivalent_sqrt_cov_tt(gsl_matrix * sqrt_cov_tt,
												   const gsl_matrix * m,
												   gsl_matrix * mat_tt,
												   gsl_vector * vect_t)
 * @param sqrt_cov_tt : in <- racine de la matrice de covariance du vecteur t, out <- racine de la matrice de covariance du vecteur t équivalent
 * @param m : matrice de passage M
 * @param mat_tt : matrice de taille (n_t, n_t) préallouée
 * @param vect_t : vecteur de taille (n_t) préalloué
 * @brief
 * Cette fonction calcule la racine de la matrice de covariance du vecteur t équivalent.
 */
void tkalman_compute_equivalent_sqrt_cov_tt(gsl_matrix * sqrt_cov_tt,
										    const gsl_matrix * m,
											gsl_matrix * mat_tt,
											gsl_vector * vect_t)
{
	gsl_matrix_memcpy(mat_tt, sqrt_cov_tt);
					  
	
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   mat_tt,
				   m,
				   0.0,
				   sqrt_cov_tt);
				   
	gsl_linalg_QR_decomp(sqrt_cov_tt,
						 vect_t);

	gsl_triangle_matrix(sqrt_cov_tt);
}

/**@fn void tkalman_compute_equivalent_sqrt_p_with_known_y(gsl_matrix * sqrt_p,
														   const gsl_matrix * m_xx,
														   gsl_matrix * mat_xx,
														   gsl_vector * vect_x)
 * @param sqrt_p : in <- racine de la matrice de covariance de x original, out <- racine de la matrice de covariance de x équivalent
 * @param m_xx : vue sur la matrice M, (0,0) à (n_x - 1, n_x - 1)
 * @param vect_x : vecteur de taille (n_x) préalloué.
 * @brief
 * Cette fonction calcule la racine de la matrice de covariance du x équivalent : sqrt_p_eq = sqrt_p Mxx' -> QR...
 * 
 */
void tkalman_compute_equivalent_sqrt_p_with_known_y(gsl_matrix * sqrt_p,
													const gsl_matrix * m_xx,
													gsl_matrix * mat_xx,
													gsl_vector * vect_x)
{
		gsl_matrix_memcpy(mat_xx, sqrt_p);
		gsl_blas_dgemm(CblasNoTrans,
					   CblasTrans,
					   1.0,
					   mat_xx,
					   m_xx,
					   0.0,
					   sqrt_p);
		gsl_linalg_QR_decomp(sqrt_p, vect_x);
		gsl_triangle_matrix(sqrt_p);
	
	
	
	
}


//Objet de prédiction
/**@fn tkalman_equivalents :: tkalman_equivalents(const gsl_matrix * f2_x,
											const gsl_matrix * sqrt_q2_xx,
											const gsl_matrix * q2_xy)
* @param f2_x : [F2xx, F2xy]
* @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
* @param q2_xy : Qxy.Qyy
* @brief
* Ce constructeur alloue les variables temp. de la prédiction du filtre de Kalman couple.
*/
tkalman_equivalents :: tkalman_equivalents (const gsl_matrix * _m,
											unsigned int size_x)
{
	tkalman_equivalents :: initialize();
	tkalman_equivalents :: setup (_m,
							      size_x);

}
/**@fn int tkalman_equivalents :: setup(const gsl_matrix * m;
										unsigned int size_x)
 * @param f2_x : [F2xx, F2xy]
 * @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
 * @param q2_xy : Qxy.Qyy
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction permet de modifier les paramètres de la prédiction (size_x et size_y)
**/
int tkalman_equivalents :: setup (const gsl_matrix * _m,
							      unsigned int size_x)
	
{
	unsigned int size_y = _m->size2 - size_x;
	
	if (size_x != _size_x || size_y != _size_y)
	{
		tkalman_equivalents :: free();
		tkalman_equivalents :: initialize();
		_size_x = size_x;
		_size_y = size_y;
		_size_t = size_x + size_y;
		if ( tkalman_equivalents :: alloc() )
		{
			tkalman_equivalents :: free();
			tkalman_equivalents :: initialize();
			return 1;
		}
		if (tkalman_equivalents :: set_params(_m, size_x))
		{
			tkalman_equivalents :: free();
			tkalman_equivalents :: initialize();
			return 1;
		}
	
		tkalman_equivalents :: create_views();
	}	

	//Calcul de m_inv
	// Mxx^-1		-Mxx Mxy
	// 0			I
	int toto;
	gsl_matrix_set_identity(m_inv);
	gsl_matrix_memcpy(&mat_xx, &m_xx);
	gsl_linalg_LU_decomp (&mat_xx, perm_x, &toto);
	gsl_linalg_LU_invert (&mat_xx, perm_x, &m_inv_xx);
	gsl_blas_dgemm(CblasNoTrans,
				   CblasNoTrans,
				   -1.0,
				   &m_inv_xx,
				   &m_xy,
				   0.0,
				   &m_inv_xy);
	
	return 0;
}

/**@fn tkalman_equivalents :: ~ tkalman_equivalents()
 * @brief
 * Destructeur
 */
tkalman_equivalents :: ~tkalman_equivalents()
{
	tkalman_equivalents :: free();
	tkalman_equivalents :: initialize();
}

/**@fn bool tkalman_equivalents :: operator!() const;
 * @return
 * - 0 si l'objet est correctement alloué
 * - 1 sinon
 * @brief
 * Check de l'objet.
 */
bool tkalman_equivalents :: operator!() const
{
	return (! (vect_t && mat_tt && m_inv && m && perm_x) );
}

/**@fn void tkalman_equivalents :: free();
 * @brief
 * Cette fonction désalloue la mémoire utilisée par les variables tmp.
 */
void tkalman_equivalents :: free()
{
	if (vect_t)
		gsl_vector_free(vect_t);
	if (mat_tt)
		gsl_matrix_free(mat_tt);
	if (m_inv)
		gsl_matrix_free(m_inv);
	if (mat_tt_bis)
		gsl_matrix_free(mat_tt_bis);
	if (perm_x)
		gsl_permutation_free(perm_x);
}

/**@fn int tkalman_equivalents :: alloc();
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction alloue la mémoire utilisée par les variables tmp.
 */
int tkalman_equivalents :: alloc()
{
	if (!vect_t)
		vect_t = gsl_vector_alloc(_size_t);
	if (!mat_tt)
		mat_tt = gsl_matrix_alloc(_size_t, _size_t);
	if (!m_inv)
		m_inv = gsl_matrix_alloc(_size_t, _size_t);
	if (!mat_tt_bis)
		mat_tt_bis = gsl_matrix_alloc(_size_t, _size_t);
	if (!perm_x)
		perm_x = gsl_permutation_alloc(_size_x);
	return (! (vect_t && mat_tt && m_inv && perm_x) );
}

/**@fn void tkalman_equivalents :: create_views();
 * @brief
 * Cette méthode crée les vues.
 */
void tkalman_equivalents :: create_views()
{
	gsl_matrix_view view;
	view = gsl_matrix_submatrix(mat_tt,
								0,
								0,
								_size_x,
								_size_x);
	mat_xx = view.matrix;
	
	view = gsl_matrix_submatrix(m_inv,
								0,
								0,
								_size_x,
								_size_x);
	m_inv_xx = view.matrix;
	
	
	view = gsl_matrix_submatrix(m_inv,
								0,
								_size_x,
								_size_x,
								_size_y);
	m_inv_xy = view.matrix;
	
	gsl_vector_view view42 = gsl_vector_subvector(vect_t, 0, _size_x);
	vect_x = view42.vector;
}

		
/**@fn void tkalman_equivalents :: initialize();
 * @brief
 * Cette fonction met les variables tmp à NULL.
 */
void tkalman_equivalents :: initialize()
{
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
	mat_tt = NULL;
	mat_tt_bis = NULL;
	vect_t = NULL;
	m = NULL;
	m_inv = NULL;
	perm_x = NULL;
}

/**@fn virtual int tkalman_equivalents :: set_params(const gsl_matrix * m;
													 unsigned int size_x)
 * @param f2_x : [F2xx, F2xy]
 * @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
 * @param q2_xy : Qxy.Qyy
 * @brief
 * Cette méthode change les paramètres de la prédiction.
 */
int tkalman_equivalents :: set_params(const gsl_matrix * _m,
									  unsigned int size_x)
{

	if (!_m || !size_x)
		return 1;
	unsigned int size_y = _m->size2 - size_x;
	m = _m;
	//Vues
	{
		gsl_matrix_const_view view = gsl_matrix_const_submatrix (m, 0, 0, size_x, size_x);
		m_xx = view.matrix;
	}
	{
		gsl_matrix_const_view view = gsl_matrix_const_submatrix (m, 0, size_x, size_x, size_y);
		m_xy = view.matrix;
	}
	return 0;
}


