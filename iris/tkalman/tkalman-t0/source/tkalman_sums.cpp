#include "tkalman_sums.hpp"
/**@fn void tkalman_get_sqrt_cov(const gsl_matrix * sqrt_p_f,
								 const gsl_matrix * sqrt_p_s_,
								 const gsl_matrix * c_s,
								 const gsl_matrix * f2_x_,
								 const gsl_matrix * sqrt_q2_xx,
								 gsl_matrix * mat_3x3x,
								 gsl_matrix * mat_3x3x_view_00,
								 gsl_matrix * mat_3x3x_view_10,
								 gsl_matrix * mat_3x3x_view_11,
								 gsl_matrix * mat_3x3x_view_21,
								 gsl_matrix * mat_3x3x_view_22,
								 gsl_vector * vect_3x);
								 
 * @param sqrt_p_f : racine de la matrice de covariance de l'état filtré (\left[Z_{n|n}\right]^{\frac{1}{2}} \left[F_2^{x} \right]^{\frac{1}{2}})
 * @param sqrt_p_s_ : racine de la matrice de covariance de l'état lissé suivant (\left[P_{n + 1 |N} K_{n|N}^T\right]^{\frac{1}{2}})
 * @param c_s : Matrice de covariance \times terme entre l'état lissé et l'état lissé suivant. (\left[P_{n + 1 |N} K_{n|N}^T\right]^{\frac{1}{2}})
 * @param f2_x_ : matrice d'évolution réduite (F_2^{x})
 * @param sqrt_q2_xx : matrice de bruit de process réduit (\left[Q_2^{x,x}\right]^{\frac{1}{2}})
 * @param mat_3x3x : matrice M allouée de taille (3n_x, 3n_x) ou (2 n_x + n_t, 2 n_x + n_t)
 * @param mat_3x3x_view_00 : 
 * vue sur M partant de (0, 0) à (n_x - 1, n_x - 1)
 * @param mat_3x3x_view_10 : 
 * vue sur M partant de (n_x, 0) à (2 n_x - 1, n_x - 1) ou (de (n_x, 0) à (n_x + n_t - 1, n_x - 1) )
 * @param mat_3x3x_view_11 :
 * vue sur M partant de (n_x, n_x) à (2 n_x - 1, 2 n_x - 1) ou (de (n_x, n_x) à (n_x + n_t - 1, n_x + n_t - 1) )
 * @param mat_3x3x_view_21 :
 * vue sur M partant de (2 n_x, n_x) à (3 n_x - 1, 2 n_x - 1) ou ( de (n_x + n_t, n_x) à (2 n_x + n_t - 1, n_x + n_t - 1)
 * @param mat_3x3x_view_21 :
 * vue sur M partant de (2 n_x, 2 n_x) à (3 n_x - 1, 3 n_x - 1) ou ( de (n_x + n_t, n_x + n_t) à (2 n_x + n_t - 1, 2 n_x + n_t - 1)
 * @param vect_3x : vecteur alloué de taille (3 n_x) ou (2 n_x + n_t)
 * @brief
 * Cette fonction calcule la racine de la matrice de covariance du vecteur aléatoire (z_{n|N}; x_{n + 1 |N}). (z = x ou t). Pour réaliser ce calcul, on réalise la décomposition QR de la matrice suivante : 
 * M = 
 * \begin{pmatrix}
 * \left[Q_2^{x,x}\right]^{\frac{1}{2}}							&	0_{n_x, n_t}	&	0_{n_x, n_x}	\newline
 * \left[Z_{n|n}\right]^{\frac{1}{2}} \left[F_2^{x} \right]^T	&	\left[Z_{n|n}\right]^{\frac{1}{2}}	&	0_{n_t, n_x}	\newline
 * 0_{n_x, n_x}													&	\left[P_{n + 1 |N} K_{n|N}^T\right]^{\frac{1}{2}}	&	\left[P_{n + 1 |N} K_{n|N}^T\right]^{\frac{1}{2}}
 * \end{pmatrix}
 * Elle s'écrit : 
 * M =
 * Q
 * \begin{pmatrix}
 * 	M_{0,0}	&	M_{0,1}	&	M_{0,2} \newline
 * 	O	    	&		M_{1,1}	&	M_{1,2} \newline 
 * 	O   	 	&		O	       &	M_{2,2}
 * \end{pmatrix}
 *  
 * La racine recherchée est la matrice
 * \left[Z_{n|N}\right]^{\frac{1}{2}} = 
 * \begin{pmatrix}
 * 	M_{1,1}	&	M_{1,2} \newline 
 * 	0    	&	M_{2,2}
 * \end{pmatrix}
 *
 * Pour pouvoir utiliser cette fonction, il est nécessaire de construire avant une vue sur la matrice M allant de (n_x + n_z, n_x + n_z) à (2 n_x + n_z - 1, 2 n_x + n_z - 1). (Le résultat est stocké dedans.). 
**/
void tkalman_get_sqrt_cov(const gsl_matrix * sqrt_p_f,
						  const gsl_matrix * sqrt_p_s_,
						  const gsl_matrix * c_s,
						  const gsl_matrix * f2_x_,
						  const gsl_matrix * sqrt_q2_xx,
						  gsl_matrix * mat_3x3x,
						  gsl_matrix * mat_3x3x_view_00,
						  gsl_matrix * mat_3x3x_view_10,
						  gsl_matrix * mat_3x3x_view_11,
						  gsl_matrix * mat_3x3x_view_21,
						  gsl_matrix * mat_3x3x_view_22,
						  gsl_vector * vect_3x)
{
	//Mise à zéro de la matrice M
	gsl_matrix_set_zero(mat_3x3x);
	
	//Construction de la matrice
	//M00
	gsl_matrix_memcpy(mat_3x3x_view_00, 
					  sqrt_q2_xx);
	//M10
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   sqrt_p_f,
				   f2_x_,
				   0.0,
				   mat_3x3x_view_10);
	//M11
	gsl_matrix_memcpy(mat_3x3x_view_11, sqrt_p_f);
	
	//M21
	gsl_matrix_memcpy(mat_3x3x_view_21, c_s);
	
	//M22
	gsl_matrix_memcpy(mat_3x3x_view_22, sqrt_p_s_);
	
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_3x3x,
						 vect_3x);
	gsl_triangle_matrix(mat_3x3x);
}

/**@fn void tkalman_get_sqrt_corr( const gsl_vector * x_s,
								   const gsl_vector * x_s_,
								   const gsl_vector * _y,
								   const gsl_vector * y,
								   const gsl_matrix * sqrt_cov_view_00,
								   const gsl_matrix * sqrt_cov_view_01,
								   const gsl_matrix * sqrt_cov_view_11,
								   gsl_matrix * mat_2tp1_2t,
								   gsl_matrix * mat_2tp1_2t_view_00, // Modif
								   gsl_matrix * mat_2tp1_2t_view_02, 
								   gsl_matrix * mat_2tp1_2t_view_22, // Modif
								   gsl_vector * mat_2tp1_2t_view_30,
								   gsl_vector * mat_2tp1_2t_view_31,
								   gsl_vector * mat_2tp1_2t_view_32,
								   gsl_vector * mat_2tp1_2t_view_33,
								   gsl_vector * vect_2t )
 * @param x_s : Espérance de l'état lissé courant (\hat{x}_{n|N}) 
 * @param x_s_ : Espérance de l'état lissé suivant (\hat{x}_{n + 1|N}) 
 * @param _y : Observation précédente ou espérance de l'observation précédente (\hat{y}_{n - 1|N}) 
 * @param y : observation courante (y_{n})
 * @param sqrt_cov_view_00 : 
 * vue sur la racine de matrice de covariance du vecteur (z_{n|N}; x_{n + 1|N}) partant de (0, 0) à (n_x - 1, n_x - 1) ou (partant de (0,0) à (n_t - 1, n_t - 1) si z = t)
 * @param sqrt_cov_view_01 : 
 * vue sur la racine de matrice de covariance du vecteur (z_{n|N}; x_{n + 1|N}) partant de (0, n_x) à (n_x - 1, 2 n_x - 1) ou (partant de (n_t, 0) à (n_t + n_x - 1 , n_t - 1) si z = t)
 * @param sqrt_cov_view_11 : 
 * vue sur la racine de matrice de covariance du vecteur (z_{n|N}; x_{n + 1|N}) partant de (n_x, n_x) à (2 n_x - 1, 2 n_x - 1) ou (partant de (n_t, n_t) à (n_t + n_x - 1 , n_t + n_x - 1) si z = t)
 * @param mat_2tp1_2t : matrice M de taille (2n_t + 1, 2 n_t) préallouée
 * @param mat_2tp1_2t_view_00 : 
 * vue sur la matrice M partant de (0, 0) à (n_x - 1, n_x - 1) ou ( de (0,0) à (n_t - 1, n_t - 1) si z = t)
 * @param mat_2tp1_2t_view_02 : 
 * vue sur la matrice M partant de (0, n_t) à (n_x - 1, n_t + n_x - 1) ou ( de (0, n_t) à (n_t - 1, n_t + n_x - 1) si z = t)
 * @param mat_2tp1_2t_view_02 :
 * vue sur la matrice M partant de (n_t, n_t) à (n_t + n_x - 1, n_t + n_x - 1)
 * @param mat_2tp1_2t_view_30 : 
 * vue vectorielle sur la matrice M allant de (2 n_t,0) à (2 n_t, n_x - 1)
 * @param mat_2tp1_2t_view_31 : 
 * vue vectorielle sur la matrice M allant de (2 n_t, n_x) à (2n_t, n_t - 1)
 * @param mat_2tp1_2t_view_32 : 
 * vue vectorielle sur la matrice M allant de (2 n_t , n_t) à (2n_t, n_t + n_x - 1)
 * @param mat_2tp1_2t_view_33 : 
 * vue vectorielle sur la matrice M allant de (2 n_t, n_t + n_x) à (2 n_t, 2 n_t - 1)
 * @param vect_2t : vecteur de taille (2n_t) préalloué.
 * @brief
 * Cette fonction calcule la racine de la matrice de corrélation du vecteur (t_{n|N}; t_{n + 1|N}). Pour réaliser ce calcul, elle effectue la décomposition QR de la matrice suivante : 
 * M = 
 * \begin{pmatrix}
 * [\left[cov]^{\frac{1}{2}}]_{0,0}				&	[\left[cov]^{\frac{1}{2}}]_{0,1}	&	0	\newline
 * 0											&	[\left[cov]^{\frac{1}{2}}]_{1,1}	&	0	\newline
 * 0											&	0									&	0	\newline
 * (\hat{x}_{n|N}^T, \hat{y}_{n - 1|N}^T)		&	\hat{x}_{n + 1|N}					&	y_{n}
 * \end{pmatrix}
 * Elle s'écrit :
 * M = 
 * \begin{pmatrix}
 * M_{0,0}	&	M_{0,1}	&	M_{0,2}	&	M_{0,3}	\newline
 * 0		&	M_{1,1}	&	M_{1,2}	&	M_{1,3}	\newline
 * 0		&	0		&	M_{2,2}	&	M_{2,3}	\newline
 * 0		&	0		&	0		&	M_{3,3}	\newline
 * 0		&	0		&	0		&	0		
 * \end{pmatrix}
 * La racine recherchée est donnée par la matrice suivante : 
 * \left[corr(t_n, t_{n+1} |N) \right]^{\frac{1}{2}} = 
 * \begin{pmatrix}
 * M_{0,0}	&	M_{0,1}	&	M_{0,2}	&	M_{0,3}	\newline
 * 0		&	M_{1,1}	&	M_{1,2}	&	M_{1,3}	\newline
 * 0		&	0		&	M_{2,2}	&	M_{2,3}	\newline
 * 0		&	0		&	0		&	M_{3,3}
 * \end{pmatrix}
 * 
 * Pour pouvoir utiliser cette fonction, il est nécessaire de construire une vue sur la matrice M partant de (0,0) à (2 n_t - 1, 2 n_t - 1).
 */
void tkalman_get_sqrt_corr( const gsl_vector * x_s,
						    const gsl_vector * x_s_,
						    const gsl_vector * _y,
						    const gsl_vector * y,
						    const gsl_matrix * sqrt_cov_view_00,
						    const gsl_matrix * sqrt_cov_view_01,
						    const gsl_matrix * sqrt_cov_view_11,
						    gsl_matrix * mat_2tp1_2t,
						    gsl_matrix * mat_2tp1_2t_view_00, // Modif
						    gsl_matrix * mat_2tp1_2t_view_02, 
						    gsl_matrix * mat_2tp1_2t_view_22, // Modif
						    gsl_vector * mat_2tp1_2t_view_30,
						    gsl_vector * mat_2tp1_2t_view_31,
						    gsl_vector * mat_2tp1_2t_view_32,
						    gsl_vector * mat_2tp1_2t_view_33,
						    gsl_vector * vect_2t )
{
	//Mise à zéro de M
	gsl_matrix_set_zero(mat_2tp1_2t);
	
	//Construction de la matrice M
	gsl_matrix_memcpy(mat_2tp1_2t_view_00, sqrt_cov_view_00);
	gsl_matrix_memcpy(mat_2tp1_2t_view_02, sqrt_cov_view_01);
	gsl_matrix_memcpy(mat_2tp1_2t_view_22, sqrt_cov_view_11);
	gsl_vector_memcpy(mat_2tp1_2t_view_30, x_s);
	gsl_vector_memcpy(mat_2tp1_2t_view_31, _y);
	gsl_vector_memcpy(mat_2tp1_2t_view_32, x_s_);
	gsl_vector_memcpy(mat_2tp1_2t_view_33, y);
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_2tp1_2t,
						 vect_2t);
	gsl_triangle_matrix(mat_2tp1_2t);
	
	
	//Suivre les instructions pour obtenir le résultat ! 
}



//Objet de prédiction
/**@fn tkalman_sums :: tkalman_sums(const gsl_matrix * f2_x_,
									const gsl_matrix * sqrt_q2_xx)
* @param f2_x : [F2xx, F2xy]
* @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
* @param q2_xy : Qxy.Qyy
* @brief
* Ce constructeur alloue les variables temp. de la prédiction du filtre de Kalman couple.
*/
tkalman_sums :: tkalman_sums(const gsl_matrix * f2_x_,
							 const gsl_matrix * sqrt_q2_xx)
{
	tkalman_sums :: initialize();
	tkalman_sums :: setup (f2_x_,
						   sqrt_q2_xx);

}
/**@fn int tkalman_sums :: setup(const gsl_matrix * f2_x_,
								 const gsl_matrix * sqrt_q2_xx)
 * @param f2_x : [F2xx, F2xy]
 * @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
 * @param q2_xy : Qxy.Qyy
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction permet de modifier les paramètres de la prédiction (size_x et size_y)
**/
int tkalman_sums :: setup (const gsl_matrix * f2_x_,
							 const gsl_matrix * sqrt_q2_xx)
	
{
	unsigned int size_x = sqrt_q2_xx->size1;
	unsigned int size_y = f2_x_->size2 - size_x;
	
	if (size_x != _size_x || size_y != _size_y)
	{
		tkalman_sums :: free();
		tkalman_sums :: initialize();
		if (tkalman_sums :: set_params(f2_x_, sqrt_q2_xx))
		{
			tkalman_sums :: free();
			tkalman_sums :: initialize();
			return 1;
		}
		
		
		
		
		_size_x = size_x;
		_size_y = size_y;
		_size_t = size_x + size_y;
		if ( tkalman_sums :: alloc() )
		{
			tkalman_sums :: free();
			tkalman_sums :: initialize();
			return 1;
		}
		tkalman_sums :: create_views();
	}	
	return 0;
}

/**@fn tkalman_sums :: ~ tkalman_sums()
 * @brief
 * Destructeur
 */
tkalman_sums :: ~tkalman_sums()
{
	tkalman_sums :: free();
	tkalman_sums :: initialize();
}

/**@fn bool tkalman_sums :: operator!() const;
 * @return
 * - 0 si l'objet est correctement alloué
 * - 1 sinon
 * @brief
 * Check de l'objet.
 */
bool tkalman_sums :: operator!() const
{
	return (! (mat_2xpt_2xpt && mat_4tp1_2t && vect_2t && vect_2xpt && _f2_x_ && _sqrt_q2_xx) );
}

/**@fn void tkalman_sums :: free();
 * @brief
 * Cette fonction désalloue la mémoire utilisée par les variables tmp.
 */
void tkalman_sums :: free()
{
	if (vect_2t)
		gsl_vector_free(vect_2t);
	if (vect_2xpt)
		gsl_vector_free(vect_2xpt);
	if (mat_2xpt_2xpt)
		gsl_matrix_free(mat_2xpt_2xpt);
	if (mat_4tp1_2t)
		gsl_matrix_free(mat_4tp1_2t);
}

/**@fn int tkalman_sums :: alloc();
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction alloue la mémoire utilisée par les variables tmp.
 */
int tkalman_sums :: alloc()
{
	if (!vect_2t)
		vect_2t = gsl_vector_alloc(2 * _size_t);
	if (!vect_2xpt)
		vect_2xpt = gsl_vector_alloc(_size_t + 2 * _size_x);
		
	if (!mat_2xpt_2xpt)
	{
		mat_2xpt_2xpt = gsl_matrix_alloc(2 * _size_x + _size_t, 
										 2 * _size_x + _size_t);
	}
	if (!mat_4tp1_2t)
	{
		mat_4tp1_2t = gsl_matrix_alloc(4 * _size_t + 1, 
									   2 * _size_t);
	}
	return (! (vect_2t && vect_2xpt && mat_2xpt_2xpt && mat_4tp1_2t   )     );
}

/**@fn void tkalman_sums :: create_views();
 * @brief
 * Cette méthode crée les vues.
 */
void tkalman_sums :: create_views()
{
	gsl_matrix_view view;
	gsl_vector_view view42;
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								0,
								0,
								_size_x,
								_size_x);
	mat_2xpt_2xpt_view_00 = view.matrix;
	
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								_size_x,
								0,
								_size_t,
								_size_x);
	mat_2xpt_2xpt_view_10 = view.matrix;
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								_size_x,
								_size_x,
								_size_t,
								_size_t);
	mat_2xpt_2xpt_view_11 = view.matrix;
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								_size_x,
								_size_x + _size_t,
								_size_t,
								_size_x);
	mat_2xpt_2xpt_view_12 = view.matrix;
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								_size_x + _size_t,
								_size_x,
								_size_x,
								_size_t);
	mat_2xpt_2xpt_view_21 = view.matrix;
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								_size_x + _size_t,
								_size_x + _size_t,
								_size_x,
								_size_x);
	mat_2xpt_2xpt_view_22 = view.matrix;
	
	
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								0,
								0,
								3 * _size_x,
								3 * _size_x);
	mat_3x3x = view.matrix;
	
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								0,
								0,
								_size_x,
								_size_x);
	mat_3x3x_view_00 = view.matrix;
	
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								_size_x,
								0,
								_size_x,
								_size_x);
	mat_3x3x_view_10 = view.matrix;
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								_size_x,
								_size_x,
								_size_x,
								_size_x);
	mat_3x3x_view_11 = view.matrix;
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								_size_x,
								_size_x + _size_x,
								_size_x,
								_size_x);
	mat_3x3x_view_12 = view.matrix;
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								_size_x + _size_x,
								_size_x,
								_size_x,
								_size_x);
	mat_3x3x_view_21 = view.matrix;
	view = gsl_matrix_submatrix(mat_2xpt_2xpt,
								_size_x + _size_x,
								_size_x + _size_x,
								_size_x,
								_size_x);
	mat_3x3x_view_22 = view.matrix;
	
	
	
	
	view = gsl_matrix_submatrix(mat_4tp1_2t,
								0,
								0,
								4 * _size_t,
								2 * _size_t);
	mat_4t2t = view.matrix;
	
	
	view = gsl_matrix_submatrix(mat_4tp1_2t,
								0,
								0,
								_size_t,
								_size_t);
	mat_4tp1_2t_view_00 = view.matrix;
	
	view = gsl_matrix_submatrix(mat_4tp1_2t,
								0,
								_size_t,
								_size_t,
								_size_t);
	mat_4tp1_2t_view_01 = view.matrix;
	
	view = gsl_matrix_submatrix(mat_4tp1_2t,
								_size_t,
								_size_t,
								_size_t,
								_size_t);
	mat_4tp1_2t_view_11 = view.matrix;
	
	
	
	
	view = gsl_matrix_submatrix(mat_4tp1_2t,
								2 * _size_t,
								0,
								2 * _size_t + 1,
								2 * _size_t);
	mat_2tp1_2t = view.matrix;
	
	view = gsl_matrix_submatrix(mat_4tp1_2t,
								2 * _size_t,
								0,
								_size_x,
								_size_x);
	mat_2tp1_2t_view_00 = view.matrix;
	
	view = gsl_matrix_submatrix(mat_4tp1_2t,
								2 * _size_t,
								0,
								_size_t,
								_size_t);
	mat_2tp1_2t_view_00_bis = view.matrix;
	
	view = gsl_matrix_submatrix(mat_4tp1_2t,
								2 * _size_t,
								_size_t,
								_size_x,
								_size_x);
	mat_2tp1_2t_view_02 = view.matrix;
	
	view = gsl_matrix_submatrix(mat_4tp1_2t,
								2 * _size_t,
								_size_t,
								_size_t,
								_size_x);
	mat_2tp1_2t_view_02_bis = view.matrix;
	
	view = gsl_matrix_submatrix(mat_4tp1_2t,
								3 * _size_t,
								_size_t,
								_size_x,
								_size_x);
	mat_2tp1_2t_view_22 = view.matrix;
	
	view42 = gsl_matrix_subrow (mat_4tp1_2t, 
								4 * _size_t, 
								0, 
								_size_x);
	mat_2tp1_2t_view_30 = view42.vector;
								
	view42 = gsl_matrix_subrow (mat_4tp1_2t, 
								4 * _size_t, 
								_size_x, 
								_size_y);
	mat_2tp1_2t_view_31 = view42.vector;
								
	view42 = gsl_matrix_subrow (mat_4tp1_2t, 
								4 * _size_t, 
								_size_t, 
								_size_x);
	mat_2tp1_2t_view_32 = view42.vector;
								
	view42 = gsl_matrix_subrow (mat_4tp1_2t, 
								4 * _size_t, 
								_size_t + _size_x, 
								_size_y);
	mat_2tp1_2t_view_33 = view42.vector;	
	

	view42 = gsl_vector_subvector(vect_2xpt, 0, 3 * _size_x);
	vect_3x = view42.vector;
}

		
/**@fn void tkalman_sums :: initialize();
 * @brief
 * Cette fonction met les variables tmp à NULL.
 */
void tkalman_sums :: initialize()
{
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
	mat_2xpt_2xpt = NULL;
	mat_4tp1_2t = NULL;
	vect_2t = NULL;
	vect_2xpt = NULL;
}

/**@fn virtual int tkalman_sums :: set_params(const gsl_matrix * f2_x_,
											  const gsl_matrix * sqrt_q2_xx);
 * @param f2_x : [F2xx, F2xy]
 * @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
 * @param q2_xy : Qxy.Qyy
 * @brief
 * Cette méthode change les paramètres de la prédiction.
 */
int tkalman_sums :: set_params(const gsl_matrix * f2_x_,
							   const gsl_matrix * sqrt_q2_xx)
{

	if (!f2_x_ || !sqrt_q2_xx)
		return 1;
	unsigned int size_x = sqrt_q2_xx->size1;
	unsigned int size_y = f2_x_->size2 - size_x;
	_f2_x_ = f2_x_;
	_sqrt_q2_xx =  sqrt_q2_xx;
	//Vues
	
	gsl_matrix_const_view view1 = gsl_matrix_const_submatrix (_f2_x_, 0, 0, size_x, size_x);
	f2_xx = view1.matrix;
	
	gsl_matrix_const_view view2 = gsl_matrix_const_submatrix (_f2_x_, 0, size_x, size_x, size_y);
	f2_xy = view2.matrix;
	return 0;
}


/**@fn void  tkalman_sums :: compute_sums(const gsl_matrix * const * sqrt_p_f,
										  const gsl_vector * const * x_s,
										  const gsl_matrix * const * sqrt_p_s,
										  const gsl_matrix * const * c_s,
										  const gsl_vector * const * y,
										  unsigned int n,
										  const gsl_vector * const * x_s_0,
										  const gsl_vector * y_m1,
										  const gsl_matrix * sqrt_q_f_0,
										  const gsl_matrix * c_s_0);
 * @param[in] sqrt_p_f : racines des matrices de covariance des états filtrés
 * @param[in] x_s : espérances des états lissés
 * @param[in] sqrt_p_s : racines des matrices de covariance des états lissés
 * @param[in] c_s :  racines des matrices de covariance des états lissés x gain
 * @param[in] y : observations
 * @param[in] n : nombre d'observation
 * @param[in] x_s_0 : espérance de l'état lissé 0.
 * @param[in] y_m1 : espérance de l'observation lissée - 1
 * @param[in] sqrt_q_f_0 : Matrice de covariance de t_{0|0}.
 * @brief
 * Cette méthode calcule la racine de la somme nécessaire à l'algorithme EM. 
 */
void tkalman_sums :: compute_sums(const gsl_matrix * const * sqrt_p_f,
								  const gsl_vector * const * x_s,
								  const gsl_matrix * const * sqrt_p_s,
								  const gsl_matrix * const * c_s,
								  const gsl_vector * const * y,
								  unsigned int n,
								  const gsl_vector * x_s_0,
								  const gsl_vector * y_m1,
								  const gsl_matrix * sqrt_q_f_0,
								  const gsl_matrix * c_s_0)
{
	
	
	//Mise à zéro
	gsl_matrix_set_zero(&mat_4t2t);

	for(unsigned i = n; i > 1; --i)
	{
		tkalman_get_sqrt_cov(sqrt_p_f[i - 1],
							 sqrt_p_s[i],
							 c_s[i - 1],
							 &f2_xx,
							 _sqrt_q2_xx,
							 &mat_3x3x,
							 &mat_3x3x_view_00,
							 &mat_3x3x_view_10,
							 &mat_3x3x_view_11,
							 &mat_3x3x_view_21,
							 &mat_3x3x_view_22,
							 &vect_3x);
		
		
		
		
		tkalman_get_sqrt_corr(x_s[i - 1],
							  x_s[i],
							  y[i - 2],
							  y[i - 1],
						      &mat_3x3x_view_11,
						      &mat_3x3x_view_12,
						      &mat_3x3x_view_22,
						      &mat_2tp1_2t,
						      &mat_2tp1_2t_view_00, // Modif
						      &mat_2tp1_2t_view_02, 
							  &mat_2tp1_2t_view_22, // Modif
						      &mat_2tp1_2t_view_30,
						      &mat_2tp1_2t_view_31,
						      &mat_2tp1_2t_view_32,
						      &mat_2tp1_2t_view_33,
							  vect_2t );
		//Décomposition QR
		gsl_linalg_QR_decomp(&mat_4t2t,
							 vect_2t);
		gsl_triangle_matrix(&mat_4t2t);
	}

	//Cas 0
	tkalman_get_sqrt_cov(sqrt_q_f_0,
						 sqrt_p_s[1],
						 c_s_0,
						 _f2_x_,
						 _sqrt_q2_xx,
						 mat_2xpt_2xpt,
						 &mat_2xpt_2xpt_view_00,
						 &mat_2xpt_2xpt_view_10,
						 &mat_2xpt_2xpt_view_11,
						 &mat_2xpt_2xpt_view_21,
						 &mat_2xpt_2xpt_view_22,
						 vect_2xpt);
	
	tkalman_get_sqrt_corr(x_s_0,
						  x_s[1],
						  y_m1,
						  y[0],
						  &mat_2xpt_2xpt_view_11,
						  &mat_2xpt_2xpt_view_12,
						  &mat_2xpt_2xpt_view_22,
						  &mat_2tp1_2t,
						  &mat_2tp1_2t_view_00_bis,
						  &mat_2tp1_2t_view_02_bis, 
						  &mat_2tp1_2t_view_22,
						  &mat_2tp1_2t_view_30,
						  &mat_2tp1_2t_view_31,
						  &mat_2tp1_2t_view_32,
						  &mat_2tp1_2t_view_33,
						  vect_2t );

		//Décomposition QR
		gsl_linalg_QR_decomp(&mat_4t2t,
							 vect_2t);
		gsl_triangle_matrix(&mat_4t2t);	
}



