/**@file tkalman_prediction.cpp
 * @author Valérian Némesin
 * @date 10/11/2011
 * @brief
 * Ce fichier contient le code source de la prédiction du filtre de Kalman triplet.
**/
#include "tkalman_prediction.hpp"




//Fonctions de prédiction
/**@fn void tkalman_get_x_p(gsl_vector * x_p,
							const gsl_vector * _x_f,
							const gsl_vector * __y,
							const gsl_vector * _y,
							const gsl_matrix * f2_xx,
							const gsl_matrix * f2_xy,
							const gsl_matrix * q2_xy)
 * @param x_p : espérance de l'état prédit
 * @param[in] _x_f : espérance de l'état filtré précédent
 * @param[in] __y : espérance de l'observation (n - 2)
 * @param[in] _y : espérance de l'observation précédente
 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
 * @param[in] f2_xy : F_2^{x,y} = F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y}
 * @param[in] q2_xy : Q_2^{x,y} = Q^{x,y} \: [Q^{y,y}]^{-1}
 * @brief
 Cette fonction calcule l'espérance de l'état prédit à partir de la formule : 
 * \hat{x_{n|n-1}} = F_2^{x,x} \: \hat{x_{n - 1 |n-1}} + Q_2^{x,y} \; y_{n - 2} + F_2^{x,y} \: y_{n - 1}
 */
void tkalman_get_x_p(gsl_vector * x_p,
					 const gsl_vector * _x_f,
					 const gsl_vector * __y,
					 const gsl_vector * _y,
					 const gsl_matrix * f2_xx,
					 const gsl_matrix * f2_xy,
					 const gsl_matrix * q2_xy)
{
	//F2xx . x_f
	gsl_blas_dgemv (CblasNoTrans, 
					1.0, 
					f2_xx, 
					_x_f, 
					0, 
					x_p);
	//F2xx . x_f + Q2xy . _y 
	gsl_blas_dgemv (CblasNoTrans, 
					1.0, 
					q2_xy, 
					_y, 
					1.0, 
					x_p);
					
	//F2xx . x_f + Q2xy . _y  + F2xy . __y
	gsl_blas_dgemv (CblasNoTrans, 
					1.0, 
					f2_xy, 
					__y, 
					1.0, 
					x_p);	
					
}

/**@fn void tkalman_get_sqrt_p_p_1(gsl_matrix * sqrt_p_p_1,
								   const gsl_matrix * sqrt_q_0_f,
								   const gsl_matrix * f2_x,
								   const gsl_matrix * sqrt_q2_xx,
								   gsl_matrix * mat_xpt_x,
								   gsl_matrix * mat_xpt_x_view_00,
								   gsl_matrix * mat_xpt_x_view_10,
								   gsl_vector * vect_x)
 * @param sqrt_p_p_1 : racine de la matrice de covariance de l'état prédit 1
 * @param sqrt_q_0_f : racine de la matrice de covariance du t_0 filtré.
 * @param f2_x : matrice dérivée de F, F2 = [F_2^{x,x}, F_2^{x,y}]
 * @param sqrt_q2_xx : racine de la matrice de covariance du bruit de process réduit Q_2^{x,x} = Q^{x,x} - Q^{x,y}[Q^{y,y}]^{-1}Q^{y,x}
 * @param mat_xpt_x : matrice de taille (x+t, x) allouée
 * @param mat_xpt_x_view_00 : vue sur la matrice mat_xpt_x (de (0,0) à (n_x - 1, n_x - 1))
 * @param mat_xpt_x_view_10 : vue sur la matrice mat_xpt_x (de (n_x,0) à (n_x + n_t - 1, n_x - 1))
 * @param vect_x : vecteur de taille (n_x) alloué
 * @brief
 * Cette fonction calcule la racine de la matrice de covariance de l'état prédit 1.
 * 
 **/
void tkalman_get_sqrt_p_p_1(gsl_matrix * sqrt_p_p_1,
						    const gsl_matrix * sqrt_q_f_0,
						    const gsl_matrix * f2_x,
						    const gsl_matrix * sqrt_q2_xx,
						    gsl_matrix * mat_xpt_x,
							gsl_matrix * mat_xpt_x_view_00,
						    gsl_matrix * mat_xpt_x_view_10,
							gsl_vector * vect_x)
{
	//Construction de la matrice
	// Q2xx
	//
	// sqrt_q_f . F2_x^T
	
	gsl_matrix_memcpy(mat_xpt_x_view_00, sqrt_q2_xx);
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   sqrt_q_f_0,
				   f2_x,
				   0.0,
				   mat_xpt_x_view_10);
				   
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_xpt_x,
						 vect_x);
	gsl_triangle_matrix(mat_xpt_x);
	
	//Recopie du résultat
	gsl_matrix_memcpy(sqrt_p_p_1, mat_xpt_x_view_00);
}


/**@fn void tkalman_get_sqrt_p_p(gsl_matrix * sqrt_p_p,
							     const gsl_matrix * _sqrt_p_f,
							     const gsl_matrix * f2_xx,
							     const gsl_matrix * sqrt_q2_xx,
							     gsl_matrix * mat_2xx,
							     gsl_matrix * mat_2xx_view_00,
							     gsl_matrix * mat_2xx_view_10,
							     gsl_vector * vect_x)
 * @param sqrt_p_p : racine de la matrice de covariance de l'état prédit courant
 * @param[in] _sqrt_p_f : racine de la matrice de covariance de l'état filtré précédent
 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
 * @param[in] sqrt_q2_xx : racine de la matrice de covariance du bruit de process réduit
 * @param mat_2xx : Matrice de taille (2x.x) préallouée.
 * @param mat_2xx_view_00 : vue sur la matrice mat_2xx allant de (0,0) à (n_x, n_x)
 * @param mat_2xx_view_10 : vue sur la matrice mat_2xx allant de (n_x, 0) à (2n_x, n_x)
 * @param vect_x : vecteur de taille (x) préalloué
 * @brief
 Cette fonction calcule la racine de la matrice de covariance de l'état prédit courant. 
 */
void tkalman_get_sqrt_p_p(gsl_matrix * sqrt_p_p,
						  const gsl_matrix * _sqrt_p_f,
						  const gsl_matrix * f2_xx,
						  const gsl_matrix * sqrt_q2_xx,
						  gsl_matrix * mat_2xx,
						  gsl_matrix * mat_2xx_view_00,
						  gsl_matrix * mat_2xx_view_10,
						  gsl_vector * vect_x)
{
	//Construction de la matrice
	// Q2xx
	//
	// sqrt_p_f . F2_xx^T
	
	gsl_matrix_memcpy(mat_2xx_view_00, sqrt_q2_xx);
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   _sqrt_p_f,
				   f2_xx,
				   0.0,
				   mat_2xx_view_10);
				   
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_2xx,
						 vect_x);
	gsl_triangle_matrix(mat_2xx);
	
	//Recopie du résultat
	gsl_matrix_memcpy(sqrt_p_p, mat_2xx_view_00);
	
}


/**@fn void tkalman_do_prediction ( gsl_vector * x_p,
								 gsl_matrix * sqrt_p_p,
								 const gsl_vector * _x_f,
								 const gsl_matrix * _sqrt_p_f,
							  	 const gsl_vector * __y,
							 	 const gsl_vector * _y,
							 	 const gsl_matrix * f2_xx,
							  	 const gsl_matrix * f2_xy,
								 const gsl_matrix * sqrt_q2_xx,
								 const gsl_matrix * q2_xy,
								 gsl_matrix * mat_2xx,
							     gsl_matrix * mat_2xx_view_00,
							     gsl_matrix * mat_2xx_view_10,
							     gsl_vector * vect_x)
 * @param x_p : espérance de l'état prédit \hat{x_{n|n-1}} = F_2^{x,x} \: \hat{x_{n - 1 |n-1}} + Q_2^{x,y} \; y_{n - 2} + F_2^{x,y} \: y_{n - 1}
 * @param p_p : Matrice de covariance de l'état prédit. P_{n|n-1} = Q_2^{x,x} + F_2^{x,x} \: P_{n-1|n-1} \: [F_2^{x,x}]^T 
 * @param[in] _x_f : espérance de l'état filtré précédent
 * @param[in] _sqrt_p_f : racine de la matrice de covariance de l'état filtré précédent
 * @param[in] __y : observation (n - 2)
 * @param[in] _y : observation précédente
 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
 * @param[in] f2_xy : F_2^{x,y} = F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y}
 * @param[in] sqrt_q2_xx : racine de la matrice de covariance du bruit de process réduit
 * @param[in] q2_xy : Q_2^{x,y} = Q^{x,y} \: [Q^{y,y}]^{-1}
 * @param mat_2xx : Matrice de taille (2x.x) préallouée.
 * @param mat_2xx_view_00 : vue sur la matrice mat_2xx allant de (0,0) à (n_x, n_x)
 * @param mat_2xx_view_10 : vue sur la matrice mat_2xx allant de (n_x, 0) à (2n_x, n_x)
 * @param vect_x : vecteur de taille (x) préalloué
 * @brief
 Cette fonction calcule l'espérance et la racine de la matrice de covariance de l'état filtré courant.

**/
void tkalman_do_prediction ( gsl_vector * x_p,
						  gsl_matrix * sqrt_p_p,
						  const gsl_vector * _x_f,
						  const gsl_matrix * _sqrt_p_f,
						  const gsl_vector * __y,
						  const gsl_vector * _y,
						  const gsl_matrix * f2_xx,
						  const gsl_matrix * f2_xy,
						  const gsl_matrix * sqrt_q2_xx,
						  const gsl_matrix * q2_xy,
						  gsl_matrix * mat_2xx,
					      gsl_matrix * mat_2xx_view_00,
					      gsl_matrix * mat_2xx_view_10,
					      gsl_vector * vect_x 
					     )
{
	//Espérance
	tkalman_get_x_p ( x_p,
					  _x_f,
					  __y,
					  _y,
					  f2_xx,
					  f2_xy,
					  q2_xy);
							   
	//Racine de la matrice de covariance				 
	tkalman_get_sqrt_p_p ( sqrt_p_p,
						   _sqrt_p_f,
						   f2_xx,
						   sqrt_q2_xx,
						   mat_2xx,
						   mat_2xx_view_00,
						   mat_2xx_view_10,
						   vect_x);
}

/**@fn void tkalman_do_prediction_1 ( gsl_vector * x_p_1,
						    gsl_matrix * sqrt_p_p_1,
						    const gsl_vector * x_f_0,
						    const gsl_vector * y_f_m1,
						    const gsl_matrix * sqrt_q_f_0,
						    const gsl_vector * y_0,
						    const gsl_matrix * f2_x,
						    const gsl_matrix * f2_xx,
						    const gsl_matrix * f2_xy,
						    const gsl_matrix * sqrt_q2_xx,
						    const gsl_matrix * q2_xy,
						    gsl_matrix * mat_xpt_x,
							gsl_matrix * mat_xpt_x_view_00,
						    gsl_matrix * mat_xpt_x_view_10,
							gsl_vector * vect_x
						   )
 * @param x_p_1 : espérance de l'état prédit 1
 * @param sqrt_p_p_1 : racine de la matrice de covariance de l'état prédit 1
 * @param[in] x_f_0 : espérance de l'état filtré 0
 * @param[in] y_f_m1 : espérance de l'observation filtrée -1
 * @param sqrt_q_0_f : racine de la matrice de covariance du t_0 filtré.
 * @param[in] y_0 : observation 0
 * @param f2_x : matrice dérivée de F, F2 = [F_2^{x,x}, F_2^{x,y}]
 * @param sqrt_q2_xx : racine de la matrice de covariance du bruit de process réduit Q_2^{x,x} = Q^{x,x} - Q^{x,y}[Q^{y,y}]^{-1}Q^{y,x}
 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
 * @param[in] f2_xy : F_2^{x,y} = F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y}
 * @param[in] q2_xy : Q_2^{x,y} = Q^{x,y} \: [Q^{y,y}]^{-1}
 * @param mat_xpt_x : matrice de taille (x+t, x) allouée
 * @param mat_xpt_x_view_00 : vue sur la matrice mat_xpt_x (de (0,0) à (n_x - 1, n_x - 1))
 * @param mat_xpt_x_view_10 : vue sur la matrice mat_xpt_x (de (n_x,0) à (n_x + n_t - 1, n_x - 1))
 * @param vect_x : vecteur de taille (n_x) alloué
 * @brief
 * Cette fonction effectue la première prédiction du filtre de Kalman couple.
 */
void tkalman_do_prediction_1 ( gsl_vector * x_p_1,
						      gsl_matrix * sqrt_p_p_1,
						      const gsl_vector * x_f_0,
						      const gsl_vector * y_f_m1,
						      const gsl_matrix * sqrt_q_f_0,
						      const gsl_vector * y_0,
						      const gsl_matrix * f2_x,
						      const gsl_matrix * f2_xx,
						      const gsl_matrix * f2_xy,
						      const gsl_matrix * sqrt_q2_xx,
						      const gsl_matrix * q2_xy,
						      gsl_matrix * mat_xpt_x,
							  gsl_matrix * mat_xpt_x_view_00,
						      gsl_matrix * mat_xpt_x_view_10,
							  gsl_vector * vect_x
						   )
{
	
	tkalman_get_x_p(x_p_1,
					x_f_0,
					y_f_m1,
					y_0,
					f2_xx,
					f2_xy,
					q2_xy);
					
					
	tkalman_get_sqrt_p_p_1(sqrt_p_p_1,
						   sqrt_q_f_0,
						   f2_x,
						   sqrt_q2_xx,
						   mat_xpt_x,
						   mat_xpt_x_view_00,
						   mat_xpt_x_view_10,
						   vect_x);			
}							

//Objet de prédiction
/**@fn tkalman_prediction :: tkalman_prediction(const gsl_matrix * f2_x,
											const gsl_matrix * sqrt_q2_xx,
											const gsl_matrix * q2_xy)
* @param f2_x : [F2xx, F2xy]
* @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
* @param q2_xy : Qxy.Qyy
* @brief
* Ce constructeur alloue les variables temp. de la prédiction du filtre de Kalman couple.
*/
tkalman_prediction :: tkalman_prediction ( const gsl_matrix * f2_x,
										   const gsl_matrix * sqrt_q2_xx,
										   const gsl_matrix * q2_xy)
{
	tkalman_prediction :: initialize();
	tkalman_prediction :: setup (f2_x,
								 sqrt_q2_xx,
								 q2_xy);
}


/**@fn int tkalman_prediction :: setup(const gsl_matrix * f2_x,
									   const gsl_matrix * sqrt_q2_xx,
									   const gsl_matrix * q2_xy)
 * @param f2_x : [F2xx, F2xy]
 * @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
 * @param q2_xy : Qxy.Qyy
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction permet de modifier les paramètres de la prédiction (size_x et size_y)
**/
int tkalman_prediction :: setup (const gsl_matrix * f2_x,
								 const gsl_matrix * sqrt_q2_xx,
								 const gsl_matrix * q2_xy)
	
{
	unsigned int size_x = f2_x->size1,
				 size_y = f2_x->size2 - size_x;
	
	if (size_x != _size_x || size_y != _size_y)
	{
		tkalman_prediction :: free();
		tkalman_prediction :: initialize();
		if (tkalman_prediction :: set_params(f2_x, sqrt_q2_xx, q2_xy))
		{
			tkalman_prediction :: free();
			tkalman_prediction :: initialize();
			return 1;
		}
	
		
		_size_x = size_x;
		_size_y = size_y;
		if ( tkalman_prediction :: alloc() )
		{
			tkalman_prediction :: free();
			tkalman_prediction :: initialize();
			return 1;
		}
		else
			tkalman_prediction :: create_views();
	}	
	
	return 0;
}

/**@fn tkalman_prediction :: ~ tkalman_prediction()
 * @brief
 * Destructeur
 */
tkalman_prediction :: ~tkalman_prediction()
{
	tkalman_prediction :: free();
	tkalman_prediction :: initialize();
}

/**@fn bool tkalman_prediction :: operator!() const;
 * @return
 * - 0 si l'objet est correctement alloué
 * - 1 sinon
 * @brief
 * Check de l'objet.
 */
bool tkalman_prediction :: operator!() const
{
	return (! (vect_x && mat_xpt_x && _f2_x && _sqrt_q2_xx && _q2_xy) );
}

/**@fn void tkalman_prediction :: free();
 * @brief
 * Cette fonction désalloue la mémoire utilisée par les variables tmp.
 */
void tkalman_prediction :: free()
{
	if (vect_x)
		gsl_vector_free(vect_x);
	if (mat_xpt_x)
		gsl_matrix_free(mat_xpt_x);
}

/**@fn int tkalman_prediction :: alloc();
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction alloue la mémoire utilisée par les variables tmp.
 */
int tkalman_prediction :: alloc()
{
	if (!vect_x)
		vect_x = gsl_vector_alloc(_size_x);
	if (!mat_xpt_x)
		mat_xpt_x = gsl_matrix_alloc(2 * _size_x + _size_y, _size_x);
	return (! (vect_x && mat_xpt_x) );
}

/**@fn void tkalman_prediction :: create_views();
 * @brief
 * Cette méthode crée les vues.
 */
void tkalman_prediction :: create_views()
{
	gsl_matrix_view view;
	view = gsl_matrix_submatrix(mat_xpt_x,
								0,
								0,
								2 * _size_x,
								_size_x);
	mat_2xx = view.matrix;
	
	view = gsl_matrix_submatrix(mat_xpt_x,
								0,
								0,
								_size_x,
								_size_x);
	mat_2xx_view_00 = view.matrix;
	mat_xpt_x_view_00 = view.matrix;
	
	view = gsl_matrix_submatrix(mat_xpt_x,
								_size_x,
								0,
								_size_x,
								_size_x);
	mat_2xx_view_10 = view.matrix;
	view = gsl_matrix_submatrix(mat_xpt_x,
								_size_x,
								0,
								_size_x + _size_y,
								_size_x);
	mat_xpt_x_view_10 = view.matrix;
}

		
/**@fn void tkalman_prediction :: initialize();
 * @brief
 * Cette fonction met les variables tmp à NULL.
 */
void tkalman_prediction :: initialize()
{
	_size_x = 0;
	_size_y = 0;
	vect_x = NULL;
	mat_xpt_x = NULL;	
	_f2_x = NULL;
	_sqrt_q2_xx = NULL;
	_q2_xy = NULL;
}

/**@fn virtual int tkalman_prediction :: set_params(const gsl_matrix * f2_x,
													const gsl_matrix * sqrt_q2_xx,
													const gsl_matrix * q2_xy)
 * @param f2_x : [F2xx, F2xy]
 * @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
 * @param q2_xy : Qxy.Qyy
 * @brief
 * Cette méthode change les paramètres de la prédiction.
 */
int tkalman_prediction :: set_params(const gsl_matrix * f2_x,
									 const gsl_matrix * sqrt_q2_xx,
								     const gsl_matrix * q2_xy)
{

	if (!f2_x || !sqrt_q2_xx || !q2_xy)
		return 1;
	unsigned int size_x = f2_x->size1,
				 size_y = f2_x->size2 - size_x;
	_f2_x = f2_x;
	_sqrt_q2_xx = sqrt_q2_xx;
	_q2_xy = q2_xy;
	//Vues
	{
		gsl_matrix_const_view view = gsl_matrix_const_submatrix (_f2_x, 0, 0, size_x, size_x);
		f2_xx = view.matrix;
	}
	{
		gsl_matrix_const_view view = gsl_matrix_const_submatrix (_f2_x, 0, size_x, size_x, size_y);
		f2_xy = view.matrix;
	}
	return 0;
}
