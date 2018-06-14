#include "tkalman_nc_prediction.hpp"
/**@fn tkalman_nc_prediction :: tkalman_nc_prediction(const gsl_matrix * f,
													  const gsl_matrix * sqrt_q,
													  unsigned int size_x) throw(exception &);
 * @param[in] f: F
 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
 * @param[in] size_x : dim. de x;
 * @brief
 * Constructeur de la classe @class tkalman_nc_prediction
 * @throw 
 * Exception (std :: bad_alloc si problème de mémoire ou invalid_argument en cas d'arguments invalides)
			 */
tkalman_nc_prediction :: tkalman_nc_prediction(const gsl_matrix * f,
											   const gsl_matrix * sqrt_q,
											   unsigned int size_x) throw(exception &)
{
	initialize(); //Tout à zéro
	try
	{
		setup(f,
			  sqrt_q,
			  size_x); //Setup
	}
	catch(exception & except) //En cas de problèmes
	{
		throw(except);
	}
}

/**@fn void tkalman_nc_prediction :: setup(const gsl_matrix * f,
										   const gsl_matrix * sqrt_q,
										   unsigned int size_x) throw(exception &);
 * @param[in] f: F
 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
 * @param[in] size_x : dim. de x;
 * @brief
 * Cette fonction libères les attributs et les réalloue.
 */
void tkalman_nc_prediction :: setup(const gsl_matrix * f,
									const gsl_matrix * sqrt_q,
								    unsigned int size_x) throw(exception &)
{
	//Vérif des arguments
	if (!f || !sqrt_q || !size_x)
		throw(invalid_argument("F, sqrt(Q) or size_x are NULL!\n"));
		
	unsigned int size_t, 
				 size_y;
	size_t = f->size1;
	size_y = size_t - size_x;
	
	//Vérif de la dim. de y;
	if (!size_y)
		throw(invalid_argument("size_y is 0!\n"));
	
	//Modif des matrices
	if (size_x != _size_x || size_y != _size_y)
	{

		free();
		initialize();
		_size_x = size_x;
		_size_y = size_y;
		_size_t = size_t;
		_f = f;
		_sqrt_q = sqrt_q;
		try
		{
			alloc();
		}
		catch (exception & except)
		{
			throw(except);
		}
	}
	else
	{
		_f = f;
		_sqrt_q = sqrt_q;
	}
	//Création des vues
	create_views();
}

/**@fn bool tkalman_nc_prediction :: operator !() const
 * @return 
 * - 0 si l'objet est normal
 * - 1 sinon.
 * @brief
 * Check de l'objet.
 **/
bool tkalman_nc_prediction :: operator !() const
{
	return ( !( mat_xpt_x && mat_2x_x && mat_2t_t && _size_x && _size_y && _size_t && _f && _sqrt_q && vect_x && vect_t) );
}

/**@fn void tkalman_nc_prediction :: compute_prediction_1(gsl_vector * x_p_1,
														  gsl_matrix * sqrt_p_p_1,
														  const gsl_vector * t0_f,
														  const gsl_matrix * sqrt_q_f_0)
  * @param x1_p : \hat{x}_{1|0}, espérance de l'état prédit 1
  * @param sqrt_p_p_1 : [P_{1|0}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit 1.
  * @param[in] t0_f : \hat{t}_{0|0}, espérance de t0 filtré
  * @param[in] sqrt_q_f_0 : [Q_{0|0}]^{\frac{1}{2}}, racine de la matrice de covariance de \hat{t}_{0|0}.
  * @brief
	Cette fonction effectue la prédiction 1 dans le filtre de Kalman non-corrélé (cf @fn void tkalman_nc_get_x1_p et @fn void tkalman_nc_get_sqrt_p1_p pour comprendre ces différentes phases).
 **/
void tkalman_nc_prediction :: compute_prediction_1(gsl_vector * x_p_1,
												   gsl_matrix * sqrt_p_p_1,
												   const gsl_vector * t0_f,
												   const gsl_matrix * sqrt_q_f_0)
{
	tkalman_nc_do_prediction_1 ( x_p_1,
								 sqrt_p_p_1,
								 t0_f,
								 sqrt_q_f_0,
								 &_f_xt,
								 &_sqrt_q_xx,
								 mat_xpt_x,
								 &mat_xpt_x_view_00,
								 &mat_xpt_x_view_10,
								 vect_x );
}

/**@fn void tkalman_nc_prediction :: compute_prediction(gsl_vector * x_p,
							   gsl_matrix * sqrt_p_p,
							   const gsl_vector * _x_f,
							   const gsl_matrix * _sqrt_p_f,
							   const gsl_vector * __y);
  * @param x_p : \hat{x}_{n + 1|n}, espérane de l'état prédit courant
  * @param sqrt_p_p : [P_{n + 1|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit courant.
  * @param[in] _x_f : \hat{x}_{n|n}, espérance de l'état filtré précédent
  * @param[in] _sqrt_p_f : [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré précédent
  * @param[in] __y : y_{n - 1}, observation
  * @brief
	Cette fonction effectue la phase de prédiction dans le filtre de Kalman non-corrélé (cf @fn void tkalman_nc_get_x_p et @fn void tkalman_nc_get_sqrt_p_p pour comprendre ces différentes phases).
 **/
void tkalman_nc_prediction :: compute_prediction(gsl_vector * x_p,
												 gsl_matrix * sqrt_p_p,
												 const gsl_vector * _x_f,
												 const gsl_matrix * _sqrt_p_f,
												 const gsl_vector * __y)
{
	tkalman_nc_do_prediction ( x_p,
							   sqrt_p_p,
							   _x_f,
							   _sqrt_p_f,
							   __y,
							   &_f_xx,
							   &_f_xy,
							   &_sqrt_q_xx,
							   mat_2x_x,
							   &mat_2x_x_view_00,
							   &mat_2x_x_view_10,
							   vect_x );
}

/**@fn void tkalman_nc_prediction :: compute_prediction_(gsl_vector * t_p,
							    gsl_matrix * sqrt_q_p,
								const gsl_vector * _t_p,
								const gsl_matrix * _sqrt_p_p);
* @param t_p : \hat{t}_{n + 1 | m}
* @param sqrt_q_p : [Q_{n + 1|m}]^{\frac{1}{2}}
* @param[in] _t_p : \hat{t}_{n | m}
* @param[in] _sqrt_q_p : [Q_{n|m}]^{\frac{1}{2}}
* @brief
* Cette fonction effectue la prédiction du filtre de Kalman triplet non corrélé sans observations (cf @fn void tkalman_nc_get_t_p et @fn void tkalman_nc_get_sqrt_q_p pour comprendre ces différentes phases).
**/
void tkalman_nc_prediction :: compute_prediction_(gsl_vector * t_p,
												  gsl_matrix * sqrt_q_p,
												  const gsl_vector * _t_p,
												  const gsl_matrix * _sqrt_q_p)
{
	
	tkalman_nc_do_prediction_ ( t_p,
								sqrt_q_p,
								_t_p,
								_sqrt_q_p,
								_f,
								_sqrt_q,
								mat_2t_t,
							    &mat_2t_t_view_00,
						        &mat_2t_t_view_10,
							    vect_t );
}

/**@fn tkalman_nc_prediction :: ~tkalman_nc_prediction();
 * @brief
 * Destructeur de la classe @class tkalman_nc_prediction
 */
tkalman_nc_prediction :: ~tkalman_nc_prediction()
{
	free();
	initialize();
}

/**@fn void tkalman_nc_prediction :: initialize();
 * @brief
 * Cette fonction met tous les attributs de l'objet à 0.
 */
void tkalman_nc_prediction :: initialize()
{
	mat_xpt_x = 0;
	mat_2x_x = 0;
	mat_2t_t = 0;
	vect_x = 0;
	vect_t = 0;
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
	_f = 0;
	_sqrt_q = 0;
}

/**@fn void tkalman_nc_prediction :: alloc() throw(exception &);
 * @brief
 * Cette fonction alloue les différents élements de la classe.
 * @throw
 * bad_alloc en cas de problème de mémoire
 */
void tkalman_nc_prediction :: alloc() throw(exception &)
{
	if (!mat_xpt_x)
	{
		try
		{
			mat_xpt_x = gsl_matrix_alloc(_size_t + _size_x, _size_x);
		}
		catch(exception & except)
		{
			throw(except);
		}	
	}
	
	if (!mat_2x_x)
	{
		try
		{
			mat_2x_x = gsl_matrix_alloc(2 * _size_x, _size_x);
		}
		catch(exception & except)
		{
			throw(except);
		}	
	}
	
	if (!mat_2t_t)
	{
		try
		{
			mat_2t_t = gsl_matrix_alloc(2 * _size_t, _size_t);
		}
		catch(exception & except)
		{
			throw(except);
		}	
	}
	
	if (!vect_x)
	{
		try
		{
			vect_x = gsl_vector_alloc(_size_x);
		}
		catch(exception & except)
		{
			throw(except);
		}	
	}
	
	if (!vect_t)
	{
		try
		{
			vect_t = gsl_vector_alloc(_size_t);
		}
		catch(exception & except)
		{
			throw(except);
		}	
	}
}

/**@fn void tkalman_nc_prediction :: free();
 * @brief
 * Cette fonction désalloue tous les attributs alloués.
 **/
void tkalman_nc_prediction :: free()
{
	if (mat_xpt_x)	
		gsl_matrix_free(mat_xpt_x);
	if (mat_2x_x)	
		gsl_matrix_free(mat_2x_x);
	if (mat_2t_t)	
		gsl_matrix_free(mat_2t_t);
	if (vect_x)	
		gsl_vector_free(vect_x);
	if (vect_t)	
		gsl_vector_free(vect_t);
}

/**@fn void tkalman_nc_prediction :: create_views();
 * @brief 
 * Cette fonction génère les différentes vues sur les matrices.
 */
void tkalman_nc_prediction :: create_views()
{
	gsl_matrix_view view;
	//Mat_2xx
	view = gsl_matrix_submatrix(mat_2x_x, 
								0, 
								0, 
								_size_x, 
								_size_x);
	mat_2x_x_view_00 = view.matrix;
	view = gsl_matrix_submatrix(mat_2x_x, 
								_size_x, 
								0, 
								_size_x, 
								_size_x);
	mat_2x_x_view_10 = view.matrix;
	//Mat_2tt
	view = gsl_matrix_submatrix(mat_2t_t, 
								0, 
								0, 
								_size_t, 
								_size_t);
	mat_2t_t_view_00 = view.matrix;
	view = gsl_matrix_submatrix(mat_2t_t, 
								_size_t, 
								0, 
								_size_t, 
								_size_t);
	mat_2t_t_view_10 = view.matrix;
	//Mat_xpt_x
	view = gsl_matrix_submatrix(mat_xpt_x, 
								0, 
								0, 
								_size_x, 
								_size_x);
	mat_xpt_x_view_00 = view.matrix;
	view = gsl_matrix_submatrix(mat_xpt_x, 
								_size_x, 
								0, 
								_size_t, 
								_size_x);
	mat_xpt_x_view_10 = view.matrix;
	
	{
		gsl_matrix_const_view view2 = gsl_matrix_const_submatrix(_f, 0, 0, _size_x, _size_t);
		_f_xt = view2.matrix;
	}
	{
		gsl_matrix_const_view view2 = gsl_matrix_const_submatrix(_f, 0, 0, _size_x, _size_x);
		_f_xx = view2.matrix;
	}
	{
		gsl_matrix_const_view view2 = gsl_matrix_const_submatrix(_f, 0, _size_x, _size_x, _size_y);
		_f_xy = view2.matrix;
	}
	{
		gsl_matrix_const_view view2 = gsl_matrix_const_submatrix(_sqrt_q, 0, 0, _size_x, _size_x);
		_sqrt_q_xx = view2.matrix;
	}
	{
		gsl_matrix_const_view view2 = gsl_matrix_const_submatrix(_sqrt_q, _size_x, _size_x, _size_y, _size_y);
		_sqrt_q_yy = view2.matrix;
	}
	
}



/**@fn void tkalman_nc_get_x1_p(gsl_vector * x1_p,
								const gsl_vector * t0_f,
								const gsl_matrix * f_xt)

 * @param x1_p :  \hat{x}_{1|0}, espérance de l'état prédit 1
 * @param[in] t0_f : \hat{t}_{0|0}, espérance de t0 filtré
 * @param[in] f_xt : F^{x,t}.
 * @brief
 * Cette fonction calcule l'espérance de l'état prédit 1 dans le filtre de Kalman triplet sans corrélation entre bruit de process et bruit de mesure. \n
 \hat{x}_{1|0} = F^{x,t}  \hat{t}_{0|0}
 **/
void tkalman_nc_get_x1_p ( gsl_vector * x1_p,
						   const gsl_vector * t0_f,
						   const gsl_matrix * f_xt )
{
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   f_xt,
				   t0_f,
				   0.0,
				   x1_p);
}



/**@fn void tkalman_nc_get_x_p(gsl_vector * x_p,
							   const gsl_vector * _x_f,
							   const gsl_vecotr * __y,
							   const gsl_matrix * f_xx,
							   const gsl_matrix * f_xy)
 * @param x_p : \hat{x}_{n + 1|n}, espérane de l'état prédit courant
 * @param[in] _x_f : \hat{x}_{n|n}, espérance de l'état filtré précédent
 * @param[in] __y : y_{n - 1}, observation
 * @param[in] f_xx : F^{x,x}
 * @param[in] f_xy : F^{x,y} 
 * @brief
 * Cette fonction calcule l'espérance de l'état prédit n dans le filtre de Kalman triplet sans corrélation entre bruit de mesure et bruit de process. Cette fonction est valable jusqu'à n = Card(Y). Après, il est nécessaire d'estimer aussi y. \n
 * \hat{x}_{n + 1|n} = F^{x,x} \hat{x}_{n|n} + F^{x,y} y_{n}
 */
void tkalman_nc_get_x_p( gsl_vector * x_p,
						 const gsl_vector * _x_f,
						 const gsl_vector * __y,
						 const gsl_matrix * f_xx,
						 const gsl_matrix * f_xy )
{
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   f_xx,
				   _x_f,
				   0.0,
				   x_p);
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   f_xy,
				   __y,
				   1.0,
				   x_p);	   
}

/**@fn void tkalman_nc_get_t_p( gsl_vector * t_p,
								const gsl_vector * _t_p,
								const gsl_matrix * f )
 * @param t_p : \hat{t}_{n + 1 | m}
 * @param[in] _t_p : \hat{t}_{n | m}
 * @param[in] f : F
 * @brief
 * Cette fonction calcule l'espérance du vecteur t prédit suivant. Il faut que m soit strictement inférieur à n. \n
 * \hat{t}_{n + 1 | m} = F \hat{t}_{n | m}
 */
void tkalman_nc_get_t_p( gsl_vector * t_p,
						 const gsl_vector * _t_p,
						 const gsl_matrix * f)
{
	gsl_blas_dgemv(CblasNoTrans,
				   1.0,
				   f,
				   _t_p,
				   0.0,
				   t_p);
	
}

/**@fn void tkalman_nc_get_sqrt_p1_p(gsl_matrix * sqrt_p_p_1,
									 const gsl_matrix * sqrt_q_f_0,
									 const gsl_matrix * f_xt,
									 const gsl_matrix * sqrt_q_xx,
									 gsl_matrix * mat_xpt_x,
									 gsl_matrix * mat_xpt_x_view_00,
									 gsl_matrix * mat_xpt_x_view_10,
									 gsl_vector * vect_x)
 * @param sqrt_p_p_1 : [P_{1|0}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit 1.
 * @param[in] sqrt_q_f_0 : [Q_{0|0}]^{\frac{1}{2}}, racine de la matrice de covariance de \hat{t}_{0|0}.
 * @param[in] f_xt : F^{x,t}.
 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
 * @param mat_xpt_x : matrice de taille (n_x+n_t, n_x) allouée.
 * @param mat_xpt_x_view_00 : vue sur la matrice mat_xpt_x (de (0,0) à (n_x - 1, n_x - 1)).
 * @param mat_xpt_x_view_10 : vue sur la matrice mat_xpt_x (de (n_x,0) à (n_x + n_t - 1, n_x - 1)).
 * @param vect_x : vecteur de taille (n_x) alloué
 * @brief
 * Cette fonction calcule la racine de la matrice de covariance de l'état prédit 1, [P_{1|0}]^{\frac{1}{2}} dans le filtre de Kalman triplet non-corrélé. \n
 * Dans un premier temps, nous construisons la matrice M : \newline
 * M = \begin{pmatrix}
 * [Q^{x,x}]^{\frac{1}{2}} \newline
 * [Q_{0|0}]^{\frac{1}{2}} [F^{x,t}]^T
 * \end{pmatrix}
 * \newline
 * Dans un second temps, nous effectuons la décomposition QR de cette matrice : 
 * M = Q
 * \begin{pmatrix}
 * [P_{1|0}]^{\frac{1}{2}}} \newline
 * 0
 * \end{pmatrix}
 * 
 */
void tkalman_nc_get_sqrt_p1_p( gsl_matrix * sqrt_p_p_1,
							   const gsl_matrix * sqrt_q_f_0,
							   const gsl_matrix * f_xt,
							   const gsl_matrix * sqrt_q_xx,
						       gsl_matrix * mat_xpt_x,
							   gsl_matrix * mat_xpt_x_view_00,
						       gsl_matrix * mat_xpt_x_view_10,
							   gsl_vector * vect_x )
{
	//Construction de la matrice
	// Q2xx
	//
	// sqrt_q_f . F2_x^T
	
	gsl_matrix_memcpy(mat_xpt_x_view_00, sqrt_q_xx);
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   sqrt_q_f_0,
				   f_xt,
				   0.0,
				   mat_xpt_x_view_10);
				   
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_xpt_x,
						 vect_x);
	gsl_triangle_matrix(mat_xpt_x);
	
	//Recopie du résultat
	gsl_matrix_memcpy(sqrt_p_p_1, mat_xpt_x_view_00);
}

/**@fn void tkalman_nc_get_sqrt_q_p( gsl_matrix * sqrt_q_p,
							  const gsl_matrix * _sqrt_q_p,
							  const gsl_matrix * f,
							  const gsl_matrix * sqrt_q,
						      gsl_matrix * mat_2t_t,
							  gsl_matrix * mat_2t_t_view_00,
						      gsl_matrix * mat_2t_t_view_10,
							  gsl_vector * vect_t )
 * @param sqrt_q_p : [Q_{n + 1|m}]^{\frac{1}{2}}
 * @param[in] _sqrt_q_p : [Q_{n|m}]^{\frac{1}{2}}
 * @param[in] f : F
 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
 * @param mat_2tt : Matrice de taille (2n_t.n_t) préallouée.
 * @param mat_2tt_view_00 : vue sur la matrice mat_2tt allant de (0,0) à (n_t - 1, n_t - 1)
 * @param mat_2tt_view_10 : vue sur la matrice mat_2tt allant de (n_t, 0) à (2n_t - 1, n_t - 1)
 * @param vect_t : vecteur de taille (n_t) préalloué
 * @brief
 * Cette fonction calcule racine de la matrice de covariance du vecteur aléatoire t_{n + 1|m}. \n
 * Dans un premier temps, nous construisons la matrice M : \newline
 * M = \begin{pmatrix}
 * [Q]^{\frac{1}{2}} \newline
 * [Q_{n|m}]^{\frac{1}{2}} F^T
 * \end{pmatrix}
 * \newline
 * Dans un second temps, nous effectuons la décomposition QR de cette matrice : 
 * M = Q
 * \begin{pmatrix}
 * [Q_{n + 1|m}]^{\frac{1}{2}}} \newline
 * 0
 * \end{pmatrix}
 **/
void tkalman_nc_get_sqrt_q_p( gsl_matrix * sqrt_q_p,
							  const gsl_matrix * _sqrt_q_p,
							  const gsl_matrix * f,
							  const gsl_matrix * sqrt_q,
						      gsl_matrix * mat_2tt,
							  gsl_matrix * mat_2tt_view_00,
						      gsl_matrix * mat_2tt_view_10,
							  gsl_vector * vect_t )
{
	//Construction de la matrice
	// Q
	//
	// sqrt_q_f . F^T
	
	gsl_matrix_memcpy(mat_2tt_view_00, sqrt_q);
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   _sqrt_q_p,
				   f,
				   0.0,
				   mat_2tt_view_10);
				   
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_2tt,
						 vect_t);
	gsl_triangle_matrix(mat_2tt);
	
	//Recopie du résultat
	gsl_matrix_memcpy(sqrt_q_p, mat_2tt_view_00);
}




/**@fn void tkalman_nc_get_sqrt_p_p(gsl_matrix * sqrt_p_p,
									const gsl_matrix * _sqrt_p_f,
									const gsl_matrix * f_xx,
									const gsl_matrix * sqrt_q_xx,
									gsl_matrix * mat_2x_x,
									gsl_matrix * mat_2x_x_view_0,
									gsl_matrix * mat_2x_x_view_1,
									gsl_vector * vect_x)
 * @param sqrt_p_p : [P_{n + 1|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit courant.
 * @param[in] _sqrt_p_f : [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré précédent
 * @param[in] f_xx : F^{x,x}
 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
 * @param mat_2xx : Matrice de taille (2n_x.n_x) préallouée.
 * @param mat_2xx_view_00 : vue sur la matrice mat_2xx allant de (0,0) à (n_x - 1, n_x - 1)
 * @param mat_2xx_view_10 : vue sur la matrice mat_2xx allant de (n_x, 0) à (2n_x - 1, n_x - 1)
 * @param vect_x : vecteur de taille (n_x) préalloué
 * @brief
 * Cette fonction calcule la racine de la matrice de covariance de l'état prédit, [P_{n + 1|n}]^{\frac{1}{2}} dans le filtre de Kalman triplet non-corrélé. \n
 * Dans un premier temps, nous construisons la matrice M : \newline
 * M = \begin{pmatrix}
 * [Q^{x,x}]^{\frac{1}{2}} \newline
 * [P_{n|n}]^{\frac{1}{2}} [F^{x,x}]^T
 * \end{pmatrix}
 * \newline
 * Dans un second temps, nous effectuons la décomposition QR de cette matrice : 
 * M = Q
 * \begin{pmatrix}
 * [P_{n + 1|n}]^{\frac{1}{2}} \newline
 * 0
 * \end{pmatrix}
 * 
 */
void tkalman_nc_get_sqrt_p_p ( gsl_matrix * sqrt_p_p,
						       const gsl_matrix * _sqrt_p_f,
						       const gsl_matrix * f_xx,
						       const gsl_matrix * sqrt_q_xx,
						       gsl_matrix * mat_2x_x,
						       gsl_matrix * mat_2x_x_view_00,
						       gsl_matrix * mat_2x_x_view_10,
						       gsl_vector * vect_x )
{
	
	
	//Construction de la matrice
	// Qxx
	//
	// sqrt_p_f . F_xx^T
	
	gsl_matrix_memcpy(mat_2x_x_view_00, sqrt_q_xx);
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   _sqrt_p_f,
				   f_xx,
				   0.0,
				   mat_2x_x_view_10);
				   
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_2x_x,
						 vect_x);
	gsl_triangle_matrix(mat_2x_x);
	
	//Recopie du résultat
	gsl_matrix_memcpy(sqrt_p_p, mat_2x_x_view_00);
	
}

/**@fn void tkalman_nc_do_prediction ( gsl_vector * x_p,
									   gsl_matrix * sqrt_p_p,
									   const gsl_vector * _x_f,
									   const gsl_matrix * _sqrt_p_f,
									   const gsl_vector * __y,
									   const gsl_matrix * f_xx,
									   const gsl_matrix * f_xy,
									   const gsl_matrix * sqrt_q_xx,
									   gsl_matrix * mat_2xx,
									   gsl_matrix * mat_2xx_view_00,
									   gsl_matrix * mat_2xx_view_10,
									   gsl_vector * vect_x )
 * @param x_p : \hat{x}_{n + 1|n}, espérane de l'état prédit courant
 * @param sqrt_p_p : [P_{n + 1|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit courant.
 * @param[in] _x_f : \hat{x}_{n|n}, espérance de l'état filtré précédent
 * @param[in] _sqrt_p_f : [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré précédent
 * @param[in] __y : y_{n - 1}, observation
 * @param[in] f_xx : F^{x,x}
 * @param[in] f_xy : F^{x,y} 
 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
 * @param mat_2xx : Matrice de taille (2x.x) préallouée.
 * @param mat_2xx_view_00 : vue sur la matrice mat_2xx allant de (0,0) à (n_x - 1, n_x - 1)
 * @param mat_2xx_view_10 : vue sur la matrice mat_2xx allant de (n_x, 0) à (2n_x - 1, n_x - 1)
 * @param vect_x : vecteur de taille (x) préalloué
 * @brief
 Cette fonction effectue la phase de prédiction dans le filtre de Kalman non-corrélé (cf @fn void tkalman_nc_get_x_p et @fn void tkalman_nc_get_sqrt_p_p pour comprendre ces différentes phases).
**/
void tkalman_nc_do_prediction ( gsl_vector * x_p,
								gsl_matrix * sqrt_p_p,
								const gsl_vector * _x_f,
								const gsl_matrix * _sqrt_p_f,
								const gsl_vector * __y,
								const gsl_matrix * f_xx,
								const gsl_matrix * f_xy,
								const gsl_matrix * sqrt_q_xx,
								gsl_matrix * mat_2xx,
								gsl_matrix * mat_2xx_view_00,
								gsl_matrix * mat_2xx_view_10,
								gsl_vector * vect_x )
{
		tkalman_nc_get_x_p( x_p,
						    _x_f,
						    __y,
						    f_xx,
						    f_xy );
						    
		tkalman_nc_get_sqrt_p_p ( sqrt_p_p,
						          _sqrt_p_f,
						          f_xx,
						          sqrt_q_xx,
						          mat_2xx,
						          mat_2xx_view_00,
						          mat_2xx_view_10,
						          vect_x );
}

/**@fn void tkalman_nc_do_prediction_ ( gsl_vector * t_p,
										gsl_matrix * _sqrt_q_p,
										const gsl_vector * _t_p,
										const gsl_matrix * _sqrt_q_p,
										const gsl_matrix * f,
										const gsl_matrix * sqrt_q,
										gsl_matrix * mat_2tt,
										gsl_matrix * mat_2tt_view_00,
										gsl_matrix * mat_2tt_view_10,
										gsl_vector * vect_t )
 * @param t_p : \hat{t}_{n + 1 | m}
 * @param sqrt_q_p : [Q_{n + 1|m}]^{\frac{1}{2}}
 * @param[in] _t_p : \hat{t}_{n | m}
 * @param[in] _sqrt_q_p : [Q_{n|m}]^{\frac{1}{2}}* 
 * @param[in] f : F
 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
 * @param mat_2tt : Matrice de taille (2n_t.n_t) préallouée.
 * @param mat_2tt_view_00 : vue sur la matrice mat_2tt allant de (0,0) à (n_t - 1, n_t - 1)
 * @param mat_2tt_view_10 : vue sur la matrice mat_2tt allant de (n_t, 0) à (2n_t - 1, n_t - 1)
 * @param vect_t : vecteur de taille (n_t) préalloué
 * @brief
 * Cette fonction effectue la prédiction du filtre de Kalman triplet non corrélé sans observations (cf @fn void tkalman_nc_get_t_p et @fn void tkalman_nc_get_sqrt_q_p pour comprendre ces différentes phases).
 */
void tkalman_nc_do_prediction_ ( gsl_vector * t_p,
								 gsl_matrix * sqrt_q_p,
								 const gsl_vector * _t_p,
								 const gsl_matrix * _sqrt_q_p,
								 const gsl_matrix * f,
								 const gsl_matrix * sqrt_q,
								 gsl_matrix * mat_2tt,
							     gsl_matrix * mat_2tt_view_00,
						         gsl_matrix * mat_2tt_view_10,
							     gsl_vector * vect_t )
{
		tkalman_nc_get_t_p( t_p,
						    _t_p,
						    f);
						    
		tkalman_nc_get_sqrt_q_p ( sqrt_q_p,
						          _sqrt_q_p,
						          f,
						          sqrt_q,
						          mat_2tt,
						          mat_2tt_view_00,
						          mat_2tt_view_10,
						          vect_t );
}							     


/**@fn void tkalman_nc_do_prediction_1 ( gsl_vector * x_p_1,
										 gsl_matrix * sqrt_p_p_1,
										 const gsl_vector * t0_f,
										 const gsl_matrix * sqrt_q_f_0,
										 const gsl_matrix * f_xt,
										 const gsl_matrix * sqrt_q_xx,
										 gsl_matrix * mat_xpt_x,
										 gsl_matrix * mat_xpt_x_view_00,
										 gsl_matrix * mat_xpt_x_view_10,
										 gsl_vector * vect_x )
 * @param x1_p : \hat{x}_{1|0}, espérance de l'état prédit 1
 * @param sqrt_p_p_1 : [P_{1|0}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit 1.
 * @param[in] t0_f : \hat{t}_{0|0}, espérance de t0 filtré
 * @param[in] sqrt_q_f_0 : [Q_{0|0}]^{\frac{1}{2}}, racine de la matrice de covariance de \hat{t}_{0|0}.
 * @param[in] f_xt : F^{x,t}.
 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
 * @param mat_xpt_x : matrice de taille (x+t, x) allouée
 * @param mat_xpt_x_view_00 : vue sur la matrice mat_xpt_x (de (0,0) à (n_x - 1, n_x - 1))
 * @param mat_xpt_x_view_10 : vue sur la matrice mat_xpt_x (de (n_x,0) à (n_x + n_t - 1, n_x - 1))
 * @param vect_x : vecteur de taille (n_x) alloué
 * @brief
 Cette fonction effectue la prédiction 1 dans le filtre de Kalman non-corrélé (cf @fn void tkalman_nc_get_x1_p et @fn void tkalman_nc_get_sqrt_p1_p pour comprendre ces différentes phases).
 */
void tkalman_nc_do_prediction_1 ( gsl_vector * x_p_1,
								  gsl_matrix * sqrt_p_p_1,
								  const gsl_vector * t0_f,
								  const gsl_matrix * sqrt_q_f_0,
								  const gsl_matrix * f_xt,
								  const gsl_matrix * sqrt_q_xx,
								  gsl_matrix * mat_xpt_x,
								  gsl_matrix * mat_xpt_x_view_00,
								  gsl_matrix * mat_xpt_x_view_10,
								  gsl_vector * vect_x )
{
	
	tkalman_nc_get_x1_p ( x_p_1,
						  t0_f,
						  f_xt );
						  
	tkalman_nc_get_sqrt_p1_p( sqrt_p_p_1,
							  sqrt_q_f_0,
							  f_xt,
							  sqrt_q_xx,
						      mat_xpt_x,
							  mat_xpt_x_view_00,
						      mat_xpt_x_view_10,
							  vect_x );			  
}



