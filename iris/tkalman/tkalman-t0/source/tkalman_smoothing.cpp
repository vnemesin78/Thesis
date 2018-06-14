/**@file tkalman_smoothing.cpp
*@author Valérian Némesin
*@brief
Ce fichier contient le code source des fonctions nécessaires au lissage dans le filtre de Kalman triplet.
**/
#include "tkalman_smoothing.hpp"
/**@fn void tkalman_get_x_s(gsl_vector * x_s,
								     const gsl_vector * x_f,
									 const gsl_vector * x_p_,
									 const gsl_vector * x_s_,
									 const gsl_matrix * gain)
 * @param x_s : Espérance de l'état lissé
 * @param[in] x_f : Espérance de l'état filtré
 * @param[in] x_p_ : Espérance de l'état prédit suivant
 * @param[in] x_s_ : Espérance de l'état lissé suivant
 * @param[in] gain : Gain de lissage
 * @brief
 * Cette fonction calcule l'espérance de l'état lissé selon la formule \n
 * \hat{x}_{n|N} = \hat{x}_{n|n} + K_{n|N} \: (\hat{x}_{n + 1|N} - \hat{x}_{n+1|n}
 * .
 */
void tkalman_get_x_s(gsl_vector * x_s,
					 const gsl_vector * x_f,
					 const gsl_vector * x_p_,
			  	     const gsl_vector * x_s_,
				     const gsl_matrix * gain) //mat xx ou tx
{
	gsl_vector_memcpy(x_s,
					  x_f);

	gsl_blas_dgemv (CblasNoTrans,
					1.0,
					gain,
					x_s_,
					1.0,
					x_s);

	gsl_blas_dgemv (CblasNoTrans,
					-1.0,
					gain,
					x_p_,
					1.0,
					x_s);
}
/**@fn void tkalman_get_smoothing_gain_0(gsl_matrix * s_gain,
									     const gsl_matrix * sqrt_p_f,
									     const gsl_matrix * sqrt_p_p_,
								         const gsl_matrix * f2_xx,
									     gsl_matrix * mat_xx,
									     gsl_permutation * perm_x)
 * @param s_gain : gain de Kalman triplet pour le lissage
 * @param[in] sqrt_p_f : racine de la covariance de l'état filtré courant
 * @param[in] sqrt_p_p_ : racine de la covariance de l'état prédit suivant
 * @param[in] f2_xx : Fxx - Qxy.Qyy?¹.Fyx
 * @param mat_xx : matrice de taille (x.x) préallouée
 * @param perm_x : permutation de taille (x) préallouée
 * @brief
 * Cette fonction calcule le gain de lissage du filtre de Kalman
 * triplet.
 */
void tkalman_get_smoothing_gain_0(gsl_matrix * s_gain,
							      const gsl_matrix * sqrt_q_f,
							      const gsl_matrix * sqrt_p_p_,
							      const gsl_matrix * f2_x_,
							      gsl_matrix * mat_xx,
							      gsl_matrix * mat_tx,
							      gsl_permutation * perm_x)
{
	//Q F2_x_^T
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   sqrt_q_f, //tt
  				   f2_x_,   //tx
				   0.0,
				   mat_tx);
	gsl_blas_dgemm(CblasTrans,
				   CblasNoTrans,
				   1.0,
				   sqrt_q_f, //tt
  				   mat_tx,   //tx
				   0.0,
				   s_gain);     
	//P_p - 1
	gsl_permutation_init(perm_x);
	gsl_linalg_LU_invert(sqrt_p_p_,
					     perm_x,
						 mat_xx);
	
	
	
	gsl_blas_dgemm(CblasNoTrans,
				   CblasNoTrans,
				   1.0,
				   s_gain, //tt
  				   mat_xx,   //tx
				   0.0,
				   mat_tx);
	
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   mat_tx, //tt
  				   mat_xx,   //tx
				   0.0,
				   s_gain);
	
}

/**@fn void tkalman_get_sqrt_p_s_and_c_s(gsl_matrix * sqrt_p_s,
							                gsl_matrix * c_s,
								     const gsl_matrix * sqrt_p_f,
								     const gsl_matrix * sqrt_p_s_,
								     const gsl_matrix * f2_xx,
								     const gsl_matrix * sqrt_q2_xx,
								     const gsl_matrix * s_gain,
								     gsl_matrix * mat_3x2x,
								     gsl_matrix * mat_3x2x_view_00,
								    gsl_matrix * mat_3x2x_view_01,
								    gsl_matrix * mat_3x2x_view_10,
								    gsl_matrix * mat_3x2x_view_11,
								    gsl_matrix * mat_3x2x_view_20,
								    gsl_matrix * mat_3x2x_view_21,
								    gsl_vector * vect_2x)
 * @param sqrt_p_s : racine de la matrice de covariance de l'état lissé courant
 * @param c_s : sqrt(P_{n+1|N}) \; K_{n|N}^T
 * @param[in] sqrt_p_f : racine de la matrice de covariance de l'état filtré courant
 * @param[in]  sqrt_p_s_ : racine de la matrice de covariance de l'état lissé suivant
 * @param[in] f2_xx : Fxx - Qxy.Qyy?¹.Fyx
 * @param[in] sqrt_q2_xx : racine de la matrice de covariance du bruit de process réduit
 * @param[in] s_gain :  gain de Kalman triplet pour le lissage
 * @param mat_3x2x : matrice de taille (3x.2x) préallouée
 * @param mat_3x2x_view_00 : vue sur la matrice mat_3x2x allant de (0,0) à (n_x, n_x)
 * @param mat_3x2x_view_01 : vue sur la matrice mat_3x2x allant de (0,n_x) à (n_x,2 n_x)
 * @param mat_3x2x_view_10 : vue sur la matrice mat_3x2x allant de (n_x,0) à (2n_x, n_x)
 * @param mat_3x2x_view_11 : vue sur la matrice mat_3x2x allant de (n_x,n_x) à (2n_x, 2n_x)
 * @param mat_3x2x_view_20 : vue sur la matrice mat_3x2x allant de (2n_x,0) à (3n_x, n_x)
 * @param mat_3x2x_view_21 : vue sur la matrice mat_3x2x allant de (2n_x,n_x) à (3n_x, 2n_x)
 * @param vect_2x : vecteur de taille (2x) préalloué
 * @brief
 Cette fonction calcule la racine de la matrice de covariance de l'état lissé courant et sqrt(P_{n+1|N}) \; K_{n|N}^T.
**/
void tkalman_get_sqrt_p_s_and_c_s(gsl_matrix * sqrt_p_s,
								  gsl_matrix * c_s,
								  const gsl_matrix * sqrt_p_f,
								  const gsl_matrix * sqrt_p_s_,
								  const gsl_matrix * f2_xx,
								  const gsl_matrix * sqrt_q2_xx,
								  const gsl_matrix * s_gain,
								  gsl_matrix * mat_3x2x,
								  gsl_matrix * mat_3x2x_view_00,
								  gsl_matrix * mat_3x2x_view_01,
								  gsl_matrix * mat_3x2x_view_10,
								  gsl_matrix * mat_3x2x_view_11,
								  gsl_matrix * mat_3x2x_view_20,
								  gsl_matrix * mat_3x2x_view_21,
								  gsl_vector * vect_2x)
{
	//Création de la matrice (2 matrices possibles)
	// sqrt_Q2xx		0
	// sqrt_p_f F2xx^T	sqrt_p_f
	// 0				sqrt(P_s_) K^T
		gsl_matrix_memcpy(mat_3x2x_view_00,
						  sqrt_q2_xx);
		gsl_matrix_set_zero(mat_3x2x_view_01);
		gsl_blas_dgemm(CblasNoTrans,
					   CblasTrans,
					   1.0,
					   sqrt_p_f,
					   f2_xx,
					   0.0,
					   mat_3x2x_view_10);
		gsl_matrix_memcpy(mat_3x2x_view_11,
						  sqrt_p_f);
		gsl_matrix_set_zero(mat_3x2x_view_20);
		//Calcul de la covariance entre les états lissés courant et suivant
		gsl_blas_dgemm(CblasNoTrans,
					   CblasTrans,
					   1.0,
					   sqrt_p_s_,
					   s_gain,
					   0.0,
					   mat_3x2x_view_21);
		//Sauvegarde du résultat
		if (c_s != NULL)
			gsl_matrix_memcpy(c_s,
							  mat_3x2x_view_21);

	//Décomposition QR
		gsl_linalg_QR_decomp(mat_3x2x,
							 vect_2x);
		gsl_triangle_matrix(mat_3x2x);

	//Recopie
		gsl_matrix_memcpy(sqrt_p_s,
						  mat_3x2x_view_11);


}

/**@fn void tkalman_do_smoothing(gsl_vector * x_s,
							  gsl_matrix * sqrt_p_s,
							  gsl_matrix * c_s,
							  const gsl_vector * x_f,
							  const gsl_matrix * sqrt_p_f,
							  const gsl_vector * x_p_,
							  const gsl_matrix * sqrt_p_p_,
							  const gsl_vector * x_s_,
							  const gsl_matrix * sqrt_p_s_,
							  const gsl_matrix * f2_xx,
							  const gsl_matrix * sqrt_q2_xx,
							  gsl_matrix * mat_xx,
							  gsl_matrix * mat_3x2x,
							  gsl_matrix * mat_3x2x_view_00,
							  gsl_matrix * mat_3x2x_view_01,
							  gsl_matrix * mat_3x2x_view_10,
							  gsl_matrix * mat_3x2x_view_11,
							  gsl_matrix * mat_3x2x_view_20,
							  gsl_matrix * mat_3x2x_view_21,
							  gsl_permutation * perm_x,
							  gsl_vector * vect_2x)
 * @param x_s : Espérance de l'état lissé
 * @param sqrt_p_s : racine de la matrice de covariance de l'état lissé courant
 * @param c_s : sqrt(P_{n+1|N}) \; K_{n|N}^T
 * @param[in] x_f : Espérance de l'état filtré
 * @param[in] sqrt_p_f : racine de la matrice de covariance de l'état filtré courant
 * @param[in] x_p_ : Espérance de l'état prédit suivant
 * @param[in] sqrt_p_p_ : racine de la covariance de l'état prédit suivant
 * @param[in] x_s_ : Espérance de l'état lissé suivant
 * @param[in]  sqrt_p_s_ : racine de la matrice de covariance de l'état lissé suivant
 * @param[in] f2_xx : Fxx - Qxy.Qyy?¹.Fyx
 * @param[in] sqrt_q2_xx : racine de la matrice de covariance du bruit de process réduit
 * @param mat_xx : matrice de taille (x.x) préallouée
 * @param mat_3x2x : matrice de taille (3x.2x) préallouée
 * @param mat_3x2x_view_00 : vue sur la matrice mat_3x2x allant de (0,0) à (n_x, n_x)
 * @param mat_3x2x_view_01 : vue sur la matrice mat_3x2x allant de (0,n_x) à (n_x,2 n_x)
 * @param mat_3x2x_view_10 : vue sur la matrice mat_3x2x allant de (n_x,0) à (2n_x, n_x)
 * @param mat_3x2x_view_11 : vue sur la matrice mat_3x2x allant de (n_x,n_x) à (2n_x, 2n_x)
 * @param mat_3x2x_view_20 : vue sur la matrice mat_3x2x allant de (2n_x,0) à (3n_x, n_x)
 * @param mat_3x2x_view_21 : vue sur la matrice mat_3x2x allant de (2n_x,n_x) à (3n_x, 2n_x)
 * @param perm_x : permutation de taille (x) préallouée
 * @param vect_2x : vecteur de taille (2x) préalloué
 * @brief
 Cette fonction effectue l'étape de lissage dans le filtre de Kalman triplet.
**/
void tkalman_do_smoothing(gsl_vector * x_s,
						  gsl_matrix * sqrt_p_s,
						  gsl_matrix * c_s,
						  const gsl_vector * x_f,
						  const gsl_matrix * sqrt_p_f,
						  const gsl_vector * x_p_,
						  const gsl_matrix * sqrt_p_p_,
						  const gsl_vector * x_s_,
						  const gsl_matrix * sqrt_p_s_,
						  const gsl_matrix * f2_xx,
						  const gsl_matrix * sqrt_q2_xx,
						  gsl_matrix * mat_xx,
						  gsl_matrix * mat_3x2x,
						  gsl_matrix * mat_3x2x_view_00,
						  gsl_matrix * mat_3x2x_view_01,
						  gsl_matrix * mat_3x2x_view_10,
						  gsl_matrix * mat_3x2x_view_11,
						  gsl_matrix * mat_3x2x_view_20,
						  gsl_matrix * mat_3x2x_view_21,
						  gsl_permutation * perm_x,
						  gsl_vector * vect_2x)
{
	//Gain
	tkalman_get_smoothing_gain_0( mat_xx,
								  sqrt_p_f,
								  sqrt_p_p_,
								  f2_xx,
							      mat_3x2x_view_00,
							      mat_3x2x_view_01,
							      perm_x);
	//Racine des matrices de covariance
	tkalman_get_sqrt_p_s_and_c_s( sqrt_p_s,
										 c_s,
										 sqrt_p_f,
										 sqrt_p_s_,
										 f2_xx,
										 sqrt_q2_xx,
										 mat_xx,
										 mat_3x2x,
										 mat_3x2x_view_00,
										 mat_3x2x_view_01,
										 mat_3x2x_view_10,
										 mat_3x2x_view_11,
										 mat_3x2x_view_20,
										 mat_3x2x_view_21,
										 vect_2x);
	//Espérance
	tkalman_get_x_s( x_s,
					 x_f,
					 x_p_,
					 x_s_,
					 mat_xx);

}

/**@fn void tkalman_do_smoothing_0(gsl_vector * x_s,
							  gsl_matrix * sqrt_p_s,
							  gsl_matrix * c_s,
							  const gsl_vector * x_f,
							  const gsl_matrix * sqrt_p_f,
							  const gsl_vector * x_p_,
							  const gsl_matrix * sqrt_p_p_,
							  const gsl_vector * x_s_,
							  const gsl_matrix * sqrt_p_s_,
							  const gsl_matrix * f2_xx,
							  const gsl_matrix * sqrt_q2_xx,
							  gsl_matrix * mat_xx,
							  gsl_matrix * mat_3x2x,
							  gsl_matrix * mat_3x2x_view_00,
							  gsl_matrix * mat_3x2x_view_01,
							  gsl_matrix * mat_3x2x_view_10,
							  gsl_matrix * mat_3x2x_view_11,
							  gsl_matrix * mat_3x2x_view_20,
							  gsl_matrix * mat_3x2x_view_21,
							  gsl_permutation * perm_x,
							  gsl_vector * vect_2x)
 * @param x_s : Espérance de l'état lissé
 * @param sqrt_p_s : racine de la matrice de covariance de l'état lissé courant
 * @param c_s : sqrt(P_{n+1|N}) \; K_{n|N}^T
 * @param[in] x_f : Espérance de l'état filtré
 * @param[in] sqrt_p_f : racine de la matrice de covariance de l'état filtré courant
 * @param[in] x_p_ : Espérance de l'état prédit suivant
 * @param[in] sqrt_p_p_ : racine de la covariance de l'état prédit suivant
 * @param[in] x_s_ : Espérance de l'état lissé suivant
 * @param[in]  sqrt_p_s_ : racine de la matrice de covariance de l'état lissé suivant
 * @param[in] f2_xx : Fxx - Qxy.Qyy?¹.Fyx
 * @param[in] sqrt_q2_xx : racine de la matrice de covariance du bruit de process réduit
 * @param mat_xx : matrice de taille (x.x) préallouée
 * @param mat_3x2x : matrice de taille (3x.2x) préallouée
 * @param mat_3x2x_view_00 : vue sur la matrice mat_3x2x allant de (0,0) à (n_x, n_x)
 * @param mat_3x2x_view_01 : vue sur la matrice mat_3x2x allant de (0,n_x) à (n_x,2 n_x)
 * @param mat_3x2x_view_10 : vue sur la matrice mat_3x2x allant de (n_x,0) à (2n_x, n_x)
 * @param mat_3x2x_view_11 : vue sur la matrice mat_3x2x allant de (n_x,n_x) à (2n_x, 2n_x)
 * @param mat_3x2x_view_20 : vue sur la matrice mat_3x2x allant de (2n_x,0) à (3n_x, n_x)
 * @param mat_3x2x_view_21 : vue sur la matrice mat_3x2x allant de (2n_x,n_x) à (3n_x, 2n_x)
 * @param perm_x : permutation de taille (x) préallouée
 * @param vect_2x : vecteur de taille (2x) préalloué
 * @brief
 Cette fonction effectue l'étape de lissage dans le filtre de Kalman triplet.
**/
void tkalman_do_smoothing_0(gsl_vector * t_s,
					        gsl_matrix * sqrt_q_s,
					        gsl_matrix * c_s,
					        const gsl_vector * t_f,
					        const gsl_matrix * sqrt_q_f,
					        const gsl_vector * x_p_,
					        const gsl_matrix * sqrt_p_p_,
					        const gsl_vector * x_s_,
					        const gsl_matrix * sqrt_p_s_,
					        const gsl_matrix * f2_x_,
					        const gsl_matrix * sqrt_q2_xx,
					        gsl_matrix * mat_tx,
					        gsl_matrix * mat_2xpt_xpt,
					        gsl_matrix * mat_2xpt_xpt_view_00,
							gsl_matrix * mat_2xpt_xpt_view_01,
							gsl_matrix * mat_2xpt_xpt_view_10,
							gsl_matrix * mat_2xpt_xpt_view_11,
							gsl_matrix * mat_2xpt_xpt_view_20,
							gsl_matrix * mat_2xpt_xpt_view_21,
							gsl_permutation * perm_x,
							gsl_vector * vect_xpt)
{

	//Gain
	tkalman_get_smoothing_gain_0( mat_tx,
								  sqrt_q_f,
								  sqrt_p_p_,
								  f2_x_,
								  mat_2xpt_xpt_view_00,
								  mat_2xpt_xpt_view_10,
								  perm_x);
	//Racine des matrices de covariance
	tkalman_get_sqrt_p_s_and_c_s( sqrt_q_s,
								  c_s,
								  sqrt_q_f,
							      sqrt_p_s_,
								  f2_x_,
								  sqrt_q2_xx,
								  mat_tx,
								  mat_2xpt_xpt,
								  mat_2xpt_xpt_view_00,
								  mat_2xpt_xpt_view_01,
							   	  mat_2xpt_xpt_view_10,
								  mat_2xpt_xpt_view_11,
								  mat_2xpt_xpt_view_20,
								  mat_2xpt_xpt_view_21,
								  vect_xpt );
	  
	//Espérance
	tkalman_get_x_s( t_s,
					 t_f,
					 x_p_,
					 x_s_,
					 mat_tx);

}




//Objet de prédiction
/**@fn tkalman_smoothing :: tkalman_smoothing(const gsl_matrix * f2_x,
											const gsl_matrix * sqrt_q2_xx)
* @param f2_x : [F2xx, F2xy]
* @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
* @param q2_xy : Qxy.Qyy
* @brief
* Ce constructeur alloue les variables temp. de la prédiction du filtre de Kalman couple.
*/
tkalman_smoothing :: tkalman_smoothing ( const gsl_matrix * f2_x,
										   const gsl_matrix * sqrt_q2_xx)
{
	tkalman_smoothing :: initialize();
	tkalman_smoothing :: setup (f2_x,
								 sqrt_q2_xx);
}


/**@fn int tkalman_smoothing :: setup(const gsl_matrix * f2_x,
									   const gsl_matrix * sqrt_q2_xx)
 * @param f2_x : [F2xx, F2xy]
 * @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
 * @param q2_xy : Qxy.Qyy
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction permet de modifier les paramètres de la prédiction (size_x et size_y)
**/
int tkalman_smoothing :: setup (const gsl_matrix * f2_x,
								 const gsl_matrix * sqrt_q2_xx)
	
{	
	unsigned int size_x = f2_x->size1,
				 size_y = f2_x->size2 - size_x;
	
	if (size_x != _size_x || size_y != _size_y)
	{
		tkalman_smoothing :: free();
		tkalman_smoothing :: initialize();
		if (tkalman_smoothing :: set_params(f2_x, sqrt_q2_xx))
		{
			tkalman_smoothing :: free();
			tkalman_smoothing :: initialize();
			return 1;
		}
	
	
		_size_x = size_x;
		_size_y = size_y;
		_size_t = size_x + size_y;
		if ( tkalman_smoothing :: alloc() )
		{
			tkalman_smoothing :: free();
			tkalman_smoothing :: initialize();
			return 1;
		}
		else
			tkalman_smoothing :: create_views();
	}	
	return 0;
}

/**@fn tkalman_smoothing :: ~ tkalman_smoothing()
 * @brief
 * Destructeur
 */
tkalman_smoothing :: ~tkalman_smoothing()
{
	tkalman_smoothing :: free();
	tkalman_smoothing :: initialize();
}

/**@fn bool tkalman_smoothing :: operator!() const;
 * @return
 * - 0 si l'objet est correctement alloué
 * - 1 sinon
 * @brief
 * Check de l'objet.
 */
bool tkalman_smoothing :: operator!() const
{
	return (! (perm_x && vect_xpt && mat_tx && mat_2xpt_xpt && _f2_x && _sqrt_q2_xx) );
}

/**@fn void tkalman_smoothing :: free();
 * @brief
 * Cette fonction désalloue la mémoire utilisée par les variables tmp.
 */
void tkalman_smoothing :: free()
{
	if (perm_x)
		gsl_permutation_free(perm_x);
	if (vect_xpt)
		gsl_vector_free(vect_xpt);
	if (mat_tx)
		gsl_matrix_free(mat_tx);
	if (mat_2xpt_xpt)
		gsl_matrix_free(mat_2xpt_xpt);
}

/**@fn int tkalman_smoothing :: alloc();
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction alloue la mémoire utilisée par les variables tmp.
 */
int tkalman_smoothing :: alloc()
{
	if (!perm_x)
		perm_x = gsl_permutation_alloc(_size_x);
	if (!vect_xpt)
		vect_xpt = gsl_vector_alloc(_size_x + _size_t);
	if (!mat_tx)
		mat_tx = gsl_matrix_alloc(_size_t, _size_x);
	if (!mat_2xpt_xpt)
		mat_2xpt_xpt = gsl_matrix_alloc(2 * _size_x + _size_t, _size_x + _size_t); 
		
	return (! (perm_x && vect_xpt && mat_tx && mat_2xpt_xpt) );
}

/**@fn void tkalman_smoothing :: create_views();
 * @brief
 * Cette méthode crée les vues.
 */
void tkalman_smoothing :: create_views()
{
	gsl_matrix_view view;
	view = gsl_matrix_submatrix(mat_tx,
								0,
								0,
								_size_x,
								_size_x);
	mat_xx = view.matrix;
	
	view = gsl_matrix_submatrix(mat_2xpt_xpt,
								0,
								0,
								3 * _size_x,
								2 * _size_x);
	mat_3x2x = view.matrix;
	//00
	view = gsl_matrix_submatrix(mat_2xpt_xpt,
								0,
								0,
								_size_x,
								_size_x);
	mat_3x2x_view_00 = view.matrix;
	mat_2xpt_xpt_view_00 = view.matrix;
	//01
	view = gsl_matrix_submatrix(mat_2xpt_xpt,
								0,
								_size_x,
								_size_x,
								_size_t);
	mat_2xpt_xpt_view_01 = view.matrix;
	
	view = gsl_matrix_submatrix(mat_2xpt_xpt,
								0,
								_size_x,
								_size_x,
								_size_x);
	mat_3x2x_view_01 = view.matrix;
	//10
	view = gsl_matrix_submatrix(mat_2xpt_xpt,
								_size_x,
								0,
								_size_t,
								_size_x);
	mat_2xpt_xpt_view_10 = view.matrix;
	
	view = gsl_matrix_submatrix(mat_2xpt_xpt,
								_size_x,
								0,
								_size_x,
								_size_x);
	mat_3x2x_view_10 = view.matrix;
	//11
	view = gsl_matrix_submatrix(mat_2xpt_xpt,
								_size_x,
								_size_x,
								_size_t,
								_size_t);
	mat_2xpt_xpt_view_11 = view.matrix;
	
	view = gsl_matrix_submatrix(mat_2xpt_xpt,
								_size_x,
								_size_x,
								_size_x,
								_size_x);
	mat_3x2x_view_11 = view.matrix;
	//21
	view = gsl_matrix_submatrix(mat_2xpt_xpt,
								_size_x + _size_t,
								_size_x,
								_size_x,
								_size_t);
	mat_2xpt_xpt_view_21 = view.matrix;
	
	view = gsl_matrix_submatrix(mat_2xpt_xpt,
								2 * _size_x,
								_size_x,
								_size_x,
								_size_x);
	mat_3x2x_view_21 = view.matrix;
	
	//Vect
	gsl_vector_view view42 = gsl_vector_subvector(vect_xpt, 0, 2 * _size_x);
	vect_2x = view42.vector;
}

		
/**@fn void tkalman_smoothing :: initialize();
 * @brief
 * Cette fonction met les variables tmp à NULL.
 */
void tkalman_smoothing :: initialize()
{
	perm_x = NULL;
	vect_xpt = NULL;
	mat_tx = NULL;
	mat_2xpt_xpt = NULL;
	_f2_x = NULL;
	_sqrt_q2_xx = NULL;
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
}

/**@fn virtual int tkalman_smoothing :: set_params(const gsl_matrix * f2_x,
													const gsl_matrix * sqrt_q2_xx)
 * @param f2_x : [F2xx, F2xy]
 * @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
 * @param q2_xy : Qxy.Qyy
 * @brief
 * Cette méthode change les paramètres de la prédiction.
 */
int tkalman_smoothing :: set_params(const gsl_matrix * f2_x,
									 const gsl_matrix * sqrt_q2_xx)
{

	if (!f2_x || !sqrt_q2_xx)
		return 1;
	unsigned int size_x = f2_x->size1,
				 size_y = f2_x->size2 - size_x;
	_f2_x = f2_x;
	_sqrt_q2_xx = sqrt_q2_xx;
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














