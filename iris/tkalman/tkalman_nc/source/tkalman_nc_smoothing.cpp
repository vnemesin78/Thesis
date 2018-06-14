#include "tkalman_nc_smoothing.hpp"
/**@fn void tkalman_nc_get_x_s(gsl_vector * x_s,
							   const gsl_vector * x_f,
							   const gsl_vector * x_p_,
							   const gsl_vector * x_s_,
							   const gsl_matrix * gain)
 * @param x_s : \hat{x}_{n|N}, espérance de l'état lissé
 * @param[in] x_f : \hat{x}_{n|n}, espérance de l'état filtré
 * @param[in] x_p_ : \hat{x}_{n + 1|n}, espérance de l'état prédit suivant
 * @param[in] x_s_ : \hat{x}_{n + 1|N}, espérance de l'état lissé suivant
 * @param[in] gain : K_{n|N}, gain de lissage
 * @brief
 * Cette fonction calcule l'espérance de l'état lissé selon la formule \n
 * \hat{x}_{n|N} = \hat{x}_{n|n} + K_{n|N} (\hat{x}_{n + 1|N} - \hat{x}_{n+1|n})
 * .
 */
void tkalman_nc_get_x_s(gsl_vector * x_s,
						const gsl_vector * x_f,
						const gsl_vector * x_p_,
						const gsl_vector * x_s_,
						const gsl_matrix * gain)
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

/**@fn void tkalman_nc_get_t_0_s(gsl_vector * t_0_s,
							     const gsl_vector * t_0_f,
								 const gsl_vector * x_1_p,
							     const gsl_vector * x_1_s,
							     const gsl_matrix * gain)
 * @param t_0_s : \hat{t}_{0|N}, espérance de l'état lissé
 * @param[in] t_0_f : \hat{t}_{0|n}, espérance de l'état filtré
 * @param[in] x_1_p : \hat{x}_{1|0}, espérance de l'état prédit suivant
 * @param[in] x_1_s : \hat{x}_{1|N}, espérance de l'état lissé suivant
 * @param[in] gain : K^t_{0|N}, gain de lissage
 * @brief
 * Cette fonction calcule l'espérance de l'état lissé selon la formule \n
 * \hat{t}_{0|N} = \hat{t}_{0|n} + K^t_{0|N} (\hat{x}_{1|N} - \hat{x}_{1|0})
 */
void tkalman_nc_get_t_0_s(gsl_vector * x_s,
						  const gsl_vector * x_f,
						  const gsl_vector * x_p_,
						  const gsl_vector * x_s_,
						  const gsl_matrix * gain)
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
/**@fn void tkalman_nc_get_smoothing_gain_0(gsl_matrix * s_gain,
										 const gsl_matrix * sqrt_q_0_f,
										 const gsl_matrix * sqrt_p_1_p,
										 const gsl_matrix * f_xt,
										 gsl_matrix * mat_xx,
										 gsl_matrix * mat_tx,
										 gsl_permutation * perm_x)
 * @param s_gain : K^t_{0|N}, gain de Kalman triplet pour le lissage
 * @param[in] sqrt_q_0_f :[Q_{0|0}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
 * @param[in] sqrt_p_1_p : [P_{1|0}]^{\frac{1}{2}}, racine de la covariance de l'état prédit suivant
 * @param[in] f_xt : F^{x,t}
 * @param mat_xx : matrice de taille (n_x.n_x) préallouée
 * @param mat_tx : matrice de taille (n_t.n_x) préallouée
 * @param perm_x : permutation de taille (n_x) préallouée
 * @brief
 * Cette fonction calcule le gain de lissage du filtre de Kalman
 * triplet : \newline
 * K^t_{0|N} = Q_{0|0} [F^{x,t}]^T P_{1|0}^{-1}  
 */
void tkalman_nc_get_smoothing_gain_0(gsl_matrix * s_gain,
									 const gsl_matrix * sqrt_q_0_f,
									 const gsl_matrix * sqrt_p_1_p,
									 const gsl_matrix * f_xt,
									 gsl_matrix * mat_xx,
									 gsl_matrix * mat_tx,
									 gsl_permutation * perm_x)
{
	//Q F2_x_^T
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   sqrt_q_0_f, //tt
  				   f_xt,   //tx
				   0.0,
				   mat_tx);
	gsl_blas_dgemm(CblasTrans,
				   CblasNoTrans,
				   1.0,
				   sqrt_q_0_f, //tt
  				   mat_tx,   //tx
				   0.0,
				   s_gain);     
	//P_p - 1
	gsl_permutation_init(perm_x);
	gsl_linalg_LU_invert(sqrt_p_1_p,
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

/**@fn void tkalman_nc_get_smoothing_gain(gsl_matrix * s_gain,
									   const gsl_matrix * sqrt_p_f,
									   const gsl_matrix * sqrt_p_p_,
									   const gsl_matrix * f_xx,
									   gsl_matrix * mat_xx,
									   gsl_permutation * perm_x)
 * @param s_gain : K_{n|N}, gain de Kalman triplet pour le lissage
 * @param[in] sqrt_p_f :[P_{n|n}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
 * @param[in] sqrt_p_p_ : [P_{n + 1|n}]^{\frac{1}{2}}, racine de la covariance de l'état prédit suivant
 * @param[in] f_xx : F^{x,x}
 * @param mat_xx : matrice de taille (n_x.n_x) préallouée
 * @param perm_x : permutation de taille (n_x) préallouée
 * @brief
 * Cette fonction calcule le gain de lissage du filtre de Kalman
 * triplet : \newline
 * K_{n|N} = P_{n|n} [F^{x,x}]^T P_{n + 1|n} 
 */
void tkalman_nc_get_smoothing_gain(gsl_matrix * s_gain,
								   const gsl_matrix * sqrt_p_f,
								   const gsl_matrix * sqrt_p_p_,
								   const gsl_matrix * f_xx,
								   gsl_matrix * mat_xx,
								   gsl_permutation * perm_x)
{
	
	//P_p - 1
	gsl_permutation_init(perm_x);
	gsl_linalg_LU_invert(sqrt_p_p_,
					     perm_x,
						 s_gain);
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   s_gain, //tt
  				   s_gain,   //tx
				   0.0,
				   mat_xx);
	
	//F2_x^T P_p^-1
	gsl_blas_dgemm(CblasTrans,
				   CblasNoTrans,
				   1.0,
				   f_xx, //tt
  				   mat_xx,   //tx
				   0.0,
				   s_gain);
	
	gsl_blas_dgemm(CblasNoTrans,
				   CblasNoTrans,
				   1.0,
				   sqrt_p_f, //tt
  				   s_gain,   //tx
				   0.0,
				   mat_xx);
				   
	gsl_blas_dgemm(CblasTrans,
				   CblasNoTrans,
				   1.0,
				   sqrt_p_f, //tt
  				   mat_xx,   //tx
				   0.0,
				   s_gain);    
				    
}

/**@fn void tkalman_nc_get_sqrt_p_s_and_c_s(gsl_matrix * sqrt_p_s,
							                gsl_matrix * c_s,
											const gsl_matrix * sqrt_p_f,
											const gsl_matrix * sqrt_p_s_,
											const gsl_matrix * f_xx,
											const gsl_matrix * sqrt_q_xx,
											const gsl_matrix * s_gain,
											gsl_matrix * mat_3x2x,
											gsl_matrix * mat_3x2x_view_00,
											gsl_matrix * mat_3x2x_view_01,
											gsl_matrix * mat_3x2x_view_10,
											gsl_matrix * mat_3x2x_view_11,
											gsl_matrix * mat_3x2x_view_20,
											gsl_matrix * mat_3x2x_view_21,
											gsl_vector * vect_2x)
 * @param sqrt_p_s : [P_{n|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant
 * @param c_s : [P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
 * @param[in] sqrt_p_f :[P_{n|n}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
 * @param[in]  sqrt_p_s_ : [P_{n + 1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant
 * @param[in] f_xx : F^{x,x}
 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
 * @param[in] s_gain :  K_{n|N}, gain de Kalman triplet pour le lissage
 * @param mat_3x2x : matrice de taille (3n_x.2n_x) préallouée
 * @param mat_3x2x_view_00 : vue sur la matrice mat_3x2x allant de (0,0) à (n_x - 1 , n_x - 1)
 * @param mat_3x2x_view_01 : vue sur la matrice mat_3x2x allant de (0,n_x) à (n_x - 1, 2 n_x - 1)
 * @param mat_3x2x_view_10 : vue sur la matrice mat_3x2x allant de (n_x,0) à (2n_x - 1 , n_x - 1)
 * @param mat_3x2x_view_11 : vue sur la matrice mat_3x2x allant de (n_x,n_x) à (2n_x - 1, 2n_x - 1)
 * @param mat_3x2x_view_20 : vue sur la matrice mat_3x2x allant de (2n_x,0) à (3n_x - 1, n_x - 1)
 * @param mat_3x2x_view_21 : vue sur la matrice mat_3x2x allant de (2n_x,n_x) à (3n_x - 1, 2n_x - 1)
 * @param vect_2x : vecteur de taille (2n_x) préalloué
 * @brief
 Cette fonction calcule la racine de la matrice de covariance de l'état lissé courant et sqrt(P_{n+1|N}) \; K_{n|N}^T.
 * On calcule la matrice
 * \begin{pmatrix}
 [Q^{x,x}]^{\frac{1}{2}}	&	0 \newline
 [P_{n|n}]^{\frac{1}{2}}	&	[ F^{x,x}]^T \newline
 0							&	[P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
 * \end{pmatrix}
 * On effectue la décomposition QR :
 * \begin{pmatrix}
[P_{n + 1|N}]^{\frac{1}{2}}	&	0 \newline
 0							&	[P_{n|N}]^{\frac{1}{2}} K_{n|N}^T
 0							&	0
 * \end{pmatrix}
**/
void tkalman_nc_get_sqrt_p_s_and_c_s(gsl_matrix * sqrt_p_s,
									 gsl_matrix * c_s,
									 const gsl_matrix * sqrt_p_f,
									 const gsl_matrix * sqrt_p_s_,
									 const gsl_matrix * f_xx,
									 const gsl_matrix * sqrt_q_xx,
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
						  sqrt_q_xx);
		gsl_matrix_set_zero(mat_3x2x_view_01);
		gsl_blas_dgemm(CblasNoTrans,
					   CblasTrans,
					   1.0,
					   sqrt_p_f,
					   f_xx,
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

/**@fn void tkalman_nc_get_sqrt_q_0_s_and_c_s(gsl_matrix * sqrt_q_0_s,
											  gsl_matrix * c_s,
											  const gsl_matrix * sqrt_q_0_f,
											  const gsl_matrix * sqrt_p_1_s,
											  const gsl_matrix * f_xt,
											  const gsl_matrix * sqrt_q_xx,
											  const gsl_matrix * s_gain,
											  gsl_matrix * mat_2xpt_xpt,
											  gsl_matrix * mat_2xpt_xpt_view_00,
											  gsl_matrix * mat_2xpt_xpt_view_01,
											  gsl_matrix * mat_2xpt_xpt_view_10,
											  gsl_matrix * mat_2xpt_xpt_view_11,
											  gsl_matrix * mat_2xpt_xpt_view_20,
											  gsl_matrix * mat_2xpt_xpt_view_21,
											  gsl_vector * vect_xpt)
 * @param sqrt_q_0_s : [Q_{0|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant
 * @param c_s : [P_{1|N}]^{\frac{1}{2}} K^t_{0|N}^T (/!\ Spécial)
 * @param[in] sqrt_q_0_f :[Q_{0|0}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
 * @param[in] sqrt_p_1_s : [P_{1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant
 * @param[in] f_xt : F^{x,t}
 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
 * @param[in] s_gain :  K_{n|N}, gain de Kalman triplet pour le lissage
 * @param mat_2xpt_xpt : matrice de taille (2n_x + n_t, n_x + n_t) préallouée
 * @param mat_2xpt_xpt_view_00 : vue sur la matrice mat_2xpt_xpt allant de (0,0) à (n_x - 1 , n_x - 1)
 * @param mat_2xpt_xpt_view_01 : vue sur la matrice mat_2xpt_xpt allant de (0,n_x) à (n_x - 1, n_x + n_t - 1)
 * @param mat_2xpt_xpt_view_10 : vue sur la matrice mat_2xpt_xpt allant de (n_x,0) à (n_x + n_t - 1 , n_x - 1)
 * @param mat_2xpt_xpt_view_11 : vue sur la matrice mat_2xpt_xpt allant de (n_x,n_x) à (n_x + n_t - 1, n_x + n_t - 1)
 * @param mat_2xpt_xpt_view_20 : vue sur la matrice mat_2xpt_xpt allant de (n_x + n_t,0) à (2n_x + n_t - 1, n_x - 1)
 * @param mat_2xpt_xpt_view_21 : vue sur la matrice mat_2xpt_xpt allant de (n_x +n_t,n_x) à (2n_x + n_t - 1, n_x + n_t - 1)
 * @param vect_xpt : vecteur de taille (n_x + n_t) préalloué
 * @brief
 Cette fonction calcule la racine de la matrice de covariance de l'état lissé courant et sqrt(P_{n+1|N}) \; K_{n|N}^T.
 * On calcule la matrice
 * \begin{pmatrix}
 [Q^{x,x}]^{\frac{1}{2}}	&	0 \newline
 [P_{n|n}]^{\frac{1}{2}}	&	[ F^{x,x}]^T \newline
 0							&	[P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
 * \end{pmatrix}
 * On effectue la décomposition QR
 * \begin{pmatrix}
[P_{n + 1|N}]^{\frac{1}{2}}	&	0 \newline
 0							&	[P_{n|N}]^{\frac{1}{2}} K_{n|N}^T
 0							&	0
 * \end{pmatrix}
**/
void tkalman_nc_get_sqrt_q_0_s_and_c_s(gsl_matrix * sqrt_p_s,
										 gsl_matrix * c_s,
										 const gsl_matrix * sqrt_p_f,
										 const gsl_matrix * sqrt_p_s_,
										 const gsl_matrix * f_xx,
										 const gsl_matrix * sqrt_q_xx,
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
						  sqrt_q_xx);
		gsl_matrix_set_zero(mat_3x2x_view_01);
		gsl_blas_dgemm(CblasNoTrans,
					   CblasTrans,
					   1.0,
					   sqrt_p_f,
					   f_xx,
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
/**@fn void tkalman_nc_do_smoothing ( gsl_vector * x_s,
									  gsl_matrix * sqrt_p_s,
									  gsl_matrix * c_s,
									  const gsl_vector * x_f,
									  const gsl_matrix * sqrt_p_f,
									  const gsl_vector * x_p_,
									  const gsl_matrix * sqrt_p_p_,
									  const gsl_vector * x_s_,
									  const gsl_matrix * sqrt_p_s_,
									  const gsl_matrix * f_xx,
									  const gsl_matrix * sqrt_q_xx,
									  gsl_matrix * mat_xx,
									  gsl_matrix * mat_3x2x,
									  gsl_matrix * mat_3x2x_view_00,
									  gsl_matrix * mat_3x2x_view_01,
									  gsl_matrix * mat_3x2x_view_10,
									  gsl_matrix * mat_3x2x_view_11,
									  gsl_matrix * mat_3x2x_view_20,
									  gsl_matrix * mat_3x2x_view_21,
									  gsl_permutation * perm_x,
									  gsl_vector * vect_2x )
 * @param x_s : \hat{x}_{n|N}, espérance de l'état lissé
 * @param sqrt_p_s : [P_{n|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant
 * @param c_s : [P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
 * @param[in] x_f : \hat{x}_{n|n}, espérance de l'état filtré
 * @param[in] sqrt_p_f :[P_{n|n}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
 * @param[in] x_p_ : \hat{x}_{n + 1|n}, espérance de l'état prédit suivant
 * @param[in] sqrt_p_p_ : [P_{n + 1|n}]^{\frac{1}{2}}, racine de la covariance de l'état prédit suivant
 * @param[in] x_s_ : \hat{x}_{n + 1|N}, espérance de l'état lissé suivant
 * @param[in] sqrt_p_s_ : [P_{n + 1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant
 * @param[in] f_xx : F^{x,x}
 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
 * @param mat_xx : matrice de taille (n_x.n_x) préallouée
 * @param mat_3x2x : matrice de taille (3n_x.2n_x) préallouée
 * @param mat_3x2x_view_00 : vue sur la matrice mat_3x2x allant de (0,0) à (n_x - 1 , n_x - 1)
 * @param mat_3x2x_view_01 : vue sur la matrice mat_3x2x allant de (0,n_x) à (n_x - 1, 2 n_x - 1)
 * @param mat_3x2x_view_10 : vue sur la matrice mat_3x2x allant de (n_x,0) à (2n_x - 1 , n_x - 1)
 * @param mat_3x2x_view_11 : vue sur la matrice mat_3x2x allant de (n_x,n_x) à (2n_x - 1, 2n_x - 1)
 * @param mat_3x2x_view_20 : vue sur la matrice mat_3x2x allant de (2n_x,0) à (3n_x - 1, n_x - 1)
 * @param mat_3x2x_view_21 : vue sur la matrice mat_3x2x allant de (2n_x,n_x) à (3n_x - 1, 2n_x - 1)
 * @param perm_x : permutation de taille (n_x) préallouée
 * @param vect_2x : vecteur de taille (2n_x) préalloué
 * @brief
 * Cette fonction effectue le lissage dans le filtre de Kalman triplet. Dans un premier temps, elle calcule le gain de lissage via la fonction @fn tkalman_nc_get_smoothing_gain. Puis dans un second temps, elle calcule la racine de la covariance de l'état lissé via @fn void tkalman_nc_get_sqrt_p_s_and_c_s. Finalement, elle calcule l'espérance de l'état lissé via @fn tkalman_nc_get_x_s .
 */
void tkalman_nc_do_smoothing ( gsl_vector * x_s,
							   gsl_matrix * sqrt_p_s,
							   gsl_matrix * c_s,
							   const gsl_vector * x_f,
							   const gsl_matrix * sqrt_p_f,
							   const gsl_vector * x_p_,
							   const gsl_matrix * sqrt_p_p_,
							   const gsl_vector * x_s_,
							   const gsl_matrix * sqrt_p_s_,
							   const gsl_matrix * f_xx,
							   const gsl_matrix * sqrt_q_xx,
							   gsl_matrix * mat_xx,
							   gsl_matrix * mat_3x2x,
							   gsl_matrix * mat_3x2x_view_00,
							   gsl_matrix * mat_3x2x_view_01,
							   gsl_matrix * mat_3x2x_view_10,
							   gsl_matrix * mat_3x2x_view_11,
							   gsl_matrix * mat_3x2x_view_20,
							   gsl_matrix * mat_3x2x_view_21,
							   gsl_permutation * perm_x,
							   gsl_vector * vect_2x )
{
	//Gain
	tkalman_nc_get_smoothing_gain( mat_xx,
								   sqrt_p_f,
								   sqrt_p_p_,
								   f_xx,
							       mat_3x2x_view_00,
							       perm_x);
	//Racine des matrices de covariance
	tkalman_nc_get_sqrt_p_s_and_c_s( sqrt_p_s,
									 c_s,
									 sqrt_p_f,
									 sqrt_p_s_,
									 f_xx,
									 sqrt_q_xx,
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
	tkalman_nc_get_x_s( x_s,
					    x_f,
					    x_p_,
					    x_s_,
					    mat_xx);
	
}

/**@fn void tkalman_nc_do_smoothing_0 ( gsl_vector * t_0_s,
										gsl_matrix * sqrt_q_0_s,
										gsl_matrix * c_s,
										const gsl_vector * t_0_f,
										const gsl_matrix * sqrt_q_0_f,
										const gsl_vector * x_1_p,
										const gsl_matrix * sqrt_p_1_p,
										const gsl_vector * x_1_s,
										const gsl_matrix * sqrt_p_1_s,
										const gsl_matrix * f_xt,
										const gsl_matrix * sqrt_q_xx,
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
 * @param t_0_s : \hat{t}_{0|N}, espérance de l'état lissé
 * @param sqrt_q_0_s : [Q_{0|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant
 * @param c_s : [P_{1|N}]^{\frac{1}{2}} K^t_{0|N}^T (/!\ Spécial)
 * @param[in] t_0_f : \hat{t}_{0|n}, espérance de l'état filtré
 * @param[in] sqrt_q_0_f :[Q_{0|0}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
 * @param[in] x_1_p : \hat{x}_{1|0}, espérance de l'état prédit suivant
 * @param[in] sqrt_p_1_p : [P_{1|0}]^{\frac{1}{2}}, racine de la covariance de l'état prédit suivant
 * @param[in] x_1_s : \hat{x}_{1|N}, espérance de l'état lissé suivant
 * @param[in] sqrt_p_1_s : [P_{1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant
 * @param[in] f_xt : F^{x,t}
 * @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
 * @param mat_tx : matrice de taille (n_t.n_x) préallouée
 * @param mat_2xpt_xpt : matrice de taille (2n_x + n_t, n_x + n_t) préallouée
 * @param mat_2xpt_xpt_view_00 : vue sur la matrice mat_2xpt_xpt allant de (0,0) à (n_x - 1 , n_x - 1)
 * @param mat_2xpt_xpt_view_01 : vue sur la matrice mat_2xpt_xpt allant de (0,n_x) à (n_x - 1, n_x + n_t - 1)
 * @param mat_2xpt_xpt_view_10 : vue sur la matrice mat_2xpt_xpt allant de (n_x,0) à (n_x + n_t - 1 , n_x - 1)
 * @param mat_2xpt_xpt_view_11 : vue sur la matrice mat_2xpt_xpt allant de (n_x,n_x) à (n_x + n_t - 1, n_x + n_t - 1)
 * @param mat_2xpt_xpt_view_20 : vue sur la matrice mat_2xpt_xpt allant de (n_x + n_t,0) à (2n_x + n_t - 1, n_x - 1)
 * @param mat_2xpt_xpt_view_21 : vue sur la matrice mat_2xpt_xpt allant de (n_x +n_t,n_x) à (2n_x + n_t - 1, n_x + n_t - 1)
 * @param vect_xpt : vecteur de taille (n_x + n_t) préalloué
 * @param perm_x : permutation de taille (n_x) préallouée
 * @brief
 * Cette fonction effectue le lissage dans le filtre de Kalman triplet. Dans un premier temps, elle calcule le gain de lissage via la fonction @fn tkalman_nc_get_smoothing_gain_0. Puis dans un second temps, elle calcule la racine de la covariance de l'état lissé via @fn void tkalman_nc_get_sqrt_q_0_s_and_c_s. Finalement, elle calcule l'espérance de l'état lissé via @fn tkalman_nc_get_t_0_s .
 */
void tkalman_nc_do_smoothing_0 ( gsl_vector * t_0_s,
							     gsl_matrix * sqrt_q_0_s,
								 gsl_matrix * c_s,
								 const gsl_vector * t_0_f,
								 const gsl_matrix * sqrt_q_0_f,
								 const gsl_vector * x_1_p,
								 const gsl_matrix * sqrt_p_1_p,
								 const gsl_vector * x_1_s,
								 const gsl_matrix * sqrt_p_1_s,
								 const gsl_matrix * f_xt,
								 const gsl_matrix * sqrt_q_xx,
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
	tkalman_nc_get_smoothing_gain_0(mat_tx,
									sqrt_q_0_f,
									sqrt_p_1_p,
									f_xt,
									mat_2xpt_xpt_view_00,
									mat_2xpt_xpt_view_10,
									perm_x);
									
	//Cov.
	tkalman_nc_get_sqrt_q_0_s_and_c_s(sqrt_q_0_s,
									  c_s,
									  sqrt_q_0_f,
									  sqrt_p_1_s,
									  f_xt,
									  sqrt_q_xx,
									  mat_tx,
									  mat_2xpt_xpt,
									  mat_2xpt_xpt_view_00,
									  mat_2xpt_xpt_view_01,
									  mat_2xpt_xpt_view_10,
									  mat_2xpt_xpt_view_11,
									  mat_2xpt_xpt_view_20,
									  mat_2xpt_xpt_view_21,
									  vect_xpt);
									  
									  
	//Esp.
	tkalman_nc_get_t_0_s(t_0_s,
						 t_0_f,
						 x_1_p,
						 x_1_s,
						 mat_tx);

}


/**@fn tkalman_nc_smoothing :: tkalman_nc_smoothing(const gsl_matrix * f_xt,
												 const gsl_matrix * sqrt_q_xx) throw(exception &);
* @param[in] f_xt : F^{x,t}
* @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
* @brief
* Constructeur de la classe @class tkalman_nc_smoothing
* @throw
* Exception en cas d'erreur.
*/
tkalman_nc_smoothing :: tkalman_nc_smoothing(const gsl_matrix * f_xt,
											 const gsl_matrix * sqrt_q_xx) throw(exception &)
{
	initialize();
	try
	{
		setup(f_xt,
			  sqrt_q_xx);
	}
	catch(exception & e)
	{
		throw(e);
	}
}



/**@fn void tkalman_nc_smoothing :: setup(const gsl_matrix * f_xt,
									   const gsl_matrix * sqrt_q_xx) throw(exception &);
* @param[in] f_xt : F^{x,t}
* @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
* @brief
* Réinitialisation de la classe @class tkalman_nc_smoothing
* @throw
* Exception en cas d'erreur.
*/
void tkalman_nc_smoothing :: setup(const gsl_matrix * f_xt,
								   const gsl_matrix * sqrt_q_xx) throw(exception &)
{
	//Vérif des arguments
	if (!f_xt || !sqrt_q_xx)
		throw(invalid_argument("Fxt, sqrt(Q) or sqrt_q_xx are NULL!\n"));
	unsigned int size_t,
				 size_x,
				 size_y;
	//Dim.
	size_x = sqrt_q_xx->size1;
	size_t = f_xt->size2;
	size_y = size_t - size_x;
	
	//Vérif.;
	if (!size_y || !size_x)
		throw(invalid_argument("size_y or size_x are 0!\n"));
	
	//Modif des matrices
	if (size_x != _size_x || size_y != _size_y)
	{

		free();
		initialize();
		_size_x = size_x;
		_size_y = size_y;
		_size_t = size_t;
		_f_xt = f_xt;
		_sqrt_q_xx = sqrt_q_xx;
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
		_f_xt = f_xt;
		_sqrt_q_xx = sqrt_q_xx;

	}
	//Création des vues
	create_views();
	
}


/**@fn tkalman_nc_smoothing :: ~ tkalman_nc_smoothing();
 * @brief
 * Destructeur de la classe @class tkalman_nc_smoothing.
 */
tkalman_nc_smoothing :: ~tkalman_nc_smoothing()
{
	free();
	initialize();
}

/**@fn void tkalman_nc_smoothing :: compute_smoothing_0( gsl_vector * t_0_s,
														 gsl_matrix * sqrt_q_0_s,
														 gsl_matrix * c_s,
														 const gsl_vector * t_0_f,
														 const gsl_matrix * sqrt_q_0_f,
														 const gsl_vector * x_1_p,
														 const gsl_matrix * sqrt_p_1_p,
														 const gsl_vector * x_1_s,
														 const gsl_matrix * sqrt_p_1_s)
 * @param t_0_s : \hat{t}_{0|N}, espérance de l'état lissé
 * @param sqrt_q_0_s : [Q_{0|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant
 * @param c_s : [P_{1|N}]^{\frac{1}{2}} K^t_{0|N}^T (/!\ Spécial)
 * @param[in] t_0_f : \hat{t}_{0|n}, espérance de l'état filtré
 * @param[in] sqrt_q_0_f :[Q_{0|0}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
 * @param[in] x_1_p : \hat{x}_{1|0}, espérance de l'état prédit suivant
 * @param[in] sqrt_p_1_p : [P_{1|0}]^{\frac{1}{2}}, racine de la covariance de l'état prédit suivant
 * @param[in] x_1_s : \hat{x}_{1|N}, espérance de l'état lissé suivant
 * @param[in] sqrt_p_1_s : [P_{1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant
 * @brief
 * Cette fonction effectue le lissage dans le filtre de Kalman triplet. Dans un premier temps, elle calcule le gain de lissage via la fonction @fn tkalman_nc_get_smoothing_gain_0. Puis dans un second temps, elle calcule la racine de la covariance de l'état lissé via @fn void tkalman_nc_get_sqrt_q_0_s_and_c_s. Finalement, elle calcule l'espérance de l'état lissé via @fn tkalman_nc_get_t_0_s .
 */
void tkalman_nc_smoothing :: compute_smoothing_0( gsl_vector * t_0_s,
												  gsl_matrix * sqrt_q_0_s,
												  gsl_matrix * c_s,
												  const gsl_vector * t_0_f,
												  const gsl_matrix * sqrt_q_0_f,
												  const gsl_vector * x_1_p,
												  const gsl_matrix * sqrt_p_1_p,
												  const gsl_vector * x_1_s,
												  const gsl_matrix * sqrt_p_1_s)
{
	tkalman_nc_do_smoothing_0 ( t_0_s,
								sqrt_q_0_s,
								c_s,
								t_0_f,
								sqrt_q_0_f,
								x_1_p,
								sqrt_p_1_p,
								x_1_s,
								sqrt_p_1_s,
								_f_xt,
								_sqrt_q_xx,
								mat_tx,
								mat_2xpt_xpt,
								&mat_2xpt_xpt_view_00,
								&mat_2xpt_xpt_view_01,
								&mat_2xpt_xpt_view_10,
								&mat_2xpt_xpt_view_11,
								&mat_2xpt_xpt_view_20,
								&mat_2xpt_xpt_view_21,
								perm_x,
								vect_xpt);	
}

/**@fn void tkalman_nc_smoothing  :: compute_smoothing( gsl_vector * x_s,
														gsl_matrix * sqrt_p_s,
														gsl_matrix * c_s,
														const gsl_vector * x_f,
														const gsl_matrix * sqrt_p_f,
														const gsl_vector * x_p_,
														const gsl_matrix * sqrt_p_p_,
														const gsl_vector * x_s_,
														const gsl_matrix * sqrt_p_s_)	
 * @param x_s : \hat{x}_{n|N}, espérance de l'état lissé
 * @param sqrt_p_s : [P_{n|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé courant
 * @param c_s : [P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
 * @param[in] x_f : \hat{x}_{n|n}, espérance de l'état filtré
 * @param[in] sqrt_p_f :[P_{n|n}]^{\frac{1}{2}}, racine de la covariance de l'état filtré courant
 * @param[in] x_p_ : \hat{x}_{n + 1|n}, espérance de l'état prédit suivant
 * @param[in] sqrt_p_p_ : [P_{n + 1|n}]^{\frac{1}{2}}, racine de la covariance de l'état prédit suivant
 * @param[in] x_s_ : \hat{x}_{n + 1|N}, espérance de l'état lissé suivant
 * @param[in] sqrt_p_s_ : [P_{n + 1|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé suivant
 * @brief
 * Cette fonction effectue le lissage dans le filtre de Kalman triplet. Dans un premier temps, elle calcule le gain de lissage via la fonction @fn tkalman_nc_get_smoothing_gain. Puis dans un second temps, elle calcule la racine de la covariance de l'état lissé via @fn void tkalman_nc_get_sqrt_p_s_and_c_s. Finalement, elle calcule l'espérance de l'état lissé via @fn tkalman_nc_get_x_s .
 */
void tkalman_nc_smoothing :: compute_smoothing( gsl_vector * x_s,
												gsl_matrix * sqrt_p_s,
												gsl_matrix * c_s,
												const gsl_vector * x_f,
												const gsl_matrix * sqrt_p_f,
												const gsl_vector * x_p_,
												const gsl_matrix * sqrt_p_p_,
												const gsl_vector * x_s_,
												const gsl_matrix * sqrt_p_s_)
{
	tkalman_nc_do_smoothing ( x_s,
							  sqrt_p_s,
							  c_s,
							  x_f,
							  sqrt_p_f,
							  x_p_,
							  sqrt_p_p_,
							  x_s_,
							  sqrt_p_s_,
							  &_f_xx,
							  _sqrt_q_xx,
							  mat_xx,
							  mat_3x2x,
							  &mat_3x2x_view_00,
							  &mat_3x2x_view_01,
							  &mat_3x2x_view_10,
							  &mat_3x2x_view_11,
							  &mat_3x2x_view_20,
							  &mat_3x2x_view_21,
							  perm_x,
							  vect_2x );	
}

/**@fn void tkalman_nc_smoothing :: free();
 * @brief
 * Cette méthode libère la mémoire occupée par les attributs de l'objet.
 */
void tkalman_nc_smoothing :: free()
{
	if (mat_xx)
		gsl_matrix_free(mat_xx);
	if (mat_tx)
		gsl_matrix_free(mat_tx);
	if (mat_2xpt_xpt)
		gsl_matrix_free(mat_2xpt_xpt);
	if (mat_3x2x)
		gsl_matrix_free(mat_3x2x);
	if (perm_x)
		gsl_permutation_free(perm_x);
	if (vect_2x)
		gsl_vector_free(vect_2x);
	if (vect_xpt)
		gsl_vector_free(vect_xpt);
}

/**@fn void tkalman_nc_smoothing :: alloc() throw(exception &);
 * @brief
 * Cette méthode alloue la mémoire nécessaire pour les attributs de l'objet
 * @throw
 * Exception en cas de prob.
 */
void tkalman_nc_smoothing :: alloc() throw(exception &)
{
	if (!mat_xx)
	{
		try
		{
			mat_xx = gsl_matrix_alloc(_size_x, _size_x);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!mat_tx)
	{
		try
		{
			mat_tx = gsl_matrix_alloc(_size_t, _size_x);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!mat_2xpt_xpt)
	{
		try
		{
			mat_2xpt_xpt = gsl_matrix_alloc(2 * _size_x + _size_t, _size_x + _size_t);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!mat_3x2x)
	{
		try
		{
			mat_3x2x = gsl_matrix_alloc(3 * _size_x, 2 * _size_x);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!perm_x)
	{
		try
		{
			perm_x = gsl_permutation_alloc(_size_x);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!vect_2x)
	{
		try
		{
			vect_2x = gsl_vector_alloc(2 * _size_x);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!vect_xpt)
	{
		try
		{
			vect_xpt = gsl_vector_alloc(_size_x + _size_t);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	
}

/**@fn void tkalman_nc_smoothing :: initialize();
 * @brief
 * Cette méthode met tous les attributs de l'objet à 0.
 * 
 */
void tkalman_nc_smoothing :: initialize()
{
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
	_f_xt = 0;
	_sqrt_q_xx = 0;
	mat_xx = 0;
	mat_tx = 0;
	mat_3x2x = 0;
	mat_2xpt_xpt = 0;
	perm_x = 0;
	vect_2x = 0;
	vect_xpt = 0;
}

/**@fn void tkalman_nc_smoothing :: create_views();
 * @brief
 * Cette méthode crée les vues...
 */
void tkalman_nc_smoothing :: create_views()
{
	//Vues const.
	{
		gsl_matrix_const_view view2 = gsl_matrix_const_submatrix(_f_xt, 0, 0, _size_x, _size_x);
		_f_xx = view2.matrix;
	}
	//Vues var.
	{
		gsl_matrix_view view;
		//mat_2xpt_xpt
		{
			view = gsl_matrix_submatrix(mat_2xpt_xpt,
										0,
										0,
										_size_x,
										_size_x);
			mat_2xpt_xpt_view_00 = view.matrix;
			
			view = gsl_matrix_submatrix(mat_2xpt_xpt,
										0,
										_size_x,
										_size_x,
										_size_t);
			mat_2xpt_xpt_view_01 = view.matrix;
			
			view = gsl_matrix_submatrix(mat_2xpt_xpt,
										_size_x,
										0,
										_size_t,
										_size_x);
			mat_2xpt_xpt_view_10 = view.matrix;
			
			view = gsl_matrix_submatrix(mat_2xpt_xpt,
										_size_x,
										_size_x,
										_size_t,
										_size_t);
			mat_2xpt_xpt_view_11 = view.matrix;
			
			view = gsl_matrix_submatrix(mat_2xpt_xpt,
										_size_x + _size_t,
										0,
										_size_x,
										_size_x);
			mat_2xpt_xpt_view_20 = view.matrix;
			
			view = gsl_matrix_submatrix(mat_2xpt_xpt,
										_size_x + _size_t,
										_size_x,
										_size_x,
										_size_t);
			mat_2xpt_xpt_view_21 = view.matrix;
		}
		
		//mat_3x2x
		{
			view = gsl_matrix_submatrix(mat_3x2x,
										0,
										0,
										_size_x,
										_size_x);
			mat_3x2x_view_00 = view.matrix;
			
			view = gsl_matrix_submatrix(mat_3x2x,
										0,
										_size_x,
										_size_x,
										_size_x);
			mat_3x2x_view_01 = view.matrix;
			
			view = gsl_matrix_submatrix(mat_3x2x,
										_size_x,
										0,
										_size_x,
										_size_x);
			mat_3x2x_view_10 = view.matrix;
			
			view = gsl_matrix_submatrix(mat_3x2x,
										_size_x,
										_size_x,
										_size_x,
										_size_x);
			mat_3x2x_view_11 = view.matrix;
			
			view = gsl_matrix_submatrix(mat_3x2x,
										_size_x + _size_x,
										0,
										_size_x,
										_size_x);
			mat_3x2x_view_20 = view.matrix;
			
			view = gsl_matrix_submatrix(mat_3x2x,
										_size_x + _size_x,
										_size_x,
										_size_x,
										_size_x);
			mat_3x2x_view_21 = view.matrix;
		}
	
	}

}

/**@fn bool tkalman_nc_smoothing :: operator!() const;
 * @return
 * - 0 si l'objet est valide
 * - 1 sinon
 * 
 */			    
bool tkalman_nc_smoothing :: operator !() const
{
	return (! (_size_x && _size_y && _f_xt && _sqrt_q_xx && mat_xx && mat_tx && mat_2xpt_xpt && mat_3x2x && perm_x && vect_2x && vect_xpt) );
}

