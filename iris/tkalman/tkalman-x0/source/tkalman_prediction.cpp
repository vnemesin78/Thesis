/**
 * @file tkalman_prediction.cpp
 * @author Valérian Némesin
 * @date 4/10/2011
 * @brief
 Ce fichier contient le code source des fonctions nécessaire à la prédiction dans l'algorithme du filtre de Kalman triplet.
**/
#include "tkalman_filtering.hpp"
/**@fn void tkalman_original_get_x_p(gsl_vector * x_p,
									 const gsl_vector * _x_f,
									 const gsl_vector * __y,
									 const gsl_vector * _y,
									 const gsl_matrix * f2_xx,
									 const gsl_matrix * f2_xy,
									 const gsl_matrix * q2_xy)
 * @param x_p : espérance de l'état prédit
 * @param[in] _x_f : espérance de l'état filtré précédent
 * @param[in] __y : observation (n - 2)
 * @param[in] _y : observation précédente
 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
 * @param[in] f2_xy : F_2^{x,y} = F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y}
 * @param[in] q2_xy : Q_2^{x,y} = Q^{x,y} \: [Q^{y,y}]^{-1}
 * @brief
 Cette fonction calcule l'espérance de l'état prédit à partir de la formule : 
 * \hat{x_{n|n-1}} = F_2^{x,x} \: \hat{x_{n - 1 |n-1}} + Q_2^{x,y} \; y_{n - 2} + F_2^{x,y} \: y_{n - 1}
 */
void tkalman_original_get_x_p(gsl_vector * x_p,
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

/**@fn void tkalman_original_get_p_p(gsl_matrix * p_p,
									 const gsl_matrix * _p_f,
									 const gsl_matrix * f2_xx,
									 const gsl_matrix * q2_xx,
									 gsl_matrix * mat_xx)
 * @param p_p : Matrice de covariance de l'état prédit.
 * @param[in] _p_f : Matrice de covariance de l'état filtré précédent
 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
 * @param[in] q2_xx : Matrice de covariance du bruit de process réduit
 * @param mat_xx : Matrice de taille (x.x) préallouée.
 * @brief
 Cette fonction calcule la matrice de covariance de l'état prédit à partir de la formule : 
 * P_{n|n-1} = Q_2^{x,x} + F_2^{x,x} \: P_{n-1|n-1} \: [F_2^{x,x}]^T 
 */
void tkalman_original_get_p_p(gsl_matrix * p_p,
							  const gsl_matrix * _p_f,
							  const gsl_matrix * f2_xx,
							  const gsl_matrix * q2_xx,
							  gsl_matrix * mat_xx)
{
	//F2_xx . P_f
	gsl_blas_dgemm (CblasNoTrans, 
					CblasNoTrans, 
					1.0, 
					f2_xx,
					_p_f, 
					0.0, 
					mat_xx);

	//Recopie de Q2xx
	gsl_matrix_memcpy(p_p, q2_xx);
	
	//Q2xx + F2_xx . P_f . F2_xx^T
	gsl_blas_dgemm (CblasNoTrans, 
					CblasTrans, 
					1.0, 
					mat_xx, 
					f2_xx,
					1.0, 
					p_p);
	
}

/**@fn void tkalman_original_prediction(gsl_vector * x_p,
										gsl_matrix * p_p,
										const gsl_vector * _x_f,
										const gsl_matrix * _p_f,
										const gsl_vector * __y,
										const gsl_vector * _y,
										const gsl_matrix * f2_xx,
										const gsl_matrix * f2_xy,
										const gsl_matrix * q2_xx,
										const gsl_matrix * q2_xy,
										gsl_matrix * mat_xx)
 * @param x_p : espérance de l'état prédit \hat{x_{n|n-1}} = F_2^{x,x} \: \hat{x_{n - 1 |n-1}} + Q_2^{x,y} \; y_{n - 2} + F_2^{x,y} \: y_{n - 1}
 * @param p_p : Matrice de covariance de l'état prédit. P_{n|n-1} = Q_2^{x,x} + F_2^{x,x} \: P_{n-1|n-1} \: [F_2^{x,x}]^T 
 * @param[in] _x_f : espérance de l'état filtré précédent
 * @param[in] _p_f : Matrice de covariance de l'état filtré précédent
 * @param[in] __y : observation (n - 2)
 * @param[in] _y : observation précédente
 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
 * @param[in] f2_xy : F_2^{x,y} = F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y}
 * @param[in] q2_xx : Matrice de covariance du bruit de process réduit
 * @param[in] q2_xy : Q_2^{x,y} = Q^{x,y} \: [Q^{y,y}]^{-1}
 * @brief
 Cette fonction effectue la prédiction dans le filtre de Kalman triplet(Voir les arguments pour les forumules.).
 */
void tkalman_original_prediction(gsl_vector * x_p,
								 gsl_matrix * p_p,
								 const gsl_vector * _x_f,
								 const gsl_matrix * _p_f,
							  	 const gsl_vector * __y,
							 	 const gsl_vector * _y,
							 	 const gsl_matrix * f2_xx,
							  	 const gsl_matrix * f2_xy,
							  	 const gsl_matrix * q2_xx,
								 const gsl_matrix * q2_xy,
								 gsl_matrix * mat_xx)
{
	tkalman_original_get_x_p(x_p,
							 _x_f,
							 __y,
							 _y,
							 f2_xx,
							 f2_xy,
							 q2_xy);
							 
	tkalman_original_get_p_p(p_p,
							 _p_f,
							 f2_xx,
							 q2_xx,
							 mat_xx);			 
}

/**@fn void tkalman_robust_get_sqrt_p_p(gsl_matrix * sqrt_p_p,
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
void tkalman_robust_get_sqrt_p_p(gsl_matrix * sqrt_p_p,
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


/**@fn void tkalman_robust_prediction ( gsl_vector * x_p,
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
void tkalman_robust_prediction ( gsl_vector * x_p,
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
{
	//Espérance
	tkalman_original_get_x_p ( x_p,
							   _x_f,
							   __y,
							   _y,
							   f2_xx,
							   f2_xy,
							   q2_xy);
							   
	//Racine de la matrice de covariance				 
	tkalman_robust_get_sqrt_p_p ( sqrt_p_p,
							      _sqrt_p_f,
							      f2_xx,
							      sqrt_q2_xx,
							      mat_2xx,
							      mat_2xx_view_00,
							      mat_2xx_view_10,
							      vect_x);


}
