/**
 * @file tkalman_filtering.cpp
 * @author Val�rian N�mesin
 * @date 4/10/2011
 * @brief
 Ce fichier contient le code source des fonctions n�cessaire au filtrage dans l'algorithme du filtre de Kalman triplet.
**/
#include "tkalman_filtering.hpp"
/**@fn void tkalman_original_get_innovation(gsl_vector * innovation,
											const gsl_vector * x_p,
											const gsl_vector * _y,
											const gsl_vector * y,
											const gsl_matrix * f_yx,
											const gsl_matrix * f_yy)
 * @param innovation : innovation
 * @param[in]�x_p : \hat{x}_{n | n - 1}, esp�rance de l'�tat pr�dit
 * @param[in] _y : y_{n - 1}, observation pr�c�dente
 * @param[in] y : y_n, observation
 * @param[in] f_yx : F^{y,x}, terme de la matrice d'�volution
 * @param[in] f_yy : F^{y,y}, terme de la matrice d'�volution
 * @brief
 Cette fonction calcule l'innovation selon la formule :
 \tilde{y}_n = y_{n} - F^{y,x} \: \hat{x}_{n | n - 1} - F^{y,y} \: y_{n - 1}
 */
void tkalman_original_get_innovation(gsl_vector * innovation,
									 const gsl_vector * x_p,
									 const gsl_vector * _y,
									 const gsl_vector * y,
									 const gsl_matrix * f_yx,
									 const gsl_matrix * f_yy)
{
	gsl_vector_memcpy(innovation,
					  y);

	gsl_blas_dgemv (CblasNoTrans,
					-1.0,
					f_yx,
					x_p,
					1.0,
					innovation);

	gsl_blas_dgemv (CblasNoTrans,
					-1.0,
					f_yy,
					_y,
					1.0,
					innovation);
}

/**@fn void tkalman_original_get_s(gsl_matrix * s,
								   const gsl_matrix * p_p,
								   const gsl_matrix * f_yx,
								   const gsl_matrix * q_yy,
								   gsl_matrix * mat_yx)
 * @param s : covariance de l'innovation
 * @param[in] p_p : P_{n|n-1}, covariance de l'�tat pr�dit
 * @param[in]�f_yx : F^{y,x}, terme de la matrice d'�volution
 * @param[in] q_yy : Q^{y,y}, matrice de covariance du bruit de mesure
 * @param mat_yx : matrice de taille (y.x) pr�allou�e.
 * @brief
 Cette fonction calcule la covariance de l'innovation d�crite par la forumle ci dessous:
 * Q^{y,y} + F^{y,x} \: P_{n | n - 1} \: [F^{y,x}]^T
 */
void tkalman_original_get_s(gsl_matrix * s,
							const gsl_matrix * p_p,
							const gsl_matrix * f_yx,
							const gsl_matrix * q_yy,
							gsl_matrix * mat_yx)
{
	//F_yx . p_p
	gsl_blas_dgemm (CblasNoTrans,
					CblasNoTrans,
					1.0,
					f_yx,
					p_p,
					0.0,
					mat_yx);
					
	//Recopie de Q2xx
	gsl_matrix_memcpy(s, q_yy);

	//Q2xx + F2_xx . P_f . F2_xx^T
	gsl_blas_dgemm (CblasNoTrans,
					CblasTrans,
					1.0,
					mat_yx,
					f_yx,
					1.0,
					s);
}
/**@fn void tkalman_original_get_gain(gsl_matrix * gain,
									  const gsl_matrix * p_p,
									  const gsl_matrix * s,
									  const gsl_matrix * f_yx,
									  gsl_matrix * mat_yx,
									  gsl_matrix * mat_yy)
 * @param gain : gain de Kalman
 * @param[in] p_p : P_{n|n-1}, covariance de l'�tat pr�dit
 * @param[in] s : covariance de l'innovation
 * @param[in] f_yx : F^{y,x}, terme de la matrice d'�volution
 * @param mat_yx : matrice de taille (y.x) pr�allou�e.
 * @param mat_yy : matrice de taille (y.y) pr�allou�e.
 * @brief
 Cette fonction calcule le gain de Kalman suivant la formule :
 * P_{n|n-1} \; [F^{y,x}]^T \; S_{n}^{-1}
 */
void tkalman_original_get_gain(gsl_matrix * gain,
							   const gsl_matrix * p_p,
							   const gsl_matrix * s,
							   const gsl_matrix * f_yx,
							   gsl_matrix * mat_yx,
							   gsl_matrix * mat_yy)
{
	//Inversion de s
	gsl_matrix_memcpy(mat_yy, s);
	gsl_linalg_cholesky_decomp (mat_yy);
	gsl_linalg_cholesky_invert (mat_yy);
	
	gsl_blas_dgemm (CblasNoTrans,
					CblasNoTrans,
					1.0,
					mat_yy,
					f_yx,
					0.0,
					mat_yx);
					
					
					
	gsl_blas_dgemm (CblasNoTrans,
					CblasTrans,
					1.0,
					p_p,
					mat_yx,
					0.0,
					gain);
}
/**@fn void tkalman_original_get_x_f(gsl_vector * x_f,
									 const gsl_vector * x_p,
									 const gsl_vector * innovation,
									 const gsl_matrix * gain)
 * @param x_f : esp�rance de l'�tat filtr�
 * @param[in] x_p : esp�rance de l'�tat pr�dit
 * @param[in] innovation : innovation
 * @param[in] gain : gain
 * @brief
 Cette fonction calcule l'esp�rance de l'�tat filtr� selon la formule :
 \hat{x}_{n|n} = \hat{x|n-1} + K_{n|n} \; \tilde{y}_n
 */
void tkalman_original_get_x_f(gsl_vector * x_f,
							  const gsl_vector * x_p,
							  const gsl_vector * innovation,
							  const gsl_matrix * gain)
{
	gsl_vector_memcpy(x_f, x_p);

	gsl_blas_dgemv (CblasNoTrans,
					1.0,
					gain,
					innovation,
					1.0,
					x_f);
}

/**@fn void tkalman_original_get_p_f(gsl_matrix * p_f,
									 const gsl_matrix * p_p,
									 const gsl_matrix * s,
									 const gsl_matrix * gain,
									 gsl_matrix * mat_yx)
 * @param p_f : Matrice de covariance de l'�tat filtr�
 * @param[in] p_p : matrice de covariance de l'�tat pr�dit
 * @param[in] s : matrice de covariance de l'innovarion
 * @param[in] gain : gain de Kalman
 * @param mat_yx : matrice de taille (y.x) pr�allou�e.
 * @brief
 Cette fonction calcule la matrice de covariance de l'�tat filtr� selon la formule :
 * P_{n|n}�= P_{n|n-1} + K_{n|n} \; S_{n} \: K_{n|n}^T
 */
void tkalman_original_get_p_f(gsl_matrix * p_f,
							  const gsl_matrix * p_p,
							  const gsl_matrix * s,
							  const gsl_matrix * gain,
							  gsl_matrix * mat_yx)
{

	//F_yx . p_p
	gsl_blas_dgemm (CblasNoTrans,
					CblasTrans,
					1.0,
					s,
					gain,
					0.0,
					mat_yx);

	//Recopie de Q2xx
	gsl_matrix_memcpy(p_f, p_p);

	//Q2xx + F2_xx . P_f . F2_xx^T
	gsl_blas_dgemm (CblasNoTrans,
					CblasNoTrans,
					-1.0,
					gain,
					mat_yx,
					1.0,
					p_f);
}
/**@fn void tkalman_original_filtering(gsl_vector * x_f,
									   gsl_matrix * p_f,
									   gsl_vector * innovation,
									   gsl_matrix * s,
									   const gsl_vector * x_p,
									   const gsl_matrix * p_p,
									   const gsl_vector * _y,
									   const gsl_vector * y,
									   const gsl_matrix * f_yx,
									   const gsl_matrix * f_yy,
									   const gsl_matrix * q_yy,
									   gsl_matrix * mat_xy, // Gain
									   gsl_matrix * mat_yx,
									   gsl_matrix * mat_yy)
 * @param x_f : esp�rance de l'�tat filtr�  \hat{x}_{n|n} = \hat{x|n-1} + K_{n|n} \; \tilde{y}_n
 * @param p_f : Matrice de covariance de l'�tat filtr� P_{n|n}�= P_{n|n-1} + K_{n|n} \; S_{n} \: K_{n|n}^T
 * @param innovation : innovation  \tilde{y}_n = y_{n} - F^{y,x} \: \hat{x}_{n | n - 1} - F^{y,y} \: y_{n - 1}
 * @param s : covariance de l'innovation Q^{y,y} + F^{y,x} \: P_{n | n - 1} \: [F^{y,x}]^T
 * @param[in] x_p : esp�rance de l'�tat pr�dit
 * @param[in] p_p : matrice de covariance de l'�tat pr�dit
 * @param[in] _y : y_{n - 1}, observation pr�c�dente
 * @param[in] y : y_n, observation
 * @param[in] f_yx : F^{y,x}, terme de la matrice d'�volution
 * @param[in] f_yy : F^{y,y}, terme de la matrice d'�volution
 * @param[in] q_yy : Q^{y,y}, matrice de covariance du bruit de mesure
 * @param mat_xy : Matrice de taille (x.y) pr�allou�e.
 * @param mat_yx : matrice de taille (y.x) pr�allou�e.
 * @param mat_yy : matrice de taille (y.y) pr�allou�e.
 * @brief
 Cette fonction effectue la partie filtrage du filtre de Kalman triplet.
 */
void tkalman_original_filtering(gsl_vector * x_f,
							    gsl_matrix * p_f,
							    gsl_vector * innovation,
								gsl_matrix * s,
								const gsl_vector * x_p,
								const gsl_matrix * p_p,
								const gsl_vector * _y,
								const gsl_vector * y,
								const gsl_matrix * f_yx,
								const gsl_matrix * f_yy,
								const gsl_matrix * q_yy,
								gsl_matrix * mat_xy, // Gain
								gsl_matrix * mat_yx,
								gsl_matrix * mat_yy)
{
	tkalman_original_get_innovation(innovation,
									x_p,
									_y,
									y,
									f_yx,
									f_yy);

	tkalman_original_get_s(s,
						   p_p,
						   f_yx,
						   q_yy,
						   mat_yx);

	tkalman_original_get_gain(mat_xy,
							  p_p,
							  s,
							  f_yx,
							  mat_yx,
							  mat_yy);

	tkalman_original_get_x_f(x_f,
							 x_p,
							 innovation,
							 mat_xy);

	tkalman_original_get_p_f(p_f,
							 p_p,
							 s,
							 mat_xy,
							 mat_yx);

}


/**@fn void tkalman_robust_get_sqrt_pf_sqrt_s_and_gain(gsl_matrix * sqrt_p_f,
									 gsl_matrix * sqrt_s,
												gsl_matrix * gain,
												const gsl_matrix * sqrt_p_p,
												const gsl_matrix * f_yx,
												const gsl_matrix * sqrt_q_yy,
												gsl_matrix * mat_tt,
												gsl_matrix * mat_tt_yy,
												gsl_matrix * mat_tt_yx,
												gsl_matrix * mat_tt_xy,
												gsl_matrix * mat_tt_xx,
												gsl_permutation * perm_y,
												gsl_vector * vect_t);
 * @param sqrt_p_f : racine de la matrice de covariance de l'�tat filtr� courant
 * @param sqrt_s : racine de la covariance de l'innovation.
 * @param gain : gain de filtrage P_{n+1|N} \; [F^{y,x}]^T \; S_{n+1}^{-1}
 * @param[in] sqrt_p_p : racine de la matrice de covariance de l'�tat pr�dit courant
 * @param[in] f_yx : terme de la matrice de transition (Fyx)
 * @param[in] sqrt_q_yy : Racine de la matrice de covariance du bruit de mesure (Qyy)
 * @param mat_tt : matrice de taille (t.t) pr�allou�e.
 * @param mat_tt_yy : vue sur la matrice mat_tt (de (0,0) � (n_y, n_y))
 * @param mat_tt_yx : vue sur la matrice mat_tt (de (0, n_y) � (n_y, n_t))
 * @param mat_tt_xy : vue sur la matrice mat_tt (de (n_y, 0) � (n_t, n_y))
 * @param mat_tt_xx : vue sur la matrice mat_tt (de (n_y, n_y) � (n_t, n_t))
 * @param perm_y : permutation de taille y pr�allou�e
 * @param vect_t : vecteur de taille t pr�allou� pour le calcul
 * @brief
 Cette fonction calcule les racines de la covariance de l'�tat filtr� courant et de la matrice de covariance de l'innovation. \n
ECe calcul s'effectue en plusieurs �tapes : nous construisons la matrice M : \n
+---------------------------+
|sqrt_q_yy      0           |\n
|                           |\n
|sqrt_pp.Fyx^T  sqrt_p_p    |\n
+---------------------------+\n
puis nous effectuons sa d�composition QR. A partir de la matrice R de la d�composition, nous obtenons :
+---------------------------------------+
|sqrt_s      sqrt_s^-1 K^T             |\n
|                           			|\n
|0           sqrt_p_f       		|\n
+---------------------------------------+\n
 */
void tkalman_robust_get_sqrt_pf_sqrt_s_and_gain(gsl_matrix * sqrt_p_f,
												gsl_matrix * sqrt_s,
												gsl_matrix * gain,
												const gsl_matrix * sqrt_p_p,
												const gsl_matrix * f_yx,
												const gsl_matrix * sqrt_q_yy,
												gsl_matrix * mat_tt,
												gsl_matrix * mat_tt_yy,
												gsl_matrix * mat_tt_yx,
												gsl_matrix * mat_tt_xy,
												gsl_matrix * mat_tt_xx,
												gsl_permutation * perm_y,
												gsl_vector * vect_t)
{
	//Construction de la matrice
	//      +---------------------------+
	//      |sqrt_q_yy      0           |
	// M =  |                           |
	//      |sqrt_pp.Fyx^T  sqrt_p_p    |
	//      +---------------------------+

	gsl_matrix_memcpy(mat_tt_yy, sqrt_q_yy);
	gsl_matrix_memcpy(mat_tt_xx, sqrt_p_p);
	gsl_matrix_set_zero(mat_tt_yx);
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   sqrt_p_p,
				   f_yx,
				   0.0,
				   mat_tt_xy);
	//D�composition QR
	gsl_linalg_QR_decomp(mat_tt,
						 vect_t);
	gsl_triangle_matrix(mat_tt);
	
	
	//Recopie des r�sultats.
	gsl_matrix_memcpy(sqrt_p_f,
					  mat_tt_xx);
	gsl_matrix_memcpy(sqrt_s,
					  mat_tt_yy);

	//Calcul du gain
		//Inversion de sqrt_s
        gsl_permutation_init(perm_y);
        gsl_linalg_LU_invert(sqrt_s, perm_y, mat_tt_yy);

		//Produit matriciel
		gsl_blas_dgemm(CblasTrans,
					   CblasTrans,
					   1.0,
					   mat_tt_yx,
					   mat_tt_yy,
					   0.0,
					   gain);
}

/**@fn void tkalman_robust_filtering(gsl_vector * x_f,
			                                   gsl_matrix * sqrt_p_f,
			                                   gsl_vector * innovation,
			                                   gsl_matrix * sqrt_s,
			                                   const gsl_vector * x_p,
			                                   const gsl_matrix * sqrt_p_p,
			                                   const gsl_vector * y,
			                                   const gsl_vector * _y,
			                                   const gsl_matrix * f_yx,
			                                   const gsl_matrix * f_yy,
			                                   const gsl_matrix * sqrt_q_yy,
			                                   gsl_matrix * mat_tt,
			                                   gsl_matrix * mat_tt_yy,
			                                   gsl_matrix * mat_tt_yx,
			                                   gsl_matrix * mat_tt_xy,
			                                   gsl_matrix * mat_tt_xx,
			                                   gsl_matrix * mat_xy,
			                                   gsl_permutation * perm_y,
			                                   gsl_vector * vect_t)
 * @param x_f : esp�rance de l'�tat filtr� courant
 * @param sqrt_p_f : racine de la matrice de covariance de l'�tat filtr� courant
 * @param innovation : esp�rance de l'innovation
 * @param sqrt_s : racine de la covariance de l'innovation.
 * @param[in]�x_p : \hat{x}_{n | n - 1}, esp�rance de l'�tat pr�dit courant
 * @param[in]�x_p : \hat{x}_{n | n - 1}, esp�rance de l'�tat pr�dit
 * @param[in] y : observation cournate
 * @param[in] _y : observation pr�c�dente
 * @param[in] f_yx : terme de la matrice de transition (Fyx)
 * @param[in] f_yy : terme de la matrice de transition (Fyy)
 * @param[in] sqrt_q_yy : Racine de la matrice de covariance du bruit de mesure (Qyy)
 * @param mat_tt : matrice de taille (t.t) pr�allou�e.
 * @param mat_tt_yy : vue sur la matrice mat_tt (de (0,0) � (n_y, n_y))
 * @param mat_tt_yx : vue sur la matrice mat_tt (de (0, n_y) � (n_y, n_t))
 * @param mat_tt_xy : vue sur la matrice mat_tt (de (n_y, 0) � (n_t, n_y))
 * @param mat_tt_xx : vue sur la matrice mat_tt (de (n_y, n_y) � (n_t, n_t))
 * @param mat_xy : matrice de taille (x.y) pr�allou�e
 * @param perm_y : permutation de taille y pr�allou�e
 * @param vect_t : vecteur de taille t pr�allou� pour le calcul
 * @brief
 * Cette fonction effectue la partie filtrage du filtre de Kalman triplet : \n
 Nous calculons dans un premier temps l'esp�rance de l'innovation. Puis dans un second temps, nous calculons les racines des matrices de covariance de l'innovarion et de l'�tat filtr� courant. Ensuite dans un troisi�me temps, nous calculons le gain de filtrage.Finalement dans un quatri�me temps, nous calculons l'esp�rance de l'�tat filtr�.
 */
void tkalman_robust_filtering(gsl_vector * x_f,
                              gsl_matrix * sqrt_p_f,
                              gsl_vector * innovation,
                              gsl_matrix * sqrt_s,
                              const gsl_vector * x_p,
                              const gsl_matrix * sqrt_p_p,
                              const gsl_vector * y,
                              const gsl_vector * _y,
                              const gsl_matrix * f_yx,
                              const gsl_matrix * f_yy,
                              const gsl_matrix * sqrt_q_yy,
                              gsl_matrix * mat_tt,
                              gsl_matrix * mat_tt_yy,
                              gsl_matrix * mat_tt_yx,
                              gsl_matrix * mat_tt_xy,
                              gsl_matrix * mat_tt_xx,
                              gsl_matrix * mat_xy,
                              gsl_permutation * perm_y,
                              gsl_vector * vect_t)
{
    //Calcul de l'esp�rance de l'innovation
    tkalman_original_get_innovation ( innovation,
									  x_p,
									  _y,
									  y,
									  f_yx,
									  f_yy);

	//Racines des matrices de covariance et du gain
	tkalman_robust_get_sqrt_pf_sqrt_s_and_gain ( sqrt_p_f,
												 sqrt_s,
												 mat_xy,
											 	 sqrt_p_p,
												 f_yx,
												 sqrt_q_yy,
												 mat_tt,
												 mat_tt_yy,
												 mat_tt_yx,
												 mat_tt_xy,
												 mat_tt_xx,
												 perm_y,
												 vect_t );

	//Esp�rance de l'�tat filtr�
	tkalman_original_get_x_f ( x_f,
							   x_p,
							   innovation,
							   mat_xy);
}
