/**@file tkalman_smoothing.hpp
*@author Val�rian N�mesin
*@brief
Ce fichier contient les prototypes des fonctions n�cessaires au lissage dans le filtre de Kalman triplet.
**/
#ifndef _TKALMAN_SMOOTHING_HPP_
	#define _TKALMAN_SMOOTHING_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"
	/**@fn void tkalman_original_get_x_s(gsl_vector * x_s,
									     const gsl_vector * x_f,
										 const gsl_vector * x_p_,
										 const gsl_vector * x_s_,
										 const gsl_matrix * gain)
	 * @param x_s : Esp�rance de l'�tat liss�
	 * @param[in] x_f : Esp�rance de l'�tat filtr�
	 * @param[in] x_p_ : Esp�rance de l'�tat pr�dit suivant
	 * @param[in] x_s_ : Esp�rance de l'�tat liss� suivant
	 * @param[in] gain : Gain de lissage
	 * @brief
	 * Cette fonction calcule l'esp�rance de l'�tat liss� selon la formule \n
	 * \hat{x}_{n|N} = \hat{x}_{n|n} + K_{n|N} \: (\hat{x}_{n + 1|N} - \hat{x}_{n+1|n}
	 * .
	 */
	void tkalman_original_get_x_s(gsl_vector * x_s,
								  const gsl_vector * x_f,
								  const gsl_vector * x_p_,
								  const gsl_vector * x_s_,
								  const gsl_matrix * gain);

	/**@fn void tkalman_original_get_p_s(gsl_matrix * p_s,
										 const gsl_matrix * p_f,
										 const gsl_matrix * p_p_,
										 const gsl_matrix * p_s_,
										 const gsl_matrix * gain,
										 gsl_matrix * mat_xx)
	 * @param p_s : Matrice de covariance de l'�tat liss�
	 * @param[in] p_f : Matrice de covariance de l'�tat filtr�
	 * @param[in] p_p_ : Matrice de covariance de l'�tat pr�dit suivant
	 * @param[in] p_s_ : Matrice de covariance de l'�tat liss� suivant
	 * @param[in] gain : gain de lissage
	 * @param mat_xx : Matrice de taille (x.x) pr�allou�e
	 * @brief
	 * Cette fonction calcule la matrice de covariance de l'�tat liss� suivant la formule
	 * P_{n|N} = P_{n|n} + K_{n|N} \: (P_{n+1|N} - P_{n+1|n}) \: K_{n|N}^T
	 */
	void tkalman_original_get_p_s(gsl_matrix * p_s,
								  const gsl_matrix * p_f,
								  const gsl_matrix * p_p_,
								  const gsl_matrix * p_s_,
								  const gsl_matrix * gain,
								  gsl_matrix * mat_xx);

	/**@fn void tkalman_original_get_smoothing_gain(gsl_matrix * gain,//mat_xx
													const gsl_matrix * p_f,
													const gsl_matrix * p_p_,
													const gsl_matrix * f2_xx,
													gsl_matrix * mat_xx)
	 * @param gain : gain de lissage
	 * @param[in] p_p_ : Matrice de covariance de l'�tat pr�dit suivant
	 * @param[in] p_f : Matrice de covariance de l'�tat filtr�
	 * @param[in] f2_xx : Matrice F_2^{x,x} = F^{x,x} - Q^{x,y}\; (Q^{y,y})^{-1} \;F^{y,x}
	 * @param mat_xx : Matrice de taille (x,x) pr�allou�e
	 * @brief
	 * Cette fonction calcule le gain du filtre de Kalman triplet suivant la formule : \n
	 * K_{n|N} = P_{n|n} \; [F_2^{x,x}]^T \; P_{n+1|n}^{-1}.
	 */
	void tkalman_original_get_smoothing_gain(gsl_matrix * gain,//mat_xx
											 const gsl_matrix * p_f,
											 const gsl_matrix * p_p_,
											 const gsl_matrix * f2_xx,
											 gsl_matrix * mat_xx);

	/**@fn void tkalman_original_smoothing(gsl_vector * x_s,
										   gsl_matrix * p_s,
										   gsl_matrix * c_s,
										   const gsl_vector * x_f,
										   const gsl_matrix * p_f,
										   const gsl_vector * x_p_,
										   const gsl_matrix * p_p_,
										   const gsl_vector * x_s_,
										   const gsl_matrix * p_s_,
										   const gsl_matrix * f2_xx,
										   gsl_matrix * mat_xx_1,
										   gsl_matrix * mat_xx_2)
	 * @param x_s : Esp�rance de l'�tat liss�
	 * @param p_s : Matrice de covariance de l'�tat liss�
	 * @param c_s : Matrice de covariance entre l'�tat liss� et l'�tat liss� suivant. (si NULL, rien n'est renvoy�.)
	 * @param[in] x_f : Esp�rance de l'�tat filtr�
	 * @param[in] p_f : Matrice de covariance de l'�tat filtr�
	 * @param[in] x_p_ : Esp�rance de l'�tat pr�dit suivant
	 * @param[in] p_p_ : Matrice de covariance de l'�tat pr�dit suivant
	 * @param[in] x_s_ : Esp�rance de l'�tat liss� suivant
	 * @param[in] p_s_ : Matrice de covariance de l'�tat liss� suivant
	 * @param[in] f2_xx : Matrice F_2^{x,x} = F^{x,x} - Q^{x,y}\; (Q^{y,y})^{-1} \;F^{y,x}
	 * @param mat_xx_1 : Matrice de taille (x,x) pr�allou�e
	 * @param mat_xx_2 : Matrice de taille (x,x) pr�allou�e
	 * @brief
	 * Cette fonction effectue le lissage du filtre de Kalman Couple.
	 *
	 */
	void tkalman_original_smoothing(gsl_vector * x_s,
									gsl_matrix * p_s,
									gsl_matrix * c_s,
									const gsl_vector * x_f,
								    const gsl_matrix * p_f,
								    const gsl_vector * x_p_,
									const gsl_matrix * p_p_,
									const gsl_vector * x_s_,
									const gsl_matrix * p_s_,
									const gsl_matrix * f2_xx,
									gsl_matrix * mat_xx_1,
									gsl_matrix * mat_xx_2);

	/**@fn void tkalman_robust_get_smoothing_gain(gsl_matrix * s_gain,
										   const gsl_matrix * sqrt_p_f,
										   const gsl_matrix * sqrt_p_p_,
									       const gsl_matrix * f2_xx,
										   gsl_matrix * mat_xx,
										   gsl_permutation * perm_x)
	 * @param s_gain : gain de Kalman triplet pour le lissage
	 * @param[in] sqrt_p_f : racine de la covariance de l'�tat filtr� courant
	 * @param[in] sqrt_p_p_ : racine de la covariance de l'�tat pr�dit suivant
	 * @param[in] f2_xx : Fxx - Qxy.Qyy?�.Fyx
	 * @param mat_xx : matrice de taille (x.x) pr�allou�e
	 * @param perm_x : permutation de taille (x) pr�allou�e
	 * @brief
	 * Cette fonction calcule le gain de lissage du filtre de Kalman
	 * triplet.
	 */
	void tkalman_robust_get_smoothing_gain(gsl_matrix * s_gain,
										   const gsl_matrix * sqrt_p_f,
										   const gsl_matrix * sqrt_p_p_,
									       const gsl_matrix * f2_xx,
										   gsl_matrix * mat_xx,
										   gsl_permutation * perm_x);

	/**@fn void tkalman_robust_get_sqrt_p_s_and_c_s(gsl_matrix * sqrt_p_s,
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
	 * @param sqrt_p_s : racine de la matrice de covariance de l'�tat liss� courant
	 * @param c_s : sqrt(P_{n+1|N}) \; K_{n|N}^T
	 * @param[in] sqrt_p_f : racine de la matrice de covariance de l'�tat filtr� courant
	 * @param[in]  sqrt_p_s_ : racine de la matrice de covariance de l'�tat liss� suivant
	 * @param[in] f2_xx : Fxx - Qxy.Qyy?�.Fyx
	 * @param[in] sqrt_q2_xx : racine de la matrice de covariance du bruit de process r�duit
	 * @param[in] s_gain :  gain de Kalman triplet pour le lissage
	 * @param mat_3x2x : matrice de taille (3x.2x) pr�allou�e
	 * @param mat_3x2x_view_00 : vue sur la matrice mat_3x2x allant de (0,0) � (n_x, n_x)
	 * @param mat_3x2x_view_01 : vue sur la matrice mat_3x2x allant de (0,n_x) � (n_x,2 n_x)
	 * @param mat_3x2x_view_10 : vue sur la matrice mat_3x2x allant de (n_x,0) � (2n_x, n_x)
	 * @param mat_3x2x_view_11 : vue sur la matrice mat_3x2x allant de (n_x,n_x) � (2n_x, 2n_x)
	 * @param mat_3x2x_view_20 : vue sur la matrice mat_3x2x allant de (2n_x,0) � (3n_x, n_x)
	 * @param mat_3x2x_view_21 : vue sur la matrice mat_3x2x allant de (2n_x,n_x) � (3n_x, 2n_x)
	 * @param vect_2x : vecteur de taille (2x) pr�allou�
	 * @brief
	 Cette fonction calcule la racine de la matrice de covariance de l'�tat liss� courant et sqrt(P_{n+1|N}) \; K_{n|N}^T.
	**/
	void tkalman_robust_get_sqrt_p_s_and_c_s(gsl_matrix * sqrt_p_s,
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
											 gsl_vector * vect_2x);

	/**@fn void tkalman_robust_smoothing(gsl_vector * x_s,
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
	 * @param x_s : Esp�rance de l'�tat liss�
	 * @param sqrt_p_s : racine de la matrice de covariance de l'�tat liss� courant
	 * @param c_s : sqrt(P_{n+1|N}) \; K_{n|N}^T
	 * @param[in] x_f : Esp�rance de l'�tat filtr�
	 * @param[in] sqrt_p_f : racine de la matrice de covariance de l'�tat filtr� courant
	 * @param[in] x_p_ : Esp�rance de l'�tat pr�dit suivant
	 * @param[in] sqrt_p_p_ : racine de la covariance de l'�tat pr�dit suivant
	 * @param[in] x_s_ : Esp�rance de l'�tat liss� suivant
	 * @param[in]  sqrt_p_s_ : racine de la matrice de covariance de l'�tat liss� suivant
	 * @param[in] f2_xx : Fxx - Qxy.Qyy?�.Fyx
	 * @param[in] sqrt_q2_xx : racine de la matrice de covariance du bruit de process r�duit
	 * @param mat_xx : matrice de taille (x.x) pr�allou�e
	 * @param mat_3x2x : matrice de taille (3x.2x) pr�allou�e
	 * @param mat_3x2x_view_00 : vue sur la matrice mat_3x2x allant de (0,0) � (n_x, n_x)
	 * @param mat_3x2x_view_01 : vue sur la matrice mat_3x2x allant de (0,n_x) � (n_x,2 n_x)
	 * @param mat_3x2x_view_10 : vue sur la matrice mat_3x2x allant de (n_x,0) � (2n_x, n_x)
	 * @param mat_3x2x_view_11 : vue sur la matrice mat_3x2x allant de (n_x,n_x) � (2n_x, 2n_x)
	 * @param mat_3x2x_view_20 : vue sur la matrice mat_3x2x allant de (2n_x,0) � (3n_x, n_x)
	 * @param mat_3x2x_view_21 : vue sur la matrice mat_3x2x allant de (2n_x,n_x) � (3n_x, 2n_x)
	 * @param perm_x : permutation de taille (x) pr�allou�e
	 * @param vect_2x : vecteur de taille (2x) pr�allou�
	 * @brief
	 Cette fonction effectue l'�tape de lissage dans le filtre de Kalman triplet.
	**/
	void tkalman_robust_smoothing(gsl_vector * x_s,
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
								  gsl_vector * vect_2x);
#endif
