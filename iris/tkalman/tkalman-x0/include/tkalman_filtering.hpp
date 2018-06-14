/**
* @file tkalman_filtering.hpp
* @author Val�rian N�mesin
* @date 04/10/2011
* @brief
Ce fichier contient les prototypes des fonctions n�cessaires � la passe-avant du filtre de Kalman triplet. \n
**/
#ifndef _TKALMAN_FILTERING_HPP_
	#define _TKALMAN_FILTERING_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"

	/**@fn void tkalman_original_get_x_p(gsl_vector * x_p,
									 const gsl_vector * _x_f,
									 const gsl_vector * __y,
									 const gsl_vector * _y,
									 const gsl_matrix * f2_xx,
									 const gsl_matrix * f2_xy,
									 const gsl_matrix * q2_xy)
	 * @param x_p : esp�rance de l'�tat pr�dit
	 * @param[in] _x_f : esp�rance de l'�tat filtr� pr�c�dent
	 * @param[in] __y : observation (n - 2)
	 * @param[in] _y : observation pr�c�dente
	 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
	 * @param[in] f2_xy : F_2^{x,y} = F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y}
	 * @param[in] q2_xy : Q_2^{x,y} = Q^{x,y} \: [Q^{y,y}]^{-1}
	 * @brief
	 Cette fonction calcule l'esp�rance de l'�tat pr�dit � partir de la formule :
	 * \hat{x_{n|n-1}} = F_2^{x,x} \: \hat{x_{n - 1 |n-1}} + Q_2^{x,y} \; y_{n - 2} + F_2^{x,y} \: y_{n - 1}
	 */
	void tkalman_original_get_x_p(gsl_vector * x_p,
								  const gsl_vector * _x_f,
								  const gsl_vector * __y,
								  const gsl_vector * _y,
								  const gsl_matrix * f2_xx,
								  const gsl_matrix * f2_xy,
								  const gsl_matrix * q2_xy);

	/**@fn void tkalman_original_get_p_p(gsl_matrix * p_p,
										 const gsl_matrix * _p_f,
										 const gsl_matrix * f2_xx,
										 const gsl_matrix * q2_xx,
										 gsl_matrix * mat_xx)
	 * @param p_p : Matrice de covariance de l'�tat pr�dit.
	 * @param[in] _p_f : Matrice de covariance de l'�tat filtr� pr�c�dent
	 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
	 * @param[in] q2_xx : Matrice de covariance du bruit de process r�duit
	 * @param mat_xx : Matrice de taille (x.x) pr�allou�e.
	 * @brief
	 Cette fonction calcule la matrice de covariance de l'�tat pr�dit � partir de la formule :
	 * P_{n|n-1} = Q_2^{x,x} + F_2^{x,x}�\: P_{n-1|n-1} \: [F_2^{x,x}]^T
	 */
	void tkalman_original_get_p_p(gsl_matrix * p_p,
								  const gsl_matrix * _p_f,
								  const gsl_matrix * f2_xx,
								  const gsl_matrix * q2_xx,
								  gsl_matrix * mat_xx);

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
	 * @param x_p : esp�rance de l'�tat pr�dit \hat{x_{n|n-1}} = F_2^{x,x} \: \hat{x_{n - 1 |n-1}} + Q_2^{x,y} \; y_{n - 2} + F_2^{x,y} \: y_{n - 1}
	 * @param p_p : Matrice de covariance de l'�tat pr�dit. P_{n|n-1} = Q_2^{x,x} + F_2^{x,x}�\: P_{n-1|n-1} \: [F_2^{x,x}]^T
	 * @param[in] _x_f : esp�rance de l'�tat filtr� pr�c�dent
	 * @param[in] _p_f : Matrice de covariance de l'�tat filtr� pr�c�dent
	 * @param[in] __y : observation (n - 2)
	 * @param[in] _y : observation pr�c�dente
	 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
	 * @param[in] f2_xy : F_2^{x,y} = F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y}
	 * @param[in] q2_xx : Matrice de covariance du bruit de process r�duit
	 * @param[in] q2_xy : Q_2^{x,y} = Q^{x,y} \: [Q^{y,y}]^{-1}
	 * @brief
	 Cette fonction effectue la pr�diction dans le filtre de Kalman triplet(Voir les arguments pour les forumules.).
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
									 gsl_matrix * mat_xx);

	/**@fn void tkalman_robust_get_sqrt_p_p(gsl_matrix * sqrt_p_p,
								     const gsl_matrix * _sqrt_p_f,
								     const gsl_matrix * f2_xx,
								     const gsl_matrix * sqrt_q2_xx,
								     gsl_matrix * mat_2xx,
								     gsl_matrix * mat_2xx_view_00,
								     gsl_matrix * mat_2xx_view_10,
								     gsl_vector * vect_x)
	 * @param sqrt_p_p : racine de la matrice de covariance de l'�tat pr�dit courant
	 * @param[in] _sqrt_p_f : racine de la matrice de covariance de l'�tat filtr� pr�c�dent
	 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
	 * @param[in] sqrt_q2_xx : racine de la matrice de covariance du bruit de process r�duit
	 * @param mat_2xx : Matrice de taille (2x.x) pr�allou�e.
	 * @param mat_2xx_view_00 : vue sur la matrice mat_2xx allant de (0,0) � (n_x, n_x)
	 * @param mat_2xx_view_10 : vue sur la matrice mat_2xx allant de (n_x, 0) � (2n_x, n_x)
	 * @param vect_x : vecteur de taille (x) pr�allou�
	 * @brief
	 Cette fonction calcule la racine de la matrice de covariance de l'�tat pr�dit courant.
	 */
	void tkalman_robust_get_sqrt_p_p(gsl_matrix * sqrt_p_p,
								     const gsl_matrix * _sqrt_p_f,
								     const gsl_matrix * f2_xx,
								     const gsl_matrix * sqrt_q2_xx,
								     gsl_matrix * mat_2xx,
								     gsl_matrix * mat_2xx_view_00,
								     gsl_matrix * mat_2xx_view_10,
								     gsl_vector * vect_x);

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
	 * @param x_p : esp�rance de l'�tat pr�dit \hat{x_{n|n-1}} = F_2^{x,x} \: \hat{x_{n - 1 |n-1}} + Q_2^{x,y} \; y_{n - 2} + F_2^{x,y} \: y_{n - 1}
	 * @param p_p : Matrice de covariance de l'�tat pr�dit. P_{n|n-1} = Q_2^{x,x} + F_2^{x,x}�\: P_{n-1|n-1} \: [F_2^{x,x}]^T
	 * @param[in] _x_f : esp�rance de l'�tat filtr� pr�c�dent
	 * @param[in] _sqrt_p_f : racine de la matrice de covariance de l'�tat filtr� pr�c�dent
	 * @param[in] __y : observation (n - 2)
	 * @param[in] _y : observation pr�c�dente
	 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
	 * @param[in] f2_xy : F_2^{x,y} = F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y}
	 * @param[in] sqrt_q2_xx : racine de la matrice de covariance du bruit de process r�duit
	 * @param[in] q2_xy : Q_2^{x,y} = Q^{x,y} \: [Q^{y,y}]^{-1}
	 * @param mat_2xx : Matrice de taille (2x.x) pr�allou�e.
	 * @param mat_2xx_view_00 : vue sur la matrice mat_2xx allant de (0,0) � (n_x, n_x)
	 * @param mat_2xx_view_10 : vue sur la matrice mat_2xx allant de (n_x, 0) � (2n_x, n_x)
	 * @param vect_x : vecteur de taille (x) pr�allou�
	 * @brief
	 Cette fonction calcule l'esp�rance et la racine de la matrice de covariance de l'�tat filtr� courant.

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
								     gsl_vector * vect_x);


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
										 const gsl_matrix * f_yy);

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
								gsl_matrix * mat_yx);

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
								   gsl_matrix * mat_yy);

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
								  const gsl_matrix * gain);

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
								  gsl_matrix * mat_yx);

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
									gsl_matrix * mat_yy);

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
													gsl_vector * vect_t);
	/**@fn void tkalman_robust_filtering(gsl_matrix * x_f,
				                                   gsl_matrix * sqrt_p_f,
				                                   gsl_vector * innovation,
				                                   gsl_matrix * sqrt_s,
				                                   const gsl_matrix * x_p,
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
	                              gsl_vector * vect_t);

#endif
