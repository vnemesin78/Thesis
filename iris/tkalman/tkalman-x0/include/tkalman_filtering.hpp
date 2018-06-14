/**
* @file tkalman_filtering.hpp
* @author Valérian Némesin
* @date 04/10/2011
* @brief
Ce fichier contient les prototypes des fonctions nécessaires à la passe-avant du filtre de Kalman triplet. \n
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
								  const gsl_matrix * q2_xy);

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
									 gsl_matrix * mat_xx);

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
								     gsl_vector * vect_x);


	/**@fn void tkalman_original_get_innovation(gsl_vector * innovation,
											const gsl_vector * x_p,
											const gsl_vector * _y,
											const gsl_vector * y,
											const gsl_matrix * f_yx,
											const gsl_matrix * f_yy)
	 * @param innovation : innovation
	 * @param[in] x_p : \hat{x}_{n | n - 1}, espérance de l'état prédit
	 * @param[in] _y : y_{n - 1}, observation précédente
	 * @param[in] y : y_n, observation
	 * @param[in] f_yx : F^{y,x}, terme de la matrice d'évolution
	 * @param[in] f_yy : F^{y,y}, terme de la matrice d'évolution
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
	 * @param[in] p_p : P_{n|n-1}, covariance de l'état prédit
	 * @param[in] f_yx : F^{y,x}, terme de la matrice d'évolution
	 * @param[in] q_yy : Q^{y,y}, matrice de covariance du bruit de mesure
	 * @param mat_yx : matrice de taille (y.x) préallouée.
	 * @brief
	 Cette fonction calcule la covariance de l'innovation décrite par la forumle ci dessous:
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
	 * @param[in] p_p : P_{n|n-1}, covariance de l'état prédit
	 * @param[in] s : covariance de l'innovation
	 * @param[in] f_yx : F^{y,x}, terme de la matrice d'évolution
	 * @param mat_yx : matrice de taille (y.x) préallouée.
	 * @param mat_yy : matrice de taille (y.y) préallouée.
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
	 * @param x_f : espérance de l'état filtré
	 * @param[in] x_p : espérance de l'état prédit
	 * @param[in] innovation : innovation
	 * @param[in] gain : gain
	 * @brief
	 Cette fonction calcule l'espérance de l'état filtré selon la formule :
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
	 * @param p_f : Matrice de covariance de l'état filtré
	 * @param[in] p_p : matrice de covariance de l'état prédit
	 * @param[in] s : matrice de covariance de l'innovarion
	 * @param[in] gain : gain de Kalman
	 * @param mat_yx : matrice de taille (y.x) préallouée.
	 * @brief
	 Cette fonction calcule la matrice de covariance de l'état filtré selon la formule :
	 * P_{n|n} = P_{n|n-1} + K_{n|n} \; S_{n} \: K_{n|n}^T
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
	 * @param x_f : espérance de l'état filtré  \hat{x}_{n|n} = \hat{x|n-1} + K_{n|n} \; \tilde{y}_n
	 * @param p_f : Matrice de covariance de l'état filtré P_{n|n} = P_{n|n-1} + K_{n|n} \; S_{n} \: K_{n|n}^T
	 * @param innovation : innovation  \tilde{y}_n = y_{n} - F^{y,x} \: \hat{x}_{n | n - 1} - F^{y,y} \: y_{n - 1}
	 * @param s : covariance de l'innovation Q^{y,y} + F^{y,x} \: P_{n | n - 1} \: [F^{y,x}]^T
	 * @param[in] x_p : espérance de l'état prédit
	 * @param[in] p_p : matrice de covariance de l'état prédit
	 * @param[in] _y : y_{n - 1}, observation précédente
	 * @param[in] y : y_n, observation
	 * @param[in] f_yx : F^{y,x}, terme de la matrice d'évolution
	 * @param[in] f_yy : F^{y,y}, terme de la matrice d'évolution
	 * @param[in] q_yy : Q^{y,y}, matrice de covariance du bruit de mesure
	 * @param mat_xy : Matrice de taille (x.y) préallouée.
	 * @param mat_yx : matrice de taille (y.x) préallouée.
	 * @param mat_yy : matrice de taille (y.y) préallouée.
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
	 * @param sqrt_p_f : racine de la matrice de covariance de l'état filtré courant
	 * @param sqrt_s : racine de la covariance de l'innovation.
	 * @param gain : gain de filtrage P_{n+1|N} \; [F^{y,x}]^T \; S_{n+1}^{-1}
	 * @param[in] sqrt_p_p : racine de la matrice de covariance de l'état prédit courant
	 * @param[in] f_yx : terme de la matrice de transition (Fyx)
	 * @param[in] sqrt_q_yy : Racine de la matrice de covariance du bruit de mesure (Qyy)
	 * @param mat_tt : matrice de taille (t.t) préallouée.
	 * @param mat_tt_yy : vue sur la matrice mat_tt (de (0,0) à (n_y, n_y))
	 * @param mat_tt_yx : vue sur la matrice mat_tt (de (0, n_y) à (n_y, n_t))
	 * @param mat_tt_xy : vue sur la matrice mat_tt (de (n_y, 0) à (n_t, n_y))
	 * @param mat_tt_xx : vue sur la matrice mat_tt (de (n_y, n_y) à (n_t, n_t))
	 * @param perm_y : permutation de taille y préallouée
	 * @param vect_t : vecteur de taille t préalloué pour le calcul
	 * @brief
	 Cette fonction calcule les racines de la covariance de l'état filtré courant et de la matrice de covariance de l'innovation. \n
	ECe calcul s'effectue en plusieurs étapes : nous construisons la matrice M : \n
	+---------------------------+
	|sqrt_q_yy      0           |\n
	|                           |\n
	|sqrt_pp.Fyx^T  sqrt_p_p    |\n
	+---------------------------+\n
	puis nous effectuons sa décomposition QR. A partir de la matrice R de la décomposition, nous obtenons :
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
	 * @param x_f : espérance de l'état filtré courant
	 * @param sqrt_p_f : racine de la matrice de covariance de l'état filtré courant
	 * @param innovation : espérance de l'innovation
	 * @param sqrt_s : racine de la covariance de l'innovation.
	 * @param[in] x_p : \hat{x}_{n | n - 1}, espérance de l'état prédit courant
	 * @param[in] x_p : \hat{x}_{n | n - 1}, espérance de l'état prédit
	 * @param[in] y : observation cournate
	 * @param[in] _y : observation précédente
	 * @param[in] f_yx : terme de la matrice de transition (Fyx)
	 * @param[in] f_yy : terme de la matrice de transition (Fyy)
	 * @param[in] sqrt_q_yy : Racine de la matrice de covariance du bruit de mesure (Qyy)
	 * @param mat_tt : matrice de taille (t.t) préallouée.
	 * @param mat_tt_yy : vue sur la matrice mat_tt (de (0,0) à (n_y, n_y))
	 * @param mat_tt_yx : vue sur la matrice mat_tt (de (0, n_y) à (n_y, n_t))
	 * @param mat_tt_xy : vue sur la matrice mat_tt (de (n_y, 0) à (n_t, n_y))
	 * @param mat_tt_xx : vue sur la matrice mat_tt (de (n_y, n_y) à (n_t, n_t))
	 * @param mat_xy : matrice de taille (x.y) préallouée
	 * @param perm_y : permutation de taille y préallouée
	 * @param vect_t : vecteur de taille t préalloué pour le calcul
	 * @brief
	 * Cette fonction effectue la partie filtrage du filtre de Kalman triplet : \n
	 Nous calculons dans un premier temps l'espérance de l'innovation. Puis dans un second temps, nous calculons les racines des matrices de covariance de l'innovarion et de l'état filtré courant. Ensuite dans un troisième temps, nous calculons le gain de filtrage.Finalement dans un quatrième temps, nous calculons l'espérance de l'état filtré.
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
