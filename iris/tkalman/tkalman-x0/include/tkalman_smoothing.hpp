/**@file tkalman_smoothing.hpp
*@author Valérian Némesin
*@brief
Ce fichier contient les prototypes des fonctions nécessaires au lissage dans le filtre de Kalman triplet.
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
	 * @param p_s : Matrice de covariance de l'état lissé
	 * @param[in] p_f : Matrice de covariance de l'état filtré
	 * @param[in] p_p_ : Matrice de covariance de l'état prédit suivant
	 * @param[in] p_s_ : Matrice de covariance de l'état lissé suivant
	 * @param[in] gain : gain de lissage
	 * @param mat_xx : Matrice de taille (x.x) préallouée
	 * @brief
	 * Cette fonction calcule la matrice de covariance de l'état lissé suivant la formule
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
	 * @param[in] p_p_ : Matrice de covariance de l'état prédit suivant
	 * @param[in] p_f : Matrice de covariance de l'état filtré
	 * @param[in] f2_xx : Matrice F_2^{x,x} = F^{x,x} - Q^{x,y}\; (Q^{y,y})^{-1} \;F^{y,x}
	 * @param mat_xx : Matrice de taille (x,x) préallouée
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
	 * @param x_s : Espérance de l'état lissé
	 * @param p_s : Matrice de covariance de l'état lissé
	 * @param c_s : Matrice de covariance entre l'état lissé et l'état lissé suivant. (si NULL, rien n'est renvoyé.)
	 * @param[in] x_f : Espérance de l'état filtré
	 * @param[in] p_f : Matrice de covariance de l'état filtré
	 * @param[in] x_p_ : Espérance de l'état prédit suivant
	 * @param[in] p_p_ : Matrice de covariance de l'état prédit suivant
	 * @param[in] x_s_ : Espérance de l'état lissé suivant
	 * @param[in] p_s_ : Matrice de covariance de l'état lissé suivant
	 * @param[in] f2_xx : Matrice F_2^{x,x} = F^{x,x} - Q^{x,y}\; (Q^{y,y})^{-1} \;F^{y,x}
	 * @param mat_xx_1 : Matrice de taille (x,x) préallouée
	 * @param mat_xx_2 : Matrice de taille (x,x) préallouée
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
	 * @param[in] sqrt_p_f : racine de la covariance de l'état filtré courant
	 * @param[in] sqrt_p_p_ : racine de la covariance de l'état prédit suivant
	 * @param[in] f2_xx : Fxx - Qxy.Qyy?¹.Fyx
	 * @param mat_xx : matrice de taille (x.x) préallouée
	 * @param perm_x : permutation de taille (x) préallouée
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
