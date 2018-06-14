/**@file tkalman_em_algorithm.hpp
 * @author Valérian Némesin
 * @brief
 Ce fichier contient les prototypes des fonctions nécessaires à l'algorithme EM du filtre de Kalman triplet.
**/
#ifndef _TKALMAN_EM_ALGORITHM_HPP_
	#define _TKALMAN_EM_ALGORITHM_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include <cmath>
	#include "gsl_triangle_matrix.hpp"
	/**@fn void tkalman_original_get_corr_t_n_t_p( gsl_matrix * corr,
												   gsl_matrix * corr_xx,
												   const gsl_vector * t_n,
												   const gsl_vector * t_p,
												   const gsl_matrix * cov);
	 * @param corr : matrice de corrélation entre t_{n|N} et t_{p|N} si corr = 0
	 sinon corr += matrice de corrélation entre t_{n|N} et t_{p|N}
	 * @param corr_xx : vue sur la matrice de corrélation (0,0) à (x,x)
	 * @param t_n : (\hat{x}_{n|N}^T, y_{n - 1})^T
	 * @param t_p : (\hat{x}_{p|N}^T, y_{p - 1})^T
	 * @param cov : matrice de covariance entre x_{n|N} et x_{p|N}.
	 * @brief
	 * Cette fonction calcule la matrice de corrélation entre t_{n|N} et t_{p|N}.
	 */
	void tkalman_original_get_corr_t_n_t_p( gsl_matrix * corr,
											gsl_matrix * corr_xx,
											const gsl_vector * t_n,
											const gsl_vector * t_p,
											const gsl_matrix * cov);
	/**@fn void tkalman_original_get_sums( gsl_matrix * c_00,
										   gsl_matrix * c_00_xx,
										   gsl_matrix * c_10,
										   gsl_matrix * c_10_xx,
										   gsl_matrix * c_11,
										   gsl_matrix * c_11_xx,
										   const gsl_vector * const * x_s,
										   const gsl_vector * const * y,
										   const gsl_matrix * const * p_s,
										   const gsl_matrix * const * c_s,
										   unsigned int n,
										   gsl_vector * vect_t_1,
										   gsl_vector * vect_t_2)
	 * @param c_00 : somme \tilde{C}_{t_n,t_n} = \sum_{n=0}^N C_{t_n,t_n}
	 * @param c_00_xx : vue sur c_00 de (0,0) à (x,x)
	 * @param c_10 : somme \tilde{C}_{t_{n + 1},t_n} = \sum_{n=0}^N C_{t_{n + 1},t_n}
	 * @param c_10_xx : vue sur c_10 de (0,0) à (x,x)
	 * @param c_11 : somme \tilde{C}_{t_{n + 1},t_{n + 1}} = \sum_{n=0}^N C_{t_{n + 1},t_{n + 1}}
	 * @param c_11_xx : vue sur c_11 de (0,0) à (x,x)
	 * @param x_s : espérances des états lissés
	 * @param y : observations
	 * @param p_s : matrice de covariance des états lissés
	 * @param n : nombre d'observations
	 * @param vect_t_1 : vecteur de taille t préalloué
	 * @param vect_t_2 : vecteur de taille t préalloué
	 * @brief
	 * Cette fonction calcule les sommes du filtre de Kalman triplet.
	 */
	void tkalman_original_get_sums( gsl_matrix * c_00,
									gsl_matrix * c_00_xx,
									gsl_matrix * c_10,
									gsl_matrix * c_10_xx,
									gsl_matrix * c_11,
									gsl_matrix * c_11_xx,
									const gsl_vector * const * x_s,
									const gsl_vector * const * y,
									const gsl_matrix * const * p_s,
									const gsl_matrix * const * c_s,
									unsigned int n,
									gsl_vector * vect_t_1,
									gsl_vector * vect_t_1_x,
									gsl_vector * vect_t_1_y,
									gsl_vector * vect_t_2,
									gsl_vector * vect_t_2_x,
									gsl_vector * vect_t_2_y);

	/**@fn void tkalman_robust_get_sqrt_cov_xx_(  const gsl_matrix * sqrt_p_f,
										  const gsl_matrix * sqrt_p_s_,
										  const gsl_matrix * c_s,
										  const gsl_matrix * f2_xx,
										  const gsl_matrix * sqrt_q2_xx,
										  gsl_matrix * mat_4x4x,
										  gsl_matrix * mat_4x4x_view_00,
										  gsl_matrix * mat_4x4x_view_11,
										  gsl_matrix * mat_4x4x_view_02,
										  gsl_matrix * mat_4x4x_view_20,
										  gsl_matrix * mat_4x4x_view_32,
										  gsl_matrix * mat_4x4x_view_33,
										  gsl_vector * vect_4x)
	 * @param[in]  sqrt_p_f : racine de la matrice de covariance de l'état filtré courant
	 * @param[in] sqrt_p_s : racine de la matrice de covariance de l'état lissé suivant
	 * @param[in] c_s : Matrice de covariance entre l'état lissé et l'état lissé suivant. (si NULL, rien n'est renvoyé.)
	 * @param[in] f2_xx : Fxx - Qxy.Qyy?¹.Fyx
	 * @param[in] sqrt_q2_xx : racine de la matrice de covariance du bruit de process réduit
	 * @param mat_4x4x : matrice de taille (4x.4x) préallouée
	 * @param mat_4x4x_view_00 : vue sur la matrice mat_4x4x  allant de (0.n_x, 0.n_x) à (1.n_x, 1.n_x).
	* @param mat_4x4x_view_11 : vue sur la matrice mat_4x4x  allant de (1.n_x, 1.n_x) à (2.n_x, 2.n_x).
	* @param mat_4x4x_view_02 : vue sur la matrice mat_4x4x  allant de (0.n_x, 2.n_x) à (1.n_x, 3.n_x).
	* @param mat_4x4x_view_20 : vue sur la matrice mat_4x4x  allant de (2.n_x, 0.n_x) à (3.n_x, 1.n_x).
	* @param mat_4x4x_view_32 : vue sur la matrice mat_4x4x  allant de (3.n_x, 2.n_x) à (4.n_x, 3.n_x).
	* @param mat_4x4x_view_33 : vue sur la matrice mat_4x4x  allant de (3.n_x, 3.n_x) à (4.n_x, 4.n_x).
	* @param vect_4x : vecteur de taille (4x) préalloué.
	* @brief
	Cette fonction calcule la racine de la matrice de covariance du vecteur (x_{n|N}; x_{n+1|N}). Le résultat est dans la vue de la matrice mat_4x4x allant de (2n_x, 2n_x) à (4n_x, 4n_x).
	**/
	void tkalman_robust_get_sqrt_cov_xx_( const gsl_matrix * sqrt_p_f,
										  const gsl_matrix * sqrt_p_s_,
										  const gsl_matrix * c_s,
										  const gsl_matrix * f2_xx,
										  const gsl_matrix * sqrt_q2_xx,
										  gsl_matrix * mat_4x4x,
										  gsl_matrix * mat_4x4x_view_00,
										  gsl_matrix * mat_4x4x_view_11,
										  gsl_matrix * mat_4x4x_view_02,
										  gsl_matrix * mat_4x4x_view_20,
										  gsl_matrix * mat_4x4x_view_32,
										  gsl_matrix * mat_4x4x_view_33,
										  gsl_vector * vect_4x);

	/**@fn void tkalman_robust_get_sqrt_corr_tt_( const gsl_vector * x_s,
								    const gsl_vector * x_s_,
								    const gsl_vector * _y,
								    const gsl_vector * y,
								    const gsl_matrix * sqrt_cov_xx_view_00,
								    const gsl_matrix * sqrt_cov_xx_view_01,
								    const gsl_matrix * sqrt_cov_xx_view_11,
								    gsl_matrix * mat_2tp12t,
								    gsl_matrix * mat_2tp12t_view_00,
								    gsl_matrix * mat_2tp12t_view_02,
								    gsl_matrix * mat_2tp12t_view_22,
								    gsl_vector * mat_2tp12t_view_30,
								    gsl_vector * mat_2tp12t_view_31,
								    gsl_vector * mat_2tp12t_view_32,
								    gsl_vector * mat_2tp12t_view_33,
								    gsl_vector * vect_2t )
	* @param[in] x_s : espérance de l'état lissé courant
	* @param[in] x_s_ : espérance de l'état lissé suivant
	* @param[in] _y : observation précédente
	* @param[in] y : observation courante
	* @param[in] sqrt_cov_xx_view_00 : vue sur la matrice sqrt_cov_xx allant de (0,0) à (n_x, n_x)
	* @param[in] sqrt_cov_xx_view_01 : vue sur la matrice sqrt_cov_xx allant de (0,n_x) à (n_x, 2n_x)
	* @param[in] sqrt_cov_xx_view_11 : vue sur la matrice sqrt_cov_xx allant de (n_x,n_x) à (2n_x, 2n_x)
	* @param mat_2tp12t : matrice de taille ((2t + 1). (2t)) préallouée
	* @param mat_2tp12t_view_00 : vue sur la matrice mat_2tp12t allant de (0,0) à (n_x, n_x)
	* @param mat_2tp12t_view_02 : vue sur la matrice mat_2tp12t allant de (0,n_t) à (n_x, n_t + n_x)
	* @param mat_2tp12t_view_22 : vue sur la matrice mat_2tp12t allant de (n_t,n_t) à (n_t + n_x, n_t + n_x)
	* @param mat_2tp12t_view_30 : vue vectorielle sur la matrice mat_2tp12t allant de (2n_t,0) à (2n_t + 1, n_x)
	* @param mat_2tp12t_view_31 : vue vectorielle sur la matrice mat_2tp12t allant de (2n_t,n_x) à (2n_t + 1, n_t)
	* @param mat_2tp12t_view_32 : vue vectorielle sur la matrice mat_2tp12t allant de (2n_t,n_t) à (2n_t + 1, n_t + n_x)
	* @param mat_2tp12t_view_33 : vue vectorielle sur la matrice mat_2tp12t allant de (2n_t, n_t + n_x) à (2n_t + 1, 2 n_t)
	* @param vect_2t : vecteur de taille (2t) préalloué
	* @brief
	Cette fonction calcule la racine de la matrice de covariance du vecteur (t_{n|N}, t_{n+1|N}). Le résultat se trouve dans la vue sur la matrice mat_2tp12t allant de (0,0) à (2n_t,2n_t).
	**/
	void tkalman_robust_get_sqrt_corr_tt_( const gsl_vector * x_s,
										   const gsl_vector * x_s_,
										   const gsl_vector * _y,
										   const gsl_vector * y,
										   const gsl_matrix * sqrt_cov_xx_view_00,
										   const gsl_matrix * sqrt_cov_xx_view_01,
										   const gsl_matrix * sqrt_cov_xx_view_11,
										   gsl_matrix * mat_2tp12t,
										   gsl_matrix * mat_2tp12t_view_00,
										   gsl_matrix * mat_2tp12t_view_02,
										   gsl_matrix * mat_2tp12t_view_22,
										   gsl_vector * mat_2tp12t_view_30,
										   gsl_vector * mat_2tp12t_view_31,
										   gsl_vector * mat_2tp12t_view_32,
										   gsl_vector * mat_2tp12t_view_33,
										   gsl_vector * vect_2t );


	/**@fn void tkalman_robust_get_sqrt_sum( const gsl_matrix * const * sqrt_p_f,
								        const gsl_vector * const * x_s,
								       const gsl_matrix * const * sqrt_p_s,
									  const gsl_matrix * const * c_s,
									  const gsl_vector * const * y,
									  unsigned int n,
									  const gsl_matrix * f2_xx,
									  const gsl_matrix * sqrt_q2_xx,
									  const gsl_vector * vect_y_zero,
									  gsl_matrix * mat_4x4x,
									  gsl_matrix * mat_4x4x_view_00,
									  gsl_matrix * mat_4x4x_view_11,
									  gsl_matrix * mat_4x4x_view_02,
									  gsl_matrix * mat_4x4x_view_20,
									  gsl_matrix * mat_4x4x_view_22,
									  gsl_matrix * mat_4x4x_view_23,
									  gsl_matrix * mat_4x4x_view_32,
									  gsl_matrix * mat_4x4x_view_33,
									  gsl_matrix * mat_4tp12t,
									  gsl_matrix * mat_4tp12t_view_mat_4t2t,
									  gsl_matrix * mat_4tp12t_view_10b,
									  gsl_matrix * mat_4tp12t_view_10b_view_00,
									  gsl_matrix * mat_4tp12t_view_10b_view_02,
									  gsl_matrix * mat_4tp12t_view_10b_view_22,
									  gsl_vector * mat_4tp12t_view_10b_view_30,
									  gsl_vector * mat_4tp12t_view_10b_view_31,
									  gsl_vector * mat_4tp12t_view_10b_view_32,
									  gsl_vector * mat_4tp12t_view_10b_view_33,
									  gsl_vector * vect_2t,
									  gsl_vector * vect_4x)
	 * @param[in] sqrt_p_f : racines des matrices de covariance des états filtrés
	 * @param[in] x_s : espérances des états lissés
	 * @param[in] sqrt_p_s : racines des matrices de covariance des états lissés
	 * @param[in] c_s :  racines des matrices de covariance des états lissés x gain
	 * @param[in] y : observations
	 * @param[in] n : nombre d'observations
	 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} [Q^{y,y}]^{-1} F^{y,x}
	 * @param[in] sqrt_q2_xx : racine de la matrice de covariance du bruit de process réduit
	 * @param[in] vect_y_zero : vecteur de taille y nul
	*@param mat_4x4x : matrice de taille (4x.4x) préallouée
	*@param mat_4x4x_view_00 : vue sur la matrice mat_4x4x allant de(0.n_x,0.n_x)  à (1.n_x, 1.n_x)
	*@param mat_4x4x_view_11 : vue sur la matrice mat_4x4x allant de(1.n_x,1.n_x)  à (2.n_x, 2.n_x)
	*@param mat_4x4x_view_02 : vue sur la matrice mat_4x4x allant de(0.n_x,2.n_x)  à (1.n_x, 3.n_x)
	*@param mat_4x4x_view_22 : vue sur la matrice mat_4x4x allant de(2.n_x,2.n_x)  à (3.n_x, 3.n_x)
	*@param mat_4x4x_view_23 : vue sur la matrice mat_4x4x allant de(2.n_x,3.n_x)  à (3.n_x, 4.n_x)
	*@param mat_4x4x_view_32 : vue sur la matrice mat_4x4x allant de(3.n_x,2.n_x)  à (4.n_x, 3.n_x)
	*@param mat_4x4x_view_33 : vue sur la matrice mat_4x4x allant de(3.n_x,3.n_x)  à (4.n_x, 4.n_x)
	*@param mat_4tp12t : matrice de taille ((4t+1).2t) préallouée
	*@param mat_4tp12t_view_mat_4t2t : vue sur la matrice mat_4tp12t allant de(0,0) à (4t, 2t)
	*@param mat_4tp12t_view_10b : vue sur la matrice mat_4tp12t allant de (2t,0) à (4t+1,2t)
	* @param mat_4tp12t_view_10b_view_00 : vue sur la matrice mat_4tp12t allant de (2t + 0,0) à (2t + n_x, n_x)
	* @param mat_4tp12t_view_10b_view_02 : vue sur la matrice mat_4tp12t allant de (2t + 0,n_t) à (2t + n_x, n_t + n_x)
	* @param mat_4tp12t_view_10b_view_22 : vue sur la matrice mat_4tp12t allant de (2t + n_t,n_t) à 2t + (n_t + n_x, n_t + n_x)
	* @param mat_4tp12t_view_10b_view_30 : vue vectorielle sur la matrice mat_4tp12t allant de (2t + 2n_t,0) à (2t + 2n_t + 1, n_x)
	* @param mmat_4tp12t_view_10b_view_31 : vue vectorielle sur la matrice mat_4tp12t allant de (2t + 2n_t,n_x) à (2t + 2n_t + 1, n_t)
	* @param mat_4tp12t_view_10b_view_32 : vue vectorielle sur la matrice mat_4tp12t allant de (2t + 2n_t,n_t) à (2t + 2n_t + 1, n_t + n_x)
	* @param mat_4tp12t_view_10b_view_33 : vue vectorielle sur la matrice mat_4tp12t allant de (2t + 2n_t, n_t + n_x) à (2t + 2n_t + 1, 2 n_t)
	* @param vect_2t : vecteur de taille (2t) préalloué
	* @param vect_4x : vecteur de taille (4x) préalloué.
	* @brief
	Cette fonction calcule la racine de la somme des corrélations des vecteurs (t_{n|N}, t_{n+1|N}). Ce résultat est stockée dans la matrice mat_4tp12t (de (0,0) à (2n_t, 2n_t)).
	**/
	void tkalman_robust_get_sqrt_sum( const gsl_matrix * const * sqrt_p_f,
							      const gsl_vector * const * x_s,
								  const gsl_matrix * const * sqrt_p_s,
								  const gsl_matrix * const * c_s,
								  const gsl_vector * const * y,
								  unsigned int n,
								  const gsl_matrix * f2_xx,
								  const gsl_matrix * sqrt_q2_xx,
								  const gsl_vector * vect_y_zero,
								  gsl_matrix * mat_4x4x,
								  gsl_matrix * mat_4x4x_view_00,
								  gsl_matrix * mat_4x4x_view_11,
								  gsl_matrix * mat_4x4x_view_02,
								  gsl_matrix * mat_4x4x_view_20,
								  gsl_matrix * mat_4x4x_view_22,
								  gsl_matrix * mat_4x4x_view_23,
								  gsl_matrix * mat_4x4x_view_32,
								  gsl_matrix * mat_4x4x_view_33,
								  gsl_matrix * mat_4tp12t,
								  gsl_matrix * mat_4tp12t_view_mat_4t2t,
								  gsl_matrix * mat_4tp12t_view_10b,
								  gsl_matrix * mat_4tp12t_view_10b_view_00,
								  gsl_matrix * mat_4tp12t_view_10b_view_02,
								  gsl_matrix * mat_4tp12t_view_10b_view_22,
								  gsl_vector * mat_4tp12t_view_10b_view_30,
								  gsl_vector * mat_4tp12t_view_10b_view_31,
								  gsl_vector * mat_4tp12t_view_10b_view_32,
								  gsl_vector * mat_4tp12t_view_10b_view_33,
								  gsl_vector * vect_2t,
								  gsl_vector * vect_4x);
	/**@fn void tkalman_original_argmax(gsl_matrix * x0,
								 gsl_matrix * p0,
								 gsl_matrix * f,
								 gsl_matrix * q,
								 gsl_matrix * c_00,
								 gsl_matrix * c_10,
								 gsl_matrix * c_11,
								  unsigned int n,
								 const gsl_vector * x_s_0,
								 const gsl_matrix * p_s_0)
	 * @param x0 : Espérance de l'état initial (NULL pour ne pas l'estimer)
	 * @param p0 : Matrice de covariance de l'état initial (NULL pour ne pas l'estimer)
	 * @param f : Matrice d'évolution
	 * @param q : Matrice de covariance du bruit
	 * @param c_00 : Somme
	 * @param c_10 : Somme
	 * @param c_11 : Somme
	 * @param n : nombre d'observations
	 * @param x_s_0 : espérance du 1er état lissé
	 * @param p_s_0 : matrice de covariance du premier état lissé
	 * @brief
	 * Cette fonction calcule les paramètres optimaux selon la fonction auxiliaire pour le filtre de Kalman triplet.
	 */
	void tkalman_original_argmax(gsl_vector * x0,
								 gsl_matrix * p0,
								 gsl_matrix * f,
								 gsl_matrix * q,
								 gsl_matrix * c_00,
								 gsl_matrix * c_10,
								 gsl_matrix * c_11,
								  unsigned int n,
								 const gsl_vector * x_s_0,
								 const gsl_matrix * p_s_0);
	/**@fn void tkalman_original_argmax(gsl_vector * x0,
								 gsl_matrix * p0,
								 gsl_matrix * f,
								 gsl_matrix * q,
								 gsl_matrix * c_00,
								 gsl_matrix * c_10,
								 gsl_matrix * c_11,
								 const gsl_vector * x_s_0,
								 const gsl_matrix * p_s_0)
	 * @param x0 : Espérance de l'état initial (NULL pour ne pas l'estimer)
	 * @param p0 : Matrice de covariance de l'état initial (NULL pour ne pas l'estimer)
	 * @param f : Matrice d'évolution
	 * @param q : Matrice de covariance du bruit
	 * @param c_00 : Somme
	 * @param c_10 : Somme
	 * @param c_11 : Somme
	 * @param x_s_0 : espérance du 1er état lissé
	 * @param p_s_0 : matrice de covariance du premier état lissé
	 * @brief
	 * Cette fonction calcule les paramètres optimaux selon la fonction auxiliaire pour le filtre de Kalman triplet.
	 */
	void tkalman_robust_argmax(gsl_vector * x0,
							   gsl_matrix * sqrt_p0,
							   gsl_matrix * f,
							   gsl_matrix * sqrt_q,
							   gsl_matrix * sqrt_c_00,
							   gsl_matrix * sqrt_c_01,
							   gsl_matrix * sqrt_c_10,
							   unsigned int n,
							   const gsl_vector * x_s_0,
							   const gsl_matrix * sqrt_p_s_0,
							   gsl_permutation * perm_t);
	
	/**@fn void kalman_get_sqrt_corr_xx(const gsl_matrix * sqrt_cov_xx,
										const gsl_vector * x_s,
										const gsl_vector * x_s_,
										gsl_matrix * mat_2xp12x,
										gsl_matrix * mat_2xp12x_view_0,
										gsl_vector * mat_2xp12x_view_10,
										gsl_vector * mat_2xp12x_view_11,
										gsl_vector * vect_2x)
	 * @param[in] sqrt_cov_xx : racine de la matrice de covariance du vecteur \left( x_{n|N}^T,  x_{n + 1|N}^T \right)^T
	 * @param[in] x_s : espérance de l'état lissé courant
	 * @param[in] x_s_ : espérance de l'état lissé suivant
	 * @param mat_2xp12x : matrice M de taille (2x+1, 2x) préallouée .
	 * @param mat_2xp12x_view_0 : vue sur la matrice M allant de (0,0) à (2 n_x - 1, 2 n_x - 1)
	 * @param  mat_2xp12x_view_10 : vue sur la matrice M allant de (2 n_x, 0) à (2 n_x, n_x - 1)
	 * @param  mat_2xp12x_view_11 : vue sur la matrice M allant de (2 n_x, n_x - 1) à (2 n_x, 2 n_x - 1)
	 * @param vect_2x : vecteur de taille (2x) préalloué
	 * @brief
	 * Cette fonction calcul la racine de la matrice de corrélation du vecteur \left( x_{n|N}^T,  x_{n + 1|N}^T \right)^T. Le résultat est stockée dans la vue sur M mat_2xp12x_view_0 .
	 */
	void kalman_get_sqrt_corr_xx(const gsl_matrix * sqrt_cov_xx,
								 const gsl_vector * x_s,
								 const gsl_vector * x_s_,
								 gsl_matrix * mat_2xp12x,
								 gsl_matrix * mat_2xp12x_view_0,
								 gsl_vector * mat_2xp12x_view_10,
								 gsl_vector * mat_2xp12x_view_11,
								 gsl_vector * vect_2x);
	/**@fn void kalman_get_sqrt_corr_xy(const gsl_matrix * sqrt_p_s,
										const gsl_vector * x_s,
										const gsl_vector * y,
										gsl_matrix * mat_tp1t,
										gsl_matrix * mat_tp1t_view_xx,
										gsl_vector * mat_tp1t_view_10,
										gsl_vector * mat_tp1t_view_11,
										gsl_vector * vect_t)
	 * @param[in] sqrt_p_s : racine carrée de la matrice de covariance de l'état lissé courant
	 * @param[in] x_s : espérance de l'état lissé courant
	 * @param[in] y : observation courante
	 * @param mat_tp1t : matrice M de taille (2n_t+1, 2n_t) préallouée .
	 * @param mat_tp1t_view_xx : vue sur la matrice M allant de (0,0) à (n_x - 1, n_x - 1)
	 * @param mat_tp1t_view_10 : vue sur la matrice M allant de (n_t, 0) à (n_t, n_x - 1)
	 * @param mat_tp1t_view_11 : vue sur la matrice M allant de (n_t, n_x - 1) à (n_t, 2 n_t - 1)
	 * @param vect_t : vecteur de taille (n_t) préalloué
	 * @brief
	 * Cette fonction calcul la racine de la matrice de corrélation du vecteur \left( x_{n|N}^T,  y_n^T \right)^T. Le résultat est stockée dans la vue sur la matrice M allant de (0,0) à (n_t - 1, n_t - 1)
	 */
	void kalman_get_sqrt_corr_xy(const gsl_matrix * sqrt_p_s,
								 const gsl_vector * x_s,
								 const gsl_vector * y,
								 gsl_matrix * mat_tp1t,
								 gsl_matrix * mat_tp1t_view_xx,
								 gsl_vector * mat_tp1t_view_10,
								 gsl_vector * mat_tp1t_view_11,
								 gsl_vector * vect_t);
	/**@fn void kalman_get_sqrt_sum_corr_xx( const gsl_matrix * const * sqrt_p_f,
											 const gsl_vector * const * x_s,
											 const gsl_matrix * const * sqrt_p_s,
											 const gsl_matrix * const * c_s,
											 unsigned int n,
											 const gsl_matrix * f2_xx,
											 const gsl_matrix * sqrt_q2_xx,
											 gsl_matrix * mat_4x4x,
											 gsl_matrix * mat_4x4x_view_00,
											 gsl_matrix * mat_4x4x_view_11,
											 gsl_matrix * mat_4x4x_view_02,
											 gsl_matrix * mat_4x4x_view_20,
											 gsl_matrix * mat_4x4x_view_32,
											 gsl_matrix * mat_4x4x_view_33,
											 gsl_matrix * mat_4x4x_view_corr,
											 gsl_matrix * mat_4xp12x,
											 gsl_matrix * mat_4xp12x_view_2xp12x,
											 gsl_matrix * mat_2xp12x_view_0,
											 gsl_vector * mat_2xp12x_view_10,
											 gsl_vector * mat_2xp12x_view_11,
											 gsl_vector * vect_2x,
											 gsl_vector * vect_4x)
	 * @param[in] sqrt_p_f : racines des matrices de covariance des états filtrés
	 * @param[in] x_s : espérances des états lissés
	 * @param[in] sqrt_p_s : racines des matrices de covariance des états lissés
	 * @param[in] c_s :  racines des matrices de covariance des états lissés x gain
	 * @param[in] n : nombre d'observations
	 * @param[in] f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} [Q^{y,y}]^{-1} F^{y,x}
	 * @param[in] sqrt_q2_xx : racine de la matrice de covariance du bruit de process réduit
	 * @param mat_4x4x : matrice de taille (4n_x,4n_x) préallouée
	 * @param mat_4x4x_view_00 : vue sur la matrice mat_4x4x allant de(0.n_x,0.n_x)  à (1.n_x - 1, 1.n_x - 1)
	 * @param mat_4x4x_view_11 : vue sur la matrice mat_4x4x allant de(1.n_x,1.n_x)  à (2.n_x - 1, 2.n_x - 1)
	 * @param mat_4x4x_view_02 : vue sur la matrice mat_4x4x allant de(0.n_x,2.n_x)  à (1.n_x - 1, 3.n_x - 1)
	 * @param mat_4x4x_view_32 : vue sur la matrice mat_4x4x allant de(3.n_x,2.n_x)  à (4.n_x - 1, 3.n_x - 1)
	 * @param mat_4x4x_view_33 : vue sur la matrice mat_4x4x allant de(3.n_x,3.n_x)  à (4.n_x - 1, 4.n_x - 1)
	 * @param mat_4x4x_view_corr : vue sur la matrice mat_4x4x allant de(2.n_x,2.n_x)  à (4.n_x - 1, 4.n_x - 1)
	 * @param mat_4xp12x : matrice de taille (4n_x + 1, 2 n_x) préallouée
	 * @param mat_4xp12x_view_2xp12x : vue sur la matrice mat_4xp12x allant de (2n_x, 0) à (4n_x + 1, 0)
	 * @param mat_2xp12x_view_0 : vue sur la matrice M allant de (0,0) à (2 n_x - 1, 2 n_x - 1)
	 * @param  mat_2xp12x_view_10 : vue sur la matrice M allant de (2 n_x, 0) à (2 n_x, n_x - 1)
	 * @param  mat_2xp12x_view_11 : vue sur la matrice M allant de (2 n_x, n_x - 1) à (2 n_x, 2 n_x - 1)
	 * @param vect_2x : vecteur de taille (2n_x) préalloué
	 * @param vect_4x : vecteur de taille (4n_x) préalloué
	 * @brief
	 * Cette fonction calcule\n
	 * \left[ \sum_{n = 0} ^ N \textrm{corr}(\left( x_{n|N}^T,  x_{n + 1|N}^T \right)^T)  \right]^{\frac{1}{2}}
	 * Le résultat est stocké dans la matrice mat_4xp12x dans la vue allant de (0,0) à (2 n_x - 1, 2 n_x - 1)
	 * 
	 */
	void kalman_get_sqrt_sum_corr_xx( const gsl_matrix * const * sqrt_p_f,
									  const gsl_vector * const * x_s,
									  const gsl_matrix * const * sqrt_p_s,
									  const gsl_matrix * const * c_s,
									  unsigned int n,
									  const gsl_matrix * f2_xx,
									  const gsl_matrix * sqrt_q2_xx,
									  gsl_matrix * mat_4x4x,
									  gsl_matrix * mat_4x4x_view_00,
									  gsl_matrix * mat_4x4x_view_11,
									  gsl_matrix * mat_4x4x_view_02,
									  gsl_matrix * mat_4x4x_view_20,
									  gsl_matrix * mat_4x4x_view_32,
									  gsl_matrix * mat_4x4x_view_33,
									  gsl_matrix * mat_4x4x_view_corr,
									  gsl_matrix * mat_4xp12x,
									  gsl_matrix * mat_4xp12x_view_2xp12x,
									  gsl_matrix * mat_2xp12x_view_0,
									  gsl_vector * mat_2xp12x_view_10,
									  gsl_vector * mat_2xp12x_view_11,
									  gsl_vector * vect_2x,
									  gsl_vector * vect_4x);
	/**@fn void kalman_get_sqrt_sum_corr_xy( const gsl_vector * const * x_s,
											 const gsl_matrix * const * sqrt_p_s,
											 const gsl_vector * const * y,
											 unsigned int n,
											 gsl_matrix * mat_2tp1t,
											 gsl_matrix * mat_2tp1t_view_tp1t,
											 gsl_matrix * mat_tp1t_view_xx,
											 gsl_matrix * mat_tp1t_view_0,
											 gsl_vector * mat_tp1t_view_10,
											 gsl_vector * mat_tp1t_view_11,
											 gsl_vector * vect_t)
	 * @param[in] x_s : espérances des états lissés
	 * @param[in] sqrt_p_s : racines des matrices de covariance des états lissés
	 * @param[in] y : observations
	 * @param[in] n : nombre d'observations
	 * @param mat_2tp1t : matrice de taille (2 n_t + 1, n_t) préallouée
	 * @param mat_2tp1t_view_tp1t : vue sur la matrice mat_2tp1t allant de (n_t, 0) à (2n_t, n_t - 1)
	 * @param mat_tp1t_view_xx : vue sur la matrice M allant de (0,0) à (n_x - 1, n_x - 1)
	 * @param mat_tp1t_view_10 : vue sur la matrice M allant de (n_t, 0) à (n_t, n_x - 1)
	 * @param mat_tp1t_view_11 : vue sur la matrice M allant de (n_t, n_x - 1) à (n_t, 2 n_t - 1)
	 * @param vect_t : vecteur de taille (n_t) préalloué
	 * @brief
	 * Cette fonction calcule la racine de la somme des corrélations des vecteurs \left( x_{n|N}^T,  y_n^T \right)^T. Le résulat est stockée dans la vue de la matrice mat_2tp1t allant de (0,0) à (n_t - 1, n_t - 1).
	 */
	void kalman_get_sqrt_sum_corr_xy( const gsl_vector * const * x_s,
									  const gsl_matrix * const * sqrt_p_s,
									  const gsl_vector * const * y,
									  unsigned int n,
									  gsl_matrix * mat_2tp1t,
									  gsl_matrix * mat_2tp1t_view_tp1t,
									  gsl_matrix * mat_tp1t_view_xx,
									  gsl_vector * mat_tp1t_view_10,
									  gsl_vector * mat_tp1t_view_11,
									  gsl_vector * vect_t);
	/**@fn void kalman_get_arg_max( gsl_vector * x0,
									 gsl_matrix * sqrt_p0,
									 gsl_matrix * f_xx,
									 gsl_matrix * f_yx,						 
									 gsl_matrix * sqrt_q_xx,
									 gsl_matrix * sqrt_q_yy,
									 gsl_matrix * sqrt_c_xx_view_00,
									 gsl_matrix * sqrt_c_xx_view_01,
									 gsl_matrix * sqrt_c_xx_view_11,
									 gsl_matrix * sqrt_c_xy,
									 gsl_matrix * sqrt_c_xy_view_00,
									 gsl_matrix * sqrt_c_xy_view_01,
									 gsl_matrix * sqrt_c_xy_view_11
									 unsigned int n,
									 const gsl_vector * x_s_0,
									 const gsl_matrix * p_s_0,
									 gsl_permutation * perm_x
								   )
	 * @param x0 : Espérance de l'état initial (NULL pour ne pas l'estimer)
	 * @param p0 : Matrice de covariance de l'état initial (NULL pour ne pas l'estimer)
	 * @param f : Matrice d'évolution
	 * @param q : Matrice de covariance du bruit
	 * @param c_00 : Somme
	 * @param c_10 : Somme
	 * @param c_11 : Somme
	 * @param x_s_0 : espérance du 1er état lissé
	 * @param p_s_0 : matrice de covariance du premier état lissé
	 * 
	 **/
	void kalman_get_arg_max( gsl_vector * x0,
							 gsl_matrix * sqrt_p0,
							 gsl_matrix * f_xx,
							 gsl_matrix * f_yx,						 
							 gsl_matrix * sqrt_q_xx,
							 gsl_matrix * sqrt_q_yy,
							 gsl_matrix * sqrt_c_xx_view_00,
							 gsl_matrix * sqrt_c_xx_view_01,
							 gsl_matrix * sqrt_c_xx_view_11,
							 gsl_matrix * sqrt_c_xy_view_00,
							 gsl_matrix * sqrt_c_xy_view_01,
							 gsl_matrix * sqrt_c_xy_view_11,
							 unsigned int n,
							 const gsl_vector * x_s_0,
							 const gsl_matrix * p_s_0,
							 gsl_permutation * perm_x
						   );


#endif
