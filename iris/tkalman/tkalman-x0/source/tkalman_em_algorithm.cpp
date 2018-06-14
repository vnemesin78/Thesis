 /**@file tkalman_em_algorithm.cpp
 * @author Valérian Némesin
 * @brief
 Ce fichier contient le code source des fonctions nécessaires à l'algorithme EM du filtre de Kalman triplet.
**/
#include "tkalman_em_algorithm.hpp"
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
										const gsl_matrix * cov)
{
	gsl_matrix_add(corr_xx,
				   cov);

	gsl_blas_dger(1.0, t_n, t_p, corr);
}

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
								gsl_vector * vect_t_2_y)
{

	gsl_matrix_set_zero(c_00);
	gsl_matrix_set_zero(c_10);
	gsl_matrix_set_zero(c_11);
	//La rang de la matrice Q est déterminé par P0 donc pour le préserver, je calcule de N à 0
	for (unsigned int i = n - 1; i >= 1; -- i)
	{
		gsl_vector_memcpy(vect_t_1_x,
						  x_s[i]);
		gsl_vector_memcpy(vect_t_1_y,
						  y[i - 1]);

		gsl_vector_memcpy(vect_t_2_x,
						  x_s[i + 1]);
		gsl_vector_memcpy(vect_t_2_y,
						  y[i]);

		tkalman_original_get_corr_t_n_t_p(c_00,
										  c_00_xx,
										  vect_t_1,
										  vect_t_1,
										  p_s[i]);

		tkalman_original_get_corr_t_n_t_p(c_10,
										  c_10_xx,
										  vect_t_2,
										  vect_t_1,
										  c_s[i]);

		tkalman_original_get_corr_t_n_t_p(c_11,
										  c_11_xx,
										  vect_t_2,
										  vect_t_2,
										  p_s[i + 1]);
	}

	gsl_vector_memcpy(vect_t_1_x,
					  x_s[0]);
	gsl_vector_set_zero(vect_t_1_y);


	gsl_vector_memcpy(vect_t_2_x,
					  x_s[1]);
	gsl_vector_memcpy(vect_t_2_y,
					  y[0]);

	tkalman_original_get_corr_t_n_t_p(c_00,
									  c_00_xx,
									  vect_t_1,
									  vect_t_1,
									  p_s[0]);

	tkalman_original_get_corr_t_n_t_p(c_10,
									  c_10_xx,
									  vect_t_2,
									  vect_t_1,
									  c_s[0]);

	tkalman_original_get_corr_t_n_t_p(c_11,
									  c_11_xx,
									  vect_t_2,
									  vect_t_2,
									  p_s[1]);

}

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
									  gsl_vector * vect_4x)
{
	//Construction de la matrice
	// sqrt(pf)F2xx^T	0	sqrt(pf)	0
	// 0				I	0		0
	// sqrt(Q2xx)		0	0		0
	// 0				0	c_s		sqrt(ps_)
		gsl_matrix_set_zero(mat_4x4x);
		gsl_blas_dgemm(CblasNoTrans,
					   CblasTrans,
					   1.0,
					   sqrt_p_f,
					   f2_xx,
					   0.0,
					   mat_4x4x_view_00);
		gsl_matrix_set_identity(mat_4x4x_view_11);
		gsl_matrix_memcpy(mat_4x4x_view_02, sqrt_p_f);
		gsl_matrix_memcpy(mat_4x4x_view_20, sqrt_q2_xx);
		gsl_matrix_memcpy(mat_4x4x_view_32, c_s);
		gsl_matrix_memcpy(mat_4x4x_view_33, sqrt_p_s_);
		

	
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_4x4x,
						 vect_4x);
	gsl_triangle_matrix(mat_4x4x);

	//Le résultat est déjà dans sqrt_cov_xx
}

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
									   gsl_vector * vect_2t )
{

	//Construction de la matrice
	// sqrt(cov_xx)00		0		sqrt(cov_xx)01		0
	// 0				0		0				0
	// 0				0		sqrt(cov_xx)11		0
	// 0				0		0				0
	//x_s				_y		x_s_				y
	gsl_matrix_set_zero(mat_2tp12t);
	gsl_matrix_memcpy(mat_2tp12t_view_00, sqrt_cov_xx_view_00);
	gsl_matrix_memcpy(mat_2tp12t_view_02, sqrt_cov_xx_view_01);
	gsl_matrix_memcpy(mat_2tp12t_view_22, sqrt_cov_xx_view_11);
	gsl_vector_memcpy(mat_2tp12t_view_30, x_s);
	gsl_vector_memcpy(mat_2tp12t_view_31, _y);
	gsl_vector_memcpy(mat_2tp12t_view_32, x_s_);
	gsl_vector_memcpy(mat_2tp12t_view_33, y);
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_2tp12t,
						 vect_2t);
	gsl_triangle_matrix(mat_2tp12t);

	//Le résultat est déjà dans la matrice sqrt_cor_tt
}


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
								  gsl_vector * vect_4x)
{
	//Mise à zéro
	gsl_matrix_set_zero(mat_4tp12t_view_mat_4t2t);

	for(unsigned i = n; i > 1; --i)
	{
		tkalman_robust_get_sqrt_cov_xx_( sqrt_p_f[i - 1],
										 sqrt_p_s[i],
										 c_s[i - 1],
										 f2_xx,
										 sqrt_q2_xx,
										 mat_4x4x,
										 mat_4x4x_view_00,
										 mat_4x4x_view_11,
										 mat_4x4x_view_02,
										 mat_4x4x_view_20,
										 mat_4x4x_view_32,
										 mat_4x4x_view_33,
										 vect_4x);

		tkalman_robust_get_sqrt_corr_tt_(x_s[i - 1],
										 x_s[i],
										 y[i - 2],
										 y[i - 1],
										 mat_4x4x_view_22,
										 mat_4x4x_view_23,
										 mat_4x4x_view_33,
										 mat_4tp12t_view_10b,
										 mat_4tp12t_view_10b_view_00,
										 mat_4tp12t_view_10b_view_02,
										 mat_4tp12t_view_10b_view_22,
										 mat_4tp12t_view_10b_view_30,
										 mat_4tp12t_view_10b_view_31,
										 mat_4tp12t_view_10b_view_32,
										 mat_4tp12t_view_10b_view_33,
										 vect_2t );

		//Décomposition QR
		gsl_linalg_QR_decomp(mat_4tp12t_view_mat_4t2t,
							 vect_2t);
		gsl_triangle_matrix(mat_4tp12t_view_mat_4t2t);
	}

	tkalman_robust_get_sqrt_cov_xx_( sqrt_p_f[0],
									 sqrt_p_s[1],
									 c_s[0],
									 f2_xx,
									 sqrt_q2_xx,
									 mat_4x4x,
									 mat_4x4x_view_00,
									 mat_4x4x_view_11,
									 mat_4x4x_view_02,
									 mat_4x4x_view_20,
									 mat_4x4x_view_32,
									 mat_4x4x_view_33,
									 vect_4x);

	tkalman_robust_get_sqrt_corr_tt_(x_s[0],
									 x_s[1],
									 vect_y_zero,
									 y[0],
									 mat_4x4x_view_22,
									 mat_4x4x_view_23,
									 mat_4x4x_view_33,
									 mat_4tp12t_view_10b,
									 mat_4tp12t_view_10b_view_00,
									 mat_4tp12t_view_10b_view_02,
									 mat_4tp12t_view_10b_view_22,
									 mat_4tp12t_view_10b_view_30,
									 mat_4tp12t_view_10b_view_31,
									 mat_4tp12t_view_10b_view_32,
									 mat_4tp12t_view_10b_view_33,
									 vect_2t );

	//Décomposition QR
	gsl_linalg_QR_decomp(mat_4tp12t_view_mat_4t2t,
						 vect_2t);
	gsl_triangle_matrix(mat_4tp12t_view_mat_4t2t);
	//Le résultat est déjà dans la matrice sqrt_sum
}

/**@fn void tkalman_original_argmax(gsl_vector * x0,
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
							 const gsl_matrix * p_s_0)
{
	//Réestimation falcutative
		//Calcul de x0, p0
			if (x0 != NULL && x_s_0 != NULL)
				gsl_vector_memcpy(x0,
								  x_s_0);


			if (p0 != NULL && p_s_0 != NULL)
				gsl_matrix_memcpy(p0,
								  p_s_0);

	//Réestimation des paramètres importants
		//Calcul de c_00?¹

			gsl_linalg_cholesky_decomp(c_00);
			gsl_linalg_cholesky_invert(c_00);
		//Calcul de F
			//F = 0.0 F + 1.0 x C_10 x C_00?¹
			gsl_blas_dgemm(CblasNoTrans,
						   CblasNoTrans,
						   1.0,
						   c_10,
						   c_00,
						   0.0,
						   f);
		//Calcul de Q
			gsl_matrix_memcpy(q,
							  c_11);
			gsl_blas_dgemm(CblasNoTrans,
						   CblasTrans,
						   -1.0,
						   f,
						   c_10,
						   1.0,
						   q);
			gsl_matrix_scale(q,
							 1.0/((double) n));

}





/**@fn void tkalman_robust_argmax(gsl_vector * x0,
						   gsl_matrix * sqrt_p0,
						   gsl_matrix * f,
						   gsl_matrix * sqrt_q,
						   gsl_matrix * sqrt_c_00,
						   gsl_matrix * sqrt_c_01,
						   gsl_matrix * sqrt_c_11,
						   unsigned int n,
						   const gsl_vector * x_s_0,
						   const gsl_matrix * sqrt_p_s_0,
						   gsl_permutation * perm_t)
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
						   gsl_matrix * sqrt_c_11,
						   unsigned int n,
						   const gsl_vector * x_s_0,
						   const gsl_matrix * sqrt_p_s_0,
						   gsl_permutation * perm_t)
{
	
	
	//Réestimation falcutative
		//Calcul de x0, p0
			if (x0 != NULL && x_s_0 != NULL)
				gsl_vector_memcpy(x0,
								  x_s_0);


			if (sqrt_p0 != NULL && sqrt_p_s_0 != NULL)
				gsl_matrix_memcpy(sqrt_p0,
								  sqrt_p_s_0);

		//Calcul de F
			//Inversion de sqrt_c_00
				gsl_permutation_init(perm_t);
				gsl_linalg_LU_invert(sqrt_c_00,
									 perm_t,
									 sqrt_q);
			//Calcul de F = sqrt_c01^T [C00^{1/2}]^{-1}
				gsl_blas_dgemm(CblasTrans,
							   CblasTrans,
							   1.0,
							   sqrt_c_01,
							   sqrt_q,
							   0.0,
							   f);

		//Calcul de Q
		gsl_matrix_memcpy(sqrt_q,
						  sqrt_c_11);

		gsl_matrix_scale(sqrt_q,
						 sqrt(1.0/n));




}

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
							 gsl_vector * vect_2x)
{
	//Construction de la matrice
	gsl_matrix_memcpy(mat_2xp12x_view_0, sqrt_cov_xx);
	gsl_vector_memcpy(mat_2xp12x_view_10, x_s);
	gsl_vector_memcpy(mat_2xp12x_view_11, x_s_);
	
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_2xp12x,
						 vect_2x);
	gsl_triangle_matrix(mat_2xp12x);

}						

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
							 gsl_vector * vect_t)
{
	//Construction de la matrice
	gsl_matrix_set_zero(mat_tp1t);
	gsl_matrix_memcpy(mat_tp1t_view_xx, sqrt_p_s);
	gsl_vector_memcpy(mat_tp1t_view_10, x_s);
	gsl_vector_memcpy(mat_tp1t_view_11, y);
	
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_tp1t,
						 vect_t);
	gsl_triangle_matrix(mat_tp1t);
}						
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
								  gsl_vector * vect_4x)
{
	//Mise à zéro de la somme
	gsl_matrix_set_zero(mat_4xp12x);
	
	
	for (unsigned int i = n - 1; i != (unsigned int) (- 1) ; -- i)
	{
		//Obtiention de la matrice de covariance
		tkalman_robust_get_sqrt_cov_xx_(sqrt_p_f[i],
										sqrt_p_s[i + 1],
										c_s[i],
										f2_xx,
										sqrt_q2_xx,
										mat_4x4x,
										mat_4x4x_view_00,
										mat_4x4x_view_11,
										mat_4x4x_view_02,
										mat_4x4x_view_20,
										mat_4x4x_view_32,
										mat_4x4x_view_33,
										vect_4x);

										
										
		//Obtention de la matrice de corrélation
		kalman_get_sqrt_corr_xx(mat_4x4x_view_corr,
								x_s[i],
								x_s[i + 1],
								mat_4xp12x_view_2xp12x,
								mat_2xp12x_view_0,
								mat_2xp12x_view_10,
								mat_2xp12x_view_11,
								vect_2x);

		//Mise à jour de la racine des sommes
		gsl_linalg_QR_decomp(mat_4xp12x,
							 vect_2x);
		gsl_triangle_matrix(mat_4xp12x);
		

	}
}
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
								  gsl_vector * vect_t)
{
	
	
	//Mise à zéro de la somme
	gsl_matrix_set_zero(mat_2tp1t);
	
	
	for (unsigned int i = n - 1; i != (unsigned int) (- 1) ; -- i)
	{
		//Obtiention de la matrice de covariance
		kalman_get_sqrt_corr_xy(sqrt_p_s[i],
								x_s[i],
								y[i],
								mat_2tp1t_view_tp1t,
								mat_tp1t_view_xx,
								mat_tp1t_view_10,
								mat_tp1t_view_11,
								vect_t);
		//Mise à jour de la racine des sommes
		gsl_linalg_QR_decomp(mat_2tp1t,
							 vect_t);
		gsl_triangle_matrix(mat_2tp1t);
	}
}



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
						 const gsl_matrix * sqrt_p_s_0,
						 gsl_permutation * perm_x
					   )
{
	//Calcul de Cxx^-1
	
	gsl_permutation_init(perm_x);
	gsl_linalg_LU_invert(sqrt_c_xx_view_00,
						 perm_x,
						 sqrt_c_xy_view_00);

	//Calcul de Q
	gsl_matrix_memcpy(sqrt_q_xx,
					  sqrt_c_xx_view_11);
	gsl_matrix_scale(sqrt_q_xx,
					 1.0/ sqrt( n));
					 
					 
	gsl_matrix_memcpy(sqrt_q_yy,
					  sqrt_c_xy_view_11);
	gsl_matrix_scale(sqrt_q_yy,
					 1.0/ sqrt( n));
	
	//Calcul de Fxx
	gsl_blas_dgemm(CblasTrans,
				   CblasTrans,
				   1.0,
				   sqrt_c_xx_view_01,
				   sqrt_c_xy_view_00,
				   0.0,
				   f_xx);
	//Calcul de Fyy
	gsl_blas_dgemm(CblasTrans,
				   CblasTrans,
				   1.0,
				   sqrt_c_xy_view_01,
				   sqrt_c_xy_view_00,
				   0.0,
				   f_yx);		   
				   
	//Réestimation falcutative
		//Calcul de x0, p0
			if (x0 != NULL && x_s_0 != NULL)
				gsl_vector_memcpy(x0,
								  x_s_0);


			if (sqrt_p0 != NULL && sqrt_p_s_0 != NULL)
				gsl_matrix_memcpy(sqrt_p0,
								  sqrt_p_s_0);	   
}

