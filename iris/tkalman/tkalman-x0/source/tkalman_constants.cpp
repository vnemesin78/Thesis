/**@file tkalman_constants.hpp
 * @author Valérian Némesin
 * @brief
 Ce fichier contient le code soure des fonctions qui calculent les différentes constantes du filtre de Kalman triplet.
 */
#include "tkalman_constants.hpp"
/**@fn void tkalman_original_get_q2_xy(gsl_matrix * q2_xy,
									   const gsl_matrix * q_xy,
									   const gsl_matrix * q_yy,
									   gsl_matrix * mat_yy)
 * @param q2_xy : q_2^{x,y} =  Q^{x,y} \: [Q^{y,y}]^{-1}
 * @param[in] q_xy : Matrice de covariance entre bruit de process et bruit de mesure
 * @param[in] q_yy : Matrice de covariance du bruit de mesure
 * @param mat_yy : Matrice de taille (y.y) préallouée.
 * @brief
 Cette fonction calcule la matrice q_2^{x,y} =  Q^{x,y} \: [Q^{y,y}]^{-1}
 */
void tkalman_original_get_q2_xy(gsl_matrix * q2_xy,
								const gsl_matrix * q_xy,
								const gsl_matrix * q_yy,
								gsl_matrix * mat_yy)
{
	//Copie

	gsl_matrix_memcpy (mat_yy, q_yy);

	//Cholesky
	gsl_linalg_cholesky_decomp (mat_yy);

	//Inversion
	gsl_linalg_cholesky_invert (mat_yy);

	//Produit
	gsl_blas_dsymm (CblasRight,
					CblasUpper,
					1.0,
					mat_yy,
					q_xy,
					0,
					q2_xy);
}

/**@fn void tkalman_original_get_q2_xx(gsl_matrix * q2_xx,
									   const gsl_matrix * q_xx,
									   const gsl_matrix * q_xy,
									   const gsl_matrix * q2_xy)
 * @param q2_xx : Matrice de covariance du bruit de process réduit
 * @param[in] q_xx : Matrice de covariance du bruit de process
 * @param[in] q_xy : Matrice de covariance entre bruit de process et bruit de mesure
 * @param[in] q2_xy : q_2^{x,y} =  Q^{x,y} \: [Q^{y,y}]^{-1}
 * @brief
 Cette fonction calcule le bruit de process réduit définit par la formule :
 * Q_2^{x,x} = Q^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: Q^{y,x}
 *
 */
void tkalman_original_get_q2_xx(gsl_matrix * q2_xx,
							 	const gsl_matrix * q_xx,
							 	const gsl_matrix * q_xy,
							 	const gsl_matrix * q2_xy)
{
	//Copie de Qxx

	gsl_matrix_memcpy(q2_xx, q_xx);

	//Qxx - Qxy . Qyy¯¹ . Qyx

	gsl_blas_dgemm (CblasNoTrans,
					CblasTrans,
					-1.0,
					q2_xy,
					q_xy,
					1.0,
					q2_xx);


}

/**@fn void tkalman_original_get_f2_x_(gsl_matrix * f2_x_,
							   		   const gsl_matrix * f_x_,
							    	   const gsl_matrix * f_y_,
							 		   const gsl_matrix * q2_xy)
 * @param f2_x_ : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
 * @param[in] f_x_ : terme de la matrice d'évolution
 * @param[in] f_y_ : terme de la matrice d'évolution
 * @param[in] q2_xy : q_2^{x,y} =  Q^{x,y} \: [Q^{y,y}]^{-1}
 * @brief
 Cette fonction calcule F_2^{x,x} définit par la formule :
 * F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
 *
 */
void tkalman_original_get_f2_x_(gsl_matrix * f2_x_,
							    const gsl_matrix * f_x_,
							    const gsl_matrix * f_y_,
							 	const gsl_matrix * q2_xy)
{
	//Copie de Qxx
	gsl_matrix_memcpy(f2_x_, f_x_);

	//Qxx - Qxy . Qyy¯¹ . Qyx

	gsl_blas_dgemm (CblasNoTrans,
					CblasNoTrans,
					-1.0,
					q2_xy,
					f_y_,
					1.0,
					f2_x_);
}

/**@fn void tkalman_original_get_constants(gsl_matrix * f2_xx,
									   gsl_matrix * f2_xy,
									   gsl_matrix * q2_xx,
									   gsl_matrix * q2_xy,
									   const gsl_matrix * f_xx,
									   const gsl_matrix * f_xy,
									   const gsl_matrix * f_yx,
									   const gsl_matrix * f_yy,
									   const gsl_matrix * q_xx,
									   const gsl_matrix * q_xy,
									   const gsl_matrix * q_yy,
									   gsl_matrix * mat_yy)
 * @param f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
 * @param f2_xy : F_2^{x,y} = F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y}
 * @param q2_xy : q_2^{x,y} =  Q^{x,y} \: [Q^{y,y}]^{-1}
 * @param q2_xx : Matrice de covariance du bruit de process réduit Q_2^{x,x} = Q^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: Q^{y,x}
 * @param[in] f_xx : terme de la matrice d'évolution
 * @param[in] f_xy : terme de la matrice d'évolution
 * @param[in] f_yx : terme de la matrice d'évolution
 * @param[in] f_yy : terme de la matrice d'évolution
 * @param[in] q_xx : Matrice de covariance du bruit de process
 * @param[in] q_xy : Matrice de covariance entre bruit de process et bruit de mesure
 * @param[in] q_yy : Matrice de covariance du bruit de mesure
 * @param mat_yy : Matrice de taille (y.y) préallouée
 * @brief
 Cette fonction calcule les constantes du filtre de Kalman triplet. (Voir les arguments pour les formules)
 */
void tkalman_original_get_constants(gsl_matrix * f2_xx,
									gsl_matrix * f2_xy,
									gsl_matrix * q2_xx,
									gsl_matrix * q2_xy,
									const gsl_matrix * f_xx,
									const gsl_matrix * f_xy,
									const gsl_matrix * f_yx,
									const gsl_matrix * f_yy,
									const gsl_matrix * q_xx,
									const gsl_matrix * q_xy,
									const gsl_matrix * q_yy,
									gsl_matrix * mat_yy)
{
	tkalman_original_get_q2_xy(q2_xy,
							   q_xy,
							   q_yy,
							   mat_yy);
	tkalman_original_get_q2_xx(q2_xx,
							   q_xx,
							   q_xy,
							   q2_xy);
	tkalman_original_get_f2_x_(f2_xx,
							   f_xx,
							   f_yx,
							   q2_xy);
	tkalman_original_get_f2_x_(f2_xy,
							   f_xy,
							   f_yy,
							   q2_xy);
}

/**@fn void tkalman_robust_get_sqrt_q2_xx_sqrt_q_yy_and_q2_xy_from_sqrt_q(gsl_matrix * sqrt_q2_xx,
													      gsl_matrix * q2_xy,
													  gsl_matrix * sqrt_q_yy,
													  const gsl_matrix * sqrt_q_view_xx,
													  const gsl_matrix * sqrt_q_view_xy,
													  const gsl_matrix * sqrt_q_view_yy,
													  gsl_matrix * mat_tt,
													  gsl_matrix * mat_tt_yy,
													  gsl_matrix * mat_tt_yx,
													  gsl_matrix * mat_tt_xy,
													  gsl_matrix * mat_tt_xx,
													  gsl_vector * vect_t,
													  gsl_permutation * perm_y)
 * @param sqrt_q2_xx : racine de la matrice de covariance réduite du bruit de process
 * @param q2_xy : terme de la matrice de passage qui permet d'annuler la corrélation entre bruit de mesure et bruit de process
 * @param sqrt_q_yy : racine de la matrice de covariance du bruit de mesure
 * @param[in] sqrt_q_view_xx : vue sur la matrice sqrt_q  (de (0,0) à (n_x, n_x))
 * @param[in] sqrt_q_view_xy : vue sur la matrice sqrt_q (de (0, n_x) à (n_x, n_t))
 * @param[in] sqrt_q_view_yy :  vue sur la matrice sqrt_q (de (n_x, n_x) à (n_t, n_t))
 * @param mat_tt : matrice de taille (t.t) préallouée.
 * @param mat_tt_yy : vue sur la matrice mat_tt (de (0,0) à (n_y, n_y))
 * @param mat_tt_yx : vue sur la matrice mat_tt (de (0, n_y) à (n_y, n_t))
 * @param mat_tt_xy : vue sur la matrice mat_tt (de (n_y, 0) à (n_t, n_y))
 * @param mat_tt_xx : vue sur la matrice mat_tt (de (n_y, n_y) à (n_t, n_t))
 * @param vect_t : vecteur de taille t préalloué pour le calcul
 * @param perm_y : permutation de taille y préallouée
 * @brief
 Cette fonction calcule les racines des matrices de covariance Q2xx et Qyy. Elle calcule aussi la matrice Q2xy = Qxy.Qyy^-1.

**/
void tkalman_robust_get_sqrt_q2_xx_sqrt_q_yy_and_q2_xy_from_sqrt_q(gsl_matrix * sqrt_q2_xx,
																   gsl_matrix * q2_xy,
																   gsl_matrix * sqrt_q_yy,
																   const gsl_matrix * sqrt_q_view_xx,
																   const gsl_matrix * sqrt_q_view_xy,
																   const gsl_matrix * sqrt_q_view_yy,
																   gsl_matrix * mat_tt,
																   gsl_matrix * mat_tt_yy,
																   gsl_matrix * mat_tt_yx,
																   gsl_matrix * mat_tt_xy,
																   gsl_matrix * mat_tt_xx,
																   gsl_vector * vect_t,
																   gsl_permutation * perm_y)
{
	//Construction de la matrice Q'
	// sqrt(Q)yy		0
	// sqrt(Q)xy		sqrt(Qxx) = sqrt(Q)xx
	gsl_matrix_memcpy(mat_tt_yy, sqrt_q_view_yy);
	gsl_matrix_set_zero(mat_tt_yx);
	gsl_matrix_memcpy(mat_tt_xy, sqrt_q_view_xy);
	gsl_matrix_memcpy(mat_tt_xx, sqrt_q_view_xx);
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_tt,
						 vect_t);
	gsl_triangle_matrix(mat_tt);

	//Récupération de sqrt(Qyy) et de sqrt(Q2xx)
	gsl_matrix_memcpy(sqrt_q2_xx, mat_tt_xx);
	gsl_matrix_memcpy(sqrt_q_yy, mat_tt_yy);

	//Calcul de Q2xy = Qyx.Qyy^-1
	gsl_permutation_init(perm_y);
    gsl_linalg_LU_invert(sqrt_q_yy, perm_y, mat_tt_yy);

	gsl_blas_dgemm(CblasTrans,
				   CblasTrans,
				   1.0,
				   mat_tt_yx,
				   mat_tt_yy,
				   0.0,
				   q2_xy);
}

/**@fn void tkalman_robust_get_constants(gsl_matrix * f2_xx,
							  gsl_matrix * f2_xy,
							  gsl_matrix * sqrt_q2_xx,
							  gsl_matrix * q2_xy,
							  gsl_matrix * sqrt_q_yy,
							  const gsl_matrix * f_xx,
							  const gsl_matrix * f_xy,
							  const gsl_matrix * f_yx,
							  const gsl_matrix * f_yy,
							  const gsl_matrix * sqrt_q_view_xx,
							  const gsl_matrix * sqrt_q_view_xy,
							  const gsl_matrix * sqrt_q_view_yy,
							  gsl_matrix * mat_tt,
							  gsl_matrix * mat_tt_yy,
							  gsl_matrix * mat_tt_yx,
							  gsl_matrix * mat_tt_xy,
							  gsl_matrix * mat_tt_xx,
							  gsl_vector * vect_t,
							  gsl_permutation * perm_y)
 * @param f2_xx : F_2^{x,x} = F^{x,x} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,x}
 * @param f2_xy : F_2^{x,y} = F^{x,y} - Q^{x,y} \: [Q^{y,y}]^{-1} \: F^{y,y}
 * @param sqrt_q2_xx : racine de la matrice de covariance réduite du bruit de process
 * @param q2_xy : terme de la matrice de passage qui permet d'annuler la corrélation entre bruit de mesure et bruit de process
 * @param sqrt_q_yy : racine de la matrice de covariance du bruit de mesure
 * @param[in] f_xx : terme de la matrice d'évolution
 * @param[in] f_xy : terme de la matrice d'évolution
 * @param[in] f_yx : terme de la matrice d'évolution
 * @param[in] f_yy : terme de la matrice d'évolution
 * @param[in] sqrt_q_view_xx : vue sur la matrice sqrt_q  (de (0,0) à (n_x, n_x))
 * @param[in] sqrt_q_view_xy : vue sur la matrice sqrt_q (de (0, n_x) à (n_x, n_t))
 * @param[in] sqrt_q_view_yy :  vue sur la matrice sqrt_q (de (n_x, n_x) à (n_t, n_t))
 * @param mat_tt : matrice de taille (t.t) préallouée.
 * @param mat_tt_yy : vue sur la matrice mat_tt (de (0,0) à (n_y, n_y))
 * @param mat_tt_yx : vue sur la matrice mat_tt (de (0, n_y) à (n_y, n_t))
 * @param mat_tt_xy : vue sur la matrice mat_tt (de (n_y, 0) à (n_t, n_y))
 * @param mat_tt_xx : vue sur la matrice mat_tt (de (n_y, n_y) à (n_t, n_t))
 * @param vect_t : vecteur de taille t préalloué pour le calcul
 * @param perm_y : permutation de taille y préallouée
 * @brief
 Cette fonction calcule de manière robuste les constantes du filtre de Kalman triplet.

**/
void tkalman_robust_get_constants(gsl_matrix * f2_xx,
								  gsl_matrix * f2_xy,
								  gsl_matrix * sqrt_q2_xx,
								  gsl_matrix * q2_xy,
								  gsl_matrix * sqrt_q_yy,
								  const gsl_matrix * f_xx,
								  const gsl_matrix * f_xy,
								  const gsl_matrix * f_yx,
								  const gsl_matrix * f_yy,
								  const gsl_matrix * sqrt_q_view_xx,
								  const gsl_matrix * sqrt_q_view_xy,
								  const gsl_matrix * sqrt_q_view_yy,
								  gsl_matrix * mat_tt,
								  gsl_matrix * mat_tt_yy,
								  gsl_matrix * mat_tt_yx,
								  gsl_matrix * mat_tt_xy,
								  gsl_matrix * mat_tt_xx,
								  gsl_vector * vect_t,
								  gsl_permutation * perm_y)
{
	//Récup. des racines
	tkalman_robust_get_sqrt_q2_xx_sqrt_q_yy_and_q2_xy_from_sqrt_q(sqrt_q2_xx,
																  q2_xy,
																  sqrt_q_yy,
																  sqrt_q_view_xx,
																  sqrt_q_view_xy,
																  sqrt_q_view_yy,
																  mat_tt,
																  mat_tt_yy,
																  mat_tt_yx,
																  mat_tt_xy,
																  mat_tt_xx,
																  vect_t,
																  perm_y);

	//Calcul des deux autres constantes
	tkalman_original_get_f2_x_(f2_xx,
							   f_xx,
							   f_yx,
							   q2_xy);

	tkalman_original_get_f2_x_(f2_xy,
							   f_xy,
							   f_yy,
							   q2_xy);

}
