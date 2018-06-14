#include "tkalman_nc_filtering.hpp"
#include <iostream>
using namespace std;
//Espérance de l'innovation
/**@fn void tkalman_nc_get_innovation(gsl_vector * innovation,
									  const gsl_vector * x_p,
									  const gsl_vector * _y,
									  const gsl_vector * y,
									  const gsl_matrix * f_yx,
									  const gsl_matrix * f_yy)
 * @param innovation : \tilde{y}_{n}, innovation
 * @param[in] x_p : \hat{x}_{n | n - 1}, espérance de l'état prédit
 * @param[in] _y : \hat{y}_{n - 1}, espérance de l'observation précédente
 * @param[in] y : \hat{y}_n, espérance de l'observation courante
 * @param[in] f_yx : F^{y,x}, terme de la matrice d'évolution
 * @param[in] f_yy : F^{y,y}, terme de la matrice d'évolution
 * @brief
 Cette fonction calcule l'innovation selon la formule :
 \tilde{y}_n = y_{n} - F^{y,x} \: \hat{x}_{n | n - 1} - F^{y,y} \: y_{n - 1}
 */
void tkalman_nc_get_innovation(gsl_vector * innovation,
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

/**@fn void tkalman_nc_get_innovation_0(gsl_vector * innovation,
									 const gsl_vector * t_0,
									 const gsl_vector * y_0,
									 const gsl_matrix * f_yt)
 * @param innovation : \tilde{y}_{0}, innovation
 * @param[in] t_p : \hat{t}_{0}
 * @param[in] y_0 : y_0, espérance de l'observation courante
 * @param[in] f_yt : F^{y,t}, terme de la matrice d'évolution
 * @brief
 Cette fonction calcule l'innovation selon la formule :
 \tilde{y}_{0} = y_{0} - F^{y,t} \hat{t}_{0}
 */
void tkalman_nc_get_innovation_0(gsl_vector * innovation,
							  const gsl_vector * t_0,
							  const gsl_vector * y_0,
							  const gsl_matrix * f_yt)
{
	gsl_vector_memcpy(innovation,
					  y_0);

	gsl_blas_dgemv (CblasNoTrans,
					-1.0,
					f_yt,
					t_0,
					1.0,
					innovation);
}

//Espérance de l'état filtré
/**@fn void tkalman_nc_get_x_f(gsl_vector * x_f,
									 const gsl_vector * x_p,
									 const gsl_vector * innovation,
									 const gsl_matrix * gain)
 * @param x_f : \hat{x}_{n|n}, espérance de l'état filtré
 * @param[in] x_p : \hat{x}_{n|n - 1}, espérance de l'état prédit
 * @param[in] innovation : \tilde{y}_{n}, innovation
 * @param[in] gain : K_{n|n}, gain de filtrage
 * @brief
 Cette fonction calcule l'espérance de l'état filtré selon la formule :
 \hat{x}_{n|n} = \hat{x|n-1} + K_{n|n} \; \tilde{y}_n
 */
void tkalman_nc_get_x_f(gsl_vector * x_f,
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


/**@fn void tkalman_nc_get_t_0_f(gsl_vector * t_0_f,
								 const gsl_vector * t_0,
								 const gsl_vector * innovation_0,
								 const gsl_matrix * gain)
 * @param t_0_f: \hat{t}_{0|0}
 * @param[in] t_0 : \hat{t}_0,
 * @param[in] innovation_0 : \tilde{y}_{0}, innovation
 * @param[in] gain : K^t_{0|0}, gain de filtrage
 * @brief
 Cette fonction calcule l'espérance de l'état filtré selon la formule :
 \hat{t}_{0|0} = \hat{t}_0 + K^t_{0|0} \; \tilde{y}_{0}
 */
void tkalman_nc_get_t_0_f(gsl_vector * t_0_f,
						  const gsl_vector * t_0,
						  const gsl_vector * innovation_0,
						  const gsl_matrix * gain)
{
	gsl_vector_memcpy(t_0_f, t_0);
	gsl_blas_dgemv (CblasNoTrans,
					1.0,
					gain,
					innovation_0,
					1.0,
					t_0_f);

}

//Racine de la matrice de cov. de l'état filtré
/**@fn void tkalman_nc_get_sqrt_pf_sqrt_s_and_gain(gsl_matrix * sqrt_p_f,
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
 * @param sqrt_p_f :  [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré
 * @param sqrt_s : [S_{n}]^{\frac{1}{2}}, racine de la covariance de l'innovation.
 * @param gain : K_{n|n} = P_{n|n - 1} \; [F^{y,x}]^T \; S_{n}^{-1},  gain de filtrage
 * @param[in] sqrt_p_p :  [P_{n|n - 1}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit courant.
 * @param[in] f_yx : F^{y,x}, terme de la matrice de transition
 * @param[in] sqrt_q_yy : [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de mesure
 * @param mat_tt : matrice de taille  (n_t, n_t) préallouée.
 * @param mat_tt_yy : vue sur la matrice mat_tt (de (0,0) à (n_y - 1, n_y - 1))
 * @param mat_tt_yx : vue sur la matrice mat_tt (de (0, n_y) à (n_y - 1, n_t - 1))
 * @param mat_tt_xy : vue sur la matrice mat_tt (de (n_y, 0) à (n_t - 1, n_y - 1))
 * @param mat_tt_xx : vue sur la matrice mat_tt (de (n_y, n_y) à (n_t - 1, n_t - 1))
 * @param perm_y : permutation de taille y préallouée
 * @param vect_t : vecteur de taille n_t préalloué pour le calcul
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
void tkalman_nc_get_sqrt_pf_sqrt_s_and_gain(gsl_matrix * sqrt_p_f,
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
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_tt,
						 vect_t);
	gsl_triangle_matrix(mat_tt);
	
	
	//Recopie des résultats.
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

//Racine de la matrice de cov. de l'état filtré
/**@fn void tkalman_nc_get_sqrt_q_0_f_sqrt_s_0_and_gain(gsl_matrix * sqrt_q_0_f,
														gsl_matrix * sqrt_s_0,
														gsl_matrix * gain,
														const gsl_matrix * sqrt_q_0,
														const gsl_matrix * f_yt,
														const gsl_matrix * sqrt_q_yy,
														gsl_matrix * mat_tpy_tpy,
														gsl_matrix * mat_tpy_tpy_view_00,
														gsl_matrix * mat_tpy_tpy_view_01,
														gsl_matrix * mat_tpy_tpy_view_10,
														gsl_matrix * mat_tpy_tpy_view_11,
														gsl_permutation * perm_y,
														gsl_vector * vect_tpy)
 * @param sqrt_q_0_f :  [Q_{0|0}]^{\frac{1}{2}}
 * @param sqrt_s_0 : [S_{0}]^{\frac{1}{2}}, racine de la covariance de l'innovation.
 * @param gain : K^t_{0|0} = Q_{0} \; [F^{y,t}]^T \; S_{0}^{-1},  gain de filtrage
 * @param[in] sqrt_q_0 :  [Q_{0}]^{\frac{1}{2}}
 * @param[in] f_yt : F^{y,t}, terme de la matrice de transition
 * @param[in] sqrt_q_yy : [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de mesure
 * @param mat_tpy_tpy : matrice de taille (n_t + n_y, n_t + n_y)  préallouée.
 * @param mat_tpy_tpy_view_00 : vue sur la matrice mat_tpy_tpy (de (0,0) à (n_y - 1, n_y - 1))
 * @param mat_tpy_tpy_view_01 : vue sur la matrice mat_tpy_tpy (de (0, n_y) à (n_y - 1, n_t + n_y - 1))
 * @param mat_tpy_tpy_view_10 : vue sur la matrice mat_tpy_tpy (n_y, 0) à (n_t + n_y - 1, n_y - 1))
 * @param mat_tpy_tpy_view_11 : vue sur la matrice mat_tpy_tpy (n_y, n_y) à (n_t + n_y - 1,  n_t + n_y - 1))
 * @param perm_y : permutation de taille y préallouée
 * @param vect_tpy : vecteur de taille (n_t + n_y) préalloué pour le calcul
 * @brief
 Cette fonction calcule les racines de la covariance de l'état filtré courant et de la matrice de covariance de l'innovation. \n
ECe calcul s'effectue en plusieurs étapes : nous construisons la matrice M : \n
+---------------------------+
|sqrt_q_yy      0           |\n
|                           |\n
|sqrt_q0.Fyt^T  sqrt_q0    |\n
+---------------------------+\n
puis nous effectuons sa décomposition QR. A partir de la matrice R de la décomposition, nous obtenons :
+---------------------------------------+
|sqrt_s0      sqrt_s0^-1 K^T             |\n
|                           			|\n
|0           sqrt_q0_f       		|\n
+---------------------------------------+\n
 */
void tkalman_nc_get_sqrt_q_0_f_sqrt_s_0_and_gain(gsl_matrix * sqrt_q_0_f,
												 gsl_matrix * sqrt_s_0,
												 gsl_matrix * gain,
												 const gsl_matrix * sqrt_q_0,
												 const gsl_matrix * f_yt,
												 const gsl_matrix * sqrt_q_yy,
												 gsl_matrix * mat_tpy_tpy,
												 gsl_matrix * mat_tpy_tpy_view_00,
												 gsl_matrix * mat_tpy_tpy_view_01,
												 gsl_matrix * mat_tpy_tpy_view_10,
												 gsl_matrix * mat_tpy_tpy_view_11,
												 gsl_permutation * perm_y,
												 gsl_vector * vect_tpy)
{
	//Construction de la matrice
	//      +---------------------------+
	//      |sqrt_q_yy      0           |
	// M =  |                           |
	//      |sqrt_pp.Fyx^T  sqrt_p_p    |
	//      +---------------------------+
	gsl_matrix_memcpy(mat_tpy_tpy_view_00, sqrt_q_yy);
	gsl_matrix_memcpy(mat_tpy_tpy_view_11, sqrt_q_0);
	gsl_matrix_set_zero(mat_tpy_tpy_view_01);
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   sqrt_q_0,
				   f_yt,
				   0.0,
				   mat_tpy_tpy_view_10);
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_tpy_tpy,
						 vect_tpy);
	gsl_triangle_matrix(mat_tpy_tpy);
	
	
	//Recopie des résultats.
	gsl_matrix_memcpy(sqrt_q_0_f,
					  mat_tpy_tpy_view_11);
	gsl_matrix_memcpy(sqrt_s_0,
					  mat_tpy_tpy_view_00);

	//Calcul du gain
		//Inversion de sqrt_s
        gsl_permutation_init(perm_y);
        gsl_linalg_LU_invert(sqrt_s_0, perm_y, mat_tpy_tpy_view_00);

		//Produit matriciel
		gsl_blas_dgemm(CblasTrans,
					   CblasTrans,
					   1.0,
					   mat_tpy_tpy_view_01,
					   mat_tpy_tpy_view_00,
					   0.0,
					   gain);
}

//Filtrage complet

/**@fn void tkalman_nc_do_filtering(gsl_vector * x_f,
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
 * @param x_f : \hat{x}_{n|n}, espérance de l'état filtré
 * @param sqrt_p_f :  [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré
 * @param innovation : \tilde{y}_{n}, innovation
 * @param sqrt_s : [S_{n}]^{\frac{1}{2}}, racine de la covariance de l'innovation.
 * @param[in] x_p : \hat{x}_{n | n - 1}, espérance de l'état prédit courant
 * @param[in] sqrt_p_p :  [P_{n|n - 1}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit courant.
 * @param[in] y : observation cournate
 * @param[in] _y : observation précédente
 * @param[in] f_yx : terme de la matrice de transition (Fyx)
 * @param[in] f_yy : terme de la matrice de transition (Fyy)
 * @param[in] sqrt_q_yy : Racine de la matrice de covariance du bruit de mesure (Qyy)
 * @param mat_tt : matrice de taille  (n_t, n_t) préallouée.
 * @param mat_tt_yy : vue sur la matrice mat_tt (de (0,0) à (n_y - 1, n_y - 1))
 * @param mat_tt_yx : vue sur la matrice mat_tt (de (0, n_y) à (n_y - 1, n_t - 1))
 * @param mat_tt_xy : vue sur la matrice mat_tt (de (n_y, 0) à (n_t - 1, n_y - 1))
 * @param mat_tt_xx : vue sur la matrice mat_tt (de (n_y, n_y) à (n_t - 1, n_t - 1))
 * @param mat_xy : matrice de taille (n_x.n_y) préallouée
 * @param perm_y : permutation de taille n_y préallouée
 * @param vect_t : vecteur de taille n_t préalloué pour le calcul
 * @brief
 * Cette fonction effectue la partie filtrage du filtre de Kalman triplet : \n
 Nous calculons dans un premier temps l'espérance de l'innovation avec la fonction @fn tkalman_nc_get_innovation. Puis dans un second temps, nous calculons les racines des matrices de covariance de l'innovarion et de l'état filtré courant avec @fn tkalman_nc_get_sqrt_pf_sqrt_s_and_gain .Finalement dans un quatrième temps, nous calculons l'espérance de l'état filtré avec @fn tkalman_nc_get_x_f.
 */
void tkalman_nc_do_filtering(gsl_vector * x_f,
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
    //Calcul de l'espérance de l'innovation
    tkalman_nc_get_innovation ( innovation,
							 x_p,
							 _y,
							 y,
							 f_yx,
							 f_yy);
	//Racines des matrices de covariance et du gain
	tkalman_nc_get_sqrt_pf_sqrt_s_and_gain ( sqrt_p_f,
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
	//Espérance de l'état filtré
	tkalman_nc_get_x_f ( x_f,
						 x_p,
						 innovation,
						 mat_xy);
}

/**@fn void tkalman_nc_do_filtering_0(gsl_vector * t_f_0,
			                       gsl_matrix * sqrt_q_f_0,
								   gsl_vector * innovation,
			                       gsl_matrix * sqrt_s_0,
			                       const gsl_vector * t_0,
			       				   const gsl_vector * x_0,
			                       const gsl_matrix * sqrt_q_0,
			                       const gsl_vector * y_0,
			                       const gsl_vector * y_m1,
			                       const gsl_matrix * f_y,
			                       const gsl_matrix * f_yx,
			                       const gsl_matrix * f_yy,
			                       const gsl_matrix * sqrt_q_yy,
			                       gsl_matrix * mat_tpy_tpy,
			                       gsl_matrix * mat_tpy_tpy_view_00,
								   gsl_matrix * mat_tpy_tpy_view_01,
			                       gsl_matrix * mat_tpy_tpy_view_10,
			                       gsl_matrix * mat_tpy_tpy_view_11,
			                       gsl_matrix * mat_ty,
			                       gsl_permutation * perm_y,
			                       gsl_vector * vect_tpy)
 * @param t_f_0 : espérance de t0 filtré
 * @param sqrt_q_f_0, : racine de la matrice de covariance de t0 filtré
 * @param innovation : espérance de l'innovation
 * @param sqrt_s_0 : racine de la covariance de l'innovation.
 * @param[in] t_0, : \hat{t}_{n | n - 1}, espérance de l'état prédit courant
 * @param[in] x_0 : \hat{x}_{n | n - 1}, espérance de l'état prédit
 * @param[in] sqrt_q_0 : racine de la matrice de covariance de t0
 * @param[in] y_0 : observation cournate
 * @param[in] y_m1 : espérance de l'observation précédente
 * @param[in] f_y : ligne des y de F
 * @param[in] f_yx : terme de la matrice de transition (Fyx)
 * @param[in] f_yy : terme de la matrice de transition (Fyy)
 * @param[in] sqrt_q_yy : Racine de la matrice de covariance du bruit de mesure (Qyy)
 * @param mat_tpy_tpy : matrice de taille (n_t + n_y.n_t + n_y) préallouée.
 * @param mat_tpy_tpy_view_00 : vue sur la matrice mat_tpy_tpy (de (0,0) à (n_y - 1, n_y - 1))
 * @param mat_tpy_tpy_view_01 : vue sur la matrice mat_tpy_tpy (de (0,n_y) à (n_y - 1, n_y + n_t - 1))
 * @param mat_tpy_tpy_view_10 : vue sur la matrice mat_tpy_tpy (de (n_y, 0) à (n_y + n_t - 1, n_y - 1))
 * @param mat_tpy_tpy_view_11 : vue sur la matrice mat_tpy_tpy (de (n_y, n_y) à (n_y + n_t - 1, n_y + n_t - 1))
 * @param mat_ty : matrice de taille (n_t.n_y) préallouée
 * @param perm_y : permutation de taille y préallouée
 * @param vect_tpy : vecteur de taille (n_t + n_y) préalloué pour le calcul
 * @brief
 * Cette fonction effectue la partie filtrage du filtre de Kalman triplet : \n
 Nous calculons dans un premier temps l'espérance de l'innovation avec la fonction @fn tkalman_nc_get_innovation_0. Puis dans un second temps, nous calculons les racines des matrices de covariance de l'innovarion et de l'état filtré courant avec @fn tkalman_nc_get_sqrt_q_0_f_sqrt_s_0_and_gain .Finalement dans un quatrième temps, nous calculons l'espérance de l'état filtré avec @fn tkalman_nc_get_t_0_f.
 */
void tkalman_nc_do_filtering_0(gsl_vector * t_0_f,
			                   gsl_matrix * sqrt_q_0_f,
							   gsl_vector * innovation,
						       gsl_matrix * sqrt_s_0,
							   const gsl_vector * t_0,
						       const gsl_matrix * sqrt_q_0,
						       const gsl_vector * y_0,
						       const gsl_matrix * f_yt,
						       const gsl_matrix * sqrt_q_yy,
						       gsl_matrix * mat_tpy_tpy,
						       gsl_matrix * mat_tpy_tpy_view_00,
						       gsl_matrix * mat_tpy_tpy_view_01,
						       gsl_matrix * mat_tpy_tpy_view_10,
						       gsl_matrix * mat_tpy_tpy_view_11,
						       gsl_matrix * mat_ty,
						       gsl_permutation * perm_y,
						       gsl_vector * vect_tpy)
{
    tkalman_nc_get_innovation_0 ( innovation,
								  t_0,
								  y_0,
								  f_yt);
								  
	
	tkalman_nc_get_sqrt_q_0_f_sqrt_s_0_and_gain ( sqrt_q_0_f,
												  sqrt_s_0,
												  mat_ty,
												  sqrt_q_0,
												  f_yt,
												  sqrt_q_yy,
											      mat_tpy_tpy,
												  mat_tpy_tpy_view_00,
												  mat_tpy_tpy_view_01,
												  mat_tpy_tpy_view_10,
												  mat_tpy_tpy_view_11,
												  perm_y,
												  vect_tpy);
												  
	tkalman_nc_get_t_0_f ( t_0_f,
						   t_0,
						   innovation,
						   mat_ty);
	
}

//Objet
/**@fn tkalman_nc_filtering :: tkalman_nc_filtering(const gsl_matrix * f_yt,
												    const gsl_matrix * sqrt_q_yy) throw(exception &);
* @param[in] f_yt : F^{y,t}
* @param[in] sqrt_q_yy : [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
* @brief
* Constructeur de la classe @class tkalman_nc_filtering
* @throw 
* Exception (std :: bad_alloc si problème de mémoire ou invalid_argument en cas d'arguments invalides)
*/
tkalman_nc_filtering :: tkalman_nc_filtering(const gsl_matrix * f_yt,
											 const gsl_matrix * sqrt_q_yy) throw(exception &)
{
	initialize();
	try
	{
		setup( f_yt,
			   sqrt_q_yy);
	}
	catch(exception & e)
	{
		throw(e);
	}
}


/**@fn void tkalman_nc_filtering :: setup(const gsl_matrix * f_yt,
										  const gsl_matrix * sqrt_q_yy) throw(exception &);
 * @param[in] f_yt : F^{y,t}
 * @param[in] sqrt_q_yy : [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
 * @brief
 * Cette fonction libère les attributs et les réalloue.
 */
void tkalman_nc_filtering :: setup(const gsl_matrix * f_yt,
								   const gsl_matrix * sqrt_q_yy) throw(exception &)
{

	//Vérif des arguments
	if (!f_yt || !sqrt_q_yy)
	{
		throw(invalid_argument("Fyt or sqrt_q_yy are NULL!\n"));
	}	
	unsigned int size_t, 
				 size_x,
				 size_y;
	size_t = f_yt->size2;
	size_y = f_yt->size1;
	size_x = size_t - size_y;
	//Vérif de la dim. de y;
	if (!size_y || !size_x)
		throw(invalid_argument("Size of x or size of y are 0!\n"));
	//Modif des matrices
	if (size_x != _size_x || size_y != _size_y)
	{

		free();
		initialize();
		_size_x = size_x;
		_size_y = size_y;
		_size_t = size_t;
		_f_yt = f_yt;
		_sqrt_q_yy = sqrt_q_yy;
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
		_f_yt = f_yt;
		_sqrt_q_yy = sqrt_q_yy;
	}

	//Création des vues
	create_views();
}

/**@fn bool tkalman_nc_filtering :: operator !() const
 * @return 
 * - 0 si l'objet est normal
 * - 1 sinon.
 * @brief
 * Check de l'objet.
 **/
bool tkalman_nc_filtering :: operator !() const
{
	return ( ! (mat_xy && mat_ty && mat_tpy_tpy && mat_tt && vect_tpy && vect_t && perm_y && _f_yt && _sqrt_q_yy ) );
}

/**@fn void tkalman_nc_filtering :: compute_filtering_0(gsl_vector * t_0_f,
														gsl_matrix * sqrt_q_0_f,
														gsl_vector * innovation,
														gsl_matrix * sqrt_s_0,
														const gsl_vector * t_0,
														const gsl_matrix * sqrt_q_0,
														const gsl_vector * y_0);
  * @param t_f_0 : espérance de t0 filtré
  * @param sqrt_q_f_0, : racine de la matrice de covariance de t0 filtré
  * @param innovation : espérance de l'innovation
  * @param sqrt_s_0 : racine de la covariance de l'innovation.
  * @param[in] t_0, : \hat{t}_{0}, espérance de l'état prédit courant
  * @param[in] sqrt_q_0 :  [Q_{0}]^{\frac{1}{2}}
  * @param[in] y_0 : y_0, espérance de l'observation courante
  * @brief
  * Cette fonction effectue la partie filtrage du filtre de Kalman triplet : \n
  Nous calculons dans un premier temps l'espérance de l'innovation avec la fonction @fn tkalman_nc_get_innovation_0. Puis dans un second temps, nous calculons les  racines des matrices de covariance de l'innovarion et de l'état filtré courant avec @fn tkalman_nc_get_sqrt_q_0_f_sqrt_s_0_and_gain .Finalement dans un quatrième temps, nous calculons l'espérance de l'état filtré avec @fn tkalman_nc_get_t_0_f.
 **/
void tkalman_nc_filtering :: compute_filtering_0(gsl_vector * t_0_f,
												 gsl_matrix * sqrt_q_0_f,
												 gsl_vector * innovation,
												 gsl_matrix * sqrt_s_0,
												 const gsl_vector * t_0,
												 const gsl_matrix * sqrt_q_0,
												 const gsl_vector * y_0)
{
	tkalman_nc_do_filtering_0(t_0_f,
							  sqrt_q_0_f,
							  innovation,
							  sqrt_s_0,
							  t_0,
							  sqrt_q_0,
							  y_0,
							  _f_yt,
							  _sqrt_q_yy,
							  mat_tpy_tpy,
							  &mat_tpy_tpy_view_00,
							  &mat_tpy_tpy_view_01,
						      &mat_tpy_tpy_view_10,
							  &mat_tpy_tpy_view_11,
							  mat_ty,
							  perm_y,
							  vect_tpy);
}

/**@fn void tkalman_nc_filtering :: compute_filtering(gsl_vector * x_f,
													  gsl_matrix * sqrt_p_f,
													  gsl_vector * innovation,
													  gsl_matrix * sqrt_s,
													  const gsl_vector * x_p,
													  const gsl_matrix * sqrt_p_p,
													  const gsl_vector * y,
													  const gsl_vector * _y);
   * @param x_f : \hat{x}_{n|n}, espérance de l'état filtré
   * @param sqrt_p_f :  [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré
   * @param innovation : \tilde{y}_{n}, innovation
   * @param sqrt_s : [S_{n}]^{\frac{1}{2}}, racine de la covariance de l'innovation.
   * @param[in] x_p : \hat{x}_{n | n - 1}, espérance de l'état prédit courant
   * @param[in] sqrt_p_p :  [P_{n|n - 1}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit courant.
   * @param[in] y : observation cournate
   * @param[in] _y : observation précédente
   * @brief
   * Cette fonction effectue la partie filtrage du filtre de Kalman triplet : \n
   Nous calculons dans un premier temps l'espérance de l'innovation avec la fonction @fn tkalman_nc_get_innovation. Puis dans un second temps, nous calculons les racines des matrices de covariance de l'innovarion et de l'état filtré courant avec @fn tkalman_nc_get_sqrt_pf_sqrt_s_and_gain .Finalement dans un quatrième temps, nous calculons l'espérance de l'état filtré avec @fn tkalman_nc_get_x_f.
 **/
void tkalman_nc_filtering :: compute_filtering(gsl_vector * x_f,
											   gsl_matrix * sqrt_p_f,
											   gsl_vector * innovation,
											   gsl_matrix * sqrt_s,
											   const gsl_vector * x_p,
											   const gsl_matrix * sqrt_p_p,
											   const gsl_vector * y,
											   const gsl_vector * _y)
{
	tkalman_nc_do_filtering( x_f,
							 sqrt_p_f,
							 innovation,
							 sqrt_s,
							 x_p,
						     sqrt_p_p,
						     y,
							 _y,
							 &_f_yx,
							 &_f_yy,
							 _sqrt_q_yy,
							 mat_tt,
							 &mat_tt_yy,
							 &mat_tt_yx,
							 &mat_tt_xy,
							 &mat_tt_xx,
						     mat_xy,
							 perm_y,
							 vect_t);	
}

/**@fn tkalman_nc_filtering :: ~tkalman_nc_filtering();
 * @brief
 * Destructeur de la classe @class tkalman_nc_filtering
 */
tkalman_nc_filtering :: ~tkalman_nc_filtering()
{
	free();
	initialize();
}

/**@fn void tkalman_nc_filtering :: initialize();
 * @brief
 * Cette fonction met tous les attributs de l'objet à 0.
 */
void tkalman_nc_filtering :: initialize()
{
	mat_xy = 0;
	mat_ty = 0;
	mat_tt = 0;
	mat_tpy_tpy = 0;
	vect_tpy = 0;
	vect_t = 0;
	perm_y = 0;
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
	_f_yt = 0;
	_sqrt_q_yy = 0;
}

/**@fn void tkalman_nc_filtering :: alloc() throw(exception &);
 * @brief
 * Cette fonction alloue les différents élements de la classe.
 * @throw
 * bad_alloc en cas de problème de mémoire
 */
void tkalman_nc_filtering :: alloc() throw(exception &)
{
	if (!mat_xy)
	{
		try
		{
			mat_xy = gsl_matrix_alloc(_size_x, _size_y);
		}
		catch(exception & e)
		{
			throw(e);
		}
		
	}
	if (!mat_ty)
	{
		try
		{
			mat_ty = gsl_matrix_alloc(_size_t, _size_y);
		}
		catch(exception & e)
		{
			throw(e);
		}
		
	}
	if (!mat_tt)
	{
		try
		{
			mat_tt = gsl_matrix_alloc(_size_t, _size_t);
		}
		catch(exception & e)
		{
			throw(e);
		}
		
	}
	if (!mat_tpy_tpy)
	{
		try
		{
			mat_tpy_tpy = gsl_matrix_alloc(_size_t + _size_y, _size_t + _size_y);
		}
		catch(exception & e)
		{
			throw(e);
		}
		
	}
	if (!vect_t)
	{
		try
		{
			vect_t = gsl_vector_alloc(_size_t);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!vect_tpy)
	{
		try
		{
			vect_tpy = gsl_vector_alloc(_size_t + _size_y);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!perm_y)
	{
		try
		{
			perm_y = gsl_permutation_alloc( _size_y);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	
}
/**@fn void tkalman_nc_filtering :: free();
 * @brief
 * Cette fonction désalloue tous les attributs alloués.
 **/
void tkalman_nc_filtering :: free()
{
	if (mat_xy)
	{
		gsl_matrix_free(mat_xy);
	}
	if (mat_ty)
	{
		gsl_matrix_free(mat_ty);
		
	}
	if (mat_tt)
	{
		gsl_matrix_free(mat_tt);
	}
	if (mat_tpy_tpy)
	{
		gsl_matrix_free(mat_tpy_tpy);
	}
	if (vect_t)
	{
		gsl_vector_free(vect_t);
	}
	if (vect_tpy)
	{
		gsl_vector_free(vect_tpy);
	}
	if (perm_y)
	{
		gsl_permutation_free(perm_y);
	}
}

/**@fn void tkalman_nc_filtering :: create_views();
 * @brief 
 * Cette fonction génère les différentes vues sur les matrices.
 */
void tkalman_nc_filtering :: create_views()
{
	gsl_matrix_view view;
	//mat_tt
		view = gsl_matrix_submatrix(mat_tt, 0, 0, _size_y, _size_y);
		mat_tt_yy = view.matrix;
		view = gsl_matrix_submatrix(mat_tt, 0, _size_y, _size_y, _size_x);
		mat_tt_yx = view.matrix;
		view = gsl_matrix_submatrix(mat_tt, _size_y, 0, _size_x, _size_y);
		mat_tt_xy = view.matrix;
		view = gsl_matrix_submatrix(mat_tt, _size_y, _size_y, _size_x, _size_x);
		mat_tt_xx = view.matrix;
		
	//mat_tpy_tpy
		view = gsl_matrix_submatrix(mat_tpy_tpy, 0, 0, _size_y, _size_y);
		mat_tpy_tpy_view_00 = view.matrix;
		view = gsl_matrix_submatrix(mat_tpy_tpy, 0, _size_y, _size_y, _size_t);
		mat_tpy_tpy_view_01 = view.matrix;
		view = gsl_matrix_submatrix(mat_tpy_tpy, _size_y, 0, _size_t, _size_y);
		mat_tpy_tpy_view_10 = view.matrix;
		view = gsl_matrix_submatrix(mat_tpy_tpy, _size_y, _size_y, _size_t, _size_t);
		mat_tpy_tpy_view_11 = view.matrix;
	
	//Fyt
		{
			gsl_matrix_const_view view2 = gsl_matrix_const_submatrix(_f_yt, 0, 0, _size_y, _size_x);
			_f_yx = view2.matrix;
		}
		{
			gsl_matrix_const_view view2 = gsl_matrix_const_submatrix(_f_yt, 0, _size_x, _size_y, _size_y);
			_f_yy = view2.matrix;
		}
	
}

