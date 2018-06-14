/** @file tkalman_filtering.cpp
 * @author Valérian Némesin
 * @date 10/11/2011
 * @brief
 Ce fichier contient le code source des fonctions nécessaire au filtrage dans l'algorithme du filtre de Kalman triplet.
**/
#include "tkalman_filtering.hpp"


/**@fn void tkalman_get_innovation(gsl_vector * innovation,
											const gsl_vector * x_p,
											const gsl_vector * _y,
											const gsl_vector * y,
											const gsl_matrix * f_yx,
											const gsl_matrix * f_yy)
 * @param innovation : innovation
 * @param[in] x_p : \hat{x}_{n | n - 1}, espérance de l'état prédit
 * @param[in] _y : y_{n - 1}, espérance de l'observation précédente
 * @param[in] y : y_n, espérance de l'observation courante
 * @param[in] f_yx : F^{y,x}, terme de la matrice d'évolution
 * @param[in] f_yy : F^{y,y}, terme de la matrice d'évolution
 * @brief
 Cette fonction calcule l'innovation selon la formule :
 \tilde{y}_n = y_{n} - F^{y,x} \: \hat{x}_{n | n - 1} - F^{y,y} \: y_{n - 1}
 */
void tkalman_get_innovation(gsl_vector * innovation,
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


/**@fn void tkalman_get_x_f(gsl_vector * x_f,
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
void tkalman_get_x_f(gsl_vector * x_f,
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

/**@fn void tkalman_get_sqrt_pf_sqrt_s_and_gain(gsl_matrix * sqrt_p_f,
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
void tkalman_get_sqrt_pf_sqrt_s_and_gain(gsl_matrix * sqrt_p_f,
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

/**@fn void tkalman_do_filtering(gsl_vector * x_f,
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
void tkalman_do_filtering(gsl_vector * x_f,
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
    tkalman_get_innovation ( innovation,
							 x_p,
							 _y,
							 y,
							 f_yx,
							 f_yy);
	//Racines des matrices de covariance et du gain
	tkalman_get_sqrt_pf_sqrt_s_and_gain ( sqrt_p_f,
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
	tkalman_get_x_f ( x_f,
			          x_p,
				      innovation,
					  mat_xy);
}

/**@fn void tkalman_do_filtering_0(gsl_vector * t_f_0,
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
 Nous calculons dans un premier temps l'espérance de l'innovation. Puis dans un second temps, nous calculons les racines des matrices de covariance de l'innovarion et de l'état filtré courant. Ensuite dans un troisième temps, nous calculons le gain de filtrage.Finalement dans un quatrième temps, nous calculons l'espérance de l'état filtré.
 */
void tkalman_do_filtering_0(gsl_vector * t_f_0,
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
{

    //Calcul de l'espérance de l'innovation
    tkalman_get_innovation ( innovation,
							 x_0,
						     y_m1,
							 y_0,
							 f_yx,
							 f_yy
							);
	//Racines des matrices de covariance et du gain
	tkalman_get_sqrt_pf_sqrt_s_and_gain ( sqrt_q_f_0,
										  sqrt_s_0,
										  mat_ty,
										  sqrt_q_0,
										  f_y,
										  sqrt_q_yy,
									      mat_tpy_tpy,
										  mat_tpy_tpy_view_00,
										  mat_tpy_tpy_view_01,
										  mat_tpy_tpy_view_10,
										  mat_tpy_tpy_view_11,
										  perm_y,
										  vect_tpy);
	//Espérance de l'état filtré
	tkalman_get_x_f ( t_f_0,
			          t_0,
				      innovation,
					  mat_ty);
}

//Objet de prédiction
/**@fn tkalman_filtering :: tkalman_filtering(const gsl_matrix * f_y,
											  const gsl_matrix * sqrt_q_yy);
* @param f_y: [Fyx, Fyy]
* @param sqrt_q_yy : racine de Qyy
* @brief
* Ce constructeur alloue les variables temp. de la prédiction du filtre de Kalman couple.
*/
tkalman_filtering :: tkalman_filtering(const gsl_matrix * f_y,
									   const gsl_matrix * sqrt_q_yy)
{
	tkalman_filtering :: initialize();
	tkalman_filtering :: setup (f_y,
								sqrt_q_yy);
}


/**@fn int tkalman_filtering :: setup(const gsl_matrix * f_y,
									   const gsl_matrix * sqrt_q_yy);
 * @param f_y: [Fyx, Fyy]
 * @param sqrt_q_yy : racine de Qyy
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction permet de modifier les paramètres de la prédiction (size_x et size_y)
**/
int tkalman_filtering :: setup (const gsl_matrix * f_y,
								const gsl_matrix * sqrt_q_yy)
	
{	
	unsigned int size_x = f_y->size2 - sqrt_q_yy->size1,
				 size_y = sqrt_q_yy->size1;
	
	if (size_x != _size_x || size_y != _size_y)
	{
		tkalman_filtering :: free();
		tkalman_filtering :: initialize();
		if (tkalman_filtering :: set_params(f_y, sqrt_q_yy))
		{
			tkalman_filtering :: free();
			tkalman_filtering :: initialize();
			return 1;
		}
		_size_x = size_x;
		_size_y = size_y;
		_size_t = size_x + size_y;
		if ( tkalman_filtering :: alloc() )
		{
			tkalman_filtering :: free();
			tkalman_filtering :: initialize();
			return 1;
		}
		else
			tkalman_filtering :: create_views();
	}	
	return 0;
}

/**@fn tkalman_filtering :: ~ tkalman_filtering()
 * @brief
 * Destructeur
 */
tkalman_filtering :: ~tkalman_filtering()
{
	tkalman_filtering :: free();
	tkalman_filtering :: initialize();
}

/**@fn bool tkalman_filtering :: operator!() const;
 * @return
 * - 0 si l'objet est correctement alloué
 * - 1 sinon
 * @brief
 * Check de l'objet.
 */
bool tkalman_filtering :: operator!() const
{
	return (! (mat_tpy_tpy && mat_ty && perm_y && vect_tpy && _f_y && _sqrt_q_yy) );
}

/**@fn void tkalman_filtering :: free();
 * @brief
 * Cette fonction désalloue la mémoire utilisée par les variables tmp.
 */
void tkalman_filtering :: free()
{
	if (mat_tpy_tpy)
		gsl_matrix_free(mat_tpy_tpy);
	if (mat_ty)
		gsl_matrix_free(mat_ty);
	if (perm_y)
		gsl_permutation_free(perm_y);
	if (vect_tpy)
		gsl_vector_free(vect_tpy);
}

/**@fn int tkalman_filtering :: alloc();
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction alloue la mémoire utilisée par les variables tmp.
 */
int tkalman_filtering :: alloc()
{
	if (!mat_tpy_tpy)
		mat_tpy_tpy = gsl_matrix_alloc(_size_t + _size_y, _size_t + _size_y);
	if (!mat_ty)
		mat_ty= gsl_matrix_alloc(_size_t, _size_y);
	if (!perm_y)
		perm_y = gsl_permutation_alloc(_size_y);
	if (!vect_tpy)
		vect_tpy = gsl_vector_alloc(_size_y + _size_t);
		
	return (! (mat_tpy_tpy && mat_ty && perm_y && vect_tpy) );
}

/**@fn void tkalman_filtering :: create_views();
 * @brief
 * Cette méthode crée les vues.
 */
void tkalman_filtering :: create_views()
{
	gsl_matrix_view view;
	view = gsl_matrix_submatrix(mat_tpy_tpy,
								0,
								0,
								_size_t,
								_size_t);
	mat_tt = view.matrix;
	
	view = gsl_matrix_submatrix(mat_tpy_tpy,
								0,
								0,
								_size_y,
								_size_y);
	mat_tt_view_00 = view.matrix;
	mat_tpy_tpy_view_00 = view.matrix;
	
	view = gsl_matrix_submatrix(mat_tpy_tpy,
								_size_y,
								0,
								_size_t,
								_size_y);
	mat_tpy_tpy_view_10 = view.matrix;
	
	
	view = gsl_matrix_submatrix(mat_tpy_tpy,
								0,
								_size_y,
								_size_y,
								_size_t);
	mat_tpy_tpy_view_01 = view.matrix;
	
	
	view = gsl_matrix_submatrix(mat_tpy_tpy,
								_size_y,
								_size_y,
								_size_t,
								_size_t);
	mat_tpy_tpy_view_11 = view.matrix;
	
	
	
	view = gsl_matrix_submatrix(mat_tpy_tpy,
								_size_y,
								0,
								_size_x,
								_size_y);
	mat_tt_view_10 = view.matrix;
	
	
	view = gsl_matrix_submatrix(mat_tpy_tpy,
								0,
								_size_y,
								_size_y,
								_size_x);
	mat_tt_view_01 = view.matrix;
	
	
	view = gsl_matrix_submatrix(mat_tpy_tpy,
								_size_y,
								_size_y,
								_size_x,
								_size_x);
	mat_tt_view_11 = view.matrix;
	
	
	view = gsl_matrix_submatrix(mat_ty,
								0,
								0,
								_size_x,
								_size_y);
	mat_xy = view.matrix;
	
	
	gsl_vector_view view2 = gsl_vector_subvector(vect_tpy, 0, _size_t);
	
	vect_t = view2.vector;
	
	
}

		
/**@fn void tkalman_filtering :: initialize();
 * @brief
 * Cette fonction met les variables tmp à NULL.
 */
void tkalman_filtering :: initialize()
{

	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
			
	mat_tpy_tpy = NULL;
	mat_ty = NULL;
	perm_y = NULL;
	vect_tpy = NULL;


	_f_y = NULL;
	_sqrt_q_yy = NULL;

}

/**@fn virtual int tkalman_filtering :: set_params(const gsl_matrix * f_y,
									   const gsl_matrix * sqrt_q_yy);
 * @param f_y: [Fyx, Fyy]
 * @param sqrt_q_yy : racine de Qyy
 * @brief
 * Cette méthode change les paramètres de la prédiction.
 */
int tkalman_filtering :: set_params(const gsl_matrix * f_y,
									const gsl_matrix * sqrt_q_yy)
{
	if (!f_y || !sqrt_q_yy)
		return 1;
		
	unsigned int size_x = f_y->size2 - sqrt_q_yy->size1,
				 size_y = sqrt_q_yy->size1;

	_f_y = f_y;
	_sqrt_q_yy = sqrt_q_yy;
	
	//Vues
	{
		gsl_matrix_const_view view = gsl_matrix_const_submatrix (_f_y, 0, 0, size_y, size_x);
		f_yx = view.matrix;
	}
	{
		gsl_matrix_const_view view = gsl_matrix_const_submatrix (_f_y, 0, size_x, size_y, size_y);
		f_yy = view.matrix;
	}
	return 0;
}
