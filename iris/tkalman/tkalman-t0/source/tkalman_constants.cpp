/**@file tkalman_constants.hpp
 * @author Valérian Némesin
 * @brief
 Ce fichier contient le code soure des fonctions qui calculent les différentes constantes du filtre de Kalman triplet.
 */
#include "tkalman_constants.hpp"

/**@fn void tkalman_get_f2_x_(gsl_matrix * f2_x_,
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
void tkalman_get_f2_x_(gsl_matrix * f2_x_,
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

/**@fn void tkalman_get_sqrt_q2_xx_sqrt_q_yy_and_q2_xy_from_sqrt_q(gsl_matrix * sqrt_q2_xx,
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
void tkalman_get_sqrt_q2_xx_sqrt_q_yy_and_q2_xy_from_sqrt_q(gsl_matrix * sqrt_q2_xx,
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
void tkalman_get_constants(gsl_matrix * f2_x_,
						   gsl_matrix * sqrt_q2_xx,
						   gsl_matrix * q2_xy,
						   gsl_matrix * sqrt_q_yy,
						   const gsl_matrix * f_x_,
						   const gsl_matrix * f_y_,
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
	tkalman_get_sqrt_q2_xx_sqrt_q_yy_and_q2_xy_from_sqrt_q(sqrt_q2_xx,
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
	tkalman_get_f2_x_(f2_x_,
					  f_x_,
					  f_y_,
					  q2_xy);


}

//Objet de prédiction
/**@fn tkalman_constants :: tkalman_constants(const gsl_matrix * f_x_,
											  const gsl_matrix * f_y_,
											  const gsl_matrix * sqrt_q);
* @param f2_x : [F2xx, F2xy]
* @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
* @param q2_xy : Qxy.Qyy
* @brief
* Ce constructeur alloue les variables temp. de la prédiction du filtre de Kalman couple.
*/
tkalman_constants :: tkalman_constants (const gsl_matrix * f_x_,
										const gsl_matrix * f_y_,
										const gsl_matrix * sqrt_q)
{
	tkalman_constants :: initialize();
	tkalman_constants :: setup (f_x_, f_y_, sqrt_q);
}


/**@fn int tkalman_constants :: setup(const gsl_matrix * f_x_,
									  const gsl_matrix * f_y_,
									  const gsl_matrix * sqrt_q);
{
 * @param f2_x : [F2xx, F2xy]
 * @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
 * @param q2_xy : Qxy.Qyy
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction permet de modifier les paramètres de la prédiction (size_x et size_y)
**/
int tkalman_constants :: setup (const gsl_matrix * f_x_,
								const gsl_matrix * f_y_,
								const gsl_matrix * sqrt_q)
{	
	unsigned int size_x = f_x_->size1,
				 size_y = f_y_->size1;

	if (size_x != _size_x || size_y != _size_y)
	{
		tkalman_constants :: free();
		tkalman_constants :: initialize();

		_size_x = size_x;
		_size_y = size_y;
		_size_t = size_x + size_y;
		if (tkalman_constants :: set_params(f_x_,
											f_y_,
											sqrt_q)
											)
		{
			tkalman_constants :: free();
			tkalman_constants :: initialize();
			return 1;
		}



		if ( tkalman_constants :: alloc() )
		{

			tkalman_constants :: free();
			tkalman_constants :: initialize();
			return 1;
		}
		else
			tkalman_constants :: create_views();
	}	

	return 0;
}

/**@fn tkalman_constants :: ~ tkalman_constants()
 * @brief
 * Destructeur
 */
tkalman_constants :: ~tkalman_constants()
{
	tkalman_constants :: free();
	tkalman_constants :: initialize();
}

/**@fn bool tkalman_constants :: operator!() const;
 * @return
 * - 0 si l'objet est correctement alloué
 * - 1 sinon
 * @brief
 * Check de l'objet.
 */
bool tkalman_constants :: operator!() const
{
	return (! (mat_tt && vect_t && perm_y && _f_x_ && _f_y_) );
}

/**@fn void tkalman_constants :: free();
 * @brief
 * Cette fonction désalloue la mémoire utilisée par les variables tmp.
 */
void tkalman_constants :: free()
{
	if (perm_y)
		gsl_permutation_free(perm_y);
	if (vect_t)
		gsl_vector_free(vect_t);
	if (mat_tt)
		gsl_matrix_free(mat_tt);
}

/**@fn int tkalman_constants :: alloc();
 * @return
 * 0 si bon déroulement de l'op.
 * @brief
 * Cette fonction alloue la mémoire utilisée par les variables tmp.
 */
int tkalman_constants :: alloc()
{
	if (!perm_y)
		perm_y = gsl_permutation_alloc(_size_y);
	if (!vect_t)
		vect_t = gsl_vector_alloc(_size_t);
	if (!mat_tt)
		mat_tt = gsl_matrix_alloc(_size_t, _size_t);
	return (! (mat_tt && vect_t && perm_y) );
}

/**@fn void tkalman_constants :: create_views();
 * @brief
 * Cette méthode crée les vues.
 */
void tkalman_constants :: create_views()
{
	gsl_matrix_view view;
	view = gsl_matrix_submatrix(mat_tt,
								0,
								0,
								_size_y,
								_size_y);
	mat_tt_yy = view.matrix;
	
	view = gsl_matrix_submatrix(mat_tt,
								0,
								_size_y,
								_size_y,
								_size_x);
	mat_tt_yx = view.matrix;
	
	view = gsl_matrix_submatrix(mat_tt,
								_size_y,
								0,
								_size_x,
								_size_y);
	mat_tt_xy = view.matrix;
	
	view = gsl_matrix_submatrix(mat_tt,
								_size_y,
								_size_y,
								_size_x,
								_size_x);
	mat_tt_xx = view.matrix;
}

		
/**@fn void tkalman_constants :: initialize();
 * @brief
 * Cette fonction met les variables tmp à NULL.
 */
void tkalman_constants :: initialize()
{
	perm_y = NULL;
	vect_t = NULL;
	mat_tt = NULL;
	_f_x_ = NULL;
	_f_y_ = NULL;
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
}

/**@fn virtual int tkalman_constants :: set_params(const gsl_matrix * f2_x,
													const gsl_matrix * sqrt_q2_xx)
 * @param f2_x : [F2xx, F2xy]
 * @param sqrt_q2_xx : racine de Qxx - Qxy Qyy Qyx
 * @param q2_xy : Qxy.Qyy
 * @brief
 * Cette méthode change les paramètres de la prédiction.
 */
int tkalman_constants :: set_params(const gsl_matrix * f_x_,
									const gsl_matrix * f_y_,
									const gsl_matrix * sqrt_q)
{

	if (!f_x_ || !f_y_ || !sqrt_q)
		return 1;
	unsigned int size_x = f_x_->size1,
				 size_y = f_y_->size1;
	_f_x_ = f_x_;
	_f_y_ = f_y_;
	//Vues
	{
		gsl_matrix_const_view view = gsl_matrix_const_submatrix (sqrt_q, 0, 0, size_x, size_x);
		sqrt_q_view_xx = view.matrix;
	}
	{
		gsl_matrix_const_view view = gsl_matrix_const_submatrix (sqrt_q, 0, size_x, size_x, size_y);
		sqrt_q_view_xy = view.matrix;
	}
	{
		gsl_matrix_const_view view = gsl_matrix_const_submatrix (sqrt_q, size_x, size_x, size_y, size_y);
		sqrt_q_view_yy = view.matrix;
	}
	return 0;
}

