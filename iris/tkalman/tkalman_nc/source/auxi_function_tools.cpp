#include "auxi_function_tools.hpp"
/**@fn unsigned int auxi_function_tools_p_vector(gsl_vector * p,
												 const gsl_vector * m)
 * @param p : vecteur contenant toutes les informations de la matrice de permutation permettant de regrouper les termes de la diagonale non nuls.
 * @param[in] m : vecteur contenant la diagonale de la matrice de contraintes M.
 * @brief
 * Cette fonction détermine la matrice de permutation permettant de regrouper tous les termes non nuls de M.
 * 
 */
unsigned int auxi_function_tools_p_vector(gsl_vector * p,
										  const gsl_vector * m)
{
	unsigned int rg = 0;
	unsigned int ker = 0;
	//Informations du vecteur M

	for (unsigned int i = 0; i < m->size; ++ i)
	{

		if ( (m->data[i * m->stride]) != 0)
		{
			p->data[rg * p->stride] = i;
			rg ++;
		}
		else
		{
			p->data[(p->size - ker - 1) * p->stride] = i;
			ker ++;
		}

	}
	//
	for (unsigned int i = m->size; i < p->size; ++ i)
	{
		p->data[rg * p->stride] = i;
		rg ++;
	}
	return rg;
}

/**@fn void auxi_function_tools_ch_var(gsl_matrix * sqrt_d,
									   gsl_matrix * y,
									   gsl_matrix * trans_mat,
									   gsl_matrix * trans_mat_view_01,
									  gsl_matrix * mat)
 * @param sqrt_d : [D]^{\frac{1}{2}}, racine du paramètre de la fonction
 * @param y : vecteur de translation du paramètre y
 * @param trans_mat : matrice de taille size(sqrt_d) allouée
 * @param trans_mat_view_01 : vue sur trans_mat
 * @param mat : matrice de taille size(sqrt_d) allouée
 * @brief
 * Cette fonction effectue un changement de variable sur le paramètre y de la fonction à maximiser.
 */
void auxi_function_tools_ch_var(gsl_matrix * sqrt_d,
								const gsl_matrix * y,
								gsl_matrix * trans_mat,
								gsl_matrix * trans_mat_view_01,
								gsl_matrix * mat)
{
	//Translation de la solution
	//Matrice de passage
	//     I    -yc
	// M = 
	//	   0     I
	gsl_matrix_set_identity(trans_mat);
	gsl_matrix_memcpy(trans_mat_view_01,
					  y);
	gsl_matrix_scale(trans_mat_view_01, -1);
	//Recopie
	
	gsl_matrix_memcpy(mat, sqrt_d);
	gsl_blas_dgemm( CblasNoTrans,
					CblasNoTrans,
					1.0,
					mat,
					trans_mat,
					0.0,
					sqrt_d);
}

/**@fn void auxi_function_tools_constrain_params(gsl_vector * p,
												 gsl_matrix * sqrt_d,
												 const gsl_vector * m,
												 gsl_vector * vect)
 * @param p : vecteur contenant l'information de la matrice de permutation (de taille size(sqrt_d))
 * @param sqrt_d : [D]^{\frac{1}{2}}, racine du paramètre de la fonction
 * @param m : vecteur de contrainte (si m(i) == 0 alors la collonne i de y ne sera pas estimée)
 * @param vect : vecteur de taille (size(sqrt_d)) alloué
 * @brief
 * Cette fonction applique la contraine sur le paramètre y et modifie les paramètres de la fonction à maximiser.
 */
void auxi_function_tools_constrain_params(gsl_matrix * sqrt_d,
										  const gsl_vector * m,
										  const gsl_vector * p,
										  unsigned int rg,
										  gsl_vector * vect)

{

	//Application des permutations à la matrice d
	for (unsigned int i = 0; i < rg; ++i)
	{
		gsl_vector_view view_dest = gsl_matrix_column (sqrt_d, i);
		gsl_vector_view view_src = gsl_matrix_column (sqrt_d, p->data[i]);
		gsl_vector_memcpy( &(view_dest.vector), &(view_src.vector) );
	}
	for (unsigned int i = rg; i < p->size; ++i)
	{
		gsl_vector_view view_dest = gsl_matrix_column (sqrt_d, i);
		gsl_vector_set_zero(&(view_dest.vector));
	}
	//Décompo QR pour obtenir une racine triangulaire
	gsl_linalg_QR_decomp (sqrt_d, vect);
	gsl_triangle_matrix(sqrt_d);
}

/**@fn void auxi_function_tools_maximize(gsl_matrix * sqrt_x,
										 gsl_matrix * y,
										 const gsl_matrix * sqrt_d,
										 const gsl_matrix * sqrt_d_view_00,
										 const gsl_matrix * sqrt_d_view_01,
										 const gsl_matrix * sqrt_d_view_11,
										 unsigned int size_x,
										 unsigned int n,
										 gsl_matrix * mat_yy,
										 gsl_permutation * perm_y)
 * @param sqrt_x : racine du premier paramètre
 * @param y : second paramètre
 * @param[in] sqrt_d : paramètre de la fonction
 * @param[in] sqrt_d_view_00 : vue sur sqrt_q
 * @param[in] sqrt_d_view_01 : vue sur sqrt_q
 * @param[in] sqrt_d_view_11 : vue sur sqrt_q
 * @param[in] size_x : dimension de x
 * @param[in] n : paramètre
 * @param mat_yy : matrice préallouée
 * @param perm_y : permutation préallouée
 * @brief
 * Cette fonction calcule les paramètres optimaux de la fonction
 */
void auxi_function_tools_maximize(gsl_matrix * sqrt_x,
								  gsl_matrix * y,
								  const gsl_matrix * sqrt_d,
								  const gsl_matrix * view_00,
								  const gsl_matrix * view_01,
								  const gsl_matrix * view_11,
								  unsigned int size_x,
								  unsigned int n,
								  gsl_matrix * mat_yy,
								  gsl_permutation * perm_y)
{
	
	//x max
	gsl_matrix_memcpy(sqrt_x, view_11);

	gsl_matrix_scale(sqrt_x, 1.0 / sqrt(n));
	
	//y max
	if (y != NULL && view_00 != NULL && view_01 != NULL)
	{
		gsl_permutation_init(perm_y);
		gsl_linalg_LU_invert(view_00,
							 perm_y,
							 mat_yy);
							 
		gsl_blas_dgemm(CblasNoTrans,
					   CblasNoTrans,
					   1.0,
					   mat_yy,
					   view_01,
					   0.0,
					   y);	
	}
}

/**@fn void auxi_function_tools_restore_y(gsl_matrix * y,
										  const gsl_matrix * y_c,
										  const gsl_matrix * z_max,
										  unsigned int rg,
										  const gsl_vector * p)
 * @param y : paramètre optimal reconstruit
 * @param y_c : valeur constante
 * @param z_max : paramètre réduit opt.
 * @param rg : rang
 * @param p : vecteur de permutation
 * @brief
 * Cette fonction reconstruit le paramètre y.
 */
void auxi_function_tools_restore_y(gsl_matrix * y,
								   const gsl_matrix * y_c,
								   const gsl_matrix * z_max,
								   unsigned int rg,
								   const gsl_vector * p)
{
	//Mise à zéro de y
	gsl_matrix_set_zero(y);
	
	//Reconstruction de y
	for (unsigned int i = 0; i < rg; ++i)
	{

		gsl_vector_const_view view_src = gsl_matrix_const_row (z_max, i);
		gsl_vector_view view_dest = gsl_matrix_row (y, p->data[i]);
		gsl_vector_memcpy( &(view_dest.vector), &(view_src.vector) );
	}
	gsl_matrix_add(y, y_c);

}
								   
/**@fn auxi_function_tools :: auxi_function_tools(const gsl_vector * m,
												  const gsl_matrix * yc,
												  unsigned int c_size) throw(exception &);
 * @param m : contrainte sur les lignes de y (si m(i) = 0 alors la ligne i de y ne sera pas estimée)
 * @param yc : partie constante de y
 * @param c_size : n_y + n_x
 * @brief
 *  Constructeur de la classe @class auxi_function_tools
 */
auxi_function_tools :: auxi_function_tools(const gsl_vector * m,
										   const gsl_matrix * yc) throw(exception &)
{
	initialize();
	try
	{
		setup(m, yc);
	}
	catch(exception & e)
	{
		throw(e);
	}
}						

/**@fn void auxi_function_tools :: setup(unsigned int y_size,
										 unsigned int c_size) throw(exception &););
 * @param y_size :
 * @param c_size :
 * @brief
 *  Setup de la classe @class auxi_function_tools
 */
void auxi_function_tools :: setup(const gsl_vector * m,
								  const gsl_matrix * yc) throw(exception &)
{
	unsigned int size_x,
				 size_y;
				 
	if (!m || !yc)
		throw(invalid_argument("m, yC or c_size are NULL!\n"));
	
	//Dim.
	size_y = yc->size1;
	size_x = yc->size2;

	
	//Vérif.;
	if (!size_y || !size_x)
		throw(invalid_argument("size_y or size_x are 0!\n"));

	//Modif des matrices
	free();
	initialize();
	
	_m = m;
	_yc = yc;
	_size_x = size_x;
	_size_y = size_y;
	_size_c = size_x + size_y;
	
	
	for (unsigned int i = 0; i < m->size; ++i)
	{
		if (m->data[i * m->stride])
			_rg ++;
	}	


	try
	{
		alloc();
	}
	catch (exception & except)
	{
		throw(except);
	}

	//Création des vues
	create_views();
	
	//Calcul du rang de M et de la perm.
	auxi_function_tools_p_vector(_p,
								 _m);

	
}

/**@fn void auxi_function_tools :: maximize(gsl_matrix * sqrt_x,
											gsl_matrix * y,
											gsl_matrix * sqrt_c,
											unsigned int n); 
 * @param sqrt_x : racine de de x
 * @param y : y
 * @param sqrt_c : racine du paramètre de la fonction
 * @param n : N
 * @brief
 * Cette fonction maximise la fonction auxiliaire.
 */
void auxi_function_tools :: maximize(gsl_matrix * sqrt_x,
									 gsl_matrix * y,
									 gsl_matrix * sqrt_d,
									 unsigned int n)
{

	auxi_function_tools_ch_var(sqrt_d,
							   _yc,
							   mat_cc_1,
							   &mat_cc_1_view_01,
							   mat_cc_2);
 
	auxi_function_tools_constrain_params(sqrt_d,
										 _m,
										 _p,
										 _rg + _size_x,
										 vect_c);

	//Création des vues
	{
		if (_rg != 0)
		{
			gsl_matrix_view view_00 = gsl_matrix_submatrix(sqrt_d, 0, 0, _rg, _rg);
			gsl_matrix_view view_01 = gsl_matrix_submatrix(sqrt_d, 0, _rg, _rg, _size_x);
			gsl_matrix_view view_11 = gsl_matrix_submatrix(sqrt_d, _rg, _rg, _size_x, _size_x);
			
			
			auxi_function_tools_maximize(sqrt_x,
										 z_max,
										 sqrt_d,
										 &(view_00.matrix),
										 &(view_01.matrix),
										 &(view_11.matrix),
										 _size_x,
										 n,
										 mat_zz,
										 perm_z);				 
		}
		else
		{
			gsl_matrix_view view_11 = gsl_matrix_submatrix(sqrt_d, _rg, _rg, _size_x, _size_x);
			auxi_function_tools_maximize(sqrt_x,
										 NULL,
										 sqrt_d,
										 NULL,
										 NULL,
										 &(view_11.matrix),
										 _size_x,
										 n,
										 NULL,
										 NULL);		
			
		}
	}					 
	auxi_function_tools_restore_y(y,
								  _yc,
								  z_max,
								  _rg,
								  _p);
}


/**@fn auxi_function_tools :: ~auxi_function_tools();
 * @brief
 * Destructeur de la classe @class auxi_function_tools
 */
auxi_function_tools :: ~auxi_function_tools()
{
	free();
	initialize();
	
}

/**@fn bool auxi_function_tools :: operator!() const;
 * @return
 * - 0 si l'objet est valide
 * - 1 sinon
 * 
 */			    
bool auxi_function_tools :: operator !() const
{
	return ( !(_size_x && _size_y && _size_c && _m && _yc && _p && mat_cc_1 && mat_cc_2 && mat_zz && z_max ) );
}


/**@fn void auxi_function_tools :: initialize();
 * @brief
 * Initialisateur de la classe @class auxi_function_tools
 */
void auxi_function_tools :: initialize()
{
	_size_x = 0;
	_size_y = 0;
	_size_c = 0;  
	_rg = 0;
	_m = 0;
	_yc = 0;
	_p = 0;		
	mat_cc_1 = 0;
	mat_cc_2 = 0;
	mat_zz = 0;
	z_max = 0;
	vect_c = 0;
	perm_z = 0;
}

/**@fn void auxi_function_tools :: free();
 * @brief
 * Cette méthode libère la mémoire occupée par les attributs de l'objet.
 */
void auxi_function_tools :: free()
{
	if (_p)
		gsl_vector_free(_p);
	if (mat_cc_1)
		gsl_matrix_free(mat_cc_1);
	if (mat_cc_2)
		gsl_matrix_free(mat_cc_2);
	if (mat_zz)
		gsl_matrix_free(mat_zz);
	if (z_max)
		gsl_matrix_free(z_max);
	if (vect_c)
		gsl_vector_free(vect_c);
	if (perm_z)
		gsl_permutation_free(perm_z);
}

/**@fn void auxi_function_tools :: alloc() throw(exception &);
 * @brief
 * Cette méthode alloue la mémoire utilisée par les attributs.
 */
void auxi_function_tools :: alloc() throw(exception &)
{
	if (!_p)
	{
		try
		{
			_p = gsl_vector_alloc(_size_c);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!vect_c)
	{
		try
		{
			vect_c = gsl_vector_alloc(_size_c);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!mat_cc_1)
	{
		try
		{
			mat_cc_1 = gsl_matrix_alloc(_size_c, _size_c);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!mat_cc_2)
	{
		try
		{
			mat_cc_2 = gsl_matrix_alloc(_size_c, _size_c);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	
	if (!mat_zz)
	{
		if (_rg)
		{
			try
			{
				mat_zz = gsl_matrix_alloc(_rg, _rg);
			}
			catch(exception & e)
			{
				throw(e);
			}
		}
	}
	if (!z_max)
	{
		if (_rg)
		{
			try
			{
				z_max = gsl_matrix_alloc(_rg, _size_x);
			}
			catch(exception & e)
			{
				throw(e);
			}
		}
	}
	if (!perm_z)
	{
		if (_rg)
		{
			try
			{
				perm_z = gsl_permutation_alloc(_rg);
			}
			catch(exception & e)
			{
				throw(e);
			}
		}
	}

}

/**@fn void auxi_function_tools :: create_views();
 * @brief
 * Cette méthode crée les vues...
 */
void auxi_function_tools :: create_views()
{
	gsl_matrix_view view;
	view = gsl_matrix_submatrix(mat_cc_1,
								0,
								_size_y,
								_size_y,
								_size_x);
	
	mat_cc_1_view_01 = view.matrix;
}

