#include "tkalman_nc_sums_base.hpp"


void tkalman_nc_get_cov_zx_ ( const gsl_matrix * sqrt_z_f,
							  const gsl_matrix * sqrt_p_s_,
							  const gsl_matrix * c_s,
							  const gsl_matrix * f_xz,
							  const gsl_matrix * sqrt_q_xx,
							  gsl_matrix * mat_2xpz_2xpz,
							  gsl_matrix * mat_2xpz_2xpz_view_00,
							  gsl_matrix * mat_2xpz_2xpz_view_10,
							  gsl_matrix * mat_2xpz_2xpz_view_11,
							  gsl_matrix * mat_2xpz_2xpz_view_21,
							  gsl_matrix * mat_2xpz_2xpz_view_22,
							  gsl_vector * vect_2xpz             )
{
	//Mise à zéro de la matrice M
	gsl_matrix_set_zero(mat_2xpz_2xpz);
	
	//Construction de la matrice
	//M00
	gsl_matrix_memcpy(mat_2xpz_2xpz_view_00, 
					  sqrt_q_xx);
	//M10
	gsl_blas_dgemm(CblasNoTrans,
				   CblasTrans,
				   1.0,
				   sqrt_z_f,
				   f_xz,
				   0.0,
				   mat_2xpz_2xpz_view_10);
	//M11
	gsl_matrix_memcpy(mat_2xpz_2xpz_view_11, sqrt_z_f);
	
	//M21
	gsl_matrix_memcpy(mat_2xpz_2xpz_view_21, c_s);
	
	//M22
	gsl_matrix_memcpy(mat_2xpz_2xpz_view_22, sqrt_p_s_);
	
	//Décomposition QR
	gsl_linalg_QR_decomp(mat_2xpz_2xpz,
						 vect_2xpz);
	gsl_triangle_matrix(mat_2xpz_2xpz);
	
}

void tkalman_nc_get_corr_tx_ ( const gsl_vector * x_s,
							   const gsl_vector * _y,
							   const gsl_vector * x_s_,
							   const gsl_matrix * sqrt_cov_xx_view_00,
							   const gsl_matrix * sqrt_cov_xx_view_01,
							   const gsl_matrix * sqrt_cov_xx_view_11,
							   gsl_matrix * mat_tpxp1_tpx,
							   gsl_matrix * mat_tpxp1_tpx_view_00,
							   gsl_matrix * mat_tpxp1_tpx_view_02, 
							   gsl_matrix * mat_tpxp1_tpx_view_22,
							   gsl_vector * mat_tpxp1_tpx_view_30,
							   gsl_vector * mat_tpxp1_tpx_view_31,
							   gsl_vector * mat_tpxp1_tpx_view_32,
							   gsl_vector * vect_tpx )
{
	
	
	//Construction de la matrice M
	gsl_matrix_set_zero ( mat_tpxp1_tpx );
	gsl_matrix_memcpy ( mat_tpxp1_tpx_view_00,
						sqrt_cov_xx_view_00);
						
		
	gsl_matrix_memcpy ( mat_tpxp1_tpx_view_02,
						sqrt_cov_xx_view_01 );
	
	gsl_matrix_memcpy ( mat_tpxp1_tpx_view_22,
						sqrt_cov_xx_view_11 );
	
	gsl_vector_memcpy( mat_tpxp1_tpx_view_30,
					   x_s);

	gsl_vector_memcpy( mat_tpxp1_tpx_view_31,
					   _y);

	gsl_vector_memcpy( mat_tpxp1_tpx_view_32,
					   x_s_);

	//Décomposition QR
	gsl_linalg_QR_decomp(mat_tpxp1_tpx, vect_tpx);
	gsl_triangle_matrix(mat_tpxp1_tpx);
	
}

void tkalman_nc_get_corr_ty ( const gsl_vector * x_s,
							  const gsl_vector * _y,
							  const gsl_vector * y,
							  const gsl_matrix * sqrt_z_s,
							  gsl_matrix * mat_tpyp1_tpy,
							  gsl_matrix * mat_tpyp1_tpy_view_00,
							  gsl_vector * mat_tpyp1_tpy_view_30,
							  gsl_vector * mat_tpyp1_tpy_view_31,
							  gsl_vector * mat_tpyp1_tpy_view_32,
							  gsl_vector * vect_tpy
							 )
{
	//Construction de la matrice M
	gsl_matrix_set_zero ( mat_tpyp1_tpy );
	gsl_matrix_memcpy ( mat_tpyp1_tpy_view_00,
						sqrt_z_s );
	
	gsl_vector_memcpy( mat_tpyp1_tpy_view_30,
					   x_s);
	
	gsl_vector_memcpy( mat_tpyp1_tpy_view_31,
					   _y);
	
	gsl_vector_memcpy( mat_tpyp1_tpy_view_32,
					   y);

	//Décomposition QR
	gsl_linalg_QR_decomp(mat_tpyp1_tpy, vect_tpy);
	gsl_triangle_matrix(mat_tpyp1_tpy);

}

tkalman_nc_sums_base :: tkalman_nc_sums_base( const gsl_matrix * f_xt,
											  const gsl_matrix * sqrt_q_xx) throw (exception &)
{
	initialize();
	//Vérif des arguments
	if (!f_xt || !sqrt_q_xx)
		throw(invalid_argument("F^{x,t}, sqrt(Q^{x,x}) are NULL!\n"));
		
	unsigned int size_t, 
				 size_x,
				 size_y;
	size_t = f_xt->size2;
	size_x = sqrt_q_xx->size1;
	size_y = size_t - size_x;
	
	//Vérif de la dim. de y;
	if (!size_y or !size_x)
		throw(invalid_argument("size_y or size_x are 0!\n"));
	
	//Modif des matrices
	_size_x = size_x;
	_size_y = size_y;
	_size_t = size_t;
	_f_xt = f_xt;
	_sqrt_q_xx = sqrt_q_xx;
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
}



void tkalman_nc_sums_base :: setup ( 	const gsl_matrix * f_xt,
										const gsl_matrix * sqrt_q_xx) throw (exception &)
{
	free();
	initialize();
	//Vérif des arguments
	if (!f_xt || !sqrt_q_xx)
		throw(invalid_argument("F^{x,t}, sqrt(Q^{x,x}) are NULL!\n"));
		
	unsigned int size_t, 
				 size_x,
				 size_y;
	size_t = f_xt->size2;
	size_x = sqrt_q_xx->size1;
	size_y = size_t - size_x;
	
	//Vérif de la dim. de y;
	if (!size_y or !size_x)
		throw(invalid_argument("size_y or size_x are 0!\n"));
	
	//Modif des matrices
	_size_x = size_x;
	_size_y = size_y;
	_size_t = size_t;
	_f_xt = f_xt;
	_sqrt_q_xx = sqrt_q_xx;
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
	
}

tkalman_nc_sums_base  :: ~tkalman_nc_sums_base ()
{
	free();
	initialize();
}

void tkalman_nc_sums_base :: initialize()
{
	_f_xt = 0;
	_sqrt_q_xx = 0;
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
	
	
	mat_2xpt_2xpt = 0;
	vect_2xpt = 0;
	
	mat_3x_3x = 0;
	vect_3x = 0;
	mat_2tp2xp1_tpx = 0;
	
	vect_tpx = 0;
	mat_2tp2yp1_tpy = 0;
	vect_tpy = 0;
	
}

void tkalman_nc_sums_base :: alloc() throw(exception &)
{
	if (!mat_2xpt_2xpt)
	{
		try
		{
			mat_2xpt_2xpt = gsl_matrix_alloc(2 * _size_x + _size_t, 2 * _size_x + _size_t);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!vect_2xpt)
	{
		try
		{
			vect_2xpt = gsl_vector_alloc(2 * _size_x + _size_t);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!mat_3x_3x)
	{
		try
		{
			mat_3x_3x = gsl_matrix_alloc(3 * _size_x, 3 * _size_x);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!vect_3x)
	{
		try
		{
			vect_3x = gsl_vector_alloc(3 * _size_x);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!mat_2tp2xp1_tpx)
	{
		try
		{
			mat_2tp2xp1_tpx = gsl_matrix_alloc(2 * _size_t + 2 * _size_x + 1,
											   _size_x + _size_t);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!vect_tpx)
	{
		try
		{
			vect_tpx = gsl_vector_alloc(_size_t + _size_x);
		}
		catch(exception & e)
		{
			throw(e);
		}
	}
	if (!mat_2tp2yp1_tpy)
	{
		try
		{
			mat_2tp2yp1_tpy = gsl_matrix_alloc(2 * _size_t + 2 * _size_y + 1,
											   _size_t + _size_y);
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
	
}

void tkalman_nc_sums_base :: free()
{
	if (mat_2xpt_2xpt)
	{
		gsl_matrix_free(mat_2xpt_2xpt);
	}
	if (vect_2xpt)
	{
		gsl_vector_free(vect_2xpt);
	}
	if (mat_3x_3x)
	{
		gsl_matrix_free(mat_3x_3x);
	}
	if (vect_3x)
	{
		gsl_vector_free(vect_3x);
	}
	if (mat_2tp2xp1_tpx)
	{
		gsl_matrix_free(mat_2tp2xp1_tpx);
	}
	if (vect_tpx)
	{
		gsl_vector_free(vect_tpx);
	}
	if (mat_2tp2yp1_tpy)
	{
		gsl_matrix_free(mat_2tp2yp1_tpy);
	}
	if (vect_tpy)
	{
		gsl_vector_free(vect_tpy);
	}
	
}

void tkalman_nc_sums_base :: create_views()
{
	
	gsl_matrix_view derdesder;
	gsl_vector_view k2000;
	//mat_2xpt_2xpt
	{
		//mat_2xpt_2xpt_view_00
		derdesder = gsl_matrix_submatrix(mat_2xpt_2xpt,
										 0,
										 0,
										 _size_x,
										 _size_x);
		mat_2xpt_2xpt_view_00 = derdesder.matrix;
		
		
		//mat_2xpt_2xpt_view_10
		derdesder = gsl_matrix_submatrix(mat_2xpt_2xpt,
										 _size_x,
										 0,
										 _size_t,
										 _size_x);
		mat_2xpt_2xpt_view_10 = derdesder.matrix;
		
		//mat_2xpt_2xpt_view_11
		derdesder = gsl_matrix_submatrix(mat_2xpt_2xpt,
										 _size_x,
										 _size_x,
										 _size_t,
										 _size_t);
		mat_2xpt_2xpt_view_11 = derdesder.matrix;
		
		//mat_2xpt_2xpt_view_12
		derdesder = gsl_matrix_submatrix(mat_2xpt_2xpt,
										 _size_x,
										 _size_x + _size_t,
										 _size_t,
										 _size_x);
		mat_2xpt_2xpt_view_12 = derdesder.matrix;
		
		
		//mat_2xpt_2xpt_view_21
		derdesder = gsl_matrix_submatrix(mat_2xpt_2xpt,
										 _size_x + _size_t,
										 _size_x,
										 _size_x,
										 _size_t);
		mat_2xpt_2xpt_view_21 = derdesder.matrix;
		
		//mat_2xpt_2xpt_view_22
		derdesder = gsl_matrix_submatrix(mat_2xpt_2xpt,
										 _size_x + _size_t,
										 _size_x + _size_t,
										 _size_x,
										 _size_x);
		mat_2xpt_2xpt_view_22 = derdesder.matrix;
	}
	//mat_3x_3x
	{
		//mat_3x_3x_view_00
		derdesder = gsl_matrix_submatrix(mat_3x_3x,
										 0,
										 0,
										 _size_x,
										 _size_x);
		mat_3x_3x_view_00 = derdesder.matrix;
		
		
		//mat_3x_3x_view_10
		derdesder = gsl_matrix_submatrix(mat_3x_3x,
										 _size_x,
										 0,
										 _size_x,
										 _size_x);
		mat_3x_3x_view_10 = derdesder.matrix;
		
		//mat_3x_3x_view_11
		derdesder = gsl_matrix_submatrix(mat_3x_3x,
										 _size_x,
										 _size_x,
										 _size_x,
										 _size_x);
		mat_3x_3x_view_11 = derdesder.matrix;
		
		//mat_3x_3x_view_12
		derdesder = gsl_matrix_submatrix(mat_3x_3x,
										 _size_x,
										 _size_x + _size_x,
										 _size_x,
										 _size_x);
		mat_3x_3x_view_12 = derdesder.matrix;
		
		
		//mat_3x_3x_view_21
		derdesder = gsl_matrix_submatrix(mat_3x_3x,
										 _size_x + _size_x,
										 _size_x,
										 _size_x,
										 _size_x);
		mat_3x_3x_view_21 = derdesder.matrix;
		
		//mat_3x_3x_view_22
		derdesder = gsl_matrix_submatrix(mat_3x_3x,
										 _size_x + _size_x,
										 _size_x + _size_x,
										 _size_x,
										 _size_x);
		mat_3x_3x_view_22 = derdesder.matrix;
	}

	//mat_2tp2xp1_tpx
	{
		//mat_2tp2xp1_tpx_view_00
		derdesder = gsl_matrix_submatrix(mat_2tp2xp1_tpx,
										 0,
										 0,
										 _size_t + _size_x,
										 _size_t + _size_x);
		
		mat_2tp2xp1_tpx_view_00 = derdesder.matrix;
		
		//mat_2tp2xp1_tpx_view_tpxp1_tpx
		derdesder = gsl_matrix_submatrix(mat_2tp2xp1_tpx,
										 _size_t + _size_x,
										 0,
										 _size_t + _size_x + 1,
										 _size_t + _size_x);
		
		mat_2tp2xp1_tpx_view_tpxp1_tpx = derdesder.matrix;
		
		//mat_tpxp1_tpx_view_00
		derdesder = gsl_matrix_submatrix(mat_2tp2xp1_tpx,
										 _size_t + _size_x,
										 0,
										 _size_x,
										 _size_x);
		
		mat_tpxp1_tpx_view_00 = derdesder.matrix;
		
		//mat_tpxp1_tpx_view_00_bis
		derdesder = gsl_matrix_submatrix(mat_2tp2xp1_tpx,
										 _size_t + _size_x,
										 0,
										 _size_t,
										 _size_t);
		
		mat_tpxp1_tpx_view_00_bis = derdesder.matrix;
		
		//mat_tpxp1_tpx_view_02
		derdesder = gsl_matrix_submatrix(mat_2tp2xp1_tpx,
										 _size_t + _size_x,
										 _size_t,
										 _size_x,
										 _size_x);
		
		mat_tpxp1_tpx_view_02 = derdesder.matrix;
		
		//mat_tpxp1_tpx_view_02_bis
		derdesder = gsl_matrix_submatrix(mat_2tp2xp1_tpx,
										 _size_t + _size_x,
										 _size_t,
										 _size_t,
										 _size_x);
		
		mat_tpxp1_tpx_view_02_bis = derdesder.matrix;
		
		//mat_tpxp1_tpx_view_22
		derdesder = gsl_matrix_submatrix(mat_2tp2xp1_tpx,
										 2 * _size_t + _size_x,
										 _size_t,
										 _size_x,
										 _size_x);
		
		mat_tpxp1_tpx_view_22 = derdesder.matrix;
		
		//mat_tpxp1_tpx_view_30
		k2000 = gsl_matrix_subrow(mat_2tp2xp1_tpx,
								  2 * _size_t + 2 * _size_x,
								  0,
								  _size_x);
		mat_tpxp1_tpx_view_30 = k2000.vector;
		
		//mat_tpxp1_tpx_view_31
		k2000 = gsl_matrix_subrow(mat_2tp2xp1_tpx,
								  2 * _size_t + 2 * _size_x,
								  _size_x,
								  _size_y);
		mat_tpxp1_tpx_view_31 = k2000.vector;
		
		//mat_tpxp1_tpx_view_32
		k2000 = gsl_matrix_subrow(mat_2tp2xp1_tpx,
								  2 * _size_t + 2 * _size_x,
								  _size_t,
								  _size_x);
		mat_tpxp1_tpx_view_32 = k2000.vector;
		
	}

	//mat_2tp2yp1_tpy
	{
		//mat_2tp2yp1_tpy_view_00
		derdesder = gsl_matrix_submatrix(mat_2tp2yp1_tpy,
										 0,
										 0,
										 _size_t + _size_y,
										 _size_t + _size_y);
		
		mat_2tp2yp1_tpy_view_00 = derdesder.matrix;
		//mat_2tp2yp1_tpy_view_tpxp1_tpx
		derdesder = gsl_matrix_submatrix(mat_2tp2yp1_tpy,
										 _size_t + _size_y,
										 0,
										 _size_t + _size_y + 1,
										 _size_t + _size_y);
		
		mat_2tp2yp1_tpy_view_tpyp1_tpy = derdesder.matrix;

		//mat_tpyp1_tpy_view_00
		derdesder = gsl_matrix_submatrix(mat_2tp2yp1_tpy,
										 _size_t + _size_y,
										 0,
										 _size_x,
										 _size_x);
		
		mat_tpyp1_tpy_view_00 = derdesder.matrix;

		//mat_tpyp1_tpy_view_00_bis
		derdesder = gsl_matrix_submatrix(mat_2tp2yp1_tpy,
										 _size_t + _size_x,
										 0,
										 _size_t,
										 _size_t);
		
		mat_tpyp1_tpy_view_00_bis = derdesder.matrix;

		//mat_tpyp1_tpy_view_30
		k2000 = gsl_matrix_subrow(mat_2tp2yp1_tpy,
								  2 * _size_t + 2 * _size_y,
								  0,
								  _size_x);
		mat_tpyp1_tpy_view_30 = k2000.vector;
		
		//mat_tpyp1_tpy_view_31
		k2000 = gsl_matrix_subrow(mat_2tp2yp1_tpy,
								  2 * _size_t + 2 * _size_y,
								  _size_x,
								  _size_y);
		mat_tpyp1_tpy_view_31 = k2000.vector;
		
		//mat_tpyp1_tpy_view_32
		k2000 = gsl_matrix_subrow(mat_2tp2yp1_tpy,
								  2 * _size_t + 2 * _size_y,
								  _size_t,
								  _size_y);
		mat_tpyp1_tpy_view_32 = k2000.vector;
		
	}
}







