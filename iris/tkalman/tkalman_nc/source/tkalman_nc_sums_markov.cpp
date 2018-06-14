#include "tkalman_nc_sums_markov.hpp"

#include <iostream>
using namespace std;
tkalman_nc_sums_markov :: tkalman_nc_sums_markov( const gsl_matrix * f_xt,
													const gsl_matrix * sqrt_q_xx,
													unsigned int p) throw (exception &)
: tkalman_nc_sums_base( f_xt,
						sqrt_q_xx)
{
	tkalman_nc_sums_markov :: initialize();
	_p = p;
	tkalman_nc_sums_markov :: alloc();
	tkalman_nc_sums_markov :: create_views();
}

void tkalman_nc_sums_markov :: setup ( const gsl_matrix * f_xt,
										 const gsl_matrix * sqrt_q_xx,
										 unsigned int p) throw (exception &)
{
	tkalman_nc_sums_base :: setup( f_xt, sqrt_q_xx);	
	tkalman_nc_sums_markov :: free();
	tkalman_nc_sums_markov :: initialize();
	_p = p;
	tkalman_nc_sums_markov :: alloc();
	tkalman_nc_sums_markov :: create_views();
}

tkalman_nc_sums_markov :: ~tkalman_nc_sums_markov()
{
	tkalman_nc_sums_markov :: free();
	tkalman_nc_sums_markov :: initialize();
}

void tkalman_nc_sums_markov :: compute_sqrt_sum_corr_tx_(	gsl_matrix * sqrt_sum_corr_tx_,
															const gsl_vector * const * z_s,
															const gsl_vector * const * y,
															const gsl_matrix * const * sqrt_z_f,
															const gsl_matrix * const * sqrt_p_s,
															const gsl_matrix * const * c_s,
															const unsigned int * r,
															unsigned int r_id,
															unsigned int n)
{
	//Mise à zéro
	gsl_matrix_set_zero(&mat_2tp2xp1_tpx_view_00);

	//Boucle
	for (unsigned int i = n; i > 1; -- i)
	{

		if ( r[i - 1] == r_id)
		{
			tkalman_nc_get_cov_zx_ ( sqrt_z_f[i - 1],
									 sqrt_p_s[i],
									 c_s[i - 1],
									 _f_xx + r_id,
									 _sqrt_q_xx + r_id,
									 mat_3x_3x,
									 &mat_3x_3x_view_00,
									 &mat_3x_3x_view_10,
									 &mat_3x_3x_view_11,
									 &mat_3x_3x_view_21,
									 &mat_3x_3x_view_22,
									 vect_3x);

			tkalman_nc_get_corr_tx_ ( z_s[i - 1],
									  y[i - 2],
									  z_s[i],
									  &mat_3x_3x_view_11,
									  &mat_3x_3x_view_12,
									  &mat_3x_3x_view_22,
									  &mat_2tp2xp1_tpx_view_tpxp1_tpx,
									  &mat_tpxp1_tpx_view_00,
									  &mat_tpxp1_tpx_view_02, 
									  &mat_tpxp1_tpx_view_22,
									  &mat_tpxp1_tpx_view_30,
									  &mat_tpxp1_tpx_view_31,
									  &mat_tpxp1_tpx_view_32,
									  vect_tpx );
					  
			gsl_linalg_QR_decomp(mat_2tp2xp1_tpx,
								 vect_tpx );
			gsl_triangle_matrix(&mat_2tp2xp1_tpx_view_00);	

		}  
								  
	}

	//0
	if ( r[0] == r_id)
	{
		tkalman_nc_get_cov_zx_ ( sqrt_z_f[0],
								 sqrt_p_s[1],
								 c_s[0],
								 _f_xt  + r_id,
								 _sqrt_q_xx  + r_id,
								 mat_2xpt_2xpt,
								 &mat_2xpt_2xpt_view_00,
								 &mat_2xpt_2xpt_view_10,
								 &mat_2xpt_2xpt_view_11,
								 &mat_2xpt_2xpt_view_21,
								 &mat_2xpt_2xpt_view_22,
								 vect_2xpt);
								 
		gsl_vector_const_view view_x = gsl_vector_const_subvector(z_s[0], 0, _size_x);
		gsl_vector_const_view view_y = gsl_vector_const_subvector(z_s[0], _size_x, _size_y);
		tkalman_nc_get_corr_tx_ ( &(view_x.vector),
								  &(view_y.vector),
								  z_s[1],
								  &mat_2xpt_2xpt_view_11,
								  &mat_2xpt_2xpt_view_12,
								  &mat_2xpt_2xpt_view_22,
								  &mat_2tp2xp1_tpx_view_tpxp1_tpx,
								  &mat_tpxp1_tpx_view_00_bis,
								  &mat_tpxp1_tpx_view_02_bis, 
								  &mat_tpxp1_tpx_view_22,
								  &mat_tpxp1_tpx_view_30,
								  &mat_tpxp1_tpx_view_31,
								  &mat_tpxp1_tpx_view_32,
								  vect_tpx );
		gsl_linalg_QR_decomp(mat_2tp2xp1_tpx,
							vect_tpx );
	gsl_triangle_matrix(&mat_2tp2xp1_tpx_view_00);				  
								  
		
	}
	



	(*sqrt_sum_corr_tx_) = mat_2tp2xp1_tpx_view_00;
	
	
}

void tkalman_nc_sums_markov :: compute_sqrt_sum_corr_ty (  	gsl_matrix * sqrt_sum_corr_ty,
															const gsl_vector * const * z_s,
															const gsl_vector * const * y,
															const gsl_matrix * const * sqrt_z_s,
															const unsigned int * r,
															unsigned int r_id,
															unsigned int n)
{
	//Mise à zéro
	gsl_matrix_set_zero(&mat_2tp2yp1_tpy_view_00);
	
	//Boucle
	for (unsigned int i = n; i > 1; -- i)
	{
		if ( r[i - 1] == r_id)
		{
			tkalman_nc_get_corr_ty ( z_s[i - 1],
									 y[i - 2],
									 y[i - 1],
									 sqrt_z_s[i - 1],
									 &mat_2tp2yp1_tpy_view_tpyp1_tpy,
									 &mat_tpyp1_tpy_view_00,
									 &mat_tpyp1_tpy_view_30,
									 &mat_tpyp1_tpy_view_31,
									 &mat_tpyp1_tpy_view_32,
									 vect_tpy
									 );
			gsl_linalg_QR_decomp(mat_2tp2yp1_tpy,
								 vect_tpy );
			gsl_triangle_matrix(&mat_2tp2yp1_tpy_view_00);	  
		}
								  
	}

	//0
	if ( r[0] == r_id)
	{
		gsl_vector_const_view view_x = gsl_vector_const_subvector(z_s[0], 0, _size_x);
		gsl_vector_const_view view_y = gsl_vector_const_subvector(z_s[0], _size_x, _size_y);
		
		tkalman_nc_get_corr_ty ( &(view_x.vector),
								 &(view_y.vector),
								 y[0],
								 sqrt_z_s[0],
								 &mat_2tp2yp1_tpy_view_tpyp1_tpy,
								 &mat_tpyp1_tpy_view_00_bis,
								 &mat_tpyp1_tpy_view_30,
								 &mat_tpyp1_tpy_view_31,
								 &mat_tpyp1_tpy_view_32,
								 vect_tpy
								 );
								 
		gsl_linalg_QR_decomp(mat_2tp2yp1_tpy,
							 vect_tpy );
		gsl_triangle_matrix(&mat_2tp2yp1_tpy_view_00);	 
	}
	
	(*sqrt_sum_corr_ty) = mat_2tp2yp1_tpy_view_00;

}


void tkalman_nc_sums_markov :: initialize()
{
	
	
	_f_xx = 0;
	_f_xy = 0;
	_p = 0;	
}

void tkalman_nc_sums_markov :: alloc() throw(exception &)
{
	_f_xx = new gsl_matrix[_p];
	_f_xy = new gsl_matrix[_p];
}

void tkalman_nc_sums_markov :: free()
{
	if (_f_xx)
		delete[] _f_xx;
	if (_f_xy)
		delete[] _f_xy;
}

void tkalman_nc_sums_markov :: create_views()
{
	for ( unsigned int i = 0; i < _p; ++ i )
	{
		{
			gsl_matrix_const_view apero = gsl_matrix_const_submatrix(_f_xt + i, 0, 0, _size_x, _size_x);
			
			_f_xx[i] = apero.matrix;
		}
		{
			gsl_matrix_const_view apero = gsl_matrix_const_submatrix(_f_xt + i, 0, _size_x, _size_x, _size_y);
			_f_xy[i] = apero.matrix;
		}
	}

}
