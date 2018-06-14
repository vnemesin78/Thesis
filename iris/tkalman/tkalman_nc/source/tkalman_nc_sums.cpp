#include "tkalman_nc_sums.hpp"

tkalman_nc_sums :: tkalman_nc_sums( const gsl_matrix * f_xt,
									const gsl_matrix * sqrt_q_xx) throw (exception &)
: tkalman_nc_sums_base( f_xt,
						sqrt_q_xx)
{
	tkalman_nc_sums :: create_views();
}

void tkalman_nc_sums :: setup ( const gsl_matrix * f_xt,
								const gsl_matrix * sqrt_q_xx) throw (exception &)
{
	tkalman_nc_sums_base :: setup( f_xt, sqrt_q_xx);
	tkalman_nc_sums :: create_views();
	
}

void tkalman_nc_sums :: compute_sqrt_sum_corr_tx_(gsl_matrix * sqrt_sum_corr_tx_,
												  const gsl_vector * const * z_s,
												  const gsl_vector * const * y,
												  const gsl_matrix * const * sqrt_z_f,
												  const gsl_matrix * const * sqrt_p_s,
												  const gsl_matrix * const * c_s,
												  unsigned int n)
{
	//Mise à zéro
	gsl_matrix_set_zero(&mat_2tp2xp1_tpx_view_00);

	//Boucle
	for (unsigned int i = n; i > 1; -- i)
	{

		tkalman_nc_get_cov_zx_ ( sqrt_z_f[i - 1],
								 sqrt_p_s[i],
								 c_s[i - 1],
								 &_f_xx,
								 _sqrt_q_xx,
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

	//0
	tkalman_nc_get_cov_zx_ ( sqrt_z_f[0],
							 sqrt_p_s[1],
							 c_s[0],
							 _f_xt,
							 _sqrt_q_xx,
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

	(*sqrt_sum_corr_tx_) = mat_2tp2xp1_tpx_view_00;
	
}

void tkalman_nc_sums :: compute_sqrt_sum_corr_ty ( gsl_matrix * sqrt_sum_corr_ty,
												   const gsl_vector * const * z_s,
												   const gsl_vector * const * y,
												   const gsl_matrix * const * sqrt_z_s,
												   unsigned int n)
{
	//Mise à zéro
	gsl_matrix_set_zero(&mat_2tp2yp1_tpy_view_00);
	
	//Boucle
	for (unsigned int i = n; i > 1; -- i)
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

	//0
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

	(*sqrt_sum_corr_ty) = mat_2tp2yp1_tpy_view_00;
}

void tkalman_nc_sums :: create_views()
{
	{
		gsl_matrix_const_view apero = gsl_matrix_const_submatrix(_f_xt, 0, 0, _size_x, _size_x);
		_f_xx = apero.matrix;
	}
	{
		gsl_matrix_const_view apero = gsl_matrix_const_submatrix(_f_xt, 0, _size_x, _size_x, _size_y);
		_f_xy = apero.matrix;
	}
}






