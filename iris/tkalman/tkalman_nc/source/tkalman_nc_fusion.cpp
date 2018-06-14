#include "tkalman_nc_fusion.hpp"
tkalman_nc_fusion :: tkalman_nc_fusion ( unsigned int size_x,
										 unsigned int size_y)
{
	initialize();
	setup( size_x, size_y );
}

void tkalman_nc_fusion :: setup ( unsigned int size_x,
								  unsigned int size_y )
{
	free();
	initialize();
	_size_x = size_x;
	_size_y = size_y;
	_size_t = size_x + size_y;
	
	//Alloc
	mat_2t_x = gsl_matrix_alloc( 2 * (_size_t + _size_x),(_size_t + _size_x) );
	mat_2t_y = gsl_matrix_alloc( 2 * (_size_t + _size_y), (_size_t + _size_y) );
	vect_x = gsl_vector_alloc( (_size_t + _size_x) );
	vect_y = gsl_vector_alloc( (_size_t + _size_y) );
	vect_t = gsl_vector_alloc( _size_t );
	mat_2t_t = gsl_matrix_alloc ( 2 * _size_t, _size_t );
	//Création des vues de merde
	gsl_matrix_view view = gsl_matrix_submatrix ( 	mat_2t_x ,
													0,
													0,
													(_size_t + _size_x),
													(_size_t + _size_x));
	mat_2t_x_view_00 = view.matrix;
	view = gsl_matrix_submatrix ( 	mat_2t_x ,
													(_size_t + _size_x),
													0,
													(_size_t + _size_x),
													(_size_t + _size_x));
	mat_2t_x_view_10 = view.matrix;
	
	view = gsl_matrix_submatrix ( 	mat_2t_y ,
													0,
													0,
													(_size_t + _size_y),
													(_size_t + _size_y));
	mat_2t_y_view_00 = view.matrix;
	view = gsl_matrix_submatrix ( 	mat_2t_y ,
													(_size_t + _size_y),
													0,
													(_size_t + _size_y),
													(_size_t + _size_y));
	mat_2t_y_view_10 = view.matrix;
	
	view = gsl_matrix_submatrix ( 	mat_2t_t ,
									0,
									0,
									_size_t,
									_size_t );
	mat_2t_t_view_0 = view.matrix;
	
	view = gsl_matrix_submatrix ( 	mat_2t_t ,
									_size_t,
									0,
									_size_t,
									_size_t );
	mat_2t_t_view_1 = view.matrix;
	
	view = gsl_matrix_submatrix ( 	mat_2t_t ,
									0,
									0,
									_size_t + 1,
									_size_t );
	mat_2t_t_view_tp1_t = view.matrix;
	
	gsl_vector_view sparadrap = gsl_matrix_row ( mat_2t_t,
												 _size_t );
	mat_2t_t_view_vector = sparadrap.vector;
	
	
	
	
	
	
	
	
}

void tkalman_nc_fusion :: fusion ( 	gsl_matrix * sqrt_sum_tx,
									gsl_matrix * sqrt_sum_ty,
									const gsl_matrix * const * sqrt_sums_tx,
									const gsl_matrix * const * sqrt_sums_ty,
									unsigned int n,
									unsigned int step )
{
	gsl_matrix_set_zero( mat_2t_x );
	gsl_matrix_set_zero( mat_2t_y );
	
	for (unsigned int i = 0; i < n; ++ i)
	{
		gsl_matrix_memcpy( &mat_2t_x_view_10, 
						   sqrt_sums_tx[i * step]);
		
		gsl_matrix_memcpy( &mat_2t_y_view_10, 
						   sqrt_sums_ty[i * step]);
						   
		gsl_linalg_QR_decomp(mat_2t_x,
							 vect_x );
		gsl_linalg_QR_decomp(mat_2t_y,
							 vect_y );
							 
		gsl_triangle_matrix(&mat_2t_x_view_00);
		gsl_triangle_matrix(&mat_2t_y_view_00);
						   
	}
						
	*sqrt_sum_tx = mat_2t_x_view_00;
	*sqrt_sum_ty = mat_2t_y_view_00;
}

void tkalman_nc_fusion :: fusion_0 (	gsl_vector * t_0,
										gsl_matrix * sqrt_q_0,
										const gsl_vector * const * t_s_0,
										const gsl_matrix * const * sqrt_q_s_0,
										unsigned int nb_signals )
{
	
	//Calcul de la moyenne
	gsl_vector_set_zero ( t_0 );
	for ( unsigned int i = 0; i < nb_signals; ++ i )
	{
		gsl_vector_add ( t_0, t_s_0[i] );
	}
	gsl_vector_scale ( t_0, 1.0/nb_signals );
	
	//Calcul de la matrice de covariance
	gsl_matrix_set_zero ( mat_2t_t );
	for ( unsigned int i = 0; i < nb_signals; ++ i )
	{
		gsl_matrix_memcpy (	&mat_2t_t_view_1,
							sqrt_q_s_0[i] );
		gsl_linalg_QR_decomp(	mat_2t_t,
								vect_t );
		gsl_triangle_matrix(	&mat_2t_t_view_0);
		
		//x0
		gsl_vector_memcpy( 	&mat_2t_t_view_vector,
							t_0 );
		gsl_vector_sub ( &mat_2t_t_view_vector,
						 t_s_0[i] );
						 
		for ( unsigned int j = _size_y; j < _size_t; ++ j )
		{
			mat_2t_t_view_vector.data[j * mat_2t_t_view_vector.stride] = 0;
		}
		
		gsl_linalg_QR_decomp(	&mat_2t_t_view_tp1_t,
								vect_t );
		gsl_triangle_matrix(	&mat_2t_t_view_0);

		//y-1
		gsl_vector_memcpy( 	&mat_2t_t_view_vector,
							t_0 );
		gsl_vector_sub ( &mat_2t_t_view_vector,
						 t_s_0[i] );
						 
		for ( unsigned int j = 0; j < _size_x; ++ j )
		{
			mat_2t_t_view_vector.data[j * mat_2t_t_view_vector.stride] = 0;
		}
		
		gsl_linalg_QR_decomp(	&mat_2t_t_view_tp1_t,
								vect_t );
		gsl_triangle_matrix(	&mat_2t_t_view_0);

	}
	
	gsl_matrix_memcpy ( sqrt_q_0, &mat_2t_t_view_0 );
	gsl_matrix_scale ( sqrt_q_0, 1.0 / sqrt( nb_signals ) );
	
	for ( unsigned int j = 0; j < _size_x; ++ j )
	{
		for ( unsigned int i = _size_x; i < _size_t; ++ i )
		{
			sqrt_q_0->data[ j * sqrt_q_0->tda + i ] = 0;
		}
	}

}

tkalman_nc_fusion :: ~tkalman_nc_fusion()
{
	free();
	initialize();
}

void tkalman_nc_fusion :: initialize()
{
	mat_2t_x = 0;
	mat_2t_y = 0;
	vect_x = 0;
	vect_y = 0;
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
	mat_2t_t = 0;
	vect_t = 0;
	
}

void tkalman_nc_fusion :: free()
{
	if ( mat_2t_x )
		gsl_matrix_free( mat_2t_x );
	if ( mat_2t_y ) 
		gsl_matrix_free ( mat_2t_y );
	if ( vect_x )
		gsl_vector_free ( vect_x );
	if ( vect_y )
		gsl_vector_free ( vect_y );
	if ( mat_2t_t )
		gsl_matrix_free ( mat_2t_t );
	if ( vect_t )
		gsl_vector_free ( vect_t );
}

