#include "c_tracking.hpp"

c_tracking :: c_tracking ()
{
	initialize();
}

c_tracking :: c_tracking ( 	const gsl_vector * t_0,
							const gsl_matrix * sqrt_q_0,
							const gsl_matrix * f,
							const gsl_matrix * sqrt_q,
							unsigned int size_x,
							const double & err_sigma,
							ostream * _err_stream )
{
	initialize();
	if ( setup ( t_0, sqrt_q_0, f, sqrt_q, size_x, err_sigma, _err_stream ) )
	{
		throw (invalid_argument("Bad args."));
	}
}


int c_tracking :: setup (  	const gsl_vector * t_0,
							const gsl_matrix * sqrt_q_0,
							const gsl_matrix * f,
							const gsl_matrix * sqrt_q,
							unsigned int size_x,
							const double & err_sigma,
							ostream * _err_stream )
{
	free();
	initialize();
	err_stream = _err_stream;
	if ( ! t_0 || ! sqrt_q_0 || !f || ! size_x || err_sigma <= 0 )
	{
		if ( err_stream)
			*err_stream << "Error : Invalid arguments for int c_tracking :: setup (  const gsl_vector *, const gsl_matrix *, const gsl_matrix *, const gsl_matrix *, unsigned int, const double & );" << endl;
		return 1;
	}

	_size_t = t_0->size;
	_size_x = size_x;
	
	//Controle des dimensions
	if ( _size_x >= _size_t )
	{
		if ( err_stream)
			*err_stream << "Error : Bad arguments for int c_tracking :: setup (  const gsl_vector *, const gsl_matrix *, const gsl_matrix *, const gsl_matrix *, unsigned int, const double & );" << endl;
		return 1;
	}
	_size_y = _size_t - _size_x;
	if ( _size_t != sqrt_q_0->size1 ||  _size_t != sqrt_q_0->size2 || _size_t != sqrt_q->size1 ||  _size_t != sqrt_q->size2 || _size_t != f->size1 ||  _size_t != f->size2 )
	{
		if ( err_stream)
			*err_stream << "Error : Bad arguments for int c_tracking :: setup (  const gsl_vector *, const gsl_matrix *, const gsl_matrix *, const gsl_matrix *, unsigned int, const double & );" << endl;
		return 1;
	}
	
	_err_sigma = err_sigma;
	
	//Alloc mémoire
	
	_x_max = gsl_vector_alloc( _size_x );
	_x_min = gsl_vector_alloc( _size_x );
	
	_t_0 = gsl_vector_alloc( _size_t );
	
	_y = gsl_vector_alloc( _size_y );
	y = gsl_vector_alloc( _size_y );

	_sqrt_q_0 = gsl_matrix_alloc( _size_t,
								  _size_t);
	_f = gsl_matrix_alloc( _size_t,
						   _size_t);
	
	_sqrt_q = gsl_matrix_alloc( _size_t,
								_size_t);
	
	x_p = gsl_vector_alloc( _size_x );
	x_f = gsl_vector_alloc( _size_x );
	innovation = gsl_vector_alloc( _size_y );
	t_f = gsl_vector_alloc( _size_t );
	
	sqrt_q_f = gsl_matrix_alloc( _size_t,
								 _size_t);
	sqrt_p_f = gsl_matrix_alloc( _size_x,
								 _size_x );
						
	sqrt_p_p = gsl_matrix_alloc( _size_x,
								 _size_x );
								
	sqrt_s = gsl_matrix_alloc( _size_y,
							   _size_y);
	
	t_p = gsl_vector_alloc( _size_t );
	vect_t = gsl_vector_alloc( _size_t );
	mat_tt = gsl_matrix_alloc( _size_t, _size_t);
	sqrt_q_p = gsl_matrix_alloc( _size_t, _size_t );
	
	
	
	
	//Recopie des données

	gsl_vector_memcpy( _t_0,
					   t_0);		   
	gsl_matrix_memcpy( _sqrt_q_0,
					   sqrt_q_0);
					   
	gsl_matrix_memcpy( _sqrt_q,
					   sqrt_q);
					   
	gsl_matrix_memcpy( _f,
					   f);
	
	
					   
	//Vues
	//Fyt	   
	{
		gsl_matrix_view view;
		view = gsl_matrix_submatrix( _f,
									 _size_x,
									 0,
									 _size_y,
									 _size_t );
		f_yt = view.matrix;
		
	}
	
	//Qyy
	{
		gsl_matrix_view view;
		view = gsl_matrix_submatrix( _sqrt_q,
									 _size_x,
									 _size_x,
									 _size_y,
									 _size_y );
		sqrt_q_yy = view.matrix;
		
	}

	//Vues
	//Fyt	   
	{
		gsl_vector_view view;
		vect_t_view_x = gsl_vector_subvector( vect_t, 0, _size_x ).vector;

	}
	//Fyt	   
	{
		gsl_matrix_view view;
		mat_tt_view_xx = gsl_matrix_submatrix( mat_tt, 0, 0, _size_x, _size_x ).matrix;

	}


	prediction = new tkalman_nc_prediction ( 	_f,
												_sqrt_q,
												_size_x
											);

	filtering = new tkalman_nc_filtering (&f_yt, &sqrt_q_yy);
	
	vect_one = gsl_vector_alloc( _size_x );
	gsl_vector_set_all ( vect_one, 1 );
	return 0;
}


int c_tracking :: setup ( 	api_parameters & params,
							ostream * _err_stream,
							const char * n_space,
							const char * t_0_name,
							const char * sqrt_q_0_name,
							const char * f_name,
							const char * sqrt_q_name,
							const char * size_x_name,
							const char * n_name )
{
	err_stream = _err_stream;
	double n;
	unsigned int size_x;
	gsl_vector t_0;
	gsl_matrix sqrt_q_0,
			   f,
			   sqrt_q;
			   
	stringstream oss;
	//T0
	oss << n_space << "::" << t_0_name;
	if ( api_get_vector(	params,
							oss.str().c_str(),
							&t_0,
							err_stream	) )
		return -1;
	oss.str("");
	//sqrt_Q0
	oss << n_space << "::" << sqrt_q_0_name;
	if ( api_get_matrix(	params,
							oss.str().c_str(),
							&sqrt_q_0,
							err_stream	) )
		return -1;
	oss.str("");
	//F
	oss << n_space << "::" << f_name;
	if ( api_get_matrix(	params,
							oss.str().c_str(),
							&f,
							err_stream	) )
		return -1;
	oss.str("");
	
	oss << n_space << "::" << sqrt_q_name;
	if ( api_get_matrix(	params,
							oss.str().c_str(),
							&sqrt_q,
							err_stream	) )
		return -1;
	oss.str("");
	
	oss << n_space << "::"  << n_name;
	if ( api_get_double(	params,
							oss.str().c_str(),
							&n,
							err_stream	) )
		return -1;
	oss.str("");
	
	oss << n_space << "::"  << size_x_name;
	if ( api_get_positive_integer(	params,
									oss.str().c_str(),
									&size_x,
									err_stream	) )
		return -1;
	oss.str("");
	
	return setup ( 	&t_0, 
					&sqrt_q_0, 
					&f, 
					&sqrt_q, 
					size_x, 
					n,
					err_stream );
}

void c_tracking :: predict_next_position( 	const double * x,
											unsigned int frame_id,
											unsigned int stride )
{
	gsl_vector_memcpy(_y, 
					   y);
	for ( unsigned int i = 0; i < _size_y; ++ i )
	{
		y->data[i] = x[i * stride];
	}
	unsigned int d_id = frame_id - _frame_id;
	if ( _previous_ok )
	{

			//Filtrage
			filtering->compute_filtering( x_f,
										 sqrt_p_f,
										 innovation,
										 sqrt_s,
										 x_p,
										 sqrt_p_p,
										 y,
										 _y);
			
			
			prediction->compute_prediction ( x_p,
											sqrt_p_p,
											x_f,
											sqrt_p_f,
											y );
			
	}
	else
	{
		//Filtrage
		filtering->compute_filtering_0( t_f,
									    sqrt_q_f,
									    innovation,
									    sqrt_s,
									    _t_0,
									    _sqrt_q_0,
									    y);
		
		//Prédiction
		prediction->compute_prediction_1( x_p,
									     sqrt_p_p,
									     t_f,
									     sqrt_q_f);
									     
									     
									     
	}
	//Construction de vect t
	gsl_vector_set_zero (vect_t);
	gsl_vector_memcpy( &vect_t_view_x, x_p );
	
	//Construction de mat_tt
	gsl_matrix_set_zero(mat_tt);
	gsl_matrix_memcpy( &mat_tt_view_xx, sqrt_p_p );
	
	for ( unsigned int i = 1; i < d_id; ++ i )
	{
		
		
		prediction->compute_prediction_(	t_p,
											sqrt_q_p,
											vect_t,
											mat_tt );
		
		
		
		gsl_matrix_memcpy( mat_tt, sqrt_q_p );
		gsl_vector_memcpy( vect_t, t_p );
	}
	//Récup. des valeurs intéressantes
	gsl_vector_memcpy( x_p, &vect_t_view_x );
	gsl_matrix_memcpy( sqrt_p_p, &mat_tt_view_xx );		     



	//Calcul de la zone de confiance
	gsl_vector_memcpy( _x_min, x_p );
	gsl_vector_memcpy( _x_max, x_p );
	
	gsl_blas_dgemv( CblasNoTrans, 
					- _err_sigma,
					sqrt_p_p,
					vect_one,
					1.0,
					_x_min );
					
	gsl_blas_dgemv( CblasNoTrans, 
					_err_sigma,
					sqrt_p_p,
					vect_one,
					1.0,
					_x_max );
	
	//Correction (sqrt_p_p n'est pas la décompostion de Cholesky)
	for ( unsigned int i = 0; i < _size_x; ++ i )
	{
		if ( _x_min->data[i] > _x_max->data[i] )
		{
			double a;
			a = _x_min->data[i];
			_x_min->data[i] = _x_max->data[i];
			_x_max->data[i] = a;
		}
	}
	_previous_ok = true;
	_frame_id = frame_id;

}

void c_tracking :: reset_tracking ( )
{
	_previous_ok = false;
	//~ _frame_id = 0;
}

void c_tracking :: reset( )
{
	reset_tracking();
	_frame_id = 0;
}

c_tracking :: ~c_tracking()
{
	free();
	initialize();
}

void c_tracking :: free()
{
	if ( vect_one )
		gsl_vector_free ( vect_one );

	if ( _x_min )
		gsl_vector_free ( _x_min );
	if ( _x_max )
		gsl_vector_free ( _x_max  );
	
	if ( _t_0 )
	{
		gsl_vector_free ( _t_0 );
	}

	if ( _sqrt_q_0 )
	{
		gsl_matrix_free ( _sqrt_q_0 );
	}
	
	if ( _f )
	{
		gsl_matrix_free ( _f );
	}

	if ( _sqrt_q )
	{
		gsl_matrix_free ( _sqrt_q );
	}
	
	if ( x_p )
	{
		gsl_vector_free ( x_p );
	}	
	if ( x_f )
	{
		gsl_vector_free ( x_f );
	}
	if ( t_f )
	{
		gsl_vector_free ( t_f );
	}
	if ( innovation )
	{
		gsl_vector_free ( innovation );
	}
	if ( y )
	{
		gsl_vector_free ( y );
	}
	if ( _y )
	{
		gsl_vector_free ( _y );
	}
	if ( sqrt_p_f )
	{
		gsl_matrix_free ( sqrt_p_f );
	}
	if ( sqrt_p_p )
	{
		gsl_matrix_free ( sqrt_p_p );
	}
	if ( sqrt_s )
	{
		gsl_matrix_free ( sqrt_s );
	}
	if ( sqrt_q_f )
	{
		gsl_matrix_free ( sqrt_q_f );
	}
	if ( filtering )
	{
		delete filtering;
	}
	if ( prediction )
	{
		delete prediction;
	}
	if ( t_p )
		gsl_vector_free( t_p );
	if ( vect_t )
		gsl_vector_free( vect_t );
	if ( sqrt_q_p )
		gsl_matrix_free( sqrt_q_p );
	if ( mat_tt )
		gsl_matrix_free( mat_tt );
	
	
}

void c_tracking :: initialize()
{
	_x_min = 0;
	_x_max = 0;
	_previous_ok = 0;
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
	_t_0 = 0;
	_sqrt_q_0 = 0;
	_f = 0;
	_sqrt_q = 0;
	x_p = 0;
	x_f = 0;
	innovation = 0;
	t_f = 0;
	_y = 0;
	y = 0;
	sqrt_q_f = 0;
	sqrt_p_f = 0;
	sqrt_p_p = 0;
	sqrt_s = 0;
	_err_sigma = 0;
	filtering = 0;
	prediction = 0;
	err_stream = NULL;
	vect_one = 0;
	_frame_id =  0;
	t_p = 0;
	vect_t = 0;
	sqrt_q_p = 0;
	mat_tt = 0;
	
}
