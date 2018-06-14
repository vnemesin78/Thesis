#include "tkalman_nc_filter_markov.hpp"
tkalman_nc_filter_markov :: tkalman_nc_filter_markov(	const gsl_vector * t0,
														const gsl_matrix * sqrt_q0,
														const gsl_matrix * const * f,
														const gsl_matrix * const * sqrt_q,
														unsigned int size_x,
														unsigned int p,
														unsigned int n_max) throw(exception &)
: tkalman_nc_filter_base( t0, sqrt_q0, size_x, n_max )
{
	tkalman_nc_filter_markov :: initialize();
	//Validité des arguments d'entrée
	if ( !(f && sqrt_q && p) )
		throw (invalid_argument("f or sqrt_q are NULL.\n"));
	//Dim
	
	//Vérif de la dim. de y;
	_p = p;
	_f = f;
	_sqrt_q = sqrt_q;
	
	//Création des vues
	tkalman_nc_filter_markov :: alloc();
	
	//Création des vues
	tkalman_nc_filter_markov :: create_views();
	//Création des objets
	tkalman_nc_filter_markov :: create_object();
}


void tkalman_nc_filter_markov :: setup(	const gsl_vector * t0,
										const gsl_matrix * sqrt_q0,
										const gsl_matrix * const * f,
										const gsl_matrix * const * sqrt_q,
										unsigned int size_x,
										unsigned int p,
										unsigned int n_max) throw(exception &)
{
	tkalman_nc_filter_base :: setup ( t0, sqrt_q0, size_x, n_max );
	tkalman_nc_filter_markov :: free();
	tkalman_nc_filter_markov :: initialize();
	//Validité des arguments d'entrée
	if ( !(f && sqrt_q && p) )
		throw (invalid_argument("f or sqrt_q are NULL.\n"));
	
	//Vérif de la dim. de y;
	_p = p;
	_f = f;
	_sqrt_q = sqrt_q;
	
	//Création des vues
	tkalman_nc_filter_markov :: alloc();
	
	//Création des vues
	tkalman_nc_filter_markov :: create_views();
	//Création des objets
	tkalman_nc_filter_markov :: create_object();
}

tkalman_nc_filter_markov :: ~tkalman_nc_filter_markov()
{
	tkalman_nc_filter_markov :: free();
	tkalman_nc_filter_markov :: initialize();
}

void tkalman_nc_filter_markov :: filter(	const gsl_vector * const * observations,
											const unsigned int * r,
											unsigned int n) throw(exception &)
{
	
	if (n > _n_max)
	{
		throw(invalid_argument("n > n_max.\n"));
	}
	
	//Prédiction 0
	gsl_vector_memcpy(_x_p[0], _t0);
	gsl_matrix_memcpy(_sqrt_p_p[0], _sqrt_q0);
	
	filtering->setup(_f_yt + r[0],
					 _sqrt_q_yy + r[0] ) ;
			
	filtering->compute_filtering_0(_x_f[0],
								   _sqrt_p_f[0],
								   _innovation[0],
								   _sqrt_s[0],
								   _t0,
								   _sqrt_q0,
								   observations[0]);
	
	prediction->setup(	_f[r[0]],
						_sqrt_q[r[0]],
						_size_x);
	prediction->compute_prediction_1(_x_p[1],
									 _sqrt_p_p[1],
									 _x_f[0],
								     _sqrt_p_f[0]);
	for (unsigned int i = 1; i < n; ++i)
	{
		
		filtering->setup(_f_yt + r[i],
						_sqrt_q_yy + r[i] ) ;
		filtering->compute_filtering(_x_f[i],
									 _sqrt_p_f[i],
									 _innovation[i],
									 _sqrt_s[i],
									 _x_p[i],
									 _sqrt_p_p[i],
									 observations[i],
									 observations[i - 1]);
		
		prediction->setup(	_f[r[i]],
							_sqrt_q[r[i]],
							_size_x);
		prediction->compute_prediction(_x_p[i + 1],
									   _sqrt_p_p[i + 1],
									   _x_f[i],
									   _sqrt_p_f[i],
									   observations[i - 1]);
	}


}


void tkalman_nc_filter_markov :: smooth(	const gsl_vector * const * observations,
											const unsigned int * r,
											unsigned int n) throw(exception &)
{
	try
	{
		filter(observations, r, n);
	}
	catch(exception & e)
	{
		throw(e);
	}

	gsl_vector_memcpy(_x_s[n], _x_p[n]);
	gsl_matrix_memcpy(_sqrt_p_s[n], _sqrt_p_p[n]);
	
	for (unsigned int i = n - 1; i > 0; --i)
	{

		smoothing->setup( _f_xt + r[i],
						  _sqrt_q_xx + r[i]);

		smoothing->compute_smoothing ( _x_s[i], 
									   _sqrt_p_s[i],
									   _c_s[i],
									   _x_f[i],
									   _sqrt_p_f[i],
									   _x_p[i + 1],
									   _sqrt_p_p[i + 1],
									   _x_s[i + 1], 
									   _sqrt_p_s[i + 1] );
	}
	smoothing->setup( _f_xt + r[0],
					  _sqrt_q_xx + r[0]);
	smoothing->compute_smoothing_0 ( _x_s[0], 
									 _sqrt_p_s[0],
									 _c_s[0],
									 _x_f[0],
									 _sqrt_p_f[0],
								     _x_p[1],
									 _sqrt_p_p[1],
									 _x_s[1], 
									 _sqrt_p_s[1] );

}

void tkalman_nc_filter_markov :: compute_sqrt_sums ( 	gsl_matrix ** sq_sum_corr_tx,
														gsl_matrix ** sq_sum_corr_ty,
														const unsigned int * r,
														const gsl_vector * const * observations,
														unsigned int n ) throw(exception &)
{
	for ( unsigned int i = 0; i < _p; ++ i)
	{
		gsl_matrix tx, ty;
		smooth(observations, r, n);
		//Calcul des sommes
		sums->compute_sqrt_sum_corr_tx_(&tx,
										x_s(),
										observations,
										sqrt_p_f(),
										sqrt_p_s(),
										c_s(),
										r,
										i,
										n);
		sums->compute_sqrt_sum_corr_ty(&ty,
									   x_s(),
									   observations,
									   sqrt_p_s(),
									   r,
									   i,
									   n);
		//Recopie des matrices sommes
		gsl_matrix_memcpy( sq_sum_corr_tx[i], &tx );
		gsl_matrix_memcpy( sq_sum_corr_ty[i], &ty );

	}
}


void tkalman_nc_filter_markov :: initialize()
{
	_f = 0;
	_f_xx = 0;
	_f_xy = 0;
	_f_xt = 0;
	_f_yx = 0;
	_f_yy = 0;
	_f_yt = 0;
	_sqrt_q = 0;
	_sqrt_q_xx = 0;
	_sqrt_q_yy = 0;
	sums = 0;
	_p = 0;
	
}

void tkalman_nc_filter_markov :: alloc() throw(exception &)
{
	_f_xx = new gsl_matrix[_p];
	_f_xy = new gsl_matrix[_p];
	_f_xt = new gsl_matrix[_p];
	_f_yx = new gsl_matrix[_p];
	_f_yy = new gsl_matrix[_p];
	_f_yt = new gsl_matrix[_p];
	
	_sqrt_q_xx = new gsl_matrix[_p];
	_sqrt_q_yy = new gsl_matrix[_p];
	
}

void tkalman_nc_filter_markov :: free()
{
	if (_f_xx)
		delete[] _f_xx;
		
	if (_f_xy)
		delete[] _f_xy;
		
	if (_f_xt)
		delete[] _f_xt;
		
	if (_f_yx)
		delete[] _f_yx;
	
	if (_f_yy)
		delete[] _f_yy;
		
	if (_f_yt)
		delete[] _f_yt;
		
	if (_sqrt_q_xx)
		delete[] _sqrt_q_xx;
		
	if (_sqrt_q_yy)
		delete[] _sqrt_q_yy;
}

void tkalman_nc_filter_markov :: create_object() throw(exception &)
{
	sums = new tkalman_nc_sums_markov(_f_xt, _sqrt_q_xx, _p);
	try
	{
		prediction = new tkalman_nc_prediction(_f[0], _sqrt_q[0], _size_x);
	}
	catch(exception & e)
	{
		throw(e);
	}
	
	try
	{
		filtering = new tkalman_nc_filtering(_f_yt, _sqrt_q_yy);
	}
	catch(exception & e)
	{
		throw(e);
	}
	
	try
	{
		smoothing = new tkalman_nc_smoothing(_f_xt, _sqrt_q_xx);
	}
	catch(exception & e)
	{
		throw(e);
	}	
	
	
}

void tkalman_nc_filter_markov :: create_views()
{
	for ( unsigned int i = 0; i < _p; ++ i)
	{
		//F
		{
			{
				gsl_matrix_const_view view = gsl_matrix_const_submatrix(_f[i], 0, 0, _size_x, _size_t);
				_f_xt[i] = view.matrix;
			}
			{
				gsl_matrix_const_view view = gsl_matrix_const_submatrix(_f[i], 0, 0, _size_x, _size_x);
				_f_xx[i] = view.matrix;
			}
			{
				gsl_matrix_const_view view = gsl_matrix_const_submatrix(_f[i], 0, _size_x, _size_x, _size_y);
				_f_xy[i] = view.matrix;
			}
			{
				gsl_matrix_const_view view = gsl_matrix_const_submatrix(_f[i], _size_x, 0, _size_y, _size_t);
				_f_yt[i] = view.matrix;
			}
			{
				gsl_matrix_const_view view = gsl_matrix_const_submatrix(_f[i], _size_x, 0, _size_y, _size_x);
				_f_yx[i] = view.matrix;
			}
			{
				gsl_matrix_const_view view = gsl_matrix_const_submatrix(_f[i], _size_x, _size_x, _size_y, _size_y);
				_f_yy[i] = view.matrix;
			}
		}
		
		//Q
		{
			{
				gsl_matrix_const_view view = gsl_matrix_const_submatrix(_sqrt_q[i], 0, 0, _size_x, _size_x);
				_sqrt_q_xx[i] = view.matrix;
			}
			{
				gsl_matrix_const_view view = gsl_matrix_const_submatrix(_sqrt_q[i], _size_x, _size_x, _size_y, _size_y);
				_sqrt_q_yy[i] = view.matrix;
			}
			
		}
	}
	//Tmp
	gsl_matrix_view view256;
	view256 = gsl_matrix_submatrix(mat_tt, 0, 0, _size_x, _size_x);
	mat_tt_view_xx = view256.matrix;
	
	gsl_vector_view view64;
	view64 = gsl_vector_subvector(vect_t, 0, _size_x);
	vect_t_view_x = view64.vector;
}


