#include "tkalman_nc_em_markov.hpp"
tkalman_nc_em_markov :: tkalman_nc_em_markov(	const gsl_vector * t0,
												const gsl_matrix * sqrt_q0,
												const gsl_matrix * f,
												const gsl_matrix * sqrt_q,
												unsigned int size_x,
												unsigned int n_max,
												unsigned int p_max,
												unsigned int nb_signal_max,
												const gsl_vector * x_mask,
												const gsl_vector * y_mask,
												bool estimate_initial_state ) throw (exception &)
: tkalman_nc_em_base(	t0,
						sqrt_q0,
						f,
						sqrt_q,
						size_x,
						n_max,
						nb_signal_max,
						x_mask,
						y_mask,
						estimate_initial_state)
{
	tkalman_nc_em_markov :: initialize();
	_p_max = p_max;
	tkalman_nc_em_markov :: alloc();
	//Recopie des paramètres
	for ( unsigned int i = 0; i < _p_max; ++ i)
	{
		gsl_matrix_memcpy(_f[i], f);
		gsl_matrix_memcpy(_sqrt_q[i], sqrt_q);
	}
	//Création des vues
		tkalman_nc_em_markov :: create_views();
	//Objets
	try
	{
		tkalman_nc_em_markov :: create_object();
	}
	catch (exception & except)
	{
		throw(except);
	}
}
void tkalman_nc_em_markov :: setup(	const gsl_vector * t0,
									const gsl_matrix * sqrt_q0,
									const gsl_matrix * f,
									const gsl_matrix * sqrt_q,
									unsigned int size_x,
									unsigned int n_max,
									unsigned int p_max,
									unsigned int nb_signal_max,
									const gsl_vector * x_mask,
									const gsl_vector * y_mask,
									bool estimate_initial_state ) throw (exception &)
{
	tkalman_nc_em_base :: setup(	t0,
									sqrt_q0,
									f,
									sqrt_q,
									size_x,
									n_max,
									nb_signal_max,
									x_mask,
									y_mask,
									estimate_initial_state);
	tkalman_nc_em_markov :: free();
	tkalman_nc_em_markov :: initialize();
	_p_max = p_max;
	tkalman_nc_em_markov :: alloc();
	//Recopie des paramètres
	for ( unsigned int i = 0; i < _p_max; ++ i)
	{
		gsl_matrix_memcpy(_f[i], f);
		gsl_matrix_memcpy(_sqrt_q[i], sqrt_q);
	}
	//Création des vues
		tkalman_nc_em_markov :: create_views();
	//Objets
	try
	{
		tkalman_nc_em_markov :: create_object();
	}
	catch (exception & except)
	{
		throw(except);
	}
}

tkalman_nc_em_markov :: ~tkalman_nc_em_markov()
{
	tkalman_nc_em_markov :: free();
	tkalman_nc_em_markov :: initialize();
}

void tkalman_nc_em_markov :: filter(	const gsl_vector * const * observations,
										const unsigned int * r,
										unsigned int n) throw(exception &)
{
		filter_obj->filter(	observations,
							r,
							n);
}

void tkalman_nc_em_markov :: smooth(	const gsl_vector * const * observations,
										const unsigned int * r,
										unsigned int n) throw(exception &)
{
		filter_obj->smooth(	observations,
							r,
							n);
}

void tkalman_nc_em_markov :: learn_parameters(	const  gsl_vector * const * const * observations,
												const unsigned int * const * r,
												const unsigned int * n,
												unsigned int p,
												unsigned int nb_iter) throw(exception &)
											  
{
	memset( sum_n, 0, sizeof(int) * _p_max);

	//Somme des n pour l'algo
	for ( unsigned int i = 0; i < p; ++ i )
	{
		for ( unsigned int j = 0; j < n[i]; ++ j )
		{
			sum_n[r[i][j]] ++;
		}
	}
	for (unsigned int i = 0; i < nb_iter; ++i)
	{

		for ( unsigned int j = 0; j < p; ++ j )
		{  
			try
			{

				filter_obj->compute_sqrt_sums ( sqrt_sums_tx + j * _p_max,
												sqrt_sums_ty + j * _p_max,
												r[j],
												observations[j],
												n[j]);
			}
			catch(exception & e)
			{
				throw (e);
			}
			if ( _estimate_initial_state )
			{
				gsl_vector_memcpy( t_0s[i], filter_obj->x_s()[0] );
				gsl_matrix_memcpy( sqrt_q_0s[i], filter_obj->sqrt_p_s()[0] );
			}
			
		}

		for ( unsigned int j = 0; j < _p_max; ++ j)
		{
			if (  sum_n[j] > 0 )
			{
			//Calcul des sommes
			gsl_matrix sq_sum_corr_tx,
					   sq_sum_corr_ty;
					   
					   
			//Fusion

			fusion->fusion ( &sq_sum_corr_tx,
							 &sq_sum_corr_ty,
							 sqrt_sums_tx + j,
							 sqrt_sums_ty + j,
							 p,
							 _p_max );

							 
			//Paramètres optimaux
			f_tools_x->maximize(_sqrt_q_xx + j,
								tmp_f_xt_t,
								&sq_sum_corr_tx,
								sum_n[j]); 
			gsl_matrix_transpose_memcpy(_f_xt + j, tmp_f_xt_t);
			
			
			f_tools_y->maximize(_sqrt_q_yy + j,
								tmp_f_yt_t,
								&sq_sum_corr_ty,
								sum_n[j]); 
			gsl_matrix_transpose_memcpy(_f_yt + j, tmp_f_yt_t);

			}
		}
		

		//Paramètres initiaux
		if ( _estimate_initial_state )
		{
			fusion->fusion_0 (_t0,
							  _sqrt_q0,
							  t_0s,
							  sqrt_q_0s,
							  p );
		}
		

		
	}
	
}

void tkalman_nc_em_markov :: initialize()
{
	filter_obj = 0;
	sqrt_sums_tx = 0;
	sqrt_sums_ty = 0;
	_f = 0;
	_f_xt = 0;
	_f_yt = 0;
	_sqrt_q = 0;
	_sqrt_q_xx = 0;
	_sqrt_q_yy = 0;
			
	sum_n = 0;
	_p_max = 0;
}


void tkalman_nc_em_markov :: alloc() throw(exception &)
{
	if (!_sqrt_q)
	{
		_sqrt_q = new gsl_matrix*[_p_max];
		for ( unsigned int i = 0; i < _p_max; ++ i)
			_sqrt_q[i] = gsl_matrix_alloc( _size_t, _size_t);
	}
	
	if (!_f)
	{
		_f = new gsl_matrix*[_p_max];
		for ( unsigned int i = 0; i < _p_max; ++ i)
			_f[i] = gsl_matrix_alloc( _size_t, _size_t);
	}
	
	if (!_f_xt)
		_f_xt = new gsl_matrix[_p_max];
	if (!_f_yt)
		_f_yt = new gsl_matrix[_p_max];

	if (!_sqrt_q_xx)
		_sqrt_q_xx = new gsl_matrix[_p_max];
	if (!_sqrt_q_yy)
		_sqrt_q_yy = new gsl_matrix[_p_max];

	if (! sum_n)
		 sum_n = new unsigned int[_p_max];

	if (!sqrt_sums_tx)
	{
		sqrt_sums_tx = new gsl_matrix*[_nb_signal_max * _p_max];
		for ( unsigned int i = 0; i < _nb_signal_max * _p_max; ++ i)
		{
			sqrt_sums_tx[i] = gsl_matrix_alloc(_size_t + _size_x, _size_t + _size_x);
		}
	}
	if (!sqrt_sums_ty)
	{
		sqrt_sums_ty = new gsl_matrix*[_nb_signal_max * _p_max];
		for ( unsigned int i = 0; i < _nb_signal_max * _p_max; ++ i)
		{
			sqrt_sums_ty[i] = gsl_matrix_alloc(_size_t + _size_y, _size_t + _size_y);
		}
	}
	
}


void tkalman_nc_em_markov :: free()
{
	if (_sqrt_q)
	{
		for ( unsigned int i = 0; i < _p_max; ++ i)
			gsl_matrix_free(_sqrt_q[i]);
		delete[] _sqrt_q;
	}
	
	if (_f)
	{
		for ( unsigned int i = 0; i < _p_max; ++ i)
			gsl_matrix_free(_f[i]);
		delete[] _f;
	}
	
	if (_f_xt)
		delete[] _f_xt;
	if (_f_yt)
		delete[] _f_yt;

	if (_sqrt_q_xx)
		delete[] _sqrt_q_xx;
	if (_sqrt_q_yy)
		delete[] _sqrt_q_yy;

	if (sum_n)
		delete[] sum_n;

	if (sqrt_sums_tx)
	{
		for ( unsigned int i = 0; i < _nb_signal_max * _p_max; ++ i)
		{
			gsl_matrix_free(sqrt_sums_tx[i]);
		}
		delete[] sqrt_sums_tx;
	}
	if (sqrt_sums_ty)
	{
		for ( unsigned int i = 0; i < _nb_signal_max * _p_max; ++ i)
		{
			gsl_matrix_free(sqrt_sums_ty[i]);
		}
		delete[] sqrt_sums_ty;
	}
	
}

void tkalman_nc_em_markov :: create_object() throw(exception &)
{
	try
	{
		filter_obj = new tkalman_nc_filter_markov(_t0, _sqrt_q0, _f, _sqrt_q, _size_x, _p_max, _n_max);
	}
	catch (exception & e)
	{
		throw(e);
	}
}

void tkalman_nc_em_markov :: create_views()
{
	for (unsigned int i = 0; i < _p_max; ++ i)
	{
		gsl_matrix_view v;
		
		v = gsl_matrix_submatrix(_f[i], 0, 0, _size_x, _size_t);
		_f_xt[i] = v.matrix;
		
		
		v = gsl_matrix_submatrix(_f[i], _size_x, 0, _size_y, _size_t);
		_f_yt[i] = v.matrix;
		
		
		v = gsl_matrix_submatrix(_sqrt_q[i], 0, 0, _size_x, _size_x);
		_sqrt_q_xx[i] = v.matrix;
		
		v = gsl_matrix_submatrix(_sqrt_q[i], _size_x, _size_x, _size_y, _size_y);
		_sqrt_q_yy[i] = v.matrix;
	}
	
	
}
