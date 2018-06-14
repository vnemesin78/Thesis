#include "tkalman_nc_filter.hpp"
#include <iostream>
using namespace std;
tkalman_nc_filter :: tkalman_nc_filter(const gsl_vector * t0,
									   const gsl_matrix * sqrt_q0,
									   const gsl_matrix * f,
									   const gsl_matrix * sqrt_q,
									   unsigned int size_x,
									   unsigned int n_max) throw(exception &)
: tkalman_nc_filter_base( t0, sqrt_q0, size_x, n_max )
{
	tkalman_nc_filter :: initialize();
	//Validité des arguments d'entrée
	if ( !(f && sqrt_q) )
		throw (invalid_argument("f or sqrt_q are NULL.\n"));
	//Dim
	if ( (sqrt_q->size1 != sqrt_q->size2) || (f->size1 != f->size2) )
		throw (invalid_argument("f or sqrt_q are not square matrix.\n"));
		
	if ( (_size_t != sqrt_q->size1) || (f->size1 != _size_t) )
		throw (invalid_argument("dim(q0), dim(Q), dim(F) or dim(t0) are different.\n"));
	
	
	//Vérif de la dim. de y;
	
	_f = f;
	_sqrt_q = sqrt_q;
	
	
	//Création des vues
	tkalman_nc_filter :: create_views();
	//Création des objets
	tkalman_nc_filter :: create_object();
}

void tkalman_nc_filter :: setup(const gsl_vector * t0,
								const gsl_matrix * sqrt_q0,
								const gsl_matrix * f,
								const gsl_matrix * sqrt_q,
								unsigned int size_x,
								unsigned int n_max) throw(exception &)
{
	
	
	tkalman_nc_filter_base :: setup ( t0, sqrt_q0, size_x, n_max);
	tkalman_nc_filter :: free();
	tkalman_nc_filter :: initialize();
	//Validité des arguments d'entrée
	if ( !(f && sqrt_q) )
		throw (invalid_argument("f or sqrt_q are NULL.\n"));
	//Dim
	if ( (sqrt_q->size1 != sqrt_q->size2) || (f->size1 != f->size2) )
		throw (invalid_argument("f or sqrt_q are not square matrix.\n"));
		
	if ( (_size_t != sqrt_q->size1) || (f->size1 != _size_t) )
		throw (invalid_argument("dim(q0), dim(Q), dim(F) or dim(t0) are different.\n"));
	
	
	//Vérif de la dim. de y;
	
	_f = f;
	_sqrt_q = sqrt_q;
	
	
	//Création des vues
	tkalman_nc_filter :: create_views();
	//Création des objets
	tkalman_nc_filter :: create_object();
	
}

tkalman_nc_filter :: ~tkalman_nc_filter()
{
	tkalman_nc_filter :: free();
	tkalman_nc_filter :: initialize();
}

void tkalman_nc_filter :: filter(const gsl_vector * const * observations,
							     unsigned int n) throw(exception &)
{
	if (n > _n_max)
	{
		throw(invalid_argument("n > n_max.\n"));
	}
	
	//Prédiction 0
	gsl_vector_memcpy(_x_p[0], _t0);
	gsl_matrix_memcpy(_sqrt_p_p[0], _sqrt_q0);
	
	filtering->compute_filtering_0(_x_f[0],
								   _sqrt_p_f[0],
								   _innovation[0],
								   _sqrt_s[0],
								   _t0,
								   _sqrt_q0,
								   observations[0]);
	
	prediction->compute_prediction_1(_x_p[1],
									 _sqrt_p_p[1],
									 _x_f[0],
								     _sqrt_p_f[0]);
	
	for (unsigned int i = 1; i < n; ++i)
	{
		filtering->compute_filtering(_x_f[i],
									 _sqrt_p_f[i],
									 _innovation[i],
									 _sqrt_s[i],
									 _x_p[i],
									 _sqrt_p_p[i],
									 observations[i],
									 observations[i - 1]);
		
		prediction->compute_prediction(_x_p[i + 1],
									   _sqrt_p_p[i + 1],
									   _x_f[i],
									   _sqrt_p_f[i],
									   observations[i - 1]);
	}
	
	
	
	
	
	
}

void tkalman_nc_filter :: smooth(const gsl_vector * const * observations,
								 unsigned int n) throw(exception &)
{
	try
	{
		filter(observations, n);
	}
	catch(exception & e)
	{
		throw(e);
	}
	gsl_vector_memcpy(_x_s[n], _x_p[n]);
	gsl_matrix_memcpy(_sqrt_p_s[n], _sqrt_p_p[n]);
	
	for (unsigned int i = n - 1; i > 0; --i)
	{
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
	
void tkalman_nc_filter :: predict(const gsl_vector * const * observations,
								  unsigned int n,
								  unsigned int p) throw(exception &)
{

	try
	{
		filter(observations, n);
	}
	catch(exception & e)
	{
		throw(e);
	}
	gsl_vector_set_zero(vect_t);
	gsl_matrix_set_zero(mat_tt);
	
	//Prédiction 0
	gsl_matrix_memcpy(_sqrt_q_pp[0], _sqrt_q0);
	gsl_vector_memcpy(_t_pp[0], _t0);
	
	//Prédiction 1
	prediction->compute_prediction_(_t_pp[1],
									_sqrt_q_pp[1],
									_x_p[0],
									_sqrt_p_p[0]);
	
	//Cas merdiques
	for (unsigned int i = 2; i <= n; ++i)
	{

		gsl_matrix_memcpy(&mat_tt_view_xx, _sqrt_p_p[i - 1]);
		gsl_vector_memcpy(&vect_t_view_x, _x_p[i - 1]);
		
		prediction->compute_prediction_(_t_pp[i],
										_sqrt_q_pp[i],
										vect_t,
										mat_tt);
	}

	for (unsigned j = 1; j < p; ++j)
	{
		for (unsigned int i = n; i > j + 1; --i)
		{

			
			prediction->compute_prediction_(_t_pp[i],
											_sqrt_q_pp[i],
											_t_pp[i - 1],
											_sqrt_q_pp[i - 1]);
		}
	}
	
	
	
	
	
}

void tkalman_nc_filter :: compute_sqrt_sums ( gsl_matrix * sq_sum_corr_tx,
											  gsl_matrix * sq_sum_corr_ty,
											  const gsl_vector * const * observations,
											  unsigned int n ) throw(exception &)
{
	gsl_matrix tx, ty;
	smooth(observations, n);
	
	//Calcul des sommes
	sums->compute_sqrt_sum_corr_tx_(&tx,
									x_s(),
									observations,
									sqrt_p_f(),
									sqrt_p_s(),
									c_s(),
									n);

	sums->compute_sqrt_sum_corr_ty(&ty,
								   x_s(),
								   observations,
								   sqrt_p_s(),
								   n);
	
	//Recopie des matrices sommes
	gsl_matrix_memcpy( sq_sum_corr_tx, &tx );
	gsl_matrix_memcpy( sq_sum_corr_ty, &ty );
	
}

void tkalman_nc_filter :: initialize()
{	
    _f = 0;
    _sqrt_q = 0;
    sums = 0;
}

void tkalman_nc_filter :: create_object() throw(exception &)
{
	try
	{
		prediction = new tkalman_nc_prediction(_f, _sqrt_q, _size_x);
	}
	catch(exception & e)
	{
		throw(e);
	}
	
	try
	{
		filtering = new tkalman_nc_filtering(&_f_yt, &_sqrt_q_yy);
	}
	catch(exception & e)
	{
		throw(e);
	}
	
	try
	{
		smoothing = new tkalman_nc_smoothing(&_f_xt, &_sqrt_q_xx);
	}
	catch(exception & e)
	{
		throw(e);
	}	
	
	try
	{
		sums = new tkalman_nc_sums( &_f_xt, &_sqrt_q_xx );
	}
	catch(exception & e)
	{
		throw(e);
	}	
	
	
	
}

void tkalman_nc_filter ::  create_views()
{
	
	//F
	{
		{
			gsl_matrix_const_view view = gsl_matrix_const_submatrix(_f, 0, 0, _size_x, _size_t);
			_f_xt = view.matrix;
		}
		{
			gsl_matrix_const_view view = gsl_matrix_const_submatrix(_f, 0, 0, _size_x, _size_x);
			_f_xx = view.matrix;
		}
		{
			gsl_matrix_const_view view = gsl_matrix_const_submatrix(_f, 0, _size_x, _size_x, _size_y);
			_f_xy = view.matrix;
		}
		{
			gsl_matrix_const_view view = gsl_matrix_const_submatrix(_f, _size_x, 0, _size_y, _size_t);
			_f_yt = view.matrix;
		}
		{
			gsl_matrix_const_view view = gsl_matrix_const_submatrix(_f, _size_x, 0, _size_y, _size_x);
			_f_yx = view.matrix;
		}
		{
			gsl_matrix_const_view view = gsl_matrix_const_submatrix(_f, _size_x, _size_x, _size_y, _size_y);
			_f_yy = view.matrix;
		}
	}
	
	//Q
	{
		{
			gsl_matrix_const_view view = gsl_matrix_const_submatrix(_sqrt_q, 0, 0, _size_x, _size_x);
			_sqrt_q_xx = view.matrix;
		}
		{
			gsl_matrix_const_view view = gsl_matrix_const_submatrix(_sqrt_q, _size_x, _size_x, _size_y, _size_y);
			_sqrt_q_yy = view.matrix;
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

void tkalman_nc_filter :: free()
{
	if ( sums )
		delete sums;
}

