#include "tkalman_em.hpp"
void tkalman_em :: initialize()
{
	tkalman_em :: initialize_params();
	tkalman_em :: initialize_objects();
	tkalman_em :: initialize_moments();
}
void tkalman_em :: initialize_params()
{
	_i = 0;
	_t_0 = NULL;
	_sqrt_q_0 = NULL;
	_f = NULL;
	_sqrt_q = NULL;
	_n = 0;
	_p = 0;
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
	f2_x_ = NULL;
	sqrt_q2_xx = NULL;
	q2_xy = NULL;
	sqrt_q_yy = NULL;
	_m = NULL;
}
void tkalman_em :: initialize_objects()
{
	constants = NULL;
	prediction = NULL;
	filtering = NULL;
	smoothing = NULL;
	equivalents = NULL;
	sums = NULL;
	arg_max = NULL;
}
void tkalman_em :: initialize_moments()
{
	_innovation = NULL;
	_sqrt_s = NULL;
	_x_p = NULL;
	_sqrt_p_p = NULL;
	_x_f = NULL;
	_sqrt_p_f = NULL;
	_x_s = NULL;
	_sqrt_p_s = NULL;
	_c_s = NULL;
}

int tkalman_em :: alloc()
{
	if (tkalman_em :: alloc_params())
	{
		return 1;
	}
	if (tkalman_em :: alloc_objects())
	{
		return 1;
	}
	if (tkalman_em :: alloc_moments())
	{
		return 1;
	}
	return 0;
}

int tkalman_em :: alloc_params()
{
	_t_0 = gsl_vector_alloc(_size_t);
	_sqrt_q_0 = gsl_matrix_alloc(_size_t, _size_t);
	_f = gsl_matrix_alloc(_size_t, _size_t);
	_sqrt_q = gsl_matrix_alloc(_size_t, _size_t);
	
	f2_x_ = gsl_matrix_alloc(_size_x, _size_t);
	
	sqrt_q2_xx = gsl_matrix_alloc(_size_x, _size_x);
	q2_xy = gsl_matrix_alloc(_size_x, _size_y);
	sqrt_q_yy = gsl_matrix_alloc(_size_y, _size_y);
	
	_m = gsl_matrix_alloc(_size_t, _size_t);
	gsl_matrix_set_identity(_m);//M = id
	return (tkalman_em :: check_params());
}

int tkalman_em :: alloc_moments()
{
	tkalman_expectation_ref_v2(_x_p, _n + 1, _size_x, _size_t);
	tkalman_expectation_ref_v2(_x_f, _n, _size_x, _size_t);
	tkalman_expectation_ref_v2(_x_s, _n + 1, _size_x, _size_t);
	tkalman_expectation_ref_v2(_innovation, _n, _size_y, _size_y);
	
	tkalman_covariance_ref_v2(_sqrt_p_p, _n + 1, _size_x, _size_x, _size_t, _size_t);
	tkalman_covariance_ref_v2(_sqrt_p_f, _n, _size_x, _size_x, _size_t, _size_t);
	tkalman_covariance_ref_v2(_sqrt_p_s, _n + 1, _size_x, _size_x, _size_t, _size_t);
	tkalman_covariance_ref_v2(_c_s, _n, _size_x, _size_x, _size_x, _size_t);
	tkalman_covariance_ref_v2(_sqrt_s, _n, _size_y, _size_y, _size_y, _size_y);
	
	
	return (tkalman_em :: check_moments());
}

int tkalman_em :: alloc_objects()
{
	create_param_views();
	//Constantes
	if (!constants)
		constants = new tkalman_constants(&f_x_,
										  &f_y_,
										  _sqrt_q);
	if (!constants)
		return 1;
	if (!(*constants))
		return 1;
	//Prédiction
	if (!prediction)
		prediction = new tkalman_prediction(f2_x_, sqrt_q2_xx, q2_xy);
	if (!prediction)
		return 1;
	if (!(*prediction))
		return 1;
	//Filtrage
	if (!filtering)
		filtering = new tkalman_filtering(&f_y_, sqrt_q_yy);
	if (!filtering)
		return 1;
	if (!(*filtering))
		return 1;
	//Lissage
	if (!smoothing)
		smoothing = new tkalman_smoothing(f2_x_,
										  sqrt_q2_xx);
	if (!smoothing)
		return 1;
	if (!(*smoothing))
		return 1;
	//Sommes
	if (!sums)
		sums = new tkalman_sums(f2_x_, sqrt_q2_xx);
	if (!sums)
		return 1;
	if (!(*sums))
		return 1;
	//Maximisation
	if (!arg_max)
		arg_max = new tkalman_argmax(sums->sqrt_c_00(),
									 sums->sqrt_c_01(),
									 sums->sqrt_c_11()
									 );
	if (!arg_max)
		return 1;
	if (!(*arg_max))
		return 1;
	//Equivalents
	if (!equivalents)
		equivalents = new tkalman_equivalents(_m, _size_x);
	if (!equivalents)
		return 1;
	if (!(*equivalents))
		return 1;
	return 0;
}

void tkalman_em :: free()
{
	tkalman_em :: free_params();
	tkalman_em :: free_moments();
	tkalman_em :: free_objects();
}

void tkalman_em :: free_params()
{
	if (_t_0)
		gsl_vector_free(_t_0);
	if (_sqrt_q_0)
		gsl_matrix_free(_sqrt_q_0);
	if (_f)
		gsl_matrix_free(_f);
	if (_sqrt_q)
		gsl_matrix_free(_sqrt_q);
	if (f2_x_)
		gsl_matrix_free(f2_x_);
	if (sqrt_q2_xx)
		gsl_matrix_free(sqrt_q2_xx);
	if (q2_xy)
		gsl_matrix_free(q2_xy);
	if (sqrt_q_yy)
		gsl_matrix_free(sqrt_q_yy);
	if (_m)
		gsl_matrix_free(_m);
}

void tkalman_em :: free_moments()
{
	tkalman_expectation_unref(_x_p, _n + 1);
	tkalman_expectation_unref(_x_f, _n);
	tkalman_expectation_unref(_x_s, _n + 1);
	tkalman_expectation_unref(_innovation, _n);
	tkalman_covariance_unref(_sqrt_p_p, _n + 1);
	tkalman_covariance_unref(_sqrt_p_f, _n);
	tkalman_covariance_unref(_sqrt_p_s, _n + 1);
	tkalman_covariance_unref(_c_s, _n);
	tkalman_covariance_unref(_sqrt_s, _n);
	
	
}

void tkalman_em :: free_objects()
{
	if (constants)
		delete constants;
	if (prediction)
		delete prediction;
	if (filtering)
		delete filtering;
	if (smoothing)
		delete smoothing;
	if (sums)
		delete sums;
	if (arg_max)
		delete arg_max;
	if (equivalents)
		delete equivalents;
}

bool tkalman_em :: check_moments() const
{
	if (! (_x_p && _x_f && _x_s && _sqrt_p_p && _sqrt_p_f && _sqrt_p_s && _c_s && _innovation && _sqrt_s) )
		return 1;
	else
	{
		for (unsigned int i = 0;  i < _n; ++ i)
		{
			if (! (_x_p[i] && _x_f[i] && _x_s[i] && _sqrt_p_p[i] && _sqrt_p_f[i] && _sqrt_p_s[i] && _c_s[i] && _innovation[i] && _sqrt_s[i]) )
				return 1;
		}
		if (! (_x_p[_n] && _x_s[_n] && _sqrt_p_p[_n] && _sqrt_p_s[_n]) )
			return 1;
	}
	
	return 0;
}

bool tkalman_em :: check_params() const
{
	return ( !(_t_0 && _sqrt_q_0 && _f && _sqrt_q && f2_x_ && sqrt_q2_xx && q2_xy && _m) );
}

bool tkalman_em :: check_objects() const
{
	//Constantes
	if (!constants)
		return 1;
	if (!(*constants))
		return 1;
	
	//Prédiction
	if (!prediction)
		return 1;
	if (!(*prediction))
		return 1;
	//Filtrage
	if (!filtering)
		return 1;
	if (!(*filtering))
		return 1;
	if (!smoothing)
		return 1;
	if (!(*smoothing))
		return 1;
	
	//Sommes
	if (!sums)
		return 1;
	if (!(*sums))
		return 1;
		
	//Maximisation
	if (!arg_max)
		return 1;
	if (!(*arg_max))
		return 1;
	//Equivalents
	if (!equivalents)
		return 1;
	if (!(*equivalents))
		return 1;
	return 0;
	
}

bool tkalman_em :: operator !() const
{
	return ( tkalman_em :: check_moments() || tkalman_em :: check_params() || tkalman_em :: check_objects());
}

tkalman_em :: ~tkalman_em()
{
	tkalman_em :: free(); // free();
	tkalman_em :: initialize(); // initialize();	
}

void tkalman_em :: create_views()
{
	tkalman_em :: create_moment_views();
}

void tkalman_em :: create_moment_views()
{
	gsl_vector_view view;
	view = gsl_vector_subvector(_x_s[0], 0, _size_x);
	x_s_0 = view.vector;
	view = gsl_vector_subvector(_x_s[0], _size_x, _size_y);
	y_s_m1 = view.vector;
	
	view = gsl_vector_subvector(_x_p[0], 0, _size_x);
	x_p_0 = view.vector;
	view = gsl_vector_subvector(_x_p[0], _size_x, _size_y);
	y_p_m1 = view.vector;
	
	view = gsl_vector_subvector(_x_f[0], 0, _size_x);
	x_f_0 = view.vector;
	view = gsl_vector_subvector(_x_f[0], _size_x, _size_y);
	y_f_m1 = view.vector;
	
}

void tkalman_em :: create_param_views()
{
	{
		gsl_matrix_view view;
		//Vues sur F
		view = gsl_matrix_submatrix(_f, 0, 0, _size_x, _size_t);
		f_x_ = view.matrix;
		view = gsl_matrix_submatrix(_f, _size_x, 0, _size_y, _size_t);
		f_y_ = view.matrix;
	
		//Vues sur M
		view = gsl_matrix_submatrix(_m, 0, 0, _size_y, _size_t);
		m_x_ = view.matrix;
	}
	{
		gsl_vector_view view;
		view = gsl_vector_subvector(_t_0, 0, _size_x);
		x_0 = view.vector;
		view = gsl_vector_subvector(_t_0, _size_x, _size_y);
		y_m1 = view.vector;
	}
}

int tkalman_em :: setup(const gsl_vector * t_0,
						const gsl_matrix * sqrt_q_0,
					    const gsl_matrix * f,
					    const gsl_matrix * sqrt_q,
					    unsigned int size_x,
						unsigned int n,
					    unsigned int p)
{
	unsigned int size_y;
				 
	size_y = sqrt_q_0->size2 - size_x;
	
	//Reset de l'objet
	if (size_x != _size_x || size_y != _size_y || _n != n || tkalman_em :: operator !())
	{
		tkalman_em :: free();
		tkalman_em :: initialize();
		_size_x = size_x;
		_size_y = size_y;
		_size_t = size_x + size_y;
		_n = n;
		if ((tkalman_em :: alloc()) )
			return 1;
		tkalman_em :: create_views();
	}
	_p = p;
	//Construction des paramètres
	gsl_vector_memcpy(_t_0, t_0);
	gsl_matrix_memcpy(_sqrt_q, sqrt_q);
	gsl_matrix_memcpy(_sqrt_q_0, sqrt_q_0);
	gsl_matrix_memcpy(_f, f);
	return 0;
}
tkalman_em :: tkalman_em(const gsl_vector * t_0,
						 const gsl_matrix * sqrt_q_0,
					     const gsl_matrix * f,
					     const gsl_matrix * sqrt_q,
					     unsigned int size_x,
						 unsigned int n,
					     unsigned int p)
{
	tkalman_em :: initialize();
	tkalman_em :: setup(t_0, sqrt_q_0, f, sqrt_q, size_x, n, p);
}

void tkalman_em :: smooth(const gsl_vector * const * observations,
						  unsigned int n)
{
	if (n > _n)
	{
		tkalman_em :: free_moments();
		_n = n;
		tkalman_em :: initialize_moments();
		tkalman_em :: alloc_moments();
	}

	for (_i = 0; _i < _p; ++ _i)
	{
		

		//Filtrage-Lissage
		tkalman_em :: filter(observations, n);

		//Sommes + Maximisation
		sums->compute_sums(_sqrt_p_f,
						   _x_s,
						   _sqrt_p_s,
						   _c_s,
						   observations,
						   n,
						   &x_s_0,
						   &y_s_m1,
						   _sqrt_p_f[0],
						   _c_s[0]);
		arg_max->compute_argmax(NULL, 
								NULL, 
								_f, 
								_sqrt_q, 
								n,
								NULL,
								NULL);				
	}
	
	//E(x) = E(y)
	tkalman_em :: constraint();
	tkalman_em :: filter(observations, n);
	
}

void tkalman_em :: filter(const gsl_vector * const * observations,
						  unsigned int n)
{
	if (n > _n)
	{
		tkalman_em :: free_moments();
		_n = n;
		tkalman_em :: initialize_moments();
		tkalman_em :: alloc_moments();
	}

	
	//Calcul des constantes
	constants->get_constants(f2_x_, sqrt_q2_xx, q2_xy, sqrt_q_yy);
	
	//Première prédiction
	gsl_vector_memcpy(_x_p[0], _t_0);
	gsl_matrix_memcpy(_sqrt_p_p[0], _sqrt_q_0);
	//Passe-avant
	//Premier filtrage
	filtering->do_filtering_0(_x_f[0],
							  _sqrt_p_f[0],
							  _innovation[0],
							  _sqrt_s[0],
							  _x_p[0],
							  &x_p_0,
							  _sqrt_p_p[0],
							  observations[0],
							  &y_p_m1);

	//Prediction 2
	prediction->do_prediction_1(_x_p[1],
								_sqrt_p_p[1],
								&x_f_0,
								&y_f_m1,
								_sqrt_p_f[0],
								observations[0]);
	for (unsigned int i = 1; i < n; ++ i)	
	{
		//filtrage i
		filtering->do_filtering(_x_f[i],
							    _sqrt_p_f[i],
								_innovation[i],
							    _sqrt_s[i],
							    _x_p[i],
								_sqrt_p_p[i],
								observations[i],
								observations[i - 1]);
		//Prédiction i + 1
		prediction->do_prediction(_x_p[i + 1],
								  _sqrt_p_p[i + 1],
								  _x_f[i],
								  _sqrt_p_f[i],
								  observations[i - 1],
								  observations[i]);
	}

	//Lissage
	//Dernier lissage
	gsl_vector_memcpy(_x_s[n], _x_p[n]);

	gsl_matrix_memcpy(_sqrt_p_s[n], _sqrt_p_p[n]);

	for (unsigned int i = n - 1; i > 0; -- i)	
	{
		
		smoothing->do_smoothing(_x_s[i],
								_sqrt_p_s[i],
								_c_s[i],
								_x_f[i],
								_sqrt_p_f[i],
								_x_p[i + 1],
								_sqrt_p_p[i + 1],
								_x_s[i + 1],
								_sqrt_p_s[i + 1]);
	}
	smoothing->do_smoothing_0(_x_s[0],
							  _sqrt_p_s[0],
							  _c_s[0],
							  _x_f[0],
							  _sqrt_p_f[0],
							  _x_p[1],
							  _sqrt_p_p[1],
							  _x_s[1],
							  _sqrt_p_s[1]);
}

void tkalman_em :: constraint()
{
	gsl_matrix_memcpy(&m_x_, &f_y_); 
	equivalents->setup(_m,
					   _size_x);
	equivalents->compute_equivalent_f(_f);
	equivalents->compute_equivalent_sqrt_q(_sqrt_q);
	equivalents->compute_equivalent_sqrt_q(_sqrt_q_0);
	equivalents->compute_equivalent_x(&x_0, &y_m1);
	
}

void tkalman_em :: get_q_0(gsl_matrix * q_0) const
{
	gsl_blas_dgemm(CblasTrans,
				   CblasNoTrans,
				   1.0,
				   _sqrt_q_0,
				   _sqrt_q_0,
				   0.0,
				   q_0);	
}

void tkalman_em :: get_q(gsl_matrix * q) const
{
	gsl_blas_dgemm(CblasTrans,
				   CblasNoTrans,
				   1.0,
				   _sqrt_q,
				   _sqrt_q,
				   0.0,
				   q);	
}
void  tkalman_em :: set_params(const gsl_vector * t_0,
							   const gsl_matrix * sqrt_q_0,
							   const gsl_matrix * f,
							   const gsl_matrix * sqrt_q)
{
	gsl_vector_memcpy(_t_0, t_0);
	gsl_matrix_memcpy(_sqrt_q, sqrt_q);
	gsl_matrix_memcpy(_sqrt_q_0, sqrt_q_0);
	gsl_matrix_memcpy(_f, f);
	
}
