#include "tkalman_nc_em.hpp"

/**@fn tkalman_nc_em :: tkalman_nc_em(const gsl_vector * t0,
									  const gsl_matrix * sqrt_q0,
									  const gsl_matrix * f,
									  const gsl_matrix * sqrt_q,
									  unsigned int size_x,
									  unsigned int n_max,
									  gsl_vector * x_mask = NULL,
									  gsl_vector * y_mask = NULL);
 * @param[in] t0 : \hat{t}_0, espérance de l'état initial.
 * @param[in] sqrt_q0: [Q_0]^{\frac{1}{2}}, racine de la matrice de covariance de l'état initial.
 * @param[in] f : F, matrice d'évolution
 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
 * @param[in] size_x : dimension de x
 * @param[in] n_max : nombre de moments alloués
 * @param[in] x_mask : Masque sur la matrice Fxt
 * si x_mask(i) == 0, alors Fxt(:, i) == 0
 * @param[in] y_mask : Masque sur la matrice Fyt
 * si y_mask(i) == 0, alors Fyt(:, i) == 0
 * @brief
 * Constructeur de la classe @class tkalman_nc_em
**/
tkalman_nc_em :: tkalman_nc_em(	const gsl_vector * t0,
								const gsl_matrix * sqrt_q0,
								const gsl_matrix * f,
								const gsl_matrix * sqrt_q,
								unsigned int size_x,
								unsigned int n_max,
								unsigned int nb_signal_max,
								const gsl_vector * x_mask,
								const gsl_vector * y_mask,
								bool estimate_initial_state) throw (exception &)
: tkalman_nc_em_base(	t0,
						sqrt_q0,
						f,
						sqrt_q,
						size_x,
						n_max,
						nb_signal_max,
						x_mask,
						y_mask,
						estimate_initial_state )
{
	tkalman_nc_em :: initialize();
	//Validité des arguments d'entrée
	if ( !(t0 && sqrt_q0 && f && sqrt_q && size_x && n_max && nb_signal_max) )
		throw (invalid_argument("t0, sqrt_q0, f, sqrt_q, size_x, p_max or n_max are NULL.\n"));
	//Dim
	if ( (sqrt_q0->size1 != sqrt_q0->size2) || (sqrt_q->size1 != sqrt_q->size2) || (f->size1 != f->size2) )
		throw (invalid_argument("sqrt_q0, f or sqrt_q are not square matrix.\n"));
		
	if ( (sqrt_q0->size1 != sqrt_q->size1) || (sqrt_q->size1 != f->size1) || (f->size1 != t0->size ) )
		throw (invalid_argument("dim(q0), dim(Q), dim(F) or dim(t0) are different.\n"));
	
	try
	{
		alloc();
	}
	catch (exception & except)
	{
		throw(except);
	}
	//Recopie des paramètres
	{
		gsl_matrix_memcpy(_f, f);
		gsl_matrix_memcpy(_sqrt_q, sqrt_q);
	}
	//Création des vues
		create_views();
	//Objets
	try
	{
		create_object();
	}
	catch (exception & except)
	{
		throw(except);
	}
	
}

void tkalman_nc_em :: set_params( 	const gsl_vector * t0,
									 const gsl_matrix * sqrt_q0,
									 const gsl_matrix * f,
									 const gsl_matrix * sqrt_q )
{
	
		gsl_vector_memcpy(_t0, t0);
		gsl_matrix_memcpy(_sqrt_q0, sqrt_q0);
		gsl_matrix_memcpy(_f, f);
		gsl_matrix_memcpy(_sqrt_q, sqrt_q);
	
	
}


void tkalman_nc_em :: setup (	const gsl_vector * t0,
								const gsl_matrix * sqrt_q0,
								const gsl_matrix * f,
								const gsl_matrix * sqrt_q,
								unsigned int size_x,
								unsigned int n_max,
								unsigned int nb_signal_max,
								const gsl_vector * x_mask,
								const gsl_vector * y_mask,
								bool estimate_initial_state) throw (exception &)
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
	tkalman_nc_em :: free();
	tkalman_nc_em :: initialize();
	
	//Validité des arguments d'entrée
	if ( !(t0 && sqrt_q0 && f && sqrt_q && size_x && n_max && nb_signal_max) )
		throw (invalid_argument("t0, sqrt_q0, f, sqrt_q, size_x, p_max or n_max are NULL.\n"));
	//Dim
	if ( (sqrt_q0->size1 != sqrt_q0->size2) || (sqrt_q->size1 != sqrt_q->size2) || (f->size1 != f->size2) )
		throw (invalid_argument("sqrt_q0, f or sqrt_q are not square matrix.\n"));
		
	if ( (sqrt_q0->size1 != sqrt_q->size1) || (sqrt_q->size1 != f->size1) || (f->size1 != t0->size ) )
		throw (invalid_argument("dim(q0), dim(Q), dim(F) or dim(t0) are different.\n"));
	
	try
	{
		alloc();
	}
	catch (exception & except)
	{
		throw(except);
	}
	//Recopie des paramètres
	{
		gsl_matrix_memcpy(_f, f);
		gsl_matrix_memcpy(_sqrt_q, sqrt_q);
	}
	//Création des vues
		create_views();
	//Objets
	try
	{
		create_object();
	}
	catch (exception & except)
	{
		throw(except);
	}
}

tkalman_nc_em :: ~tkalman_nc_em()
{
	tkalman_nc_em :: free();
	tkalman_nc_em :: initialize();
}


void tkalman_nc_em :: filter(	const gsl_vector * const * observations,
								unsigned int n) throw(exception &)
{
	try
	{
		filter_obj->filter(observations, n);
	}
	catch(exception & e)
	{
		throw (e);
	}
}

void tkalman_nc_em :: smooth (	const gsl_vector * const * observations,
								unsigned int n) throw(exception &)
{
	try
	{
		filter_obj->smooth(observations, n);
	}
	catch(exception & e)
	{
		throw (e);
	}
}



void tkalman_nc_em :: predict(	const gsl_vector * const * observations,
								unsigned int n,
								unsigned int p) throw(exception &)
{
	try
	{
		filter_obj->predict(observations, n, p);
	}
	catch(exception & e)
	{
		throw (e);
	}
}

void  tkalman_nc_em :: learn_parameters(	const  gsl_vector * const * const * observations,
												const unsigned int * n,
												unsigned int p,
												unsigned int nb_iter) throw(exception &)
{	
	unsigned int sum_n = 0;
	
	//Somme des n pour l'algo
	for ( unsigned int i = 0; i < p; ++ i )
	{
		sum_n += n[i];
	}
	
	
	for (unsigned int i = 0; i < nb_iter; ++i)
	{
		for ( unsigned int j = 0; j < p; ++ j )
		{  
			try
			{
				filter_obj->compute_sqrt_sums ( sqrt_sums_tx[j],
												sqrt_sums_ty[j],
												observations[j],
												n[j]);
			}
			catch(exception & e)
			{
				throw (e);
			}
			//Recopie des éléments initiaux
			if ( _estimate_initial_state )
			{
				gsl_vector_memcpy( t_0s[j], filter_obj->x_s()[0] );
				gsl_matrix_memcpy( sqrt_q_0s[j], filter_obj->sqrt_p_s()[0] );

				
				
			}
		}
		//Calcul des sommes
		gsl_matrix sq_sum_corr_tx,
				   sq_sum_corr_ty;
		//Fusion
		fusion->fusion ( &sq_sum_corr_tx,
						 &sq_sum_corr_ty,
						 sqrt_sums_tx,
						 sqrt_sums_ty,
						 p );
		//Paramètres optimaux
		
		f_tools_x->maximize(&_sqrt_q_xx,
							tmp_f_xt_t,
							&sq_sum_corr_tx,
							sum_n); 
		gsl_matrix_transpose_memcpy(&_f_xt, tmp_f_xt_t);
		
		f_tools_y->maximize(&_sqrt_q_yy,
							tmp_f_yt_t,
							&sq_sum_corr_ty,
							sum_n); 
		gsl_matrix_transpose_memcpy(&_f_yt, tmp_f_yt_t);
		
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

void tkalman_nc_em :: initialize()
{
	//Objets
	filter_obj = 0;
				
	//Sommes
	sqrt_sums_tx = 0;
	sqrt_sums_ty = 0;
			
	_f = 0;
	_sqrt_q = 0;
	
}

void tkalman_nc_em :: alloc() throw(exception &)
{
	if (!_sqrt_q)
	{
		try
		{
			_sqrt_q = gsl_matrix_alloc(_size_t, _size_t);
		}
		catch(exception & e)
		{
			throw (e);
		}
	}
	
	if (!_f)
	{
		try
		{
			_f = gsl_matrix_alloc(_size_t, _size_t);
		}
		catch(exception & e)
		{
			throw (e);
		}
	}
	


	if (!sqrt_sums_tx)
	{
		sqrt_sums_tx = new gsl_matrix*[_nb_signal_max];
		for ( unsigned int i = 0; i < _nb_signal_max; ++ i)
		{
			sqrt_sums_tx[i] = gsl_matrix_alloc(_size_t + _size_x, _size_t + _size_x);
		}
	}
	if (!sqrt_sums_ty)
	{
		sqrt_sums_ty = new gsl_matrix*[_nb_signal_max];
		for ( unsigned int i = 0; i < _nb_signal_max; ++ i)
		{
			sqrt_sums_ty[i] = gsl_matrix_alloc(_size_t + _size_y, _size_t + _size_y);
		}
	}
	
}

void tkalman_nc_em :: free()
{
	if (filter_obj)
		delete filter_obj;

	if (_sqrt_q)
	{
		gsl_matrix_free(_sqrt_q);
	}
	
	if (_f)
	{
		gsl_matrix_free(_f);
	}
	
	if (sqrt_sums_tx)
	{
		for ( unsigned int i = 0; i < _nb_signal_max; ++ i)
		{
			gsl_matrix_free(sqrt_sums_tx[i]);
		}
		delete[] sqrt_sums_tx;
	}
	if (sqrt_sums_ty)
	{
		for ( unsigned int i = 0; i < _nb_signal_max; ++ i)
		{
			gsl_matrix_free(sqrt_sums_ty[i]);
		}
		delete[] sqrt_sums_ty;
	}
}

void tkalman_nc_em :: create_views()
{
	gsl_matrix_view v;
	
	v = gsl_matrix_submatrix(_f, 0, 0, _size_x, _size_t);
	_f_xt = v.matrix;
	
	
	v = gsl_matrix_submatrix(_f, _size_x, 0, _size_y, _size_t);
	_f_yt = v.matrix;
	
	
	v = gsl_matrix_submatrix(_sqrt_q, 0, 0, _size_x, _size_x);
	_sqrt_q_xx = v.matrix;
	
	v = gsl_matrix_submatrix(_sqrt_q, _size_x, _size_x, _size_y, _size_y);
	_sqrt_q_yy = v.matrix;
}

void tkalman_nc_em :: create_object() throw(exception &)
{
	try
	{
		filter_obj = new tkalman_nc_filter(_t0, _sqrt_q0, _f, _sqrt_q, _size_x, _n_max);
	}
	catch (exception & e)
	{
		throw(e);
	}
}

