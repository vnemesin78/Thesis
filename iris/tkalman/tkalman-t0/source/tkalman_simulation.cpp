#include "tkalman_simulation.hpp"
/**@fn void do_tkalman_simulation_2(gsl_vector ** x,
							     gsl_vector ** y,
							     const gsl_vector * t0,
							     const gsl_matrix * sqrt_q0,
							     const gsl_matrix * f,
							     const gsl_matrix * sqrt_q,
							     const unsigned int n,
							     gsl_vector * vect_t_1,
							     gsl_vector * vect_t_1_view_x,
							     gsl_vector * vect_t_1_view_y,
							     gsl_vector * vect_t_2,
							     gsl_rng * r)
 * @param x : états cachés simulés (Doit être alloué pour (n + 1) éléments)
 * @param y : observations simulées (Doit être alloué pour n éléments)
 * @param t0 : espérance de t_0
 * @param sqrt_q0 : racine de la matrice de covariance de t_0
 * @param f : matrice d'évolution
 * @param sqrt_q : racine de la matrice de covariance du bruit
 * @param n : nombre d'observations
 * @param vect_t_1 : vecteur de taille t prélloué
 * @param vect_t_1_view_x : vue sur le vecteur vect_t_1 (0 à n_x - 1)
 * @param vect_t_1_view_y : vue sur le vecteur vect_t_1 (n_x à n_t)
 * @param vect_t_2 : vecteur de taille t alloué
 * @param r : générateur de nombre (préalloué)
 * @brief
 * Cette fonction simule des données selon le modèle du filtre de Kalman couple.
**/
void do_tkalman_simulation_2(gsl_vector ** x,
						     gsl_vector ** y,
						     const gsl_vector * t0,
						     const gsl_matrix * sqrt_q0,
						     const gsl_matrix * f,
						     const gsl_matrix * sqrt_q,
						     const unsigned int n,
						     gsl_vector * vect_t_1,
						     gsl_vector * vect_t_1_view_x,
						     gsl_vector * vect_t_1_view_y,
						     gsl_vector * vect_t_2,
						     gsl_rng * r)
{
	unsigned int size_y, 
				 size_x,
				 size_t;
	size_x = x[0]->size;
	size_y = f->size1 - size_x;
	size_t = size_x + size_y;
	//Tirage de t0
	for ( unsigned int i = 0; i < size_t; ++ i )
	{
		vect_t_2->data[i * vect_t_2->stride] = gsl_ran_ugaussian (r);
	}
	gsl_vector_memcpy(vect_t_1, t0);
	gsl_blas_dgemv(CblasTrans,
				   1.0,
				   sqrt_q0,
				   vect_t_2,
				   1.0,
				   vect_t_1);
	 		   
	gsl_vector_memcpy(x[0],
					  vect_t_1_view_x);
	//Tirage des tn
	for ( unsigned int i = 1; i <= n; ++ i)
	{
		//Calcul de F T
		gsl_blas_dgemv(CblasNoTrans,
					   1.0,
					   f,
					   vect_t_1,
					   0.0,
					   vect_t_2);
		
		//Bruit
		for (unsigned int j = 0; j < size_t; ++ j)
		{
			vect_t_1->data[j * vect_t_1->stride] = gsl_ran_ugaussian (r);	
		}
		gsl_blas_dgemv(CblasTrans,
					   1.0,
					   sqrt_q,
					   vect_t_1,
					   1.0,
					   vect_t_2);
		gsl_vector_memcpy(vect_t_1, vect_t_2);
		gsl_vector_memcpy(x[i], vect_t_1_view_x);
		gsl_vector_memcpy(y[i - 1],  vect_t_1_view_y);
	}
}

/**@fn void do_no_tkalman_simulation_2(gsl_vector ** y,
									   const gsl_vector * const * x,
									   const gsl_matrix * sqrt_q_yy,
									   const unsigned int n,
									   gsl_vector * vect,
									   gsl_rng * r);
 * @param y : observations (x + b)
 * @param x : é"tats cachés
 * @param sqrt_q_yy : cov(b)
 * @param n : nombre d'observations
 * @param vect : ecteur de taille x alloué
 * @param r : gen. de nombre aléatoire
 * @brief
 * Simulateur de signaux bruité selon un bruit gaussien.
 */
void do_no_tkalman_simulation_2(gsl_vector ** y,
								const gsl_vector * const * x,
								const gsl_matrix * sqrt_q_yy,
								const unsigned int n,
								gsl_vector * vect,
						        gsl_rng * r)
{
	unsigned int size;
	size = x[0]->size;
	//Tirage des tn
	for ( unsigned int i = 0; i < n; ++ i)
	{
		gsl_vector_memcpy(y[i], x[i]);
		//Bruit
		for (unsigned int j = 0; j < size; ++ j)
		{
			vect->data[j * vect->stride] = gsl_ran_ugaussian (r);	
		}
		gsl_blas_dgemv(CblasTrans,
					   1.0,
					   sqrt_q_yy,
					   vect,
					   1.0,
					   y[i]);
	}
}















tkalman_simulation_2 :: tkalman_simulation_2(const gsl_vector * t0,
											 const gsl_matrix * sqrt_q0,
											 const gsl_matrix * f,
											 const gsl_matrix * sqrt_q,
											 gsl_rng * r,
											 unsigned int size_x,
											 unsigned int n)
{
	tkalman_simulation_2 :: initialize();
	tkalman_simulation_2 :: setup(t0,
								  sqrt_q0,
								  f,
								  sqrt_q,
								  r,
								  size_x,
								  n);
}

int tkalman_simulation_2 :: setup(const gsl_vector * t0,
								  const gsl_matrix * sqrt_q0,
								  const gsl_matrix * f,
								  const gsl_matrix * sqrt_q,
								  gsl_rng * r,
								  unsigned int size_x,
								  unsigned int n)
{
	tkalman_simulation_2 :: free();
	tkalman_simulation_2 :: initialize();
	
	_t0 = t0;
	_sqrt_q0 = sqrt_q0;
	_sqrt_q = sqrt_q;
	_f = f;
	rng = r;
	_n = n;
	_size_x = size_x;
	_size_t = t0->size;
	_size_y = _size_t - _size_x;
	
	
	if (tkalman_simulation_2 :: alloc())
	{
		tkalman_simulation_2 :: free();
		tkalman_simulation_2 :: initialize();
		return 1;
	}
	tkalman_simulation_2 :: create_tmp_view();
	return 0;
}



void tkalman_simulation_2 :: do_simulation()
{
	do_tkalman_simulation_2(_x,
						    _y,
						    _t0,
						    _sqrt_q0,
						    _f,
						    _sqrt_q,
						    _n,
						    vect_t_1,
						    &vect_t_1_view_x,
						    &vect_t_1_view_y,
						    vect_t_2,
						    rng);
}

void tkalman_simulation_2 :: get_instant_mean_square_error(gsl_vector ** i_mse,
														   const gsl_vector * const * x_est,
														   unsigned int n,
														   double scale)
{
	if ( ( n - 1 ) > _n)
		n = _n + 1;
	tkalman_get_instant_mean_square_error(i_mse,
										  _x,
										  x_est,
										  n,
										  scale,
										  &vect_t_1_view_y);
}

tkalman_simulation_2 :: ~tkalman_simulation_2()
{
	tkalman_simulation_2 :: free();
	tkalman_simulation_2 :: initialize();
}

bool tkalman_simulation_2 :: operator!() const
{
	return ((tkalman_simulation_2 :: check_moments() || tkalman_simulation_2 :: check_params() || tkalman_simulation_2 :: check_tmp() ) ); 
}


void tkalman_simulation_2 :: initialize()
{
	tkalman_simulation_2 :: initialize_moments();
	tkalman_simulation_2 :: initialize_params();
	tkalman_simulation_2 :: initialize_tmp();
}

void tkalman_simulation_2 :: initialize_params()
{
	_n = 0;
	_sqrt_q = 0;
	_sqrt_q0 = 0;
	_f = 0;
	_t0 = 0;
	_size_x = 0;
	_size_y = 0;
	_size_t = 0;
	rng = 0;
}
void tkalman_simulation_2 :: initialize_moments()
{
	_x = NULL;
	_y = NULL;
}
void tkalman_simulation_2 :: initialize_tmp()
{
	vect_t_1 = NULL;
	vect_t_2 = NULL;
}

void tkalman_simulation_2 :: free()
{
	 tkalman_simulation_2 :: free_moments();
	 tkalman_simulation_2 :: free_tmp(); 
}

void tkalman_simulation_2 :: free_moments()
{
	tkalman_expectation_unref(_x, _n + 1);
	tkalman_expectation_unref(_y, _n);
}

void tkalman_simulation_2 :: free_tmp()
{
	if (vect_t_1)
		gsl_vector_free(vect_t_1);
	if (vect_t_2)
		gsl_vector_free(vect_t_2);
	
}

int tkalman_simulation_2 :: alloc()
{
	return (tkalman_simulation_2 :: alloc_tmp() || tkalman_simulation_2 :: alloc_moments() );
}
int tkalman_simulation_2 :: alloc_tmp()
{
	if (!vect_t_1)
		vect_t_1 = gsl_vector_alloc(_size_t);
	
	if (!vect_t_2)
		vect_t_2 = gsl_vector_alloc(_size_t);
	return tkalman_simulation_2 :: check_tmp();
}
int tkalman_simulation_2 :: alloc_moments()
{
	tkalman_expectation_ref_v2(_x, _n + 1, _size_x, _size_x);
	tkalman_expectation_ref_v2(_y, _n, _size_y, _size_y);
	return tkalman_simulation_2 :: check_moments();
}

int tkalman_simulation_2 :: check_params() const
{
	return (! (_size_x && _size_y && _size_t && _t0 && _sqrt_q0 && _sqrt_q && _f && rng));
	
}
int tkalman_simulation_2 :: check_tmp() const
{
	return (! (vect_t_1 && vect_t_2) );
}
int tkalman_simulation_2 :: check_moments() const
{
	if (!_x || !_y)
		return 1;
	for (unsigned int i = 0; i <_n ; ++i)
	{
		if (! _x[i] || !_y[i])
			return 1;
	}
	if (!(_x)[_n])
		return 1;
	return 0;
}

void tkalman_simulation_2 :: create_tmp_view()
{
	gsl_vector_view view;
	view = gsl_vector_subvector(vect_t_1, 0, _size_x);
	vect_t_1_view_x = view.vector;
	
	view = gsl_vector_subvector(vect_t_1, _size_x, _size_y);
	vect_t_1_view_y = view.vector;
}
