#include "tkalman_simulation.hpp"
#include <iostream>
using namespace std;
/**@fn void tkalman_simulation(gsl_vector ** x,
							   gsl_vector ** y,
							   const gsl_vector * x0,
							   const gsl_matrix * sqrt_p0,
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
 * @param x0 : espérance de l'état initial
 * @param sqrt_p0 : racine de la matrice de covariance de l'état initial
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
void tkalman_simulation(gsl_vector ** x,
						gsl_vector ** y,
						const gsl_vector * x0,
						const gsl_matrix * sqrt_p0,
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
	size_x = x0->size;
	size_y = f->size1 - size_x;
	size_t = size_x + size_y;
	//Tirage de x0
	for ( unsigned int i = 0; i < size_x; ++ i )
	{
		vect_t_1_view_x->data[i * vect_t_1_view_x->stride] = gsl_ran_ugaussian (r);
	}

	gsl_vector_memcpy(x[0], x0);
	
	gsl_blas_dgemv(CblasTrans,
				   1.0,
				   sqrt_p0,
				   vect_t_1_view_x,
				   1.0,
				   x[0]);

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


/**@fn void no_tkalman_simulation(gsl_vector ** y,
								   const gsl_vector * const * x,
								   const gsl_matrix * fyx,
								   const gsl_matrix * sqrt_q_yy,
								   const unsigned int n,
								   gsl_vector * vect_y,
								   gsl_rng * r)

 * @param y : observations simulées (Doit être alloué pour n éléments)
 * @param x : états cachés simulés (Doit être alloué pour (n + 1) éléments)
 * @param Fyx : relation entre x et y
 * @param sqrt_q_yy : racine de la matrice de covariance du bruit de mesure
 * @param n : nombre d'observations
 * @param vect_y : vecteur de taille y préalloué
 * @param r : générateur de nombres aléatoires
 * @brief
 * Cette fonction génère des observations à partir d'observations connuees.
**/
void no_tkalman_simulation(gsl_vector ** y,
						   const gsl_vector * const * x,
						   const gsl_matrix * fyx,
						   const gsl_matrix * sqrt_q_yy,
						   const unsigned int n,
						   gsl_vector * vect_y,
						   gsl_rng * r)
{
	unsigned int size_y;
	
	size_y = y[0]->size;
	
	for ( unsigned int i = 0; i < n; ++ i)
	{
		
		//Calcul de F T
		gsl_blas_dgemv(CblasNoTrans,
					   1.0,
					   fyx,
					   x[i],
					   0.0,
					   y[i]);
		
		for (unsigned int j = 0; j < size_y; ++ j)
		{
			vect_y->data[j * vect_y->stride] = gsl_ran_ugaussian (r)	;
		}		
		
		gsl_blas_dgemv(CblasTrans,
					   1.0,
					   sqrt_q_yy,
					   vect_y,
					   1.0,
					   y[i]);	
	}
}





