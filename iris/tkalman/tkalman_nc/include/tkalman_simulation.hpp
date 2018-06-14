
#ifndef _TKALMAN_SIMULATION_HPP_
#define _TKALMAN_SIMULATION_HPP_
	#include "tkalman_nc_em.hpp"
	#include <gsl/gsl_rng.h>
	#include <gsl/gsl_randist.h>

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
	void do_tkalman_simulation(	gsl_vector ** x,
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
								gsl_rng * r);
								 
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
	void do_no_tkalman_simulation(	gsl_vector ** y,
									const gsl_vector * const * x,
									const gsl_matrix * sqrt_q_yy,
									const unsigned int n,
									gsl_vector * vect,
									gsl_rng * r);
									
									
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
	void do_tkalman_simulation_markov(	gsl_vector ** x,
										gsl_vector ** y,
										const gsl_vector * t0,
										const gsl_matrix * sqrt_q0,
										const gsl_matrix * const * f,
										const gsl_matrix * const * sqrt_q,
										const unsigned int * process,
										const unsigned int n,
										gsl_vector * vect_t_1,
										gsl_vector * vect_t_1_view_x,
										gsl_vector * vect_t_1_view_y,
										gsl_vector * vect_t_2,
										gsl_rng * r);
									

#endif
