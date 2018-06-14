/**@file tkalman_simulation.hpp
 * @author Valérian Némesin
 * @brief
 * Ce fichier contient le prototype des fonctions de simulation. Une fonction simule des données FKC et une des données non-FKC.
 * 
 */
#ifndef _TKALMAN_SIMULATION_HPP_
	#define _TKALMAN_SIMULATION_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_blas.h>
	#include <gsl/gsl_rng.h>
	#include <gsl/gsl_randist.h>
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
							gsl_rng * r);
							
							
							
							
							
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
							   gsl_rng * r);
#endif
