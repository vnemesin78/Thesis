/**@file tkalman_em_initialization.hpp
 * @author Valérian Némesin
 * @brief
 * Ce fichier contient les prototypes des fonctions nécessaires à l'initialisation automatique de l'algorithme EM.
 */
#ifndef _TKALMAN_EM_INITIALIZATION_HPP_
	#define _TKALMAN_EM_INITIALIZATION_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	
	/**@fn void tkalman_EM_initialization(gsl_vector * x_0,
										  gsl_matrix * sqrt_p_0,
										  gsl_matrix * sqrt_q_yy,
										  const gsl_vector * const * observations,
										  const unsigned int n,
										  gsl_vector * vect_x)
	 * @param x_0 : espérance de l'état initial pour l'EM
	 * @param sqrt_p_0 : racine de la matrice de covariance de l'état initial pour l'EM
	 * @param sqrt_q_xx : racine de la matrice de covariance du bruit de process
	 * @param sqrt_q_yy : racine de la matrice de covariance du bruit de mesure
	 * @param[in] observation : observations
	 * @param[in] n : nombre d'observations
	 * @param vect_x : vecteur de taille x alloué
	 * @brief
	 * Cette fonction estime l'espérance de l'état initial, sa matrice de covariance et la matrice de covariance du bruit de mesure.
	 * @warning x et y doivent avoir même dimension
	 */
	void tkalman_EM_initialization(gsl_vector * x_0,
								   gsl_matrix * sqrt_p_0,
								   gsl_matrix * sqrt_q_yy,
								   const gsl_vector * const * observations,
								   const unsigned int n,
								   gsl_vector * vect_x);
	
	
#endif
