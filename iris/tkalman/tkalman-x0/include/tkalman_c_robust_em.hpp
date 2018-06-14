#ifndef TKALMAN_C_ROBUST_EM_HPP_INCLUDED
    #define TKALMAN_C_ROBUST_EM_HPP_INCLUDED
  #include "tkalman_robust_em.hpp"
    /**@class tkalman_c_robust_em
    *@brief
    Cette classe gère l'algorithme EM du filtre de Kalman couple. Elle ajoute une contrainte sur la matrice F :
    F = Fxx Fxy
        I   0
    *@warning
    * Dans cette classe, il est nécessaire que les dimensions de x et y soient les mêmes.
    **/
    class tkalman_c_robust_em : public tkalman_robust_em
    {
        public:
			/**@fn tkalman_c_robust_em :: tkalman_x_original_em(const gsl_vector * x0,
															  const gsl_matrix * p0,
															  const gsl_matrix * f,
															  const gsl_matrix * q,
															  unsigned int n = 0,
															  unsigned int p = 0,
															  bool data = false);
			 * @param[in] x0 : Espérance de l'état initial
			 * @param[in] p0 : Matrice de covariance de l'état initial (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
			 * @param[in] f : Matrice d'évolution
			 * @param[in] q : Matrice de covariance (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
			 * @param[in] n : Nombre d'observations (0 par défaut)
			 * @param[in] p : Nombre d'itérations de l'EM (0 par défaut et dans ce cas, cela équivaut à un filtrage simple)
			 * @param[in] data : Booléen (True = suivi de l'EM et stockage des paramètres et de la vraisemblance à chaque itération, False = Pas de suivi de l'EM)
			 * @brief
			 * Constructeur
			*/
			tkalman_c_robust_em(const gsl_vector * x0,
							      const gsl_matrix * sqrt_p0,
								  const gsl_matrix * f,
								  const gsl_matrix * sqrt_q,
								  unsigned int n = 0,
								  unsigned int p = 0,
								  bool data = false);


			/**@fn virtual void tkalman_c_orignal_em  :: filter(const gsl_vector * const * observations,
															  unsigned int nb_observations) = 0
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette méthode effectue le filtrage des données par le filtre de Kalman Triple non supervisé.
			 */
			virtual void filter(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn virtual void tkalman_c_orignal_em :: smooth(const gsl_vector * const * observations,
															  unsigned int nb_observations)
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette méthode effectue le lissage des données par le filtre de Kalman Triple non supervisé.
			 */
			virtual void smooth(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn bool tkalman_c_orignal_em :: operator!()
			 * @return
			 - 0 si l'objet est valide
			 - 1 sinon
			 * @brief
			 Cette méthode teste la validité de chaque attribut.
			**/
			virtual bool operator!() const;

        protected :
            /**@fn void tkalman_c_orignal_em  :: do_em_algorithm(const gsl_vector * const * observations);
			 * @param[in] observations : observations
			 * @brief
			 * Cette méthode estime les paramètres du filtre de Kalman couple à partir d'un jeu d'observations.
			 */
			void do_em_algorithm(const gsl_vector * const * observations);

    };



#endif // TKALMAN_C_ROBUST_EM_HPP_INCLUDED
