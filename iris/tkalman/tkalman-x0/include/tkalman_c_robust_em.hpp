#ifndef TKALMAN_C_ROBUST_EM_HPP_INCLUDED
    #define TKALMAN_C_ROBUST_EM_HPP_INCLUDED
  #include "tkalman_robust_em.hpp"
    /**@class tkalman_c_robust_em
    *@brief
    Cette classe g�re l'algorithme EM du filtre de Kalman couple. Elle ajoute une contrainte sur la matrice F :
    F = Fxx Fxy
        I   0
    *@warning
    * Dans cette classe, il est n�cessaire que les dimensions de x et y soient les m�mes.
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
			 * @param[in] x0 : Esp�rance de l'�tat initial
			 * @param[in] p0 : Matrice de covariance de l'�tat initial (remplac�e par sa d�composition de Cholesky dans certaines des classes filles)
			 * @param[in] f : Matrice d'�volution
			 * @param[in] q : Matrice de covariance (remplac�e par sa d�composition de Cholesky dans certaines des classes filles)
			 * @param[in] n : Nombre d'observations (0 par d�faut)
			 * @param[in] p : Nombre d'it�rations de l'EM (0 par d�faut et dans ce cas, cela �quivaut � un filtrage simple)
			 * @param[in] data : Bool�en (True = suivi de l'EM et stockage des param�tres et de la vraisemblance � chaque it�ration, False = Pas de suivi de l'EM)
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
			 Cette m�thode effectue le filtrage des donn�es par le filtre de Kalman Triple non supervis�.
			 */
			virtual void filter(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn virtual void tkalman_c_orignal_em :: smooth(const gsl_vector * const * observations,
															  unsigned int nb_observations)
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette m�thode effectue le lissage des donn�es par le filtre de Kalman Triple non supervis�.
			 */
			virtual void smooth(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn bool tkalman_c_orignal_em :: operator!()
			 * @return
			 - 0 si l'objet est valide
			 - 1 sinon
			 * @brief
			 Cette m�thode teste la validit� de chaque attribut.
			**/
			virtual bool operator!() const;

        protected :
            /**@fn void tkalman_c_orignal_em  :: do_em_algorithm(const gsl_vector * const * observations);
			 * @param[in] observations : observations
			 * @brief
			 * Cette m�thode estime les param�tres du filtre de Kalman couple � partir d'un jeu d'observations.
			 */
			void do_em_algorithm(const gsl_vector * const * observations);

    };



#endif // TKALMAN_C_ROBUST_EM_HPP_INCLUDED
