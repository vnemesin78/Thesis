#ifndef KALMAN_ROBUST_EM_HPP_INCLUDED
#define KALMAN_ROBUST_EM_HPP_INCLUDED
    #include "tkalman_robust_filter.hpp"
	#include "tkalman_em_algorithm.hpp"
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"
	class kalman_robust_em : public tkalman_robust_filter
	{
		public:
			/**@fn kalman_robust_em :: kalman_robust_em (const gsl_vector * x0,
															  const gsl_matrix * sqrt_p0,
															  const gsl_matrix * f,
															  const gsl_matrix * sqrt_q,
															  unsigned int n = 0,
															  unsigned int p = 0,
															  bool data = false);
			 * @param[in] x0 : Espérance de l'état initial
			 * @param[in] sqrt_p0 : Matrice de covariance de l'état initial (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
			 * @param[in] f : Matrice d'évolution
			 * @param[in] sqrt_q : Matrice de covariance (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
			 * @param[in] n : Nombre d'observations (0 par défaut)
			 * @param[in] p : Nombre d'itérations de l'EM (0 par défaut et dans ce cas, cela équivaut à un filtrage simple)
			 * @param[in] data : Booléen (True = suivi de l'EM et stockage des paramètres et de la vraisemblance à chaque itération, False = Pas de suivi de l'EM)
			 * @brief
			 * Constructeur
			*/
			kalman_robust_em (const gsl_vector * x0,
                              const gsl_matrix * sqrt_p0,
                              const gsl_matrix * f,
                              const gsl_matrix * sqrt_q,
                              unsigned int n = 0,
                              unsigned int p = 0,
                              bool data = false);

			/**@fn virtual int kalman_robust_em  :: setup(const gsl_vector * x0,
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
			 * Cette méthode réinitialise l'objet.
			 */
			virtual int setup(const gsl_vector * x0,
							  const gsl_matrix * p0,
							  const gsl_matrix * f,
							  const gsl_matrix * q,
							  unsigned int n = 0,
							  unsigned int p = 0,
							  bool data = false);


			/**@fn virtual void kalman_robust_em   :: filter(const gsl_vector * const * observations,
															  unsigned int nb_observations) = 0
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette méthode effectue le filtrage des données par le filtre de Kalman Triple non supervisé.
			 */
			virtual void filter(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn virtual void kalman_robust_em :: smooth(const gsl_vector * const * observations,
															  unsigned int nb_observations)
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette méthode effectue le lissage des données par le filtre de Kalman Triple non supervisé.
			 */
			virtual void smooth(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn kalman_robust_em  :: ~kalman_robust_em()
			 * @brief
			 * Destructeur
			 *
			 */
			~kalman_robust_em(void);

			/**@fn bool kalman_robust_em :: operator!() const
			 * @return
			 - 0 si l'objet est valide
			 - 1 sinon
			 * @brief
			 Cette méthode teste la validité de chaque attribut.
			**/
			virtual bool operator!() const;

			//Accesseurs
			/**@fn const gsl_vector * const * kalman_robust_em  :: get_x0_est()
			 * @return
			 Liste des x0 estimés
			 * @brief
			 Cette fonction renvoie les différents x0 estimés par l'algorithme EM.
			**/
			inline const gsl_vector * const * get_x0_est() const
			{
                return _x0_est;
			}

			/**@fn const gsl_matrix * const * kalman_robust_em :: get_sqrt_p0_est()
			 * @return
			 Liste des p0 estimés
			 * @brief
			 Cette fonction renvoie les différents p0 estimés par l'algorithme EM.
			**/
			inline  const gsl_matrix * const * get_sqrt_p0_est() const
			{
                return _p0_est;
			}

			/**@fn const gsl_matrix * const * kalman_robust_em :: get_f_est()
			 * @return
			 Liste des F estimés
			 * @brief
			 Cette fonction renvoie les différents F estimés par l'algorithme EM.
			**/
			inline const gsl_matrix * const * get_f_est() const
			{
                return _f_est;
			}

			/**@fn const gsl_matrix * const * kalman_robust_em :: get_sqrt_q_est()
			 * @return
			 Liste des Q estimés
			 * @brief
			 Cette fonction renvoie les différents Q estimés par l'algorithme EM.
			**/
			inline const gsl_matrix * const * get_sqrt_q_est() const
			{
                return _q_est;
			}
			/**@fn inline const double * kalman_robust_em :: get_log_likelihood()
			 * @return
			 Valeurs de la log-vraisemblance au cours des itérations de l'EM
			 * @brief
			 Cette fonction renvoie les valeurs de la log vraisemblance au cours des itérations de l'EM.
			**/
			inline const double * get_log_likelihood() const
			{
				return _log_likelihood;
			}

			/**@fn inline unsigned int kalman_robust_em :: get_nb_iter()
			 * @return
			 Nombre d'itérations de l'algorithme EM
			 @brief
			 Cette fonction renvoie le nombre d'itérations de l'algorithme EM.
			**/
			inline unsigned int get_nb_iter() const
			{
				return nb_iter;
			}
		protected:
			/**@fn void kalman_robust_em :: do_em_algorithm(const gsl_vector * const * observations);
			 * @param[in] observations : observations
			 * @brief
			 * Cette méthode estime les paramètres du filtre de Kalman couple à partir d'un jeu d'observations.
			 */
			void do_em_algorithm(const gsl_vector * const * observations);

			/**@fn void kalman_robust_em :: estimate_parameters(const gsl_vector * const * observations,
																   unsigned int i);
			 * @param[in] observations : observations
			 * @param[in] i : itération de l'EM
			 * @brief
			 * Cette méthode ré-estime les paramètres du filtre de Kalman couple à partir d'un jeu d'observations.
			 */
			void estimate_parameters(const gsl_vector * const * observations,
									 unsigned int i);

			/**@fn void kalman_robust_em :: follow(unsigned int i);
			 * @param i : itération de l'EM
			 * @brief
			 * Cette méthode enregistre les paramètres de l'itération i de l'algorithme EM.
			 **/
			void follow(unsigned int i);


			/**@fn void kalman_robust_em :: initialize();
			 * @brief
			 Cette méthode met tous les attributs à zéro.
			**/
			void initialize();

			/**@fn void kalman_robust_em :: initialize_tmp();
			 * @brief
			 Cette méthode met tous les temporaires à zéros
			**/
			void initialize_tmp();

			/**@fn int kalman_robust_em :: alloc(unsigned int p, bool data);
			 * @param p : nombre d'itérations de l'EM
			 * @param data : booléen
			 * @return
			  - 0 si l'allocation des attributs de l'objet s'est bien déroulée
			  - 1 sinon
			  *@brief
			  Cette méthode alloue les attributs de l'objet.
			 */
			int alloc(unsigned int p, bool data);

			/**@fn int kalman_robust_em :: alloc_tmp();
			 * @return
			  - 0 si l'allocation des temporaires de l'objet s'est bien déroulée
			  - 1 sinon
			  *@brief
			  Cette méthode alloue les temporaires de l'objet.
			 */
			int alloc_tmp();

			/**@fn int kalman_robust_em :: alloc_data()
			 * @brief
			 * Cette méthode alloue les données
			 *
			 */
			int alloc_data();

			/**@fn void kalman_robust_em :: free();
			 * @brief
			 Cette méthode libère la mémoire utilisée par les attributs de l'objet.
			**/
			void free();

			/**@fn void kalman_robust_em :: free_tmp();
			 * @brief
			 Cette méthode libère la mémoire utilisée par les temporaires de l'objet.
			**/
			void free_tmp();


			/**@fn void kalman_robust_em :: free_data();
			 *  @brief
			 Cette méthode libère la mémoire utilisée par les données de suivi de l'EM.
			**/
			void free_data();

			/**@fn bool kalman_robust_em :: check_tmp()
			 * @return
			 - 0 si les tmp ont été bien alloués
			 - 1 sinon
			 * @brief
			 Cette méthode teste la bonne allocation des variables temporaires utilisées par l'objet.
			**/
			bool check_tmp() const;

			/**@fn bool kalman_robust_em :: check_data()
			 * @return
			 - 0 si les données de suivi de l'EM ont été bien alloués
			 - 1 sinon
			 * @brief
			 Cette méthode teste la bonne allocation des données de suivi de l'EM.
			**/
			bool check_data() const;

			//Nombre d'itération de l'EM
				unsigned int nb_iter;


			//Suivi de l'EM
				double * _log_likelihood;
				gsl_vector ** _x0_est;
				gsl_matrix ** _p0_est;
				gsl_matrix ** _f_est;
				gsl_matrix ** _q_est;
            //tmp du suivi
                    gsl_permutation * perm_y_1;

			//Tmp de l'EM
				//Sommes
				gsl_matrix * mat_4x4x_1;
					gsl_matrix mat_4x4x_1_view_00;
					gsl_matrix mat_4x4x_1_view_11;
					gsl_matrix mat_4x4x_1_view_02;
					gsl_matrix mat_4x4x_1_view_20;
					gsl_matrix mat_4x4x_1_view_32;
					gsl_matrix mat_4x4x_1_view_33;
					gsl_matrix mat_4x4x_1_view_corr;
				gsl_matrix * mat_4xp12x_1;
					gsl_matrix mat_4xp12x_1_view_00;
					gsl_matrix mat_4xp12x_1_view_01;
					gsl_matrix mat_4xp12x_1_view_11;
					gsl_matrix mat_4xp12x_1_view_2xp12x;
					gsl_matrix mat_2xp12x_view_0;
					gsl_vector mat_2xp12x_view_10;
					gsl_vector mat_2xp12x_view_11;
				gsl_matrix * mat_2tp1t_1;
					gsl_matrix mat_2tp1t_1_view_00;
					gsl_matrix mat_2tp1t_1_view_01;
					gsl_matrix mat_2tp1t_1_view_11;
					gsl_matrix mat_2tp1t_1_view_tp1t;
					gsl_matrix mat_tp1t_view_xx;
					gsl_vector mat_tp1t_view_10;
					gsl_vector mat_tp1t_view_11;
				gsl_vector * vect_4x_1;
	};


#endif // kalman_robust_em_HPP_INCLUDED
