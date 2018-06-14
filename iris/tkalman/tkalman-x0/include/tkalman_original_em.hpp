/**@file tkalman_original_em.hpp
 * @author Valérian Némesin
 * @brief
 * Ce fichier contient le prototype de la classe qui permet de créer des objets pour effectuer l'algorithme EM du filtre de Kalman triplet.
 */
#ifndef _TKALMAN_ORIGINAL_EM_HPP_
	#define _TKALMAN_ORIGINAL_EM_HPP_
	#include "tkalman_original_filter.hpp"
	#include "tkalman_em_algorithm.hpp"
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"
	
	/**@class
	 * @brief
	 * Cette classe permet de gérer un algorithme EM du filtre de Kalman triplet de Desbouvries.
	 * 
	 * 
	 **/
	class tkalman_original_em : public tkalman_original_filter
	{
		public:
			/**@fn tkalman_original_em :: tkalman_original_em(const gsl_vector * x0,
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
			tkalman_original_em(const gsl_vector * x0,
							    const gsl_matrix * p0,
								const gsl_matrix * f,
								const gsl_matrix * q,
								unsigned int n = 0,
								unsigned int p = 0,
								bool data = false);

			/**@fn virtual int tkalman_original_em :: setup(const gsl_vector * x0,
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


			/**@fn virtual void tkalman_original_em :: filter(const gsl_vector * const * observations,
															  unsigned int nb_observations) = 0
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette méthode effectue le filtrage des données par le filtre de Kalman Triple non supervisé.
			 */
			virtual void filter(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn virtual void tkalman_original_em :: smooth(const gsl_vector * const * observations,
															  unsigned int nb_observations)
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette méthode effectue le lissage des données par le filtre de Kalman Triple non supervisé.
			 */
			virtual void smooth(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn tkalman_original_em :: ~tkalman_original_em
			 * @brief
			 * Destructeur
			 *
			 */
			~tkalman_original_em(void);

			/**@fn bool tkalman_original_em :: operator!()
			 * @return
			 - 0 si l'objet est valide
			 - 1 sinon
			 * @brief
			 Cette méthode teste la validité de chaque attribut.
			**/
			virtual bool operator!() const;

			//Accesseurs
			/**@fn const gsl_vector * const * tkalman_original_em :: get_x0_est()
			 * @return
			 Liste des x0 estimés
			 * @brief
			 Cette fonction renvoie les différents x0 estimés par l'algorithme EM.
			**/
			inline const gsl_vector * const * get_x0_est() const
			{
                return _x0_est;
			}

			/**@fn const gsl_matrix * const * tkalman_original_em :: get_p0_est()
			 * @return
			 Liste des p0 estimés
			 * @brief
			 Cette fonction renvoie les différents p0 estimés par l'algorithme EM.
			**/
			inline const gsl_matrix * const * get_p0_est() const
			{
                return _p0_est;
			}

			/**@fn const gsl_matrix * const * tkalman_original_em :: get_f_est()
			 * @return
			 Liste des F estimés
			 * @brief
			 Cette fonction renvoie les différents F estimés par l'algorithme EM.
			**/
			inline const gsl_matrix * const * get_f_est() const
			{
                return _f_est;
			}

			/**@fn const gsl_matrix * const * tkalman_original_em :: get_q_est()
			 * @return
			 Liste des Q estimés
			 * @brief
			 Cette fonction renvoie les différents Q estimés par l'algorithme EM.
			**/
			inline const gsl_matrix * const * get_q_est() const
			{
                return _q_est;
			}

			/**@fn inline const double * tkalman_original_em :: get_log_likehehood()
			 * @return
			 Valeurs de la log-vraisemblance au cours des itérations de l'EM
			 * @brief
			 Cette fonction renvoie les valeurs de la log vraisemblance au cours des itérations de l'EM.
			**/
			inline const double * get_log_likelihood() const
			{
				return _log_likelihood;
			}

			/**@fn inline unsigned int  tkalman_original_em :: get_nb_iter()
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
			/**@fn void tkalman_original_em :: do_em_algorithm(const gsl_vector * const * observations);
			 * @param[in] observations : observations
			 * @brief
			 * Cette méthode estime les paramètres du filtre de Kalman couple à partir d'un jeu d'observations.
			 */
			void do_em_algorithm(const gsl_vector * const * observations);

			/**@fn void tkalman_original_em :: estimate_parameters(const gsl_vector * const * observations,
																   unsigned int i);
			 * @param[in] observations : observations
			 * @param[in] i : itération de l'EM
			 * @brief
			 * Cette méthode ré-estime les paramètres du filtre de Kalman couple à partir d'un jeu d'observations.
			 */
			void estimate_parameters(const gsl_vector * const * observations,
									 unsigned int i);



			/**@fn void tkalman_original_em :: follow(const gsl_vector * const * observations, unsigned int i);
			 * @param i : itération de l'EM
			 * @brief
			 * Cette méthode enregistre les paramètres de l'itération i de l'algorithme EM.
			 **/
			void follow(const gsl_vector * const * observations,
						unsigned int i);


			/**@fn void tkalman_original_em :: initialize();
			 * @brief
			 Cette méthode met tous les attributs à zéro.
			**/
			void initialize();

			/**@fn void tkalman_original_em :: initialize_tmp();
			 * @brief
			 Cette méthode met tous les temporaires à zéros
			**/
			void initialize_tmp();

			/**@fn int tkalman_original_em :: alloc(unsigned int p, bool data);
			 * @param p : nombre d'itérations de l'EM
			 * @param data : booléen
			 * @return
			  - 0 si l'allocation des attributs de l'objet s'est bien déroulée
			  - 1 sinon
			  *@brief
			  Cette méthode alloue les attributs de l'objet.
			 */
			int alloc(unsigned int p, bool data);

			/**@fn int tkalman_original_em :: alloc_tmp();
			 * @return
			  - 0 si l'allocation des temporaires de l'objet s'est bien déroulée
			  - 1 sinon
			  *@brief
			  Cette méthode alloue les temporaires de l'objet.
			 */
			int alloc_tmp();

			/**@fn int tkalman_original_em :: alloc_data()
			 * @brief
			 * Cette méthode alloue les données
			 *
			 */
			int alloc_data();

			/**@fn void tkalman_original_em :: free();
			 * @brief
			 Cette méthode libère la mémoire utilisée par les attributs de l'objet.
			**/
			void free();

			/**@fn void tkalman_original_em :: free_tmp();
			 * @brief
			 Cette méthode libère la mémoire utilisée par les temporaires de l'objet.
			**/
			void free_tmp();


			/**@fn void tkalman_original_em :: free_data();
			 *  @brief
			 Cette méthode libère la mémoire utilisée par les données de suivi de l'EM.
			**/
			void free_data();

			/**@fn bool tkalman_original_em :: check_tmp()
			 * @return
			 - 0 si les tmp ont été bien alloués
			 - 1 sinon
			 * @brief
			 Cette méthode teste la bonne allocation des variables temporaires utilisées par l'objet.
			**/
			bool check_tmp() const;

			/**@fn bool tkalman_original_em :: check_data()
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


			//Tmp de l'EM
				//Sommes
					gsl_matrix * c_00;
						gsl_matrix c_00_xx;
					gsl_matrix * c_10;
						gsl_matrix c_10_xx;
					gsl_matrix * c_11;
						gsl_matrix c_11_xx;
				//Tmp tmp
					gsl_vector * vect_t_1;
						gsl_vector vect_t_1_x;
						gsl_vector vect_t_1_y;
					gsl_vector * vect_t_2;
						gsl_vector vect_t_2_x;
						gsl_vector vect_t_2_y;
                    gsl_matrix * mat_yy_1;
                    gsl_permutation * perm_y_1;


	};

#endif
