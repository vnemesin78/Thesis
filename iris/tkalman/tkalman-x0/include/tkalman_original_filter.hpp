/**@file tkalman_orignal_filter.hpp
 * @author Valérian Némesin
 * @brief
 Ce fichier contient le prototype de la classe qui permet de créer des objets qui serviront de filtre de Kalman triplet.
 */

#ifndef _TKALMAN_ORIGINAL_FILTER_HPP_
	#define _TKALMAN_ORIGINAL_FILTER_HPP_
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "tkalman_base.hpp"
	#include "tkalman_constants.hpp"
	#include "tkalman_filtering.hpp"
	#include "tkalman_smoothing.hpp"
	#include "gsl_triangle_matrix.hpp"
	/**@class tkalman_orignal_filter
	 *
	 */
	class tkalman_original_filter : public tkalman_base
	{
		public :
			/**@fn  tkalman_original_filter :: tkalman_orignal_filter(const gsl_vector * x0,
																	 const gsl_matrix * p0,
																	 const gsl_matrix * f,
																	 const gsl_matrix * q,
																	 unsigned int n = 0);
			 * @param[in] x0 : Espérance de l'état initial
			 * @param[in] p0 : Matrice de covariance de l'état initial (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
			 * @param[in] f : Matrice d'évolution
			 * @param[in] q : Matrice de covariance (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
			 * @param[in] n : Nombre d'observations (0 par défaut)
			 * @brief
			 constructeur de l'objet.
			 */
			tkalman_original_filter(const gsl_vector * x0,
								    const gsl_matrix * p0,
								    const gsl_matrix * f,
								    const gsl_matrix * q,
								    unsigned int n = 0);

			/**@fn virtual int tkalman_original_filter :: setup(const gsl_vector * x0,
															   const gsl_matrix * p0,
															   const gsl_matrix * f,
															   const gsl_matrix * q,
															   unsigned int n = 0);
			 * @param[in] x0 : Espérance de l'état initial
			 * @param[in] p0 : Matrice de covariance de l'état initial (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
			 * @param[in] f : Matrice d'évolution
			 * @param[in] q : Matrice de covariance (remplacée par sa décomposition de Cholesky dans certaines des classes filles)
			 * @param[in] n : Nombre d'observations (0 par défaut)
			 * @return
			 * - 0 si l'objet est valide
			 * - 1 en cas de problème
			 * @brief
			 * Cette méthode initialise l'objet.
			 */
			virtual int setup(const gsl_vector * x0,
							  const gsl_matrix * p0,
							  const gsl_matrix * f,
							  const gsl_matrix * q,
							  unsigned int n = 0);

			/**@fn virtual void tkalman_original_filter :: filter(const gsl_vector * const * observations,
													   unsigned int nb_observations) = 0
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette méthode effectue le filtrage des données par le filtre de Kalman Triple.
			 */
			virtual void filter(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn virtual bool tkalman_original_filter :: operator!() const;
			 * @return
			 * - 0 si l'objet est bien initialisé.
			 * - 1 sinon
			 * @brief
			 * Cette fonction controle les allocations mémoires des différents attributs de l'objet.
			 */
			virtual bool operator!() const;


			/**@fn virtual void tkalman_original_filter :: smooth(const gsl_vector * const * observations,
																unsigned int nb_observations)
			 * @param[in] observations : observations
			 * @param[in] nb_observations : nombre d'observations
			 * @brief
			 Cette méthode effectue le lissage des données par le filtre de Kalman Triple.
			 */
			virtual void smooth(const gsl_vector * const * observations,
								unsigned int nb_observations);

			/**@fn tkalman_original_filter :: ~tkalman_original_filter();
			 * @brief
			 destructeur de l'objet.
			 */
			virtual ~tkalman_original_filter();

			/**@fn int tkalman_original_filter :: check_positivity() const;
			 * @return
			 * - 0 si Q et P0 sont positives
			 * - 1 sinon.
			 */
			int check_positivity() const;
		protected :

            /**@fn void tkalman_original_filter :: filter_without_equivalents(const gsl_vector * const * observations);
             * @param observations : observations
             * @brief
             * Cette méthode effectue le filtrage avec les paramètres réduits.
             */
            void filter_without_equivalents(const gsl_vector * const * observations);



            /**@fn void tkalman_original_filter :: smooth_without_equivalents(const gsl_vector * const * observations);
            * @param observations : observations
            * @brief
            * Cette méthode effectue le lissage avec les paramètres réduits.
            */
            void smooth_without_equivalents(const gsl_vector * const * observations);




			/**@fn void tkalman_original_filter :: compute_equivalents_p_p_and_p_f();
			 * @brief
			 * Cette méthode calcule les p_p et p_f équivalents.
			 */
			void compute_equivalents_p_p_and_p_f();

			/**@fn void tkalman_original_filter :: compute_equivalents_p_s()
			 * @brief
			 * Cette méthode calcule les p_s équivalents.
			 */
			void compute_equivalents_p_s();

			/**@fn void tkalman_original_filter :: compute_constants()
			 * @brief
			 * Cette méthode calcule les constantes associées aux paramètres du filtre de Kalman couple.
			 *
			 */
			void compute_constants();


			/**@fn void tkalman_original_filter :: initialize();
			 * @brief
			 Cette méthode met tous les attributs à zéro.
			**/
			void initialize();

			/**@fn void tkalman_original_filter :: initialize_tmp();
			 * @brief
			 Cette méthode met tous les temporaires à zéros
			**/
			void initialize_tmp();

			/**@fn int tkalman_original_filter :: alloc();
			 * @return
			  - 0 si l'allocation des attributs de l'objet s'est bien déroulée
			  - 1 sinon
			  *@brief
			  Cette méthode alloue les attributs de l'objet.
			 */
			int alloc();

			/**@fn int tkalman_original_filter :: alloc_tmp();
			 * @return
			  - 0 si l'allocation des temporaires de l'objet s'est bien déroulée
			  - 1 sinon
			  *@brief
			  Cette méthode alloue les temporaires de l'objet.
			 */
			int alloc_tmp();

			/**@fn void tkalman_original_filter :: free();
			 * @brief
			 Cette méthode libère la mémoire utilisée par les attributs de l'objet.
			**/
			void free();

			/**@fn void tkalman_original_filter :: free_tmp();
			 * @brief
			 Cette méthode libère la mémoire utilisée par les temporaires de l'objet.
			**/
			void free_tmp();

			/**@fn bool tkalman_original_filter :: check_tmp() const;
			 * @return
			 * - 0 si les tmps sont bien alloués
			 * - 1 sinon
			 * @brief
			 * Cette fonction controle l'allocation mémoire des tmp.
			 */
			bool check_tmp() const;


			//Tmp
			gsl_matrix * mat_tt_1;
				gsl_matrix mat_tt_1_xx;
				gsl_matrix mat_tt_1_xy;
				gsl_matrix mat_tt_1_yx;
				gsl_matrix mat_tt_1_yy;


	};

#endif
