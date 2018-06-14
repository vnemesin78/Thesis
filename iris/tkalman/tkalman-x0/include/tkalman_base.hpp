/**@file tkalman_base.hpp
 * @author Valérian Némesin
 * @brief
 Ce fichier contient le prototype de la classe mère à tous les types de filtre de Kalman triplet.
 *
 *
 */
#ifndef _TKALMAN_BASE_H_
	#define _TKALMAN_BASE_H_

	#include <cmath>
	#include <gsl/gsl_matrix.h>
	#include <gsl/gsl_linalg.h>
	#include <gsl/gsl_blas.h>
	#include "gsl_triangle_matrix.hpp"

	#ifndef tkalman_esperance_unref
		/**@def tkalman_esperance_unref(moment, nb)
		 * @param moment : moment
		 * @param nb : nombre
		 * @brief
		 Cette macro désalloue un moment statistique vectoriel.
		 *
		 */
		#define tkalman_esperance_unref(moment, nb) if (moment) \
													{\
														for( unsigned int i = 0; i < (nb) ; ++ i) \
														{\
															if ( (moment)[i] )\
															{\
																gsl_vector_free((moment)[i]);\
															}\
														}\
														delete [] (moment);\
													}
	#endif	
	
	#ifndef tkalman_covariance_unref
		/**@def tkalman_covariance_unref(moment, nb)
		 * @param moment : moment
		 * @param nb : nombre
		 * @brief
		 Cette macro désalloue un moment statistique matriciel.
		 *
		 */
		#define tkalman_covariance_unref(moment, nb) if (moment) \
													 {\
														for( unsigned int i = 0; i < (nb) ; ++ i) \
														{\
															if ( (moment)[i] )\
															{\
																gsl_matrix_free((moment)[i]);\
															}\
														}\
														delete [] (moment);\
													 }
	#endif
	
	#ifndef tkalman_esperance_ref
		/**@def tkalman_esperance_ref(moment, nb, size)
		 * @param moment : moment
		 * @param nb : nombre
		 * @param size : dimension
		 * @brief
		 Cette macro alloue un moment statistique vectoriel.
		 */
		#define tkalman_esperance_ref(moment, nb, size) if ( !(moment) ) \
														{\
															(moment) = new gsl_vector*[(nb)];\
															if ( (moment) )\
															{\
																for( unsigned int i = 0; i < (nb) ; ++ i) \
																{\
																	(moment)[i] = gsl_vector_calloc(size);\
																}\
															}\
														}
	#endif
	
	#ifndef tkalman_covariance_ref
		/**@def tkalman_covariance_ref(moment, nb, size1, size2)
		 * @param moment : moment
		 * @param nb : nombre
		 * @param size1 : dimension
		 * @param size2 : dimension
		 * @brief
		 Cette macro alloue un moment statistique vectoriel.
		 */
		#define tkalman_covariance_ref(moment, nb, size1, size2) 	if ( !(moment) ) \
																	{\
																		(moment) = new gsl_matrix*[(nb)];\
																		if ( (moment) )\
																		{\
																			for( unsigned int i = 0; i < (nb) ; ++ i) \
																			{\
																				(moment)[i] = gsl_matrix_calloc(size1, size2);\
																			}\
																		}\
																	}
	#endif


	/**@class tkalman_base
	 * @brief
	 Cette classe est la classe mère de tous les filtre de Kalman triplet programmés.
	 */
	class tkalman_base
	{
		public:
		//Constructeur
		/**@fn  tkalman_base :: tkalman_base(const gsl_vector * x0,
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
		 tkalman_base(const gsl_vector * x0,
					  const gsl_matrix * p0,
					  const gsl_matrix * f,
					  const gsl_matrix * q,
					  unsigned int n = 0);

		//setup
		/**@fn virtual int tkalman_base :: setup(const gsl_vector * x0,
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

		//Filtrage
		/**@fn virtual void tkalman_base :: filter(const gsl_vector * const * observations,
						   		   				   unsigned int nb_observations) = 0
		 * @param[in] observations : observations
		 * @param[in] nb_observations : nombre d'observations
		 * @brief
		 Cette méthode effectue le filtrage des données par le filtre de Kalman Triple.
		 */
		virtual void filter(const gsl_vector * const * observations,
						    unsigned int nb_observations) = 0;

		//Lissage
		/**@fn virtual void tkalman_base :: smooth(const gsl_vector * const * observations,
						   		   				   unsigned int nb_observations) = 0
		 * @param[in] observations : observations
		 * @param[in] nb_observations : nombre d'observations
		 * @brief
		 Cette méthode effectue le lissage des données par le filtre de Kalman Triple.
		 */
		virtual void smooth(const gsl_vector * const * observations,
						    unsigned int nb_observations) = 0;

		//Changement de repère
		/**@fn virtual void tkalman_base :: get_equivalent(const gsl_matrix * p);
		 * @param[in] p : matrice de passage entre les système
		 * @brief
		 Cette fonction modifie la matrice de passage du filtre de Kalman triplet. Le filtre utilisé pour débruiter sera le filtre original. Cependant la méthode restaura les moments du filtre souhaité à partir du filtre original.
		 */
		virtual void get_equivalent(const gsl_matrix * p);

		/**@fn virtual void tkalman_base :: get_equivalent(const gsl_matrix * p,
														   const gsl_matrix * p_inv);
		 * @param[in] p : matrice de passage entre les système
		 * @param[in] p_inv : inverse de la matrice de passage
		 * @brief
		 Cette fonction modifie la matrice de passage du filtre de Kalman triplet. Le filtre utilisé pour débruiter sera le filtre original. Cependant la méthode restaura les moments du filtre souhaité à partir du filtre original.
		 */
		virtual void get_equivalent(const gsl_matrix * p,
									const gsl_matrix * p_inv);

		//Accesseurs
			//Paramètres
			/**@fn void tkalman_base :: get_x0(gsl_vector * x0) const;
			 * @param x0 : Espérance de l'état initial(Préalloué)
			 * @brief
			 Cette fonction calcule l'espérance de l'état initial. (P^{x,x}\;\hat{x}_0)
			 */
			void get_x0(gsl_vector * x0) const;


			inline const gsl_vector * x0() const
			{
				return _x0;
			}


			/**@fn void tkalman_base :: get_p0(gsl_vector * p0) const;
			 * @param p0 : Matrice de covariance de l'état initial (Préallouée)
			 * @param mat_xx : matrice temporaire de taille (x.x) (Préallouée)
			 * @brief
			 Cette fonction calcule la matrice de covariance de l'état initial. (P^{x,x}\;P_0\;[P^{x,x}]^{-1})
			 */
			virtual void get_p0(gsl_matrix * p0,
								gsl_matrix * mat_xx) const;
								
			inline const gsl_matrix * p0() const
			{
					return _p0;
			}

			/**@fn void tkalman_base :: get_f(gsl_matrix * f,
											  gsl_matrix * mat_tt) const;
			 * @param f : Matrice d'évolution (Préallouée)
			 * @param mat_tt : matrice de taille (t.t) (Préallouée)
			 * @brief
			 Cette fonction calcule la matrice d'évolution. ((P\;F\;P^{-1})
			 */
			void get_f(gsl_matrix * f,
					   gsl_matrix * mat_tt) const;
					   
					   
			inline const gsl_matrix * f() const
			{
					return _f;
			}
					   

			/**@fn void tkalman_base :: get_q(gsl_matrix * f,
											  gsl_matrix * mat_tt) const;
			 * @param q : Matrice de covariance du bruit (Préallouée)
			 * @param mat_tt : matrice de taille (t.t) (Préallouée)
			 * @brief
			 Cette fonction calcule la matrice de covariance. ((P\;Q\;P^{T})
			 */
			virtual void get_q(gsl_matrix * q,
					   		   gsl_matrix * mat_tt) const;

			inline const gsl_matrix * q() const
			{
					return _q;
			}

			/**@fn inline unsigned int tkalman_base :: size_x() const
			 * @return Dim. de x
			 */
			inline unsigned int size_x() const
			{
				return _size_x;
			}

			/**@fn inline unsigned int tkalman_base :: size_y() const
			 * @return Dim. de y
			 */
			inline unsigned int size_y() const
			{
				return _size_y;
			}

			/**@fn inline unsigned int tkalman_base :: size_t() const
			 * @return Dim. de t = Dim. x + Dim. y
			 */
			inline unsigned int size_t() const
			{
				return _size_t;
			}

			//Moments statistiques
			/**@fn inline const gsl_vector * const * tkalman_base :: x_p() const
			 * @return Espérances des états prédits.(n + 1 élements)
			 */
			inline const gsl_vector * const * x_p() const
			{
				return _x_p;
			}

			/**@fn inline const gsl_vector * const * tkalman_base :: x_f() const
			 * @return Espérances des états filtrés. (n éléments)
			 */
			inline const gsl_vector * const * x_f() const
			{
				return _x_f;
			}

			/**@fn inline const gsl_vector * const * tkalman_base :: x_s() const
			 * @return Espérances des états lissés. (n + 1 éléments)
			 */
			inline const gsl_vector * const * x_s() const
			{
				return _x_s;
			}

			/**@fn inline const gsl_matrix * const * tkalman_base :: p_p() const
			 * @return Matrice de covariance des états prédits (n + 1 éléments)
			 * @warning
			 Cet accesseur renvoie la décomposition de Cholesky dans certaines des classes filles!
			 */
			inline const gsl_matrix * const * p_p() const
			{
				return _p_p;
			}

			/**@fn inline const gsl_matrix * const * tkalman_base :: p_s() const
			 * @return Matrice de covariance des états lissés (n + 1 éléments)
			 * @warning
			 Cet accesseur renvoie la décomposition de Cholesky dans certaines des classes filles!
			 */
			inline const gsl_matrix * const * p_s() const
			{
				return _p_s;
			}

			/**@fn inline const gsl_matrix * const * tkalman_base :: p_f() const
			 * @return Matrice de covariance des états filtrés (n éléments)
			 * @warning
			 Cet accesseur renvoie la décomposition de Cholesky dans certaines des classes filles!
			 */
			inline const gsl_matrix * const * p_f() const
			{
				return _p_f;
			}

			//Nombre d'obs
			/**@fn inline int tkalman_base :: n() const
			 * @return Nombre d'observations
			 */
			inline unsigned int n() const
			{
				return _n;
			}

            /**@fn virtual double tkalman_base :: log_likelihood(const gsl_vector * const * observations,
             * 													 gsl_matrix * mat_yy_1,
																 gsl_vector * vect_y,
																 gsl_permutation * perm_y,
																 gsl_matrix * mat_yy_2) const
             * @param mat_yy_1 : matrice temporaire de taille (y.y) (préallouée)
             * @param vect_y : vecteur temporaire de taille (y) (préalloué)
             * @param perm_y : permutation de taille (y) préallouée
             * @param mat_yy_2 : matrice temporaire de taille (y.y) (préallouée)
             * @return
             Valeur du log-vraisemblance.
             */
            virtual double log_likelihood(const gsl_vector * const * observations,
										  gsl_matrix * mat_yy_1,
										  gsl_vector * vect_y,
										  gsl_permutation * perm_y,
										  gsl_matrix * mat_yy_2) const;

            //Check
            /**@fn bool tkalman_base :: operator ! () const;
             * @return
             * - 1 si l'objet est invalide
             * - 0 sinon
             * @brief
             Cette méthode contrôle la validité de l'objet.
             */
            bool operator ! () const;

            //Destructeur
            /**@fn tkalman_base :: ~tkalman_base();
             * @brief
             destructeur de l'objet.
             */
            virtual ~tkalman_base();

		protected:
		//Checks
			/**@fn bool tkalman_base :: check_params() const;
			 * @return
			 * - 0 si les paramètres sont valides
			 * - 1 sinon
			 */
			bool check_params() const;

			/**@fn bool tkalman_base :: check_moments() const;
			 * @return
			 * - 0 si les moments sont valides.
			 * - 1 sinon
			 */
			bool check_moments() const;

			/**@fn bool tkalman_base :: check_tmp() const;
			 * @return
			 * - 0 si les tmp sont valides.
			 * - 1 sinon
			 */
			bool check_tmp() const;

		//Réajustement des estimations en fonction de P
			/**@fn void tkalman_base :: compute_equivalents_x_f_and_x_p(const gsl_vector * const * observations);
			 * @param observations : observations
			 * @brief
			 * Cette méthode calcule les équivalents avant le lissage.
			**/
			void compute_equivalents_x_f_and_x_p(const gsl_vector * const * observations);

			/**@fn void tkalman_base :: compute_equivalents_x_s();
			 * @param observations : observations
			 * @brief
			 * Cette méthode calcule les équivalents après le lissage.
			**/
			void compute_equivalents_x_s(const gsl_vector * const * observations);

		//Initialisation
			/**@fn void tkalman_base :: initialize();
			 * @brief
			 Cette méthode met tous les attributs à zéro.
			**/
			void initialize();

			/**@fn void tkalman_base :: initialize_params();
			 * @brief
			 Cette méthode met tous les paramètres à zéros
			**/
			void initialize_params();

			/**@fn void tkalman_base :: initialize_moments();
			 * @brief
			 Cette méthode met tous les moments à zéros
			**/
			void initialize_moments();

			/**@fn void tkalman_base :: initialize_tmp();
			 * @brief
			 Cette méthode met tous les temporaires à zéros
			**/
			void initialize_tmp();

		//Lib. mémoire
			/**@fn void tkalman_base :: free();
			 * @brief
			 Cette méthode libère la mémoire utilisée par les attributs de l'objet.
			**/
			void free();

			/**@fn void tkalman_base :: free_params();
			 * @brief
			 Cette  méthode libère la mémoire utilisée par les paramètres de l'objet.
			**/
			void free_params();

			/**@fn void tkalman_base :: free_moments();
			 * @brief
			 Cette méthode libère la mémoire utilisée par les moments de l'objet.
			**/
			void free_moments();

			/**@fn void tkalman_base :: free_tmp();
			 * @brief
			 Cette méthode libère la mémoire utilisée par les temporaires de l'objet.
			**/
			void free_tmp();

		//Alloc. mémoire
			/**@fn int tkalman_base :: alloc();
			 * @return
			  - 0 si l'allocation des attributs de l'objet s'est bien déroulée
			  - 1 sinon
			  *@brief
			  Cette méthode alloue les attributs de l'objet.
			 */
			int alloc();

			/**@fn int tkalman_base :: alloc_params();
			 * @return
			  - 0 si l'allocation des paramètres de l'objet s'est bien déroulée
			  - 1 sinon
			  *@brief
			  Cette méthode alloue les paramètres de l'objet.
			 */
			int alloc_params();

			/**@fn int tkalman_base :: alloc_moments();
			 * @return
			  - 0 si l'allocation des moments de l'objet s'est bien déroulée
			  - 1 sinon
			  *@brief
			  Cette méthode alloue les moments de l'objet.
			 */
			int alloc_moments();

			/**@fn int tkalman_base :: alloc_tmp();
			 * @return
			  - 0 si l'allocation des temporaires de l'objet s'est bien déroulée
			  - 1 sinon
			  *@brief
			  Cette méthode alloue les temporaires de l'objet.
			 */
			int alloc_tmp();

		//Paramètres

			//Matrices de passage entre système
			gsl_matrix * _p;
				gsl_matrix p_xx;
				gsl_matrix p_xy;
				gsl_matrix p_yx;
				gsl_matrix p_yy;
				//Termes intéressants pour le cas où (size_x != size_y)
				gsl_matrix p_ux; 
				gsl_matrix p_uy;
				
				
				
				
			gsl_matrix * _p_inv;
				gsl_matrix p_inv_xx;
				gsl_matrix p_inv_xy;
				gsl_matrix p_inv_yx;
				gsl_matrix p_inv_yy;

			//Espérance de l'état initial
			gsl_vector * _x0;

			//Matrice de covariance de l'état initial (Remplacée par la décomposition de Chol.)
			gsl_matrix * _p0;

			//Matrice d'évolution
			gsl_matrix * _f;
				gsl_matrix f_xx;
				gsl_matrix f_xy;
				gsl_matrix f_yx;
				gsl_matrix f_yy;
			//Matrice de covariance (Remplacée par la décomposition de Chol.)
			gsl_matrix * _q;
				gsl_matrix q_xx;
				gsl_matrix q_xy;
				gsl_matrix q_yx;
				gsl_matrix q_yy;
			//Constantes dérivées des paramètres
				gsl_matrix * q2_xx; // Qxx - Qxy.Qyy⁻¹.Qyx // Remplacé après par chol.(...)
				gsl_matrix * q2_xy; // Qxy.Qyy⁻¹

				gsl_matrix * f2_xx; // Fxx - Qxy.Qyy⁻¹.Fyx
				gsl_matrix * f2_xy; // Fxx - Qxy.Qyy⁻¹.Fyy
			//Dim.
			unsigned int _size_x;
			unsigned int _size_y;
			unsigned int _size_t;
		//Constante
			gsl_vector * vect_y_zero;

		//Moments

			//Nombre d'observations
			unsigned int _n;

			//Espérances
			gsl_vector ** _x_p;
			gsl_vector ** _x_f;
			gsl_vector ** _x_s;
			gsl_vector ** _innovation;

			//Matrices de covariances (Remplacée par la décomposition de Chol.)
			gsl_matrix ** _p_p;
			gsl_matrix ** _p_f;
			gsl_matrix ** _p_s;
			gsl_matrix ** _s;

            //Matrices de covariance entre x_n et x_n+1
            gsl_matrix ** _c_s;

		//tmp (chg. de repère)
			gsl_matrix * mat_xx_1;
			gsl_permutation * perm_x_1;
			gsl_vector * vect_x_1;

	};
#endif
