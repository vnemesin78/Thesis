	 
#ifndef _TKALMAN_NC_FILTER_BASE_HPP_
	#define _TKALMAN_NC_FILTER_BASE_HPP_
	#include "gsl_triangle_matrix.hpp"
	#include "tkalman_nc_prediction.hpp"
	#include "tkalman_nc_filtering.hpp"
	#include "tkalman_nc_smoothing.hpp"
	#include "tkalman_nc_sums.hpp"
	#include <gsl/gsl_matrix.h> //Matrices - Vecteurs
	#include <gsl/gsl_linalg.h> //Décompo. QR + Inversions LU
	#include <gsl/gsl_blas.h> //Produits matriciels.
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	using namespace std;
	
	/**@def tkalman_esperance_unref(moment, nb)
	 * @param moment : moment
	 * @param nb : nombre
	 * @brief
	 Cette macro désalloue un moment statistique vectoriel.
	 *
	**/
	#ifndef tkalman_expectation_unref

		#define tkalman_expectation_unref(moment, nb) if (moment) \
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
		
	#ifndef tkalman_expectation_ref_v2

		#define tkalman_expectation_ref_v2(moment, nb, size_x, size_t) if (! (moment) ) \
		{\
			try \
			{\
				(moment) = (gsl_vector **) malloc( sizeof(gsl_vector * ) * (nb) );\
			}\
			catch(exception & e)\
			{\
				throw(e);\
			}\
			if (moment)\
			{\
				try \
				{\
					(moment)[0] = gsl_vector_calloc(size_t);\
				}\
				catch(exception & e)\
				{\
					throw(e);\
				}\
				for (unsigned j = 1 ; j < (nb); ++ j)\
				{\
					try \
					{\
						(moment)[j] = gsl_vector_calloc(size_x);\
					}\
					catch(exception & e)\
					{\
						throw(e);\
					}\
				}\
			}\
		}
	#endif

	#ifndef tkalman_covariance_ref_v2

		#define tkalman_covariance_ref_v2(moment, nb, size_1, size_2, size_1_0, size_2_0) if (! (moment) ) \
		{\
			try \
			{\
				(moment) = (gsl_matrix **) malloc( sizeof(gsl_matrix * ) * (nb) );\
			}\
			catch(exception & e)\
			{\
				throw(e);\
			}\
			if (moment)\
			{\
				try \
				{\
					(moment)[0] = gsl_matrix_calloc(size_1_0, size_2_0);\
				}\
				catch(exception & e)\
				{\
					throw(e);\
				}\
				for (unsigned j = 1 ; j < (nb) ; ++ j)\
				{\
					try \
					{\
						(moment)[j] = gsl_matrix_calloc(size_1, size_2);\
					}\
					catch(exception & e)\
					{\
						throw(e);\
					}\
				}\
			}\
		}
	#endif


	
	 /**@class tkalman_nc_filter_base
	  * @brief
	  * Classe mère pour les filtre de Kalman
	  * 
	  */
	class tkalman_nc_filter_base
	{
		
		public:
			/**@fn
			 * @param[in] t0 : \hat{t}_0, espérance de l'état initial.
			 * @param[in] sqrt_q0: [Q_0]^{\frac{1}{2}}, racine de la matrice de covariance de l'état initial.
			 * @param[in] size_x : dimension de x
			 * @param[in] n_max : nombre de moments alloués
			 * @brief
			 * Constructeur
			 */
			tkalman_nc_filter_base ( const gsl_vector * t0,
									 const gsl_matrix * sqrt_q0,
									 unsigned int size_x,
									 unsigned int n_max ) throw(exception &);
		
			/**@fn
			 * @param[in] n_max : nombre de moments max
			 * @brief
			 * Modifie le nombre de moment alloué
			 */
			virtual void set_n_max( unsigned int n_max ) throw(exception &);
		
		
			/**@fn
			 * @param[in] t0 : \hat{t}_0, espérance de l'état initial.
			 * @param[in] sqrt_q0: [Q_0]^{\frac{1}{2}}, racine de la matrice de covariance de l'état initial.
			 * @param[in] size_x : dimension de x
			 * @param[in] n_max : nombre de moments alloués
			 * @brief
			 * Setup
			 */
			virtual void setup (	 const gsl_vector * t0,
									const gsl_matrix * sqrt_q0,
									unsigned int size_x,
									unsigned int n_max ) throw(exception &);
		
			/**@fn 
			 * @brief
			 * Destructeur
			 **/
			~tkalman_nc_filter_base();
			
			
		//Méthodes internes
		protected:
			/**@fn
			 * @brief
			 * Lib. mémoire
			 * 
			 **/
			void free();
			/**@fn 
			 * @brief
			 * Initialisation
			 **/
			void initialize();
			/**@fn
			 * @brief
			 * Alloc. mémoire
			 **/
			void alloc() throw(exception &);
			
		//Attributs
		protected:
			//Paramètres initiaux
			const gsl_vector * _t0;
			const gsl_matrix * _sqrt_q0;
			
			//Dimensions
			unsigned int _size_x;
			unsigned int _size_y;
			unsigned int _size_t;
			
			//Moments
			unsigned int _n_max;
			gsl_vector ** _t_pp;
			gsl_vector ** _x_p;
			gsl_vector ** _x_f;
			gsl_vector ** _x_s;
			gsl_vector ** _innovation;
			
			gsl_matrix ** _sqrt_q_pp;
			gsl_matrix ** _sqrt_p_p;
			gsl_matrix ** _sqrt_p_s;
			gsl_matrix ** _sqrt_p_f;
			gsl_matrix ** _sqrt_s;
			gsl_matrix ** _c_s;
			
			//Objets
			tkalman_nc_prediction * prediction;
			tkalman_nc_filtering * filtering;
			tkalman_nc_smoothing * smoothing;
			
			//Tmp
			gsl_vector * vect_t;
				gsl_vector vect_t_view_x;
				
			gsl_matrix * mat_tt;
				gsl_matrix mat_tt_view_xx;
		//Accesseurs	
		public:
			/**@fn
			 * @return t0 : \hat{t}_0, espérance de l'état initial.
			 * 
			 */
			inline const gsl_vector * t0() const
			{
				return _t0;
			}
			
			/**@fn
			 * @return sqrt_q0: [Q_0]^{\frac{1}{2}}, racine de la matrice de covariance de l'état initial.
			 * 
			 */
			inline const gsl_matrix * sqrt_q0() const
			{
				return _sqrt_q0;
			}
			
			/**@fn
			 * @return 
			 * Dim. de x
			 */
			inline unsigned int size_x() const
			{
				return _size_x;
			}
			
			/**@fn
			 * @return 
			 * Dim. de y
			 */
			inline unsigned int size_y() const
			{
				return _size_y;
			}
			
			/**@fn
			 * @return 
			 * Dim. de t
			 */
			inline unsigned int size_t() const
			{
				return _size_t;
			}
		 
			/**@fn
			 * @return
			 * \hat{x}_{n|n - 1}, espérance des états prédits.
			 */
			inline const gsl_vector * const * x_p() const
			{
				return _x_p;
			}

			/**@fn
			 * @return
			 * \hat{x}_{n|n}, espérance des états filtrés.
			 */
			inline const gsl_vector * const * x_f() const
			{
				return _x_f;
			}

			/**@fn
			 * @return
			 * \hat{x}_{n|N}, espérance des états lissés.
			 */
			inline const gsl_vector * const * x_s() const
			{
				return _x_s;
			}

			/**@fn
			 * @return
			 * \hat{t}_{n|N}, Prédiction à l'ordre p
			 */
			inline const gsl_vector * const * t_p() const
			{
				return _t_pp;
			}

			/**@fn
			 * @return
			 * [P_{n|n - 1}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit 
			 **/
			inline const gsl_matrix * const * sqrt_p_p() const
			{
				return _sqrt_p_p;
			}
			
			/**@fn
			 * @return
			 * [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré 
			 **/
			inline const gsl_matrix * const * sqrt_p_f() const
			{
				return _sqrt_p_f;
			}
			
			/**@fn
			 * @return
			 * [P_{n|N}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état lissé 
			 **/
			inline const gsl_matrix * const * sqrt_p_s() const
			{
				return _sqrt_p_s;
			}
			
			/**@fn
			 * @return
			 * [P_{n + 1|N}]^{\frac{1}{2}} K_{n|N}^T
			 */
			inline const gsl_matrix * const * c_s() const
			{
				return _c_s;
			}
			
			/**@fn
			 * @return
			 * [Q_{n|n - p}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit à l'ordre p
			 **/
			inline const gsl_matrix * const * sqrt_q_p() const
			{
				return _sqrt_q_pp;
			}
			
			/**@fn
			 * @return 
			 * Nombre d'états supportés - 1.
			 */
			inline unsigned int n_max() const
			{
				return _n_max;
			}
			
			
		
	};
#endif
