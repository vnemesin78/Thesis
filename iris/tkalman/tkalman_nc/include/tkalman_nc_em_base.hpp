
#ifndef _TKALMAN_NC_EM_BASE_HPP_
	#define _TKALMAN_NC_EM_BASE_HPP_
	#include "gsl_triangle_matrix.hpp"
	#include "auxi_function_tools.hpp"
	#include "tkalman_nc_fusion.hpp"
	#include <gsl/gsl_matrix.h> //Matrices - Vecteurs
	#include <gsl/gsl_linalg.h> //Décompo. QR + Inversions LU
	#include <gsl/gsl_blas.h> //Produits matriciels.
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	#include <iostream> //A virer
	using namespace std;
	/**@class tkalman_nc_em_base
	 * 
	 * 
	 */
	class tkalman_nc_em_base
	{
		public:
			/**@fn
			 * @param[in] t0 : \hat{t}_0, espérance de l'état initial.
			 * @param[in] sqrt_q0: [Q_0]^{\frac{1}{2}}, racine de la matrice de covariance de l'état initial.
			 * @param[in] size_x : dimension de x
			 * @param[in] n_max : nombre de moments alloués
			 * @param[in] p_max : nombre maximal de signaux
			 * @param[in] x_mask : Masque sur la matrice Fxt
			 * si x_mask(i) == 0, alors Fxt(:, i) == 0
			 * @param[in] y_mask : Masque sur la matrice Fyt
			 * si y_mask(i) == 0, alors Fyt(:, i) == 0
			 * @brief
			 * Constructeur
			**/
			tkalman_nc_em_base(	const gsl_vector * t0,
								const gsl_matrix * sqrt_q0,
								const gsl_matrix * f,
								const gsl_matrix * sqrt_q,
								unsigned int size_x,
								unsigned int n_max,
								unsigned int nb_signal_max,
								const gsl_vector * x_mask = NULL,
								const gsl_vector * y_mask = NULL,
								bool estimate_initial_state = false) throw (exception &);
	
			/**@fn
			 * @param[in] t0 : \hat{t}_0, espérance de l'état initial.
			 * @param[in] sqrt_q0: [Q_0]^{\frac{1}{2}}, racine de la matrice de covariance de l'état initial.
			 * @param[in] size_x : dimension de x
			 * @param[in] n_max : nombre de moments alloués
			 * @param[in] p_max : nombre maximal de signaux
			 * @param[in] x_mask : Masque sur la matrice Fxt
			 * si x_mask(i) == 0, alors Fxt(:, i) == 0
			 * @param[in] y_mask : Masque sur la matrice Fyt
			 * si y_mask(i) == 0, alors Fyt(:, i) == 0
			 * @brief
			 * Setup
			 **/
			virtual void setup (	const gsl_vector * t0,
									const gsl_matrix * sqrt_q0,
									const gsl_matrix * f,
									const gsl_matrix * sqrt_q,
									unsigned int x,
									unsigned int n_max,
									unsigned int p_max,
									const gsl_vector * x_mask = NULL,
									const gsl_vector * y_mask = NULL,
									bool estimate_initial_state = false) throw (exception &);
									




			/**@fn
			 * @brief
			 * Destructeur de la classe @class tkalman_nc_em.
			 */
			~tkalman_nc_em_base();
									
	
	
	protected:
			/**@fn
			 * @brief
			 * Cette fonction met tous les attributs de l'objet à 0.
			 */
			void initialize();
			
			/**@fn
			 * @brief
			 * Cette fonction alloue les différents élements de la classe.
			 * @throw
			 * bad_alloc en cas de problème de mémoire
			 */
			void alloc() throw(exception &);
			
			/**@fn
			 * @brief
			 * Cette fonction désalloue tous les attributs alloués.
			 **/
			void free();
		
			/**@fn
			 * @brief
			 * Cette fonction crée les attributs
			 */
			void create_object() throw(exception &);
		
	
	
		//Fonction auxiliaire
		auxi_function_tools * f_tools_x;
		auxi_function_tools * f_tools_y;
		//Fusion des résultats
		tkalman_nc_fusion * fusion;
		
		//Sommes
		unsigned int _nb_signal_max;
		
		//Paramètres
		unsigned int _size_x,
					 _size_y,
					 _size_t;
		unsigned int _n_max;
		gsl_vector * _t0;
		gsl_matrix * _sqrt_q0;
		
		gsl_vector * _x_mask;
		gsl_vector * _y_mask;
		bool _estimate_initial_state;
		//Tmp
		gsl_matrix * tmp_f_xt_t;
		gsl_matrix * tmp_f_yt_t;
		gsl_matrix * _f_xt_t;
		gsl_matrix * _f_yt_t;
		
		gsl_vector ** t_0s;
		gsl_matrix ** sqrt_q_0s;
		
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
		 * Nombre d'états supportés - 1.
		 */
		inline unsigned int n_max() const
		{
			return _n_max;
		}

	};
#endif
