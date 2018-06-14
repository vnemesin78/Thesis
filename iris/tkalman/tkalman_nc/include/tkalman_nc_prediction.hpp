/**@fn tkalman_nc_prediction.hpp
 * @author Valérian Némesin
 **/
#ifndef _TKALMAN_NC_PREDICTION_HPP_
	#define _TKALMAN_NC_PREDICTION_HPP_
	#include "gsl_triangle_matrix.hpp"
	#include <gsl/gsl_matrix.h> //Matrices - Vecteurs
	#include <gsl/gsl_linalg.h> //Décompo. QR + Inversions LU
	#include <gsl/gsl_blas.h> //Produits matriciels.
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	using namespace std;
	/**@class tkalman_nc_prediction
	 * @brief
	 * Cette classe permet de réaliser les différents types de prédiction du filtre de Kalman couple sans corrélation.  Les différentes méthodes à utiliser sont :
	 * - @fn void tkalman_nc_prediction :: compute_prediction_1
	 * - @fn void tkalman_nc_prediction :: compute_prediction
	 * - @fn void tkalman_nc_prediction :: compute_prediction_
	 */
	class tkalman_nc_prediction
	{
		public: 
			/**@fn tkalman_nc_prediction :: tkalman_nc_prediction(const gsl_matrix * f,
																  const gsl_matrix * sqrt_q,
																  unsigned int size_x) throw(exception &);
			 * @param[in] f: F
			 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
			 * @param[in] size_x : dim. de x;
			 * @brief
			 * Constructeur de la classe @class tkalman_nc_prediction
			 * @throw 
			 * Exception (std :: bad_alloc si problème de mémoire ou invalid_argument en cas d'arguments invalides)
			 */
			tkalman_nc_prediction(const gsl_matrix * f,
								  const gsl_matrix * sqrt_q,
								  unsigned int size_x) throw(exception &);
			
			/**@fn void tkalman_nc_prediction :: setup(const gsl_matrix * f,
													   const gsl_matrix * sqrt_q,
													   unsigned int size_x) throw(exception &);
			 * @param[in] f: F
			 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
			 * @param[in] size_x : dim. de x;
			 * @brief
			 * Cette fonction libères les attributs et les réalloue.
			 */
			void setup(const gsl_matrix * f,
					   const gsl_matrix * sqrt_q,
					   unsigned int size_x) throw(exception &);
			
			/**@fn bool tkalman_nc_prediction :: operator !() const
			 * @return 
			 * - 0 si l'objet est normal
			 * - 1 sinon.
			 * @brief
			 * Check de l'objet.
			 **/
			virtual bool operator !() const;
			
			/**@fn void tkalman_nc_prediction :: compute_prediction_1(gsl_vector * x_p_1,
																	  gsl_matrix * sqrt_p_p_1,
																	  const gsl_vector * t0_f,
																	  const gsl_matrix * sqrt_q_f_0)
			  * @param x1_p : \hat{x}_{1|0}, espérance de l'état prédit 1
			  * @param sqrt_p_p_1 : [P_{1|0}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit 1.
			  * @param[in] t0_f : \hat{t}_{0|0}, espérance de t0 filtré
			  * @param[in] sqrt_q_f_0 : [Q_{0|0}]^{\frac{1}{2}}, racine de la matrice de covariance de \hat{t}_{0|0}.
			  * @brief
				Cette fonction effectue la prédiction 1 dans le filtre de Kalman non-corrélé (cf @fn void tkalman_nc_get_x1_p et @fn void tkalman_nc_get_sqrt_p1_p pour comprendre ces différentes phases).
			 **/
			void compute_prediction_1(gsl_vector * x_p_1,
									  gsl_matrix * sqrt_p_p_1,
									  const gsl_vector * t0_f,
									  const gsl_matrix * sqrt_q_f_0);
			
			/**@fn void tkalman_nc_prediction :: compute_prediction(gsl_vector * x_p,
										   gsl_matrix * sqrt_p_p,
										   const gsl_vector * _x_f,
										   const gsl_matrix * _sqrt_p_f,
										   const gsl_vector * __y);
			  * @param x_p : \hat{x}_{n + 1|n}, espérane de l'état prédit courant
			  * @param sqrt_p_p : [P_{n + 1|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit courant.
			  * @param[in] _x_f : \hat{x}_{n|n}, espérance de l'état filtré précédent
			  * @param[in] _sqrt_p_f : [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré précédent
			  * @param[in] __y : y_{n - 1}, observation
			  * @brief
				Cette fonction effectue la phase de prédiction dans le filtre de Kalman non-corrélé (cf @fn void tkalman_nc_get_x_p et @fn void tkalman_nc_get_sqrt_p_p pour comprendre ces différentes phases).
			 **/
			void compute_prediction(gsl_vector * x_p,
								    gsl_matrix * sqrt_p_p,
								    const gsl_vector * _x_f,
								    const gsl_matrix * _sqrt_p_f,
									const gsl_vector * __y);
			
			/**@fn void tkalman_nc_prediction :: compute_prediction_(gsl_vector * t_p,
										    gsl_matrix * sqrt_q_p,
										    const gsl_vector * _t_p,
										    const gsl_matrix * _sqrt_p_p);
			  * @param t_p : \hat{t}_{n + 1 | m}
			  * @param sqrt_q_p : [Q_{n + 1|m}]^{\frac{1}{2}}
			  * @param[in] _t_p : \hat{t}_{n | m}
			  * @param[in] _sqrt_q_p : [Q_{n|m}]^{\frac{1}{2}}
			  * @brief
			  * Cette fonction effectue la prédiction du filtre de Kalman triplet non corrélé sans observations (cf @fn void tkalman_nc_get_t_p et @fn void tkalman_nc_get_sqrt_q_p pour comprendre ces différentes phases).
			 **/
			void compute_prediction_(gsl_vector * t_p,
								     gsl_matrix * sqrt_q_p,
								     const gsl_vector * _t_p,
								     const gsl_matrix * _sqrt_q_p);
			
			/**@fn tkalman_nc_prediction :: ~tkalman_nc_prediction();
			 * @brief
			 * Destructeur de la classe @class tkalman_nc_prediction
			 */
			~tkalman_nc_prediction();
			
			
			//Accéseurs
			/**@fn inline unsigned int  tkalman_nc_prediction :: size_x() const
			 * @return 
			 * Dim. de x
			 */
			inline unsigned int size_x() const
			{
				return _size_x;
			}
			
			/**@fn inline unsigned int  tkalman_nc_prediction :: size_y() const
			 * @return 
			 * Dim. de y
			 */
			inline unsigned int size_y() const
			{
				return _size_y;
			}
			
			/**@fn inline unsigned int  tkalman_nc_prediction :: size_t() const
			 * @return 
			 * Dim. de t
			 */
			inline unsigned int size_t() const
			{
				return _size_t;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_prediction :: f_xt() const
			 * @return 
			 * F^{x,t}
			 **/
			inline const gsl_matrix * f_xt() const
			{
				return &_f_xt;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_prediction :: f_xx() const
			 * @return 
			 * F^{x,x}
			 **/
			inline const gsl_matrix * f_xx() const
			{
				return &_f_xx;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_prediction :: f_xy() const
			 * @return 
			 * F^{x,y}
			 **/
			inline const gsl_matrix * f_xy() const
			{
				return &_f_xy;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_prediction :: f() const
			 * @return 
			 * F
			 **/
			inline const gsl_matrix * f() const
			{
				return _f;
			}
			
			
			
			/**@fn inline const gsl_matrix * tkalman_nc_prediction :: sqrt_q_xx() const
			 * @return 
			 * [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
			 **/
			inline const gsl_matrix * sqrt_q_xx() const
			{
				return &_sqrt_q_xx;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_prediction :: sqrt_q_yy() const
			 * @return 
			 * [Q^{y,y}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de mesure.
			 **/
			inline const gsl_matrix * sqrt_q_yy() const
			{
				return &_sqrt_q_yy;
			}
			
			/**@fn inline const gsl_matrix * tkalman_nc_prediction :: sqrt_q()  const
			 * @return 
			 * Racine de la matrice de covariance du bruit
			 **/
			inline const gsl_matrix * sqrt_q() const
			{
				return _sqrt_q;
			}
			
			
		
		
		protected: 
			
			/**@fn void tkalman_nc_prediction :: initialize();
			 * @brief
			 * Cette fonction met tous les attributs de l'objet à 0.
			 */
			void initialize();
			
			/**@fn void tkalman_nc_prediction :: alloc() throw(exception &);
			 * @brief
			 * Cette fonction alloue les différents élements de la classe.
			 * @throw
			 * bad_alloc en cas de problème de mémoire
			 */
			void alloc() throw(exception &);
			
			/**@fn void tkalman_nc_prediction :: free();
			 * @brief
			 * Cette fonction désalloue tous les attributs alloués.
			 **/
			void free();
			
			/**@fn void tkalman_nc_prediction :: create_views();
			 * @brief 
			 * Cette fonction génère les différentes vues sur les matrices.
			 */
			void create_views();
			
		//Tmps
			//Matrice de taille (x + t, x)
			gsl_matrix * mat_xpt_x,
						 mat_xpt_x_view_00,
						 mat_xpt_x_view_10;
			//Matrice de taille (2x, x)		 
			gsl_matrix * mat_2x_x,
						 mat_2x_x_view_00,
						 mat_2x_x_view_10;
			//Matrice de taille (2t, t)
			gsl_matrix * mat_2t_t,
						 mat_2t_t_view_00,
						 mat_2t_t_view_10;
			//Vecteurs
			gsl_vector * vect_x,
					   * vect_t;
			
		//Paramètres
			unsigned int _size_x;
			unsigned int _size_y;	
			unsigned int _size_t;
			const gsl_matrix * _f;
				gsl_matrix _f_xt;
				gsl_matrix _f_xx;
				gsl_matrix _f_xy;
			const gsl_matrix * _sqrt_q;
				gsl_matrix _sqrt_q_xx;
				gsl_matrix _sqrt_q_yy;
	};
	
	
	
	
	
	/**@fn void tkalman_nc_get_x1_p(gsl_vector * x1_p,
									const gsl_vector * t0_f,
									const gsl_matrix * f_xt)
	* @param x1_p :  \hat{x}_{1|0}, espérance de l'état prédit 1
	* @param[in] t0_f : \hat{t}_{0|0}, espérance de t0 filtré
	* @param[in] f_xt : F^{x,t}.
	* @brief
	* Cette fonction calcule l'espérance de l'état prédit 1 dans le filtre de Kalman triplet sans corrélation entre bruit de process et bruit de mesure. \n
	\hat{x}_{1|0} = F^{x,t}  \hat{t}_{0|0}
	**/
	void tkalman_nc_get_x1_p ( gsl_vector * x1_p,
							   const gsl_vector * t0_f,
							   const gsl_matrix * f_xt );
	/**@fn void tkalman_nc_get_x_p(gsl_vector * x_p,
							   const gsl_vector * _x_f,
							   const gsl_vecotr * __y,
							   const gsl_matrix * f_xx,
							   const gsl_matrix * f_xy)
	* @param x_p : \hat{x}_{n + 1|n}, espérane de l'état prédit courant
	* @param[in] _x_f : \hat{x}_{n|n}, espérance de l'état filtré précédent
	* @param[in] __y : y_{n - 1}, observation
	* @param[in] f_xx : F^{x,x}
	* @param[in] f_xy : F^{x,y} 
	* @brief
	* Cette fonction calcule l'espérance de l'état prédit n dans le filtre de Kalman triplet sans corrélation entre bruit de mesure et bruit de process. Cette fonction est valable jusqu'à n = Card(Y). Après, il est nécessaire d'estimer aussi y. \n
	* \hat{x}_{n + 1|n} = F^{x,x} \hat{x}_{n|n} + F^{x,y} y_{n}
	*/
	void tkalman_nc_get_x_p( gsl_vector * x_p,
							 const gsl_vector * _x_f,
							 const gsl_vector * __y,
							 const gsl_matrix * f_xx,
							 const gsl_matrix * f_xy );
							 
	/**@fn void tkalman_nc_get_t_p( gsl_vector * t_p,
									const gsl_vector * _t_p,
									const gsl_matrix * f )
	 * @param t_p : \hat{t}_{n + 1 | m}
	 * @param _t_p : \hat{t}_{n | m}
	 * @param f : F
	 * @brief
	 * Cette fonction calcule l'espérance du vecteur t prédit suivant. Il faut que m soit strictement inférieur à n. \n
	 * \hat{t}_{n + 1 | m} = F \hat{t}_{n | m}
	 */
	void tkalman_nc_get_t_p( gsl_vector * t_p,
							 const gsl_vector * _t_p,
							 const gsl_matrix * f);
							 
							 
							 
	/**@fn void tkalman_nc_get_sqrt_p1_p(gsl_matrix * sqrt_p_p_1,
									 const gsl_matrix * sqrt_q_f_0,
									 const gsl_matrix * f_xt,
									 const gsl_matrix * sqrt_q_xx,
									 gsl_matrix * mat_xpt_x,
									 gsl_matrix * mat_xpt_x_view_00,
									 gsl_matrix * mat_xpt_x_view_10,
									 gsl_vector * vect_x)
	* @param sqrt_p_p_1 : [P_{1|0}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit 1.
	* @param[in] sqrt_q_f_0 : [Q_{0|0}]^{\frac{1}{2}}, racine de la matrice de covariance de \hat{t}_{0|0}.
	* @param[in] f_xt : F^{x,t}.
	* @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
	* @param mat_xpt_x : matrice de taille (x+t, x) allouée.
	* @param mat_xpt_x_view_00 : vue sur la matrice mat_xpt_x (de (0,0) à (n_x - 1, n_x - 1)).
	* @param mat_xpt_x_view_10 : vue sur la matrice mat_xpt_x (de (n_x,0) à (n_x + n_t - 1, n_x - 1)).
	* @param vect_x : vecteur de taille (n_x) alloué
	* @brief
	* Cette fonction calcule la racine de la matrice de covariance de l'état prédit 1, [P_{1|0}]^{\frac{1}{2}} dans le filtre de Kalman triplet non-corrélé. \n
	* Dans un premier temps, nous construisons la matrice M : \newline
	* M = \begin{pmatrix}
	* [Q^{x,x}]^{\frac{1}{2}} \newline
	* [Q_{0|0}]^{\frac{1}{2}} [F^{x,t}]^T
	* \end{pmatrix}
	* \newline
	* Dans un second temps, nous effectuons la décomposition QR de cette matrice : 
	* M = Q
	* \begin{pmatrix}
	* [P_{1|0}]^{\frac{1}{2}}} \newline
	* 0
	* \end{pmatrix}
	* 
	*/
	void tkalman_nc_get_sqrt_p1_p( gsl_matrix * sqrt_p_p_1,
								   const gsl_matrix * sqrt_q_f_0,
								   const gsl_matrix * f_xt,
								   const gsl_matrix * sqrt_q_xx,
								   gsl_matrix * mat_xpt_x,
								   gsl_matrix * mat_xpt_x_view_00,
								   gsl_matrix * mat_xpt_x_view_10,
								   gsl_vector * vect_x );
	/**@fn void tkalman_nc_get_sqrt_p_p(gsl_matrix * sqrt_p_p,
									const gsl_matrix * _sqrt_p_f,
									const gsl_matrix * f_xx,
									const gsl_matrix * sqrt_q_xx,
									gsl_matrix * mat_2x_x,
									gsl_matrix * mat_2x_x_view_0,
									gsl_matrix * mat_2x_x_view_1,
									gsl_vector * vect_x)
	* @param sqrt_p_p : [P_{n + 1|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit courant.
	* @param[in] _sqrt_p_f : [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré précédent
	* @param[in] f_xx : F^{x,x}
	* @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
	* @param mat_2xx : Matrice de taille (2x.x) préallouée.
	* @param mat_2xx_view_00 : vue sur la matrice mat_2xx allant de (0,0) à (n_x - 1, n_x - 1)
	* @param mat_2xx_view_10 : vue sur la matrice mat_2xx allant de (n_x, 0) à (2n_x - 1, n_x - 1)
	* @param vect_x : vecteur de taille (x) préalloué
	* @brief
	* Cette fonction calcule la racine de la matrice de covariance de l'état prédit, [P_{n + 1|n}]^{\frac{1}{2}} dans le filtre de Kalman triplet non-corrélé. \n
	* Dans un premier temps, nous construisons la matrice M : \newline
	* M = \begin{pmatrix}
	* [Q^{x,x}]^{\frac{1}{2}} \newline
	* [P_{n|n}]^{\frac{1}{2}} [F^{x,x}]^T
	* \end{pmatrix}
	* \newline
	* Dans un second temps, nous effectuons la décomposition QR de cette matrice : 
	* M = Q
	* \begin{pmatrix}
	* [P_{n + 1|n}]^{\frac{1}{2}} \newline
	* 0
	* \end{pmatrix}
	* 
	*/
	void tkalman_nc_get_sqrt_p_p ( gsl_matrix * sqrt_p_p,
								   const gsl_matrix * _sqrt_p_f,
								   const gsl_matrix * f_xx,
								   const gsl_matrix * sqrt_q_xx,
								   gsl_matrix * mat_2x_x,
								   gsl_matrix * mat_2x_x_view_0,
								   gsl_matrix * mat_2x_x_view_1,
								   gsl_vector * vect_x );
	/**@fn void tkalman_nc_get_sqrt_q_p( gsl_matrix * sqrt_q_p,
								  const gsl_matrix * _sqrt_q_p,
								  const gsl_matrix * f,
								  const gsl_matrix * sqrt_q,
								  gsl_matrix * mat_2t_t,
								  gsl_matrix * mat_2t_t_view_00,
								  gsl_matrix * mat_2t_t_view_10,
								  gsl_vector * vect_t )
	 * @param sqrt_p_p : [Q_{n + 1|m}]^{\frac{1}{2}}
	 * @param[in] _sqrt_p_p : [Q_{n|m}]^{\frac{1}{2}}
	 * @param[in] f : F
	 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
	 * @param mat_2tt : Matrice de taille (2n_t.n_t) préallouée.
	 * @param mat_2tt_view_00 : vue sur la matrice mat_2xx allant de (0,0) à (n_t - 1, n_t - 1)
	 * @param mat_2tt_view_10 : vue sur la matrice mat_2xx allant de (n_t, 0) à (2n_t - 1, n_t - 1)
	 * @param vect_x : vecteur de taille (x) préalloué
	 * @brief
	 * Cette fonction calcule racine de la matrice de covariance du vecteur aléatoire t_{n + 1|m}. \n
	 * Dans un premier temps, nous construisons la matrice M : \newline
	 * M = \begin{pmatrix}
	 * [Q]^{\frac{1}{2}} \newline
	 * [Q_{n|m}]^{\frac{1}{2}} F^T
	 * \end{pmatrix}
	 * \newline
	 * Dans un second temps, nous effectuons la décomposition QR de cette matrice : 
	 * M = Q
	 * \begin{pmatrix}
	 * [Q_{n + 1|m}]^{\frac{1}{2}}} \newline
	 * 0
	 * \end{pmatrix}
	 **/
	void tkalman_nc_get_sqrt_q_p( gsl_matrix * sqrt_q_p,
								  const gsl_matrix * _sqrt_q_p,
								  const gsl_matrix * f,
								  const gsl_matrix * sqrt_q,
								  gsl_matrix * mat_2tt,
								  gsl_matrix * mat_2tt_view_00,
								  gsl_matrix * mat_2tt_view_10,
								  gsl_vector * vect_t );
	/**@fn void tkalman_nc_do_prediction ( gsl_vector * x_p,
									   gsl_matrix * sqrt_p_p,
									   const gsl_vector * _x_f,
									   const gsl_matrix * _sqrt_p_f,
									   const gsl_vector * __y,
									   const gsl_matrix * f_xx,
									   const gsl_matrix * f_xy,
									   const gsl_matrix * sqrt_q_xx,
									   gsl_matrix * mat_2xx,
									   gsl_matrix * mat_2xx_view_00,
									   gsl_matrix * mat_2xx_view_10,
									   gsl_vector * vect_x )
	* @param x_p : \hat{x}_{n + 1|n}, espérane de l'état prédit courant
	* @param sqrt_p_p : [P_{n + 1|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit courant.
	* @param[in] _x_f : \hat{x}_{n|n}, espérance de l'état filtré précédent
	* @param[in] _sqrt_p_f : [P_{n|n}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état filtré précédent
	* @param[in] __y : y_{n - 1}, observation
	* @param[in] f_xx : F^{x,x}
	* @param[in] f_xy : F^{x,y} 
	* @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
	* @param mat_2xx : Matrice de taille (2x.x) préallouée.
	* @param mat_2xx_view_00 : vue sur la matrice mat_2xx allant de (0,0) à (n_x - 1, n_x - 1)
	* @param mat_2xx_view_10 : vue sur la matrice mat_2xx allant de (n_x, 0) à (2n_x - 1, n_x - 1)
	* @param vect_x : vecteur de taille (x) préalloué
	* @brief
	Cette fonction effectue la phase de prédiction dans le filtre de Kalman non-corrélé (cf @fn void tkalman_nc_get_x_p et @fn void tkalman_nc_get_sqrt_p_p pour comprendre ces différentes phases).
	**/
	void tkalman_nc_do_prediction ( gsl_vector * x_p,
									gsl_matrix * sqrt_p_p,
									const gsl_vector * _x_f,
									const gsl_matrix * _sqrt_p_f,
									const gsl_vector * __y,
									const gsl_matrix * f_xx,
									const gsl_matrix * f_xy,
									const gsl_matrix * sqrt_q_xx,
									gsl_matrix * mat_2xx,
									gsl_matrix * mat_2xx_view_00,
									gsl_matrix * mat_2xx_view_10,
									gsl_vector * vect_x );
	/**@fn void tkalman_nc_do_prediction_1 ( gsl_vector * x_p_1,
										 gsl_matrix * sqrt_p_p_1,
										 const gsl_vector * t0_f,
										 const gsl_matrix * sqrt_q_f_0,
										 const gsl_matrix * f_xt,
										 const gsl_matrix * sqrt_q_xx,
										 gsl_matrix * mat_xpt_x,
										 gsl_matrix * mat_xpt_x_view_00,
										 gsl_matrix * mat_xpt_x_view_10,
										 gsl_vector * vect_x )
	* @param x1_p : \hat{x}_{1|0}, espérance de l'état prédit 1
	* @param sqrt_p_p_1 : [P_{1|0}]^{\frac{1}{2}}, racine de la matrice de covariance de l'état prédit 1.
	* @param[in] t0_f : \hat{t}_{0|0}, espérance de t0 filtré
	* @param[in] sqrt_q_f_0 : [Q_{0|0}]^{\frac{1}{2}}, racine de la matrice de covariance de \hat{t}_{0|0}.
	* @param[in] f_xt : F^{x,t}.
	* @param[in] sqrt_q_xx : [Q^{x,x}]^{\frac{1}{2}}, racine de la matrice de covariance du bruit de process.
	* @param mat_xpt_x : matrice de taille (n_x+n_t, x) allouée
	* @param mat_xpt_x_view_00 : vue sur la matrice mat_xpt_x (de (0,0) à (n_x - 1, n_x - 1))
	* @param mat_xpt_x_view_10 : vue sur la matrice mat_xpt_x (de (n_x,0) à (n_x + n_t - 1, n_x - 1))
	* @param vect_x : vecteur de taille (n_x) alloué
	* @brief
	Cette fonction effectue la prédiction 1 dans le filtre de Kalman non-corrélé (cf @fn void tkalman_nc_get_x1_p et @fn void tkalman_nc_get_sqrt_p1_p pour comprendre ces différentes phases).
	*/
	void tkalman_nc_do_prediction_1 ( gsl_vector * x_p_1,
									  gsl_matrix * sqrt_p_p_1,
									  const gsl_vector * t0_f,
									  const gsl_matrix * sqrt_q_f_0,
									  const gsl_matrix * f_xt,
									  const gsl_matrix * sqrt_q_xx,
									  gsl_matrix * mat_xpt_x,
									  gsl_matrix * mat_xpt_x_view_00,
									  gsl_matrix * mat_xpt_x_view_10,
									  gsl_vector * vect_x );

	
	/**@fn void tkalman_nc_do_prediction_ ( gsl_vector * t_p,
											const gsl_matrix * sqrt_q_p,
											const gsl_vector * _t_p,
											const gsl_matrix * _sqrt_q_p,
											const gsl_matrix * f,
											const gsl_matrix * sqrt_q,
											gsl_matrix * mat_2tt,
											gsl_matrix * mat_2tt_view_00,
											gsl_matrix * mat_2tt_view_10,
											gsl_vector * vect_t )
	 * @param t_p : \hat{t}_{n + 1 | m}
	 * @param sqrt_q_p : [Q_{n + 1|m}]^{\frac{1}{2}}
	 * @param[in] _t_p : \hat{t}_{n | m}
	 * @param[in] _sqrt_q_p : [Q_{n|m}]^{\frac{1}{2}}* 
	 * @param[in] f : F
	 * @param[in] sqrt_q : [Q]^{\frac{1}{2}}, racine de la matrice de covariance du bruit
	 * @param mat_2tt : Matrice de taille (2n_t.n_t) préallouée.
	 * @param mat_2tt_view_00 : vue sur la matrice mat_2tt allant de (0,0) à (n_t - 1, n_t - 1)
	 * @param mat_2tt_view_10 : vue sur la matrice mat_2tt allant de (n_t, 0) à (2n_t - 1, n_t - 1)
	 * @param vect_t : vecteur de taille (n_t) préalloué
	 * @brief
	 * Cette fonction effectue la prédiction du filtre de Kalman triplet non corrélé sans observations (cf @fn void tkalman_nc_get_t_p et @fn void tkalman_nc_get_sqrt_q_p pour comprendre ces différentes phases).
	 */
	void tkalman_nc_do_prediction_ ( gsl_vector * t_p,
									 gsl_matrix * sqrt_q_p,
									 const gsl_vector * _t_p,
									 const gsl_matrix * _sqrt_q_p,
									 const gsl_matrix * f,
									 const gsl_matrix * sqrt_q,
									 gsl_matrix * mat_2tt,
									 gsl_matrix * mat_2tt_view_00,
									 gsl_matrix * mat_2tt_view_10,
									 gsl_vector * vect_t );
	
	
#endif
