/**@file auxi_function_tools.hpp
 * @brief
 * Ce fichier contient les fonctions qui permettent de maximiser la fonction : 
 * f_C(x,y) = -N log|X| - 1/2 tr(X^{-1}(\begin{pmatrix} - y \\ I \end{pmatrix} C ( - y \;   I )))
 */
#ifndef _AUXI_FUNCTION_TOOLS_HPP_
	#define _AUXI_FUNCTION_TOOLS_HPP_
	#include "gsl_triangle_matrix.hpp"
	#include <gsl/gsl_matrix.h> //Matrices - Vecteurs
	#include <gsl/gsl_linalg.h> //Décompo. QR + Inversions LU
	#include <gsl/gsl_blas.h> //Produits matriciels.
	#include <exception> //Exceptions (pour ne pas faire planter le programme en cas de problèmes de mémoire ou autre)
	#include <stdexcept>
	#include <iostream> //A virer
	#include <cmath>
	using namespace std;
	/**@class auxi_function_tools
	 * @brief
	 * Cette classe permet de maximiser la fonction \f$ f_C(x,y) = -N log|X| - 1/2 tr(X^{-1}(\begin{pmatrix} - y \\ I \end{pmatrix} C ( - y \;   I ))) \f$ nécessaire à l'algorithme EM du filtre de Kalman triplet.
	 * 
	 */
	class auxi_function_tools
	{
		public :
			/**@fn auxi_function_tools :: auxi_function_tools(const gsl_vector * m,
															  const gsl_matrix * yc,
															  unsigned int c_size) throw(exception &);
			 * @param m : contrainte sur les lignes de y (si m(i) = 0 alors la ligne i de y ne sera pas estimée)
			 * @param yc : partie constante de y
			 * @param c_size : n_y + n_x
			 * @brief
			 *  Constructeur de la classe @class auxi_function_tools
			 */
			auxi_function_tools(const gsl_vector * m,
								const gsl_matrix * yc) throw(exception &);
			
			
			/**@fn void auxi_function_tools :: setup(unsigned int y_size,
													 unsigned int c_size) throw(exception &););
			 * @param y_size :
			 * @param c_size :
			 * @brief
			 *  Setup de la classe @class auxi_function_tools
			 */
			void setup(const gsl_vector * m,
					   const gsl_matrix * yc) throw(exception &);
					   
			/**@fn void auxi_function_tools :: maximize(gsl_matrix * sqrt_x,
													    gsl_matrix * y,
													    gsl_matrix * sqrt_c,
													    unsigned int n); 
			 * @param sqrt_x : racine de de x
			 * @param y : y
			 * @param sqrt_c : racine du paramètre de la fonction
			 * @param n : N
			 * @brief
			 * Cette fonction maximise la fonction auxiliaire.
			 */
			void maximize(gsl_matrix * sqrt_x,
						  gsl_matrix * y,
						  gsl_matrix * sqrt_d,
						  unsigned int n); 
			/**@fn auxi_function_tools :: ~auxi_function_tools();
			 * @brief
			 * Destructeur de la classe @class auxi_function_tools
			 * 
			 */
			~auxi_function_tools();
			
			/**@fn inline unsigned int auxi_function_tools :: size_x() const
			 * @return 
			 * Dim. de x
			 */
			inline unsigned int size_x() const
			{
				return _size_x;
			}
			
			/**@fn inline unsigned int auxi_function_tools :: size_y() const
			 * @return 
			 * Dim. de y
			 */
			inline unsigned int size_y() const
			{
				return _size_y;
			}
			
			/**@fn inline unsigned int auxi_function_tools :: size_c() const
			 * @return 
			 * Dim. de c
			 */
			inline unsigned int size_c() const
			{
				return _size_c;
			}
			
			/**@fn inline unsigned int auxi_function_tools :: rg() const
			 * @return 
			 * Rang de M
			 */
			inline unsigned int rg() const
			{
				return _rg;
			}
			
			/**@fn inline const gsl_vector * auxi_function_tools :: m() const
			 * @return 
			 * Données de la contraintes
			 **/
			inline const gsl_vector * m() const
			{
				return _m;
			}	
			
			/**@fn inline const gsl_vector * auxi_function_tools :: p() const
			 * @return 
			 * Données de la permutation
			 **/
			inline const gsl_vector * p() const
			{
				return _p;
			}	
			
			/**@fn inline const gsl_matrix * auxi_function_tools :: yc() const
			 * @return 
			 * Partie constante de F
			 **/
			inline const gsl_matrix * yc() const
			{
				return _yc;
			}	
			
			/**@fn bool auxi_function_tools :: operator!() const;
			 * @return
			 * - 0 si l'objet est valide
			 * - 1 sinon
			 * 
			 */			    
			virtual bool operator !() const;	
			
			
			
		protected:
			/**@fn void auxi_function_tools :: initialize();
			 * @brief
			 * Initialisateur de la classe @class auxi_function_tools
			 */
			void initialize();
			
			/**@fn void auxi_function_tools :: free();
			 * @brief
			 * Cette méthode libère la mémoire occupée par les attributs de l'objet.
			 */
			void free();
			
			/**@fn void auxi_function_tools :: alloc() throw(exception &);
			 * @brief
			 * Cette méthode alloue la mémoire utilisée par les attributs.
			 */
			void alloc() throw(exception &);
			
			/**@fn void auxi_function_tools :: create_views();
			 * @brief
			 * Cette méthode crée les vues...
			 */
			void create_views();
			
			
			unsigned int _size_x;
			unsigned int _size_y;
			unsigned int _size_c;
			unsigned int _rg;
			const gsl_vector * _m;
			const gsl_matrix * _yc;
			gsl_vector * _p;
			
			
			gsl_matrix * mat_cc_1;
				gsl_matrix mat_cc_1_view_01;
			gsl_matrix * mat_cc_2;
			gsl_matrix * mat_zz;
			gsl_matrix * z_max;
			
			gsl_vector * vect_c;
			gsl_permutation * perm_z;

	};
	
	/**@fn unsigned int auxi_function_tools_p_vector(gsl_vector * p,
													 const gsl_vector * m)
	 * @param p : vecteur contenant toutes les informations de la matrice de permutation permettant de regrouper les termes de la diagonale non nuls.
	 * @param[in] m : vecteur contenant la diagonale de la matrice de contraintes M.
	 * @brief
	 * Cette fonction détermine la matrice de permutation permettant de regrouper tous les termes non nuls de M.
	 * 
	 */
	unsigned int auxi_function_tools_p_vector(gsl_vector * p,
											  const gsl_vector * m);

	/**@fn void auxi_function_tools_ch_var(gsl_matrix * sqrt_d,
										   gsl_matrix * y,
										   gsl_matrix * trans_mat,
										   gsl_matrix * trans_mat_view_01,
										  gsl_matrix * mat)
	 * @param sqrt_d : [D]^{\frac{1}{2}}, racine du paramètre de la fonction
	 * @param y : vecteur de translation du paramètre y
	 * @param trans_mat : matrice de taille size(sqrt_d) allouée
	 * @param trans_mat_view_01 : vue sur trans_mat
	 * @param mat : matrice de taille size(sqrt_d) allouée
	 * @brief
	 * Cette fonction effectue un changement de variable sur le paramètre y de la fonction à maximiser.
	 */
	void auxi_function_tools_ch_var(gsl_matrix * sqrt_d,
									const gsl_matrix * y,
									gsl_matrix * trans_mat,
									gsl_matrix * trans_mat_view_01,
									gsl_matrix * mat);
	/**@fn void auxi_function_tools_constrain_params(gsl_vector * p,
													 gsl_matrix * sqrt_d,
													 const gsl_vector * m,
													 gsl_vector * vect)
	 * @param p : vecteur contenant l'information de la matrice de permutation (de taille size(sqrt_d))
	 * @param sqrt_d : [D]^{\frac{1}{2}}, racine du paramètre de la fonction
	 * @param m : vecteur de contrainte (si m(i) == 0 alors la collonne i de y ne sera pas estimée)
	 * @param vect : vecteur de taille (size(sqrt_d)) alloué
	 * @brief
	 * Cette fonction applique la contraine sur le paramètre y et modifie les paramètres de la fonction à maximiser.
	 */
	void auxi_function_tools_constrain_params(gsl_matrix * sqrt_d,
											  const gsl_vector * m,
											  const gsl_vector * p,
											  unsigned int rg,
											  gsl_vector * vect);

	/**@fn void auxi_function_tools_maximize(gsl_matrix * sqrt_x,
											 gsl_matrix * y,
											 const gsl_matrix * sqrt_d,
											 const gsl_matrix * sqrt_d_view_00,
											 const gsl_matrix * sqrt_d_view_01,
											 const gsl_matrix * sqrt_d_view_11,
											 unsigned int size_x,
											 unsigned int n,
											 gsl_matrix * mat_yy,
											 gsl_permutation * perm_y)
	 * @param sqrt_x : racine du premier paramètre
	 * @param y : second paramètre
	 * @param[in] sqrt_d : paramètre de la fonction
	 * @param[in] sqrt_d_view_00 : vue sur sqrt_q
	 * @param[in] sqrt_d_view_01 : vue sur sqrt_q
	 * @param[in] sqrt_d_view_11 : vue sur sqrt_q
	 * @param[in] size_x : dimension de x
	 * @param[in] n : paramètre
	 * @param mat_yy : matrice préallouée
	 * @param perm_y : permutation préallouée
	 * @brief
	 * Cette fonction calcule les paramètres optimaux de la fonction
	 */
	void auxi_function_tools_maximize(gsl_matrix * sqrt_x,
									  gsl_matrix * y,
									  const gsl_matrix * sqrt_d,
									  const gsl_matrix * sqrt_d_view_00,
									  const gsl_matrix * sqrt_d_view_01,
									  const gsl_matrix * sqrt_d_view_11,
									  unsigned int size_x,
									  unsigned int n,
									  gsl_matrix * mat_yy,
									  gsl_permutation * perm_y);
	/**@fn void auxi_function_tools_restore_y(gsl_matrix * y,
											  const gsl_matrix * y_c,
											  const gsl_matrix * z_max,
											  unsigned int rg,
											  const gsl_vector * p)
	 * @param y : paramètre optimal reconstruit
	 * @param y_c : valeur constante
	 * @param z_max : paramètre réduit opt.
	 * @param rg : rang
	 * @param p : vecteur de permutation
	 * @brief
	 * Cette fonction reconstruit le paramètre y.
	 */
	void auxi_function_tools_restore_y(gsl_matrix * y,
									   const gsl_matrix * y_c,
									   const gsl_matrix * z_max,
									   unsigned int rg,
									   const gsl_vector * p);



#endif
