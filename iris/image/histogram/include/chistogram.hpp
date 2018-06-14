/**@file chistogram.hpp
 * 
 */
#ifndef _C_HISTOGRAM_HPP_
	#define _C_HISTOGRAM_HPP_
	#define DEFAULT_NB_BINS 256
	#define EPSILON_DEF 0.05

	
	#define C_HISTOGRAM_STRETCH_1( type )\
	template int c_histogram :: stretch( 	type * data,\
											unsigned int width,\
											unsigned int height,\
											unsigned int width_step,\
											double eps_low,\
											double eps_up );
	
	#define C_HISTOGRAM_STRETCH_2( type1, type2 )\
	template int c_histogram :: stretch( 	type1 * data,\
											const type2 * data_in,\
											unsigned int width,\
											unsigned int height,\
											unsigned int width_step_in,\
											unsigned int width_step_out,\
											double eps_low,\
											double eps_up );
											
	#define C_HISTOGRAM_EQ_2( type1, type2 )\
	template void c_histogram :: equalize( 	type1 * data,\
											const type2 * data_in,\
											unsigned int width,\
											unsigned int height,\
											unsigned int width_step_in,\
											unsigned int width_step_out );		
						
	#define C_HISTOGRAM_C_ALPHA_2( type1, type2 )\
	template int c_histogram :: correct_alpha(	type2 * data_out,\
												const type1 * data_in,\
												unsigned int width,\
												unsigned int height,\
												unsigned int width_step_in,\
												unsigned int width_step_out,\
												double new_mean );
	
	#define C_HISTOGRAM_C_ALPHA( type )\
	template int c_histogram :: correct_alpha(	type * data,\
												unsigned int width,\
												unsigned int height,\
												unsigned int width_step,\
												double new_mean );
	
	
		
	#define C_HISTOGRAM_EQ( type )\
	template void c_histogram :: equalize( 	type * data,\
											unsigned int width,\
											unsigned int height,\
											unsigned int width_step );											
																
											
	#define C_HISTOGRAM_COMPUTE( type )\
	template void c_histogram :: compute(	const type * data,\
											unsigned int width,\
											unsigned int height,\
											unsigned int width_step,\
											const double & sigma,\
											unsigned int c);\
	template void  c_histogram :: compute(	const type * data,\
											const unsigned char * mask_data,\
											unsigned int width,\
											unsigned int height,\
											unsigned int width_step,\
											unsigned int mask_width_step,\
											const double & sigma,\
											unsigned int c);
							
	#define C_HISTOGRAM_COMPUTE_HIST( type )\
	template void c_histogram :: compute_histogram (	const type * data,\
														unsigned int width,\
														unsigned int height,\
														unsigned int width_step);
							
											
											
	
	/**@class c_histogram
	 * 
	 **/
	class c_histogram
	{
		public:
			/**@fn c_histogram :: c_histogram (	unsigned int nb_bins );
			 * @param[in] nb_bins : nombre d'états
			 * @brief
			 * Constructeur
			 */
			c_histogram (	unsigned int nb_bins = 0 );
		
			/**@fn void setup (	unsigned int nb_bins );
			 * @param[in] nb_bins : nombre d'états
			 * @brief
			 * Setup
			 * 
			 **/
			void setup (	unsigned int nb_bins = 0 );
		
			/**@fn void c_histogram :: compute(	const unsigned char * data,
												unsigned int width,
												unsigned int height,
												unsigned int width_step = 0,
												const double & sigma = 1,
												unsigned int c = 0);
			 * @param[in] data : données de l'image
			 * @param[in] width : largeur de l'image
			 * @param[in] height : hauteur de l'image
			 * @param[in] width_step : taille réelle d'une ligne de l'image (si == 0 alors remplacée par width)
			 * @param[in] nb_bins : nombre d'états possibles pour un pixel (si == 0 alors remplacé par _nb_bins)
			 * @param[in] sigma : écart type de la gaussienne
			 * @param[in] c : taille de l'intervalle : [-M, M] avec M = c sigma + 1 (Il faut choisir entre 3 et 6)
			 */
			template <class type> void compute(	const type * data,
												unsigned int width,
												unsigned int height,
												unsigned int width_step = 0,
												const double & sigma = 1,
												unsigned int c = 0);

			/**@fn void c_histogram :: compute(	const unsigned char * data,
												unsigned int width,
												unsigned int height,
												unsigned int width_step = 0,
												const double & sigma = 1,
												unsigned int c = 0);
			 * @param[in] data : données de l'image
			 * @param[in] width : largeur de l'image
			 * @param[in] height : hauteur de l'image
			 * @param[in] width_step : taille réelle d'une ligne de l'image (si == 0 alors remplacée par width)
			 * @param[in] nb_bins : nombre d'états possibles pour un pixel (si == 0 alors remplacé par _nb_bins)
			 * @param[in] sigma : écart type de la gaussienne
			 * @param[in] c : taille de l'intervalle : [-M, M] avec M = c sigma + 1 (Il faut choisir entre 3 et 6)
			 */
			template <class type> void compute(	const type * data,
												const unsigned char * mask_data,
												unsigned int width,
												unsigned int height,
												unsigned int width_step = 0,
												unsigned int mask_width_step = 0,
												const double & sigma = 1,
												unsigned int c = 0);
			
			/**@fn void template <class type> int c_histogram :: stretch( 	type * data,
																			unsigned int width,
																			unsigned int height,
																			unsigned int width_step,
																			double eps );
			 * @param[in] data : données de l'image
			 * @param[in] width : largeur de l'image
			 * @param[in] height : hauteur de l'image
			 * @param[in] width_step : taille réelle d'une ligne de l'image (si == 0 alors remplacée par width)
			 * @param[in] eps : tolérance
			 * @brief
			 * Réhausse l'histogramme entre 0 et 255.
			 * 
			 **/
			template <class type> int stretch( 	type * data,
												unsigned int width,
												unsigned int height,
												unsigned int width_step,
												double eps_low,
												double eps_up );
			
			/**@fn
			 * @param[in] data : données de l'image
			 * @param[in] width : largeur de l'image
			 * @param[in] height : hauteur de l'image
			 * @param[in] width_step : taille réelle d'une ligne de l'image (si == 0 alors remplacée par width)
			 * @param[in] eps : tolérance
			 * @brief
			 * Réhausse l'histogramme entre 0 et 255.
			 */
			template<class type1, class type2> int stretch( 	type1 * data,
																const type2 * data_in,
																unsigned int width,
																unsigned int height,
																unsigned int width_step_in,
																unsigned int width_step_out,
																double eps_low,
																double eps_up );
								

			template<class type1, class type2> void equalize(	type2 * data_out,
																const type1 * data_in,
																unsigned int width,
																unsigned int height,
																unsigned int width_step_in,
																unsigned int width_step_out);

			template<class type> void equalize(	type * data,
												unsigned int width,
												unsigned int height,
												unsigned int width_step);
													
			template<class type1, class type2> int correct_alpha(	type2 * data_out,
																	const type1 * data_in,
																	unsigned int width,
																	unsigned int height,
																	unsigned int width_step_in,
																	unsigned int width_step_out,
																	double new_mean );

			template<class type> int correct_alpha(	type * data,
														unsigned int width,
														unsigned int height,
														unsigned int width_step,
														double new_mean );
			
			
													
			/**@fn void c_histogram :: search_maxima (	double * maxima,
														unsigned int & nb_max,
														unsigned int nb_max_max )
			 * @param[out] maxima : maxima
			 * @param[out] nb_max : nombre de maxima détectés
			 * @param nb_max_max : nombre maximum de maxima
			 * @brief
			 * Cette méthode cherche les maxima de l'histogramme.
			 */
			void search_maxima (	unsigned int * maxima,
									unsigned int & nb_max,
									unsigned int nb_max_max ) const;
			
			/**@fn void c_histogram :: search_minima (	double * minima,
														unsigned int & nb_min,
														unsigned int nb_min_max )
			 * @param[out]  minima :  minima
			 * @param[out] nb_min : nombre de minima détectés
			 * @param nb_min_max : nombre maximum de minima
			 * @brief
			 * Cette méthode cherche les maxima de l'histogramme.
			 */
			void search_minima (	unsigned int * minima,
									unsigned int & nb_min,
									unsigned int nb_min_max ) const;
			
			/**@fn c_histogram :: ~c_histogram();
			 * @brief
			 * Destructeur
			 */
			~c_histogram();
		
			/**@fn
			 * @brief
			 * Retourne la moyenne
			 * 
			 **/
			double get_mean( ) const;
			
			
			/**@fn
			 * @brief
			 * Retourne la moyenne
			 * 
			 **/
			double get_variance( ) const;
		
		
			/**@fn
			 * @brief
			 * Retourne la moyenne
			 * 
			 **/
			double get_moment( unsigned int p ) const;
		
		
			double get_percentile( unsigned int n ) const;
		
		
		
			//Accesseurs
			/**@fn inline const unsigned int * c_histogram :: data() const
			 * @return données de l'histogramme
			 */
			inline const double * data() const
			{
				return _data;
			}
		
			/**@fn inline const unsigned int * c_histogram :: s_data() const
			 * @return données de l'histogramme filtré
			 */
			inline const double * s_data() const
			{
				return _s_data;
			}
		
			/**@fn inline const unsigned int * c_histogram :: d_data() const
			 * @return Dérivée de l'histogramme
			 */
			inline const double * d_data() const
			{
				return _d_data;
			}
		
			/**@fn inline const unsigned int * c_histogram :: d2_data() const
			 * @return Dérivée de l'histogramme
			 */
			inline const double * d2_data() const
			{
				return _d_d_data;
			}
		
			/**@fn inline unsigned int c_histogram :: nb_bins() const
			 * @return
			 * nombre d'états
			 */
			inline unsigned int nb_bins() const
			{
				return _nb_bins;
			}
		
			/**@fn inline unsigned int c_histogram :: nb_pixels() const
			 * @return
			 * Nombre de pixels
			**/
			inline unsigned int nb_pixels() const
			{
				return _nb_pixels;
			}
		protected:
			/**@fn void c_histogram :: initialize();
			 * @brief
			 * TOut à 0.
			 */
			void initialize();
			
			/**@fn void c_histogram :: free();
			 * @brief
			 * Lib. mémoire
			 **/
			void free();
		
			/**@fn void c_histogram :: compute_histogram (	const unsigned char * data,
															unsigned int width,
															unsigned int height,
															unsigned int width_step = 0,
															unsigned int nb_bins = 0 );
			 * @param[in] data : données de l'image
			 * @param[in] width : largeur de l'image
			 * @param[in] height : hauteur de l'image
			 * @param[in] width_step : taille réelle d'une ligne de l'image (si == 0 alors remplacée par width)
			 * @param[in] nb_bins : nombre d'états possibles pour un pixel (si == 0 alors remplacé par _nb_bins)
			 * @brief
			 * Cette méthode calcule l'histogramme d'une image.
			 */
			template <class type> void compute_histogram (	const type * data,
															unsigned int width,
															unsigned int height,
															unsigned int width_step = 0);
			
			
			
			
			
			/**@fn void c_histogram :: filter (	unsigned int c,
												const double & sigma );
			 * @param[in] c : taille de l'intervalle : [-M, M] avec M = c sigma + 1 (Il faut choisir entre 3 et 6)
			 * @param[in] sigma : écart type de la gaussienne
			 * @brief
			 * Cette méthode filtre l'histogramme par une gaussienne
			 */
			void filter (	unsigned int c,
							const double & sigma );
			
			/**@fn void c_histogram :: compute_derivates();
			 * @brief
			 * Calcul des dérivées.
			 */
			void compute_derivates();
			
			
		
		
		
			unsigned int _nb_pixels;
			unsigned int _nb_bins;
			
			double * _data; // nb_bins array
			double * _s_data;
			double * _d_data; // nb_bins - 1 array
			double * _d_d_data; // nb_bins - 2 array
	};
	
#endif
