/**@file clabel.hpp 
 */
#ifndef _CLABEL_HPP_
	#define _CLABEL_HPP_
	#define LABEL_FUNCTION(type) \
	template int c_label :: label( 	const type * data,\
										unsigned int width_step,\
										unsigned char bg_color,\
										unsigned int width,\
										unsigned int height);\
	template int c_label :: copy_region ( type * data, \
											unsigned int region_id,\
											unsigned int width_step,\
											type color,\
											unsigned int width,\
											unsigned int height);
	
	/**@class c_label
	 * @brief
	 * Cette classe permet l'étiquettage des régions connexes d'un masque binaire.
	 * 
	 **/
	class c_label
	{
		public:
			// Constructeurs
			/**@fn c_label :: c_label(unsigned int width, 
									  unsigned int height);
			 * @param width : largeur
			 * @param height : hauteur
			 * @brief
			 * Constructeur 
			**/
			c_label(unsigned int width, 
				    unsigned int height);
				    	
			/**@fn c_label :: ~c_label();
			 * @brief
			 * Destructeur
			 */
			~c_label();
		
			//Setup
			/**@fn void c_label :: setup(unsigned int width, 
										 unsigned int height);
			 * @param width : largeur
			 * @param height : hauteur
			 * @brief
			 * Setup.
			 **/
			void setup(unsigned int width, 
					   unsigned int height);
		
			//Label
			/**@fn template <class type> c_label :: label( const type * data, 
														    unsigned int width_step = 0, 
															unsigned char bg_color = 0,
															unsigned int width = 0,
															unsigned int height = 0);
			 * @param data : données à étiquetter
			 * Doit être un masque binaire (0, a)
			 * @param width_step : taille des lignes des données
			 * @param bg_color : couleur du fond
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @brief
			 * Cette méthode étiquete les différentes régions connexes du masque
			 */
			template <class type> int label(	const type * data, 
												unsigned int width_step = 0, 
												unsigned char bg_color = 0,
												unsigned int width = 0,
												unsigned int height = 0);
			
			//Recopy
			
			
			/**@fn template <class type> int c_label :: copy_region ( 	type * data, 
																		unsigned int region_id,
																		unsigned int width_step = 0, 
																		unsigned char color = 0,
																		unsigned int width = 0,
																		unsigned int height = 0);
			 * @param data : données à étiquetter
			 * Doit être un masque binaire (0, a)
			 * @param width_step : taille des lignes des données
			 * @param bg_color : couleur du fond
			 * @param width : largeur de l'image
			 * @param height : hauteur de l'image
			 * @brief
			 * Cette méthode étiquete les différentes régions connexes du masque
			 */
			template <class type> int copy_region ( 	type * data, 
														unsigned int region_id,
														unsigned int width_step = 0, 
														type color = 255,
														unsigned int width = 0,
														unsigned int height = 0);
			
			
			
			//Accesseurs
			/**@fn inline unsigned int c_label :: width_step() const
			 * @return Taille réelle d'une ligne
			 */
			inline unsigned int width_step() const
			{
				return _width_max;
			}
			
			/**@fn inline unsigned int c_label :: height_max() const
			 * @return Nombre de ligne max géré par l'objet.
			 */
			inline unsigned int height_max() const
			{
				return _height_max;
			}
			
			
			/**@fn inline unsigned int c_label :: width() const
			 * @return
			 * Largeur
			 * 
			 **/
			inline unsigned int width() const
			{
				return _width;
			}
			/**@fn inline unsigned int c_label :: height() const
			 * @return
			 * Hauteur
			 **/
			inline unsigned int height() const
			{
				return _height;
			}
			
			/**@fn inline unsigned int c_label :: nb_pixels() const
			 * @return
			 * Nombre de pixels
			 **/
			inline unsigned int nb_pixels() const
			{
				return _nb_pixels;
			}
			
			/**@fn inline const unsigned int * c_label :: label_map() const
			 * @return
			 * Carte des régions connexes
			 */
			inline const unsigned int * label_map() const
			{
				return _label_map;
			}
		
			/**@fn inline const unsigned int * c_label :: surfaces() const
			 * @return
			 * Aire des différentes régions connexes.
			 */
			inline const unsigned int * surfaces() const
			{
				return _surfaces;
			}
		
			/**@fn inline unsigned int c_label :: nb_regions() const
			 * @return
			 * Nombre de régions connexes
			 */
			inline unsigned int nb_regions() const
			{
				return _nb_regions;
			}
		
			/**@fn int get_bounding_box(unsigned int & x,
										unsigned int & y,
										unsigned int & w,
										unsigned int & h
										unsigned int region_id) const;
			 * @param[out] x : abscisse du coin HG
			 * @param[out] y : ordonnée du coin HG
			 * @param[out] w : largeur de la boite englobante
			 * @param[out] h : hauteur de la boite englobante
			 * @param[in] region_id : numéro de la région connexe
			 * @return
			 * - 0 si region_id < nb_regions
			 * - 1 sinon
			 * @brief
			 * Cette fonction renvoie la taille de la boite englobante de la région connexe n° id
			 */
			int get_bounding_box(unsigned int & x,
								 unsigned int & y,
								 unsigned int & w,
								 unsigned int & h,
								 unsigned int region_id) const;
			
			//Methodes
			/**@fn int c_label :: erase_region(unsigned int region_id);
			 * @param[in] region_id : numéro de la région connexe
			 * @return 
			 * - 0 si region_id < nb_regions
			 * - 1 sinon
			 * @brief
			 * Cette fonction supprime la région connexe region_id.
			 */
			int erase_region(unsigned int region_id);
			
			/**@fn int c_label :: erase_regions(const unsigned int * region_ids,
												unsigned int nb_regions);
			 * @param[in] region_ids : numéros des région connexe
			 * @param[in] nb_regions : nombre de régions à supprimers
			 * @return 
			 * - 0 si region_ids < nb_regions
			 * - 1 sinon
			 * @brief
			 * Cette fonction supprime les régions connexe region_ids.
			 */
			int erase_regions(const unsigned int * region_ids,
							  unsigned int nb_regions);
			
			/**@fn unsigned int c_label :: erase_small_component(unsigned int threshold_surface);
			 * @param : threshold_surface : seuil
			 * @return
			 * - Nombre de régions supprimées
			 * @brief
			 * Cette méthode supprime les régions connexes de surface inférieure à threshold_surface
			 */
			unsigned int erase_small_component(unsigned int threshold_surface);
			
			/**@fn unsigned int c_label :: keep_biggest_component();
			 * @return
			 * - Nombre de régions supprimées
			 * @brief
			 * Cette méthode garde la région connexe la plus grande en surface et supprime tous les autres.
			**/
			unsigned int keep_biggest_component();
			
			
			/**@fn unsigned int c_label :: get_biggest_region_id() const;
			 * @return
			 * Plus grande région.
			 */
			unsigned int get_biggest_region_id() const;
		protected:
			/**@fn void c_label :: free();
			 * @brief
			 * Libère la mémoire utilisée par les attributs de l'objet
			 */
			void free();
			
			/**@fn void c_label :: initialize();
			 * @brief
			 * Mets tous les attributs de l'objet à 0.
			 */
			void initialize();
		
			// variables
			unsigned int _width_max;
			unsigned int _height_max;
			unsigned int _width;
			unsigned int _height;
			unsigned int _nb_pixels;
			unsigned int _nb_pixels_max;
			unsigned int * _label_map;

			unsigned int _nb_regions;
			unsigned int * _surfaces; 
			
			unsigned int * x_min,
						 * y_min,
						 * x_max,
						 * y_max;
			//tmp
			unsigned * labels;
			unsigned * labels_bis;
};
	
#endif
