/**@file c_buffer_data.hpp
 * 
 */
#ifndef _C_BUFFER_DATA_HPP_
	#define _C_BUFFER_DATA_HPP_
	/**@class
	 * @brief
	 * Cette structure permet de gérer un buffer de données
	 */
	class c_buffer_data
	{
		public:
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 */
			c_buffer_data( );
			
			
			/**@fn
			 * @brief
			 * Constructeur
			 * 
			 */
			c_buffer_data( 	const void * params,
							void * (*alloc_function) ( const void * params ),
							void (*copy_function) ( void * data,
													const void * data_src ),
							void (*free_function) ( void * data ) );
			
			/**@fn
			 * @brief
			 * Setup.
			 * 
			 **/
			int setup ( 	const void * params,
							void * (*alloc_function) ( const void * params ),
							void (*copy_function) ( void * data,
													const void * data_src ),
							void (*free_function) ( void * data ) );
			
			/**@fn
			 * @brief
			 * Destructeur
			 * 
			 */
			~c_buffer_data();
		
			/**@fn
			 * @brief
			 * Set the attributes of the object
			 * 
			 **/
			void set( unsigned int id,
					 const double & score,
					 const void * data );
					 
			void set(	const c_buffer_data & data );
			
			
			/**@fn
			 * @brief
			 * Get the attributes.
			 * 
			 **/
			void get( 	unsigned int & id,
						double & score,
						void * data ) const;
			
			void get(	c_buffer_data & data ) const;
			
			/**@fn
			 * @brief
			 * Op. d'ordre
			 * 
			 **/
			bool operator< ( c_buffer_data & buffer_data ) const;
		
			/**@fn
			 * @brief
			 * Id de l'objet
			 * 
			 **/
			inline unsigned int id() const
			{
				return _id;
			}
			
			/**@fn
			 * @brief
			 * Efface le score et l'id.
			 * 
			 **/
			void erase();
		
			/**@fn
			 * @brief
			 * Score de l'objet
			 * 
			 **/
			inline double score() const
			{
				return _score;
			}
			
			/**@fn
			 * @brief
			 * Données de l'objet
			 * 
			 **/
			inline void * data()
			{
				return _data;
			}
			
			/**@fn
			 * @brief
			 * Données de l'objet
			 * 
			 **/
			inline void * data_const() const
			{
				return _data;
			}
		
		protected:
			/**@fn
			 */
			void initialize();
			
			/**@fn
			 */
			void free();
			
			void (*free_function) ( void * data );
			void (*copy_function) ( void * data_tg,
									const void * data_src );
			
			void * _data;
			double _score;
			unsigned int _id;
	};
	
#endif
