#include "c_buffer.hpp"
#include <unistd.h>
c_buffer :: c_buffer ( )
{
	initialize();
}

void c_buffer ::reset()
{
	start_score_editting();
	for ( unsigned int i = 0; i < _nb_objects; ++ i )
		erase_object( i );
	end_score_editting();
}



c_buffer :: c_buffer ( unsigned int nb_obj,
						void * (*alloc_function) ( const void * params ),
						void (*copy_function) ( void * fait_chier,
												const void * data ),
						void (*free_function) ( void * data ),
						const void * params,
						ostream * _err_stream )
{
	initialize();
	err_stream = _err_stream;
	if ( setup ( nb_obj,
				 alloc_function,
				 copy_function,
				 free_function,
				 params ) )
		throw ( invalid_argument( "Argument(s) of c_buffer") );
}

int c_buffer :: setup ( 	unsigned int nb_obj,
							void * (*alloc_function) ( const void * params ),
							void (*copy_function) ( void * fait_chier,
													const void * data ),
							void (*free_function) ( void * data ),
							const void * params,
							ostream * _err_stream )
{
	free();
	initialize();
	_err_stream = err_stream;
	if ( 	nb_obj == 0 		||
			alloc_function == 0	||
			free_function == 0	||
			copy_function == 0	 )
	{
		if ( err_stream )
			*err_stream << "Error : Argument(s) of c_buffer :: setup!" << endl;
		return 1;
	}
	
	_nb_objects = nb_obj;
	buffer_data = new c_buffer_data[nb_obj];
	
	for ( unsigned int i = 0; i < nb_obj; ++ i )
	{
		if ( buffer_data[i].setup( 	params, 
									alloc_function, 
									copy_function, 
									free_function ) )
		{
			if ( err_stream )
				*err_stream << "Error: Invalid parameters for buffer_data!" << endl;
			return 1;
		}
	}
	
	sorted_buffer_ids = new unsigned int[nb_obj];
	for ( unsigned int i = 0; i < nb_obj; ++ i )
	{
		sorted_buffer_ids[i] = i;
	}
	
	_object_in_reading = new int[nb_obj];
	_object_in_writting = new int[nb_obj];
	memset(	_object_in_reading, 
			0, 
			sizeof(int) * nb_obj );
	memset(	_object_in_writting, 
			0, 
			sizeof(int) * nb_obj );
			
	return 0;
}


int c_buffer :: add_object(	const c_buffer_data * data )
{
	return add_object( 	data->id(), 
						data->score(),
						data->data_const() );
}

int c_buffer :: add_object( 	unsigned int _id,
								double _score,
								const void * _data )
{
	unsigned int id, i;
	
	//Indique au buffer que le score va être éditer
	start_score_editting();
	//Recherche la position de l'objet
	for ( i = 0; i < _nb_objects; ++ i )
	{
		id = sorted_buffer_ids[i];
		double score = buffer_data[ id ].score();
		if ( _score > score )
		{

			id = sorted_buffer_ids[_nb_objects - 1];
			//Tri des scores
			for ( unsigned int j = ( _nb_objects - 2 ); j != (i - 1) ; -- j )
			{
				sorted_buffer_ids[ j + 1 ] = sorted_buffer_ids[ j ];
			}
			sorted_buffer_ids[ i ] = id;
			//Ecriture de l'objet
			write_object( 	_id, 
							_score, 
							_data, 
							id );
			break;
		}	
	}	
	
	//Déblocage du score
	end_score_editting();
	if ( i == _nb_objects )
		return 1;
	return 0;
}

int c_buffer :: delete_object( const c_buffer_data * data )
{

	//cout << "delete" << endl;
	start_score_editting();
	for ( unsigned int i = 0; i < _nb_objects; ++ i )
	{
		unsigned int k;
		unsigned int id;
			id = sorted_buffer_ids[i];
		//Lecture du score de l'objet
		k = buffer_data[ id ].id();
		if ( data->id() == k )
		{
			//Delete

			erase_object( id );

			//Tri des scores

			for ( unsigned int j = i; j < ( _nb_objects - 1 ); ++ j )
			{
				sorted_buffer_ids[ j ] = sorted_buffer_ids[ j + 1 ];
			}
			sorted_buffer_ids[_nb_objects - 1] = id;

			break;
			
		}
	}

	end_score_editting();
	//cout << "delete - ok" << endl;
	return 0;
}

int c_buffer :: get_object (	c_buffer_data * data,
								unsigned int score_id )
{
	read_object(	data, 
					sorted_buffer_ids[score_id] );
	return 0;
}
int c_buffer :: get_object_with_id (	c_buffer_data * data,
										unsigned int obj_id )
{
	bool find = false;
	for ( unsigned int i = 0; ( i < _nb_objects ) && ( !find ); ++ i )
	{
		start_object_reading( i );
			if ( obj_id == buffer_data[i].id() )
			{
				find = true;
				read_object(	data, 
								i );
			}
		end_object_reading( i );
	}
	return find;
}
int c_buffer :: inclued ( const c_buffer_data * data )
{
	bool find = false;
	for ( unsigned int i = 0; i < _nb_objects && !find; ++ i )
	{
		start_object_reading( i );
			if ( data->id() == buffer_data[i].id() && buffer_data[i].score() > 0 )
				find = true;
		end_object_reading( i );
	}
	return find;
}


c_buffer :: ~c_buffer()
{
	free();
	initialize();
}

void c_buffer :: free()
{
	if ( buffer_data )
		delete[] buffer_data;

	if ( sorted_buffer_ids )
		delete[] sorted_buffer_ids;
	if ( _object_in_reading )
		delete[] _object_in_reading;
	if ( _object_in_writting )
		delete[] _object_in_writting;

}

void c_buffer :: initialize()
{
	_object_in_reading = 0;
	_object_in_writting = 0;
	_score_editting = 0;
		
	_nb_objects = 0;
			
	buffer_data = 0;
	sorted_buffer_ids = 0;

	err_stream = NULL;
}

double c_buffer :: score_min()
{
	double score;
	start_score_editting();
		score = buffer_data[sorted_buffer_ids[_nb_objects - 1]].score();
	end_score_editting();
	return score;
}

void c_buffer :: erase_object( unsigned int id )
{
	start_object_writting( id );
		buffer_data[id].erase();
	end_object_writting( id );
}

void c_buffer :: write_object(	unsigned int id,
								double score,
								const void * data,
								unsigned int obj_id )
{
	start_object_writting( obj_id );
		buffer_data[obj_id].set( id, score, data);
	end_object_writting( obj_id );
}

void c_buffer :: read_object(	c_buffer_data * data,
								unsigned int id )
{
	start_object_reading( id );
		buffer_data[id].get( *data );
	end_object_reading( id );
}


void c_buffer :: start_object_writting ( unsigned int obj_id )
{	
	while ( _object_in_reading[obj_id] 	||
			_object_in_writting[obj_id] )
	{
		usleep( 1000 );
	}
	
	_object_in_writting[obj_id] = 1;
}

void c_buffer :: end_object_writting ( unsigned int obj_id )
{	
	_object_in_writting[obj_id] = 0;
}

void c_buffer :: start_object_reading ( unsigned int obj_id )
{
	while ( _object_in_writting[obj_id] )
	{
		usleep( 1000 );
	}
	_object_in_reading[obj_id] ++;	
}

void c_buffer :: end_object_reading ( unsigned int obj_id )
{	
	_object_in_reading[obj_id] --;
			
}

void c_buffer :: start_score_editting ()
{
	while ( _score_editting )
	{
		usleep( 1000 );
	}

	_score_editting = 1;
}

void c_buffer :: end_score_editting ()
{
	_score_editting = 0;
}











