#include "c_buffer_data.hpp"
#include <iostream>
using namespace std;
c_buffer_data :: c_buffer_data( )
{
	initialize();
}

c_buffer_data :: c_buffer_data( 	const void * params,
									void * (*_alloc_function) ( const void * params ),
									void (*_copy_function) ( 	void * sdq,
																const void * data ),
									void (*_free_function) ( void * data ) )
{
	initialize();
	setup (	params, 
			_alloc_function,
			_copy_function,
			_free_function );
}

int c_buffer_data :: setup ( 	const void * params,
								void * (*_alloc_function) ( const void * params ),
								void (*_copy_function)  ( 	void * dsdqsdqa,
															const void * data ),
								void (*_free_function) ( void * data ) )
{
	free();
	initialize();
	
	if ( ! params || ! _alloc_function || ! _copy_function || !_free_function )
		return 1;

	
	//Alloc
	_data = _alloc_function ( params );
	copy_function = _copy_function;
	free_function = _free_function;
	
	return 0;
}

c_buffer_data :: ~c_buffer_data()
{
	free();
	initialize();
}

void c_buffer_data :: set( 	unsigned int id,
							const double & score,
							const void * data )
{
	_id = id;
	_score = score;
	copy_function ( _data, data );
}
void c_buffer_data :: set( const c_buffer_data & data )
{
	set ( 	data._id, 
			data._score, 
			data._data );
}
void c_buffer_data :: get( 	unsigned int & id,
							double & score,
							void * data ) const
{
	id = _id;
	score = _score;
	copy_function ( data, _data );
}

void c_buffer_data :: get( c_buffer_data & data ) const
{
	get ( 	data._id, 
			data._score, 
			data._data );
}



void c_buffer_data :: initialize()
{
	free_function = 0;
	copy_function = 0;
	_data = 0;
	_score = -1;
	_id = 0;
}

void c_buffer_data :: free()
{
	if ( _data )
		free_function( _data );
}

bool c_buffer_data :: operator< ( c_buffer_data & buffer_data ) const
{
	if ( _score < buffer_data._score )
		return true;
	else
		return false;
}

void c_buffer_data :: erase()
{
	_id = 0;
	_score = -1;
}
