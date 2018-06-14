#include "lib_api.hpp"
#include <cstring>
#include <fstream>

using namespace std;
void normalize_name( char * name )
{
	for ( unsigned int i = 0; name[i] != '\0'; ++ i )
	{
		if ( name[i] >= 'a' && name[i] <= 'z' )
			name[i] = name[i] + 'A' - 'a';
		
		if ( name[i] == '.' || name[i] == ':' )
			name[i] = '_';
	}
}

void write_matrix(  ofstream & file,
					const gsl_matrix * matrix )
{
	file << "var = gsl_matrix_calloc ( " << matrix->size1 << "," << matrix->size2 << ");\\" << endl;
	for ( unsigned int i = 0; i < matrix->size1; ++ i )
	{
		for ( unsigned int j = 0; j < matrix->size2; ++ j )
		{
			file << "var->data[" << i << " * var->tda + " << j << "] = " << matrix->data[i * matrix->tda + j] << ";";
		}
	}
	file << endl;
}

void write_nb ( ofstream & file,
				double nb )
{
	file << "var = "<< nb << ";" << endl;
}

void write_string ( ofstream & file,
					const char * string )
{
	file << "var = new char["<< strlen(string)+1 << "];\\" << endl;
	file << "memcpy( var," << "\"" << string << "\"," << ( strlen(string) + 1 ) * sizeof(char) <<");" << endl;
}






int main(  int argc, char ** argv ) 
{
	if  ( argc > 1 )
	{
		if ( ! strcmp( argv[1], "--help" ) || ! strcmp( argv[1], "-H" ) )
		{
cout << 
"                              HELP\n"
"\n"
" Arg[1] : hpp file\n"
" Arg[2 - N] : cfg files\n"
" Brief\n"
" Convert \"cfg\" files to hpp file\n";
return 0;
		}
	}
	
	if ( argc < 3 ) 
	{
		cout << "Error : Missing arguments" << endl;
		return 1;
	}
	
	//Parametres
	api_parameters params;
	for ( int i = 2 ; i < argc; ++ i )
	{
		params.load ( argv[i] );
	}
	
	ofstream file ( argv[1] );
	if ( !file )
	{
		cout << "Error : Can't open " << argv[1] << " !" << endl;
		return 1;
	}
	
	normalize_name( argv[1] );
	file << "#ifndef " << argv[1] << endl;
	file << "#define " << argv[1] << endl;
	cout <<  params.nb_variables() << endl;
	for ( unsigned int i = 0; i < params.nb_variables(); ++ i )
	{
		char * buffer;
		buffer = new char[ strlen( params.get_variables()[i]->get_name() ) + 1];
		memcpy ( buffer, params.get_variables()[i]->get_name(), strlen( params.get_variables()[i]->get_name() ) + 1 );
		normalize_name( buffer );
		
		file << "#define " << buffer << "(var)\\" << endl;
		if ( params.get_variables()[i]->is_matrix() )
			write_matrix( file, params.get_variables()[i]->get_matrix_const() );
		else if ( params.get_variables()[i]->is_float() )
			write_nb(file, params.get_variables()[i]->get_float() );
		else if ( params.get_variables()[i]->is_string() )
			write_string(file, params.get_variables()[i]->get_string_const() );
		
		
		delete[] buffer;
	}
	
	
	
	
	
	file << "#endif" << endl;
	file.close();
	return 0;
}
