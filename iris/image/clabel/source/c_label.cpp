#include "c_label.hpp"
#include <iostream>
#include <cstring>
using namespace std;
c_label :: c_label(unsigned int width, 
				   unsigned int height)
{
	initialize();
	setup(width, height);
}

c_label :: ~c_label()
{
	free();
	initialize();
}

void c_label :: setup(unsigned int width, 
					  unsigned int height)
{

	free();

	initialize();
	_width_max = width;
	_height_max = height;
	_nb_pixels_max = _width_max * _height_max;
	_nb_pixels = 0;
	//Alloc.
	_label_map = new unsigned int[_nb_pixels_max];
	memset ( 	_label_map,
				0,
				_nb_pixels_max * sizeof( unsigned int ) ); //Pour valgrin tout puissant!
	
	
	_nb_regions = 0;
	memset(_label_map, 0, _nb_pixels * sizeof(unsigned int));
	_surfaces = new unsigned int[_nb_pixels_max];
	memset ( 	_surfaces,
				0,
				_nb_pixels_max * sizeof( unsigned int ) ); //Pour valgrin tout puissant!
	
	
	
	
	
	
	_surfaces[0] = _nb_pixels;
	x_min = new unsigned int[_nb_pixels_max];
	memset ( 	x_min,
				0,
				_nb_pixels_max * sizeof( unsigned int ) ); //Pour valgrin tout puissant!

	y_min = new unsigned int[_nb_pixels_max];
	memset ( 	y_min,
				0,
				_nb_pixels_max * sizeof( unsigned int ) ); //Pour valgrin tout puissant!
	
	x_max = new unsigned int[_nb_pixels_max];
	memset ( 	x_max,
				0,
				_nb_pixels_max * sizeof( unsigned int ) ); //Pour valgrin tout puissant!
	y_max = new unsigned int[_nb_pixels_max];
	memset ( 	y_max,
				0,
				_nb_pixels_max * sizeof( unsigned int ) ); //Pour valgrin tout puissant!
	labels = new unsigned int[_nb_pixels_max];
	memset ( 	labels,
				0,
				_nb_pixels_max * sizeof( unsigned int ) ); //Pour valgrin tout puissant!
	
	labels_bis = new unsigned int[_nb_pixels_max];
	memset ( 	labels_bis,
				0,
				_nb_pixels_max * sizeof( unsigned int ) ); //Pour valgrin tout puissant!
	
	
}

int c_label :: get_bounding_box(unsigned int & x,
								unsigned int & y,
								unsigned int & w,
								unsigned int & h,
								unsigned int region_id) const
{
	if (region_id >= _nb_regions)
	{
		return 1;
	}
	else if (_surfaces[region_id] == 0)
	{
		return 1;
	}
	else
	{
		x = x_min[region_id];
		y = y_min[region_id];
		w = x_max[region_id] - x_min[region_id] + 1;
		h = y_max[region_id] - y_min[region_id] + 1;
		return 0;
	}
}

void c_label :: free()
{
	if (_surfaces)
		delete[] _surfaces;
	if (_label_map)
		delete[] _label_map;
	if (x_min)
		delete[] x_min;
	if (y_min)
		delete[] y_min;
	if (x_max)
		delete[] x_max;
	if (y_max)
		delete[] y_max;

	if (labels)
		delete[] labels;
	if (labels_bis)
		delete[] labels_bis;
}

void c_label :: initialize()
{
	_width = 0;
	_width_max = 0;
	_height = 0;
	_height_max = 0;
	_nb_pixels = 0;
	_nb_pixels_max = 0;
	_label_map = 0;
	_nb_regions = 0;
	_surfaces = 0;
	x_min = 0;
	y_min = 0;
	x_max = 0;
	y_max = 0;
	labels = 0;
	labels_bis = 0;
}

template <class type> int c_label :: label( const type * data, 
											unsigned int width_step, 
											unsigned char bg_color,
											unsigned int width,
											unsigned int height)
{

	if ( ! width )
		width = _width_max;
	if (! height)
		height = _height_max;
	if ( !width_step )
		width_step = width;
	
	if (width > _width_max)
		return -1;
	
	if (height > _height_max)
		return -1;
	
	//Nombre de pixels
	_nb_pixels = width * height;	
	_width = width;
	_height = height;
	_surfaces[0] = 0;
	//Mise à zéro
	for (unsigned int i = 0; i < height; ++ i)
	{
		memset(_label_map + i * _width_max, 0, sizeof(unsigned int) * width);
	}
	labels[0] = 0;
	_nb_regions = 0;
	//Segmentation
	//Pour chaque pixel
	unsigned int i, 
				 j, 
				 n_pixel = width_step,
				 n_pixel_2 = _width_max,
				 u_label,
				 l_label;
	//Première ligne
	{
		//Première collonne
		if (data[0] != bg_color) //On teste si le pixel appartient au fond
		{
			_nb_regions ++;
			labels[_nb_regions] = _nb_regions;
			//Etiquetage
			_label_map[0] = 1;
		}
		
		//Autres collones
		for(j = 1; j < _width; ++ j) 
		{
			if (data[j] != bg_color) //On teste si le pixel appartient au fond
			{
				l_label = _label_map[j - 1];
				if (l_label == 0) //Si le pixel n'est pas lié à une région
				{
					_nb_regions ++;
					labels[_nb_regions] = _nb_regions;
					//Etiquetage
					_label_map[j] = _nb_regions;
					
				}
				else
				{
					_label_map[j] = l_label;
				}	
			}
		}
	}
	//Autres lignes
	for (i = 1; i < _height; ++i)
	{
		n_pixel = i * width_step;
		n_pixel_2 = i * _width_max;
		//Première colonne
		if (data[n_pixel] != bg_color) //On teste si le pixel appartient au fond
		{
			u_label = _label_map[n_pixel_2 - _width_max];
			if (u_label == 0) //Si le pixel n'est pas lié à une région
			{
				_nb_regions ++;
				labels[_nb_regions] = _nb_regions;
				//Etiquetage
				_label_map[n_pixel_2] = _nb_regions;
				
			}
			else
			{
				_label_map[n_pixel_2] = u_label;
			}	
		}
		
		
		//Autres colonnes
		for(j = 1; j < _width; ++ j) 
		{
			n_pixel ++;
			n_pixel_2 ++;
			if (data[n_pixel] != bg_color) //On teste si le pixel appartient au fond
			{
				l_label = _label_map[n_pixel_2 - 1];
				u_label = _label_map[n_pixel_2 - _width_max];
				if (l_label == 0 && u_label == 0) //Si le pixel n'est pas lié à une région
				{

					_nb_regions ++;
					labels[_nb_regions] = _nb_regions;
					//Etiquetage
					_label_map[n_pixel_2] = _nb_regions;
					
				}
				else if (l_label == 0) //u_label !=0 
				{
					_label_map[n_pixel_2] = u_label;
				}
				else if (u_label == 0) // l_label != 0
				{
					_label_map[n_pixel_2] = l_label;	
				}
				else //Cas merdique
				{
					_label_map[n_pixel_2] = u_label;
						
					//liaison des deux régions
					unsigned uu = u_label,
							 ll = l_label;
					while(labels[ll] != ll)
					{
						ll = labels[ll];
					}
					
					while(labels[uu] != uu) //Recherche de la région de base
					{
						uu = labels[uu];
					}
					
					if (uu > ll) //Cela signifie que les régions ne sont pas liées
					{
						labels[uu] = ll; //On lie les deux régions
						labels[u_label] = ll;
						labels[l_label] = ll;
					}
					else if (uu < ll)
					{
						labels[ll] = uu; //On lie les deux régions
						labels[u_label] = uu;
						labels[l_label] = uu;
					}
					else
					{
						labels[u_label] = uu;
						labels[l_label] = uu;
					}
				}	
			}
		}
	}
	_nb_regions ++;
    //Permet de compacter le nombre de régions
    unsigned nb_regions_2 = 0;
    for( i = 0; i < _nb_regions; ++i )
    {
        if (labels[i] == i)
        {
            labels_bis[i] = (nb_regions_2 ++);
            _surfaces[nb_regions_2] = 0;
            x_min[nb_regions_2] = _width - 1;
            y_min[nb_regions_2] = _height - 1;
            x_max[nb_regions_2] = 0;
            y_max[nb_regions_2] = 0;
        }
        else
        {
			unsigned int _label = labels[i];
			while (labels[_label] != _label)
			{
				_label = labels[ _label ];
			}
			labels_bis[i] = labels_bis[_label];
		}
    }

    //Compactage des régions
    _nb_regions = nb_regions_2;
	n_pixel_2 = 0;
    for(i = 0; i < _height; i ++)
    {
		n_pixel_2 = i * _width_max;
        for(j = 0; j < _width; j ++) //Pour chaque pixel
        {
			unsigned int _label = labels_bis[_label_map[ n_pixel_2 ] ];
			_label_map[n_pixel_2] = _label;
			_surfaces[ _label ] ++;
			if ( x_min[ _label ] > j)
				x_min[ _label ] = j;
			if ( x_max[ _label ] < j)
				x_max[ _label ] = j;
			
			if ( y_min[ _label ] > i)
				y_min[ _label ] = i;
			if ( y_max[ _label ] < i)
				y_max[ _label ] = i;
			n_pixel_2 ++;
			
        }

    }	

	return 0;
}

int c_label :: erase_region(unsigned int region_id)
{
	if (region_id >= _nb_regions)
	{
		return 1;
	}
	else if (_surfaces[region_id] == 0)
	{
		return 1;
	}
	else
	{
		unsigned int n_pixel = 0;
		for (unsigned int i = 0; i < _height; ++ i)
		{
			for (unsigned int j = 0; j < _width; ++ j)
			{
				n_pixel = i * _width_max + j;
				if (_label_map[n_pixel] == region_id)
				{
					_label_map[n_pixel] = 0;
				}
			}
		}
		_surfaces[0] += _surfaces[region_id];
		_surfaces[region_id] = 0;
		
		return 0;
	}
}

int c_label :: erase_regions(const unsigned int * region_ids,
							 unsigned int nb_regions)
{
	unsigned int count = 0;
	for (unsigned int i = 0; i < nb_regions; ++ i)
	{
		count += (! erase_region(i) );
	}
	return count;
}

unsigned int c_label :: erase_small_component(unsigned int threshold_surface)
{
	unsigned int count = 0;
	for (unsigned int i = 1; i < _nb_regions; ++ i)
	{
		if (_surfaces[i] != 0)
		{
			if (_surfaces[i] < threshold_surface)
				count += (! erase_region(i) );
		}
	}
	return count;
}

unsigned int c_label :: keep_biggest_component()
{
	unsigned int count = 0;
	if (_nb_regions == 1)
		return 0;
		
	unsigned int max = _surfaces[1];
	for (unsigned int i = 2; i < _nb_regions; ++ i)
	{
		if (_surfaces[i] >= max)
			max = _surfaces[i];
		else
			erase_region(i);
	}
	for (unsigned int i = 1; i < _nb_regions; ++ i)
	{
		if (_surfaces[i] < max)
			erase_region(i);
		else
			count = i;
	}
	
	return count;
}
unsigned int c_label :: get_biggest_region_id() const
{
	if (_nb_regions == 1)
		return 0;
		
	unsigned int max = _surfaces[1];
	unsigned int id = 1;
	for (unsigned int i = 2; i < _nb_regions; ++ i)
	{
		if (_surfaces[i] > max)
		{
			max = _surfaces[i];
			id = i;
		}
	}
	return id;
}



template <class type> int c_label :: copy_region ( 	type * data, 
														unsigned int region_id,
														unsigned int width_step, 
														type color,
														unsigned int width,
														unsigned int height)
{
	if ( width == 0 )
		width = _width;
	if ( height == 0 )
		height = _height;
	if ( width_step == 0 )
		width_step = width;
	
	if ( 	width != _width 	|| 
			height != _height )
		return 1;
		
	if ( region_id > _nb_regions )
		return 1;
		
	memset( data, 0, sizeof( type ) * width_step * height );	
	for ( unsigned int row = y_min[region_id]; row <= y_max[region_id]; ++ row )
	{
		for ( unsigned int col = x_min[region_id]; col < x_max[region_id]; ++ col )
		{
			if ( _label_map[ row * _width_max + col] == region_id )
				data[ row * width_step + col ] = color;
			
		}
	}
	
	
	return 0;
}







LABEL_FUNCTION(unsigned char)
LABEL_FUNCTION(char)

LABEL_FUNCTION(unsigned int)
