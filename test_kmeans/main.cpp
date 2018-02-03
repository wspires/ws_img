/**
	@file   main.cpp
	@author Wade Spires
	@date   2006/10/20
	@brief  Program kmeans.

	Copyright (c) 2006 Wade Spires. All rights reserved.
 */

// c++ headers
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

// tools headers
#include "ws_tools.hpp"

// image headers
#include "ws_img.hpp"

#include "K_Means.hpp"

using std::cerr;  using std::cout;  using std::endl;
using std::sort;
using std::string;
using std::vector;

using namespace ws_tools;
using namespace ws_img;
typedef Image::size_type  size_type;
typedef Image::value_type value_type;

/**
	@brief Structure for setting command-line options.

	The following is a list of mandatory program options:
		'-i img_dir' Directory to read input images from.

	The following is a list of optional program options:
		'-h' Display usage information.
		'-o out_dir' Directory to write output images to.
 */
struct CL_Options
{
	CL_Options( int argc, char** argv )
	: prog_name( argv[0] ), out_dir( "out" ), K(5)
	{
		// set each option from command-line
		while( --argc != 0 )
		{
			++argv;
			if( argv == NULL )
			{
				print_usage();
			}
			string arg = *argv;

			// by convention, options start with a dash
			if( arg[0] == '-' )
			{
				if( arg.size() == 1 )
				{
					print_usage();
				}

				// handle each option type
				switch( arg[1] )
				{
					// set '-h' option
					case 'h':
						print_usage();
						break;

					// set '-i img_dir' option
					case 'i':
						++argv;
						if( --argc == 0 )
						{
							print_usage();
						}
						img_dir = *argv;
						break;

					// set '-k num_means' option
					case 'k':
						++argv;
						if( --argc == 0 )
						{
							print_usage();
						}
						K = atoi( *argv );
						break;

					// set '-o out_dir' option
					case 'o':
						++argv;
						if( --argc == 0 )
						{
							print_usage();
						}
						out_dir = *argv;
						break;

					default:
						print_usage();
						break;
				}
			}
		}
		if( img_dir == "" )
		{
			print_usage();
		}
		if( K <= 0 )
		{
			print_usage();
		}
	}

	/**
		Print error message showing program usage.
	 */
	void print_usage( )
	{
		cerr << "usage: " << prog_name
			<< " -i img_dir [ -o out_dir -K num_means ]"
			<< endl;
		exit( EXIT_FAILURE );
	}

	string prog_name;  //< Name of program
	string img_dir;    //< Input directory
	string out_dir;    //< Output directory
	int    K;          //< Number of means
};

int main( int argc, char** argv )
{
	// process command-line arguments
	CL_Options options( argc, argv );
	string img_dir = options.img_dir;
	string out_dir = options.out_dir;
	unsigned K     = options.K;

	// read list of images to process
	vector<string> image_names = dir_traverse( img_dir, Image::is_image );
	sort( image_names.begin(), image_names.end() );

	// ensure output directory exists (and create it if it doesn't)
	check_dir( out_dir );

	// allocate space for the data
	Matrix data;
	if( !image_names.empty() )
	{
		Image img( image_names[0] );
		data = Matrix( image_names.size(),
				img.sub_row_size() * img.sub_col_size() );
	}
fprintf( stderr, "%s\n", data.get_region().to_string().c_str() );

	// don't try to form more clusters than we have data items
	K = std::min<unsigned>( K, image_names.size() );

	// copy each image n into row n of the data matrix
	for( unsigned n = 0; n != image_names.size(); ++n )
	{
		// read in image
		Image img( image_names[n] );

		// form a row of data from the image
		size_type dj = data.col_begin();
		for( size_type i = img.row_begin(); i != img.row_end(); ++i )
		{
			for( size_type j = img.col_begin(); j != img.col_end(); ++j )
			{
				data(n,dj) = img(i,j);
				++dj;
			}
		}
	}

	// perform K-means clustering
	Matrix means = K_Means::k_means( K, data );

	// write result to output directory
	string out_name = out_dir + format_string( "%u_means.m", K );
	means.write( out_name );

	return( EXIT_SUCCESS );
}
