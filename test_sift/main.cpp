/**
	@file   test_sift.cpp
	@author Wade Spires
	@date   2006/1/18
	@brief  Scale-Invariant Feature Transform (SIFT)
 */

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include <algorithm>
#include <iostream>
#include <list>
#include <string>
#include <vector>

// tools headers
#include "Progress_Bar.hpp"
#include "ws_tools.hpp"

// Image headers
#include "Edge_Detector.hpp"
#include "Gauss.hpp"
#include "Image.hpp"
#include "Image_Region.hpp"
#include "Image_Coordinate.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"
#include "Video_IO.hpp"

// local headers
#include "sift/Sift.hpp"

using std::cerr;  using std::cout;  using std::endl;
using std::list;
using std::max;
using std::sort;
using std::string;
using std::vector;

// my namespace
using namespace ws_img;

typedef Image::size_type   size_type;
typedef Image::value_type  value_type;
typedef Image::pixel_type  pixel_type;
typedef Image_Coordinate   Img_Coord;

// if > 0, then detailed diagnostics are written to the screen (defined in
// log_mesg.h)
extern unsigned verbose;

/**
	@brief Structure for setting command-line options.

	The following is a list of mandatory program options:
		-i img_dir Directory containing input images
		-o out_dir Directory to write output images

	The following is a list of optional program options:
 */
struct CL_Options
{
	CL_Options( int argc, char** argv )
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
					// set '-i img_dir' option
					case 'i':
						++argv;
						if( --argc == 0 )
						{
							print_usage();
						}

						img_dir = *argv;
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

					// set '-v' option
					case 'v':
						verbose = 1;  // set global verbose flag
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
		if( out_dir == "" )
		{
			print_usage();
		}
	}

	/**
		Print error message showing program usage.
	 */
	void print_usage( )
	{
		cerr << "usage: sift -i img_dir -o out_dir [ -v ]" << endl;
		exit( EXIT_FAILURE );
	}

	string img_dir;
	string out_dir;
};

int main( int argc, char** argv )
{
	// process command-line options
	CL_Options options( argc, argv );
	string img_dir = options.img_dir;
	string out_dir = options.out_dir;

	// read list of files (images and videos) to process
	vector<string> file_names = ws_tools::dir_traverse( img_dir );
	sort( file_names.begin(), file_names.end() );

	// ensure output directory exists (and create it if it doesn't)
	ws_tools::check_dir( out_dir );

	Progress_Bar bar( file_names.size() );
	for( unsigned k = 0; k != file_names.size(); ++k )
	{
		// read in image and convert to gray-scale
		if( verbose > 0 )
		{
			fprintf( stdout, "Processing '%s'\n", file_names[k].c_str() );
		}

		// create base name for outputting results
		string base_name = out_dir
				+ ws_tools::get_base_name( ws_tools::get_file_name( file_names[k] ) );

		// process single image
		if( Image::is_image( file_names[k] ) )
		{
			Image img( file_names[k] );

			// run SIFT on image
			Keypoint_List keypoints = sift( img );

			keypoints.write( base_name + ".key" );

			// write image with keypoints drawn on it
			draw_keypoints( img, keypoints );
			string img_name = base_name + ".jpg";
			if( verbose > 0 )
			{
				fprintf( stdout, "Writing '%s'\n", img_name.c_str() );
				fprintf( stderr, "%s", k == (file_names.size() - 1) ? "" : "\n" );
			}
			else
			{
				bar.update( k );
			}
			img.write( img_name );
		}

		// process video
		else if( Image::is_video( file_names[k] ) )
		{
			Image frame;  // single video frame

			// process each frame of gray-scale video
			Video_IO video( file_names[k] );
			while( video.next_frame( frame ) )
			{
				string frame_name = base_name
						+ ws_tools::format_string( "_%05d", video.get_frame_num() );

				// run SIFT on frame
				Keypoint_List keypoints = sift( frame );

				keypoints.write( frame_name + ".key" );

				// write image with keypoints drawn on it
				draw_keypoints( frame, keypoints );
				string img_name = frame_name + ".jpg";
				if( verbose > 0 )
				{
					fprintf( stdout, "Writing '%s'\n", img_name.c_str() );
					fprintf( stderr, "%s", k == (file_names.size() - 1) ?
							"" : "\n" );
				}
				frame.write( img_name );
			}
		}
	}

	if( verbose == 0 )
	{
		fprintf( stderr, "\n" );
	}

	return( EXIT_SUCCESS );
}
