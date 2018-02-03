/**
	@file   main.cpp
	@author Wade Spires
	@date   2005/11/21
	@brief  Test program for using matrix class.
 */

// c++ headers
#include <iostream>
#include <string>
#include <vector>

// my headers
#include "Wavelet.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"
#include "Linear_Algebra.hpp"
#include "Image.hpp"
#include "ws_tools.hpp"
#include "err_mesg.h"

using namespace ws_img;

using std::cerr;  using std::cout;  using std::endl;
using std::string;
using std::vector;

/**
	@brief Structure for setting command-line options.
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
					// set '-a' option
					case 'a':
						break;

					default:
						print_usage();
						break;
				}
			}
		}
	}

	/**
		Print error message showing program usage.
	 */
	void print_usage( )
	{
		cerr << "usage: test_wavelet" << endl;
		exit( EXIT_FAILURE );
	}
};

void test_matrix( );
void test_image( );

int main( int argc, char** argv )
{
	CL_Options options( argc, argv );

	// test_matrix();
	test_image();

	return( EXIT_SUCCESS );
}

void
test_matrix( )
{
	fprintf( stderr, "test_matrix\n" );
	fprintf( stderr, "---------------\n" );
	Matrix::size_type size = 8;

	Matrix m(size);
	for( Matrix::size_type i = m.row_begin(); i != m.row_end(); ++i )
	{
		for( Matrix::size_type j = m.col_begin(); j != m.col_end(); ++j )
		{
			m(i,j) = j * j * j;
			m(i,j) = i + j;
		}
	}
	// m = Matrix::random(size, 100);
	// m.set_region( 3, 3, 4, 4 );
	m.write();
	fprintf( stderr, "\n" );

	// Lin_Alg::Wavelet( m );

	// Wavelet wavelet( Wavelet::Haar );
	// Wavelet wavelet( Wavelet::Daub_4 );
	Wavelet wavelet( Wavelet::Bspline_103 );

	// wavelet.forward_transform( m );
	wavelet.forward_transform( m, Wavelet::Standard );
	// wavelet.forward_transform( m );

	string band_dir = "band_dir/";
	ws_tools::check_dir( band_dir );

	// get and print subbands
	unsigned num_levels = 1000;
	Wavelet_Subbands subbands = wavelet.get_subbands( m, num_levels );
	fprintf( stderr, "Done getting subbands\n\n" );
	Wavelet_Subbands::iterator iter     = subbands.begin();
	Wavelet_Subbands::iterator end_iter = subbands.end();
	for( ; iter != end_iter; ++iter )
	{
		const Wavelet_Level_Band& level_band = iter->first;
		const Wavelet_Subband&    subband    = iter->second;

		Image band_img = subband;

		string band_name = band_dir + level_band.to_string() + ".m";

		fprintf( stderr, "%s\n", level_band.to_string().c_str() );
		band_img.write( band_name );
	}

	fprintf( stderr, "Lifted signal:\n" );
	//signal.write( "lift.m", 1 );
	// img.octave_form( lift_scheme, num_levels );
	Image img = m;
	img.write( band_dir + "lift.m" );

	m.write();
	fprintf( stderr, "\n" );

	// wavelet.inverse_transform( m, Wavelet::Standard );

	// m.reset_region();
	// m.write();
}

void
test_image( )
{
	fprintf( stderr, "test_image\n" );
	fprintf( stderr, "---------------\n" );

	// string img_name = "images/Germany.pgm";
	string img_name = "Germany256.pgm";
	Image img( img_name );

	for( Image::size_type i = img.row_begin(); i != img.row_end(); ++i )
	{
		for( Image::size_type j = img.col_begin(); j != img.col_end(); ++j )
		{
			// img(i,j) = j * j * j;
			// img(i,j) = i + j;
		}
	}
	// img.filter_gauss();
	Image_Region region;
	region = Image_Region( 0, 0, 256, 256 );
	region = img.get_region();
	img.set_region( region );
	// img.write();

	// Wavelet wavelet( Wavelet::Haar );
	// Wavelet wavelet( Wavelet::Daub_4 );
	Wavelet wavelet( Wavelet::Bspline_202 );

	// img.wavelet_transform( wavelet, Wavelet::Standard );
	img.wavelet_transform( wavelet, Wavelet::Pyramid );
	img.reset_region();
	fprintf( stderr, "Writing 'wavelet.pgm'\n" );
	img.write( "wavelet.pgm" );

	/*
	// get and print subbands
	Wavelet_Subbands subbands = img.get_subbands( lift_scheme, num_levels );
	fprintf( stderr, "Done getting subbands\n\n" );
	Wavelet_Subbands::iterator iter     = subbands.begin();
	Wavelet_Subbands::iterator end_iter = subbands.end();
	for( ; iter != end_iter; ++iter )
	{
		const Wavelet_Level_Band& level_band = iter->first;
		const Wavelet_Subband&    subband    = iter->second;

		Image band_img = subband;

		string band_name = band_dir + level_band.to_string() + ".pgm";

		fprintf( stderr, "%s\n", level_band.to_string().c_str() );
		band_img.write( band_name );
	}

	fprintf( stderr, "Lifted signal:\n" );
	//signal.write( "lift.m", 1 );
	// img.octave_form( lift_scheme, num_levels );
	img.write( band_dir + "lift.pgm" );
	 */


	img.set_region( region );
	// img.inverse_wavelet_transform( wavelet, Wavelet::Standard );
	img.inverse_wavelet_transform( wavelet, Wavelet::Pyramid );
	img.reset_region();
	img.write( "inverse.pgm" );
}
