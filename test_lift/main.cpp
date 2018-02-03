/**
	@file   test_lift.cpp
	@author Wade Spires
	@date   2005/10/10
	@brief  Test program for using lifting scheme.
 */

#include <iostream>
#include <algorithm>
#include <map>
#include <string>
#include <vector>

#include <cmath>
#include <ctime>
#include <cstdlib>

#include "Image.hpp"
#include "Lift_Scheme.hpp"
#include "Wavelet_Subband.hpp"
#include "ws_tools.hpp"
#include "Progress_Bar.hpp"

using namespace ws_img;

using std::cerr;  using std::cout;  using std::endl;
using std::map;
using std::string;
using std::vector;

void print_usage( );

void   lift_img( const string&, unsigned, unsigned, unsigned );
void   lift_img_dir( const string&, unsigned, unsigned, unsigned );
void   quantize( Image& );
double get_entropy( const Image& );
void   check_for_equality( const Image&, const Image&, const string& );

void test_1D_lift( unsigned, unsigned, unsigned );
void test_2D_lift( );
void test_1D_lift_int( );

static const double EPS = .000001;

// output files
const string out_lift_name   = "lift.pgm";
const string out_unlift_name = "unlift.pgm";

int main( int argc, char** argv )
{
	if( argc != 5 )
	{
		print_usage();
	}

	unsigned N          = atoi( argv[1] );
	unsigned N_tilde    = atoi( argv[2] );
	unsigned num_levels = atoi( argv[3] );
	string img          = argv[4];

	// lift_img( img, N, N_tilde, num_levels );
	lift_img_dir( img, N, N_tilde, num_levels );

	// test_1D_lift( N, N_tilde, num_levels );

	// test_2D_lift();
	// test_1D_lift_int();

	return( EXIT_SUCCESS );
}

void print_usage( )
{
	cerr << "usage: test_lift "
		<< "num_dual_moments num_real_moments num_levels "
		<< "in_img"
		<< endl;
	exit( EXIT_FAILURE );
}

/**
	Apply 2D lifting scheme on image.
	@param[in] img_name Name of image
	@param[in] N Number of dual vanishing moments
	@param[in] N_tilde Number of real vanishing moments
	@param[in] num_levels Number of levels in the transformation
 */
void lift_img( const string& img_name, unsigned N, unsigned N_tilde,
	unsigned num_levels )
{
	time_t start_time = time( NULL );

	fprintf( stderr, "Processing image '%s'\n", img_name.c_str() );
	Image img( img_name );

	// create copy so we can compare after the inverse for equality
	Image copy = img;

	// test lifting scheme on partial image region 
	Image_Region region;
	region = Image_Region( 50, 150, 125, 125 );  // use sub-region
	region = img.get_region();  // use entire region
	// img.set_region( region );

	double entropy = get_entropy( img );
	fprintf( stderr, "   Entropy of original image: %lf\n", entropy );

	// perform lifting
	Lift_Scheme lift_scheme( N, N_tilde );
	img.lift( lift_scheme, num_levels );

	// quantize( img );
	// img.normalize();
	// img.octave_form( lift_scheme, num_levels );
	// img.write( "lift.pgm", "a" );
	// img.reset_region();
	img.write( out_lift_name );

	entropy = get_entropy( img );
	fprintf( stderr, "   Entropy of lifted image:   %lf\n", entropy );

	// perform inverse lifting
	// img.set_region( region );
	img.inverse_lift( lift_scheme, num_levels );
	// img.reset_region();
	img.write( out_unlift_name );

	// verify the inverse image is the same as the orignal image
	// check_for_equality( copy, img, img_name );

	time_t end_time = time( NULL );
	fprintf( stdout, "Finished after %u second(s)\n",
			unsigned(end_time - start_time) );
}

/**
	Apply 2D lifting scheme on all images in directory.
	@param[in] dir_name Name of image directory
	@param[in] N Number of dual vanishing moments
	@param[in] N_tilde Number of real vanishing moments
	@param[in] num_levels Number of levels in the transformation
 */
void lift_img_dir( const string& dir_name, unsigned N, unsigned N_tilde,
	unsigned num_levels )
{
	time_t start_time = time( NULL );  // time the entire run
	time_t lift_time  = 0;             // time taken for lifting only

	vector<string> image_files = ws_tools::dir_traverse( dir_name, Image::is_image );

	// Progress_Bar bar( image_files.size() );

	Lift_Scheme lift_scheme( N, N_tilde );

	// entropies for original and lifted images
	Vector orig_entropies( image_files.size() );
	Vector lift_entropies( image_files.size() );

	for( unsigned i = 0; i != image_files.size(); ++i )
	{
		string img_name = image_files[i].c_str();
		fprintf( stderr, "Processing '%s'\n", img_name.c_str() );
		Image img( img_name );

		// create copy so we can compare after the inverse for equality
		Image copy = img;

		double entropy = get_entropy( img );
		orig_entropies(i) = entropy;
		fprintf( stderr, "   Entropy of original image: %lf\n", entropy );

		// perform lifting
		time_t lift_start = time( NULL );
		img.lift( lift_scheme, num_levels );
		lift_time += time(NULL) - lift_start;


	string band_dir = "band_dir/";
	ws_tools::check_dir( band_dir );

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


		// quantize( img );
		// img.normalize();
		// img.octave_form( lift_scheme, num_levels );
		img.write( "lift.pgm" );
		img.write( out_lift_name );
		img.inverse_octave_form( lift_scheme, num_levels );

		entropy = get_entropy( img );
		lift_entropies(i) = entropy;
		fprintf( stderr, "   Entropy of lifted image:   %lf\n", entropy );

		// perform inverse lifting
		lift_start = time( NULL );
		img.inverse_lift( lift_scheme, num_levels );
		img.write( out_unlift_name );
		lift_time += time(NULL) - lift_start;

		// verify the inverse image is the same as the orignal image
		check_for_equality( copy, img, img_name );

		// bar.update( i );
		fprintf( stderr, "\n" );
	}

	time_t end_time = time( NULL );
	fprintf( stdout, "Finished %u image(s) after %u second(s)\n",
			unsigned(image_files.size()), unsigned(end_time - start_time) );
	fprintf( stdout, "Time for lifting: %u second(s)\n", unsigned(lift_time) );

	fprintf( stderr, "\n" );

	fprintf( stderr, "Entropy Summary\n" );
	fprintf( stderr, "---------------\n" );

	fprintf( stderr, "\n" );
	fprintf( stderr, "Original Image(s):\n" );
	double min, max;
	orig_entropies.get_min_max( &min, &max );
	fprintf( stderr, "(min, max) = (%lf, %lf)\n", min, max );

	double avg, var, std;
	orig_entropies.stat( &avg, &var, &std );
	fprintf( stderr, "(avg, var, std) = (%lf, %lf, %lf)\n", avg, var, std );

	fprintf( stderr, "\n" );
	fprintf( stderr, "Lifted Image(s):\n" );
	lift_entropies.get_min_max( &min, &max );
	fprintf( stderr, "(min, max) = (%lf, %lf)\n", min, max );

	lift_entropies.stat( &avg, &var, &std );
	fprintf( stderr, "(avg, var, std) = (%lf, %lf, %lf)\n", avg, var, std );
}

/**
	Quantize image by setting pixel to 0 if its absolute value is less than some
	threshold.
 */
void quantize( Image& img )
{
	Image::value_type thresh = 1.0;

	for( unsigned i = img.row_begin(); i != img.row_end(); ++i )
	{
		for( unsigned j = img.col_begin(); j != img.col_end(); ++j )
		{
			// 1. round to nearest integer
			img(i,j) = round( img(i,j) );

			// 2. eliminate all pixels within range of the threshold
			if( fabs( img(i,j) ) < thresh )
			{
				img(i,j) = 0;
			}
		}
	}
}

/**
	Get first-order entropy of image.
	@param[in] img Image to find entropy for
 */
double get_entropy( const Image& img )
{
	double entropy = 0;

	typedef map<Image::value_type, unsigned> Histogram;
	Histogram histogram;

	Image::const_iterator img_iter = img.const_begin();
	Image::const_iterator img_end  = img.const_end();

	// accumulate the frequency of each value in the image
	for( ; img_iter != img_end; ++img_iter )
	{
		++histogram[ *img_iter ];
	}

	// calculate entropy
	double pixel_count = img.row_size() * img.col_size();
	Histogram::const_iterator hist_iter = histogram.begin();
	Histogram::const_iterator hist_end  = histogram.end();
	for( ; hist_iter != hist_end; ++hist_iter )
	{
		double prob = hist_iter->second / pixel_count;
		entropy += prob * log2( prob );
	}

	return( -entropy );
}

/**
	Compare original image with reconstructed image.
	@param[in] img Original image
	@param[in] img_lift Lifted image
 */
void
check_for_equality( const Image& img, const Image& img_lift,
		const string& img_name )
{
	Image diff = img - img_lift;
	for( Image::size_type i = img.row_begin(); i != img.row_end(); ++i )
	{
		for( Image::size_type j = img.row_begin(); j != img.col_end(); ++j )
		{
			// if( diff(i,j) != 0.0 )
			if( diff(i,j) > EPS )
			{
				// fprintf( stderr, "%s: (%u, %u) (%lf vs %lf)\n",
						// img_name.c_str(), i, j, img(i,j), img_lift(i,j) );

				err_quit( "Equality check failed for %s: (%u, %u) (%lf vs %lf)\n",
						img_name.c_str(),
						i, j,
						img(i,j), img_lift(i,j) );
			}
		}
	}
}

/**
	Test 1D lifting scheme.
 */
void test_1D_lift( unsigned N, unsigned N_tilde, unsigned num_levels )
{
	fprintf( stderr, "1D Lifting\n" );

 	unsigned size = 18;
 	size = 8;

	Vector signal( size );
	for( unsigned i = 0; i != size; ++i )
	{
		signal(i) = i;
	}

	fprintf( stderr, "Original signal:\n" );
	signal.write();
	fprintf( stderr, "\n" );

	Lift_Scheme lift_scheme( N, N_tilde );
	fprintf( stderr, "Calling lift\n" );
	lift_scheme.lift( signal, num_levels );

	fprintf( stderr, "Lifted signal:\n" );
	signal.write();
	fprintf( stderr, "\n" );

	lift_scheme.inverse_lift( signal, num_levels );

	fprintf( stderr, "Inverted signal:\n" );
	signal.write();
	fprintf( stderr, "\n" );
}

/**
	Test 2D lifting scheme.
 */
void test_2D_lift( )
{
	fprintf( stderr, "2D Lifting\n" );

	unsigned N       = 4;
	unsigned N_tilde = 4;
	unsigned num_levels = 10000;

 	unsigned size = 32;

	Matrix signal( size, size / 4 );
	for( unsigned i = 0; i != signal.row_size(); ++i )
	{
		for( unsigned j = 0; j != signal.col_size(); ++j )
		{
			signal(i, j) = (i + 1) * (j + 1);
			// signal(i, j) = i + j;
			// signal(i, j) = j;
		}
	}
	// signal = Matrix::random( size, size, 10 );
	// signal = Matrix::random( size, size / 4, 10 );

	//signal.write( "signal.m", 1 );

	fprintf( stderr, "Original signal:\n" );
	// signal.write();
	fprintf( stderr, "\n" );

	Lift_Scheme lift_scheme( N, N_tilde );
	lift_scheme.lift( signal, num_levels );
	// signal.write();

	signal.write();

	string band_dir = "band_dir/";
	ws_tools::check_dir( band_dir );

	// get and print subbands
	Wavelet_Subbands subbands = lift_scheme.get_subbands( signal, num_levels );
	Wavelet_Subbands::iterator iter     = subbands.begin();
	Wavelet_Subbands::iterator end_iter = subbands.end();
	for( ; iter != end_iter; ++iter )
	{
		const Wavelet_Level_Band& level_band = iter->first;
		const Wavelet_Subband&    subband    = iter->second;

		string band_name = band_dir + level_band.to_string() + ".m";

		fprintf( stderr, "%s\n", level_band.to_string().c_str() );
		// subband.write( band_name );
		subband.write();
	}

	fprintf( stderr, "Lifted signal:\n" );
	//signal.write( "lift.m", 1 );
	lift_scheme.octave_form( signal, num_levels );
	// signal.write( band_dir + "lift.m" );
	signal.write();
	fprintf( stderr, "\n" );

	lift_scheme.inverse_lift( signal, num_levels );

	//signal.write( "inv.m", 1 );

	fprintf( stderr, "Inverted signal:\n" );
	// signal.write();
	fprintf( stderr, "\n" );
}

/**
	Test 1D integer lifting scheme.
 */
void test_1D_lift_int( )
{
	fprintf( stderr, "1D Lifting\n" );

	unsigned N       = 4;
	unsigned N_tilde = 4;
	unsigned num_levels = 10000;

 	unsigned size = 23;

	Vector signal( size );
	for( unsigned i = 0; i != size; ++i )
	{
		signal(i) = i * i * i * i / 5;
		signal(i) = i * i * i * 50;
		signal(i) = i;
	}
	Vector copy = signal;
	signal = Vector::random( size, 50 );

	fprintf( stderr, "Original signal:\n" );
	signal.write( "", 2 );
	fprintf( stderr, "\n" );

	Lift_Scheme lift_scheme( N, N_tilde );
	lift_scheme.lift_int( signal, num_levels );

	fprintf( stderr, "Lifted signal:\n" );
	signal.write( "", 2 );
	fprintf( stderr, "\n" );

	lift_scheme.inverse_lift_int( signal, num_levels );

	fprintf( stderr, "Inverted signal:\n" );
	signal.write( "", 2 );
	fprintf( stderr, "\n" );

	// make sure tranform was inverted correctly
	Vector diff = copy - signal;
	for( Vector::size_type i = diff.vec_begin(); i != diff.vec_end(); ++i )
	{
		if( diff(i) != 0 )
		{
			err_quit( "Equality check failed for vec at %u (%lf vs %lf)\n",
					i, signal(i), copy(i) );
		}
	}
	fprintf( stderr, "Equality check: Pass!\n" );
}
