/**
	@file   example.cpp
	@author Wade Spires
	@date   2005/08/02
	@brief  Example program for using Image class.
 */

#include <iostream>
#include <string>
#include <vector>

#include <cstdio>
#include <ctime>

#include "Edge_Detector.hpp"
#include "Gauss.hpp"
#include "Image.hpp"
#include "Image_Region.hpp"
#include "Image_Coordinate.hpp"
#include "Linear_Algebra.hpp"
#include "ws_tools.hpp"

using namespace ws_img;

using std::cerr;  using std::cout;  using std::endl;
using std::string;
using std::vector;

Image modify( const Image & );

// image function tests
void test_draw( Image& );
void test_rotate( Image& );
void test_rotate2( Image& );
void test_stat( Image& );
void test_normalize( Image& );
void test_gauss( Image& );
void test_matrix( );
void comp_matrix( );

int main( int argc, char** argv )
{
	if( argc != 3 )
	{
		cerr << "usage: example input_img output_img" << endl;
		exit( EXIT_FAILURE );
	}

	char* in_img_name  = argv[ 1 ];
	char* out_img_name = argv[ 2 ];

	// read in the input image from the specified file
	Image in_img( in_img_name );

	Image out_img = in_img;
	test_rotate2( out_img );
	// out_img.to_color();
	out_img.write( out_img_name );

	/*
	// make a modified copy of the input image
	Image out_img = modify( in_img );
	out_img.write( out_img_name );
	 */

	/* perform additional tests */

	//vector<double> v = in_img.make_gauss_filter( unsigned(3) );
	//in_img.make_gauss_filter( 1.5 );

	// additional tests
	//test_stat( img );
	//test_draw( img );
	//test_rotate( in_img );

	//test_normalize( in_img );

	test_gauss( in_img );

	// test_matrix();
	// comp_matrix();

	return( 0 );
}

/**
	Modify the image to show some of the Image class functionality.
	@param[in] in_img Input image
	@retval[out] out_img Output image
 */
Image modify( const Image & in_img )
{
	// generate an empty output image with half the dimension
	Image out_img =
		Image( in_img.row_size() / 2, in_img.col_size() / 2 );

	// copy every other pixel into the output image with a higher intensity
	for( unsigned i = 0; i != out_img.row_size(); i++ )
	{
		for( unsigned j = 0; j != out_img.col_size(); j++ )
		{
			// get new pixel intensity value
			Image::pixel_type new_pixel = in_img( 2 * i, 2 * j );

			// if new value exceeds the largest pixel value an image can have
			if( new_pixel > out_img.get_max_value() )
			{
				// set it to maximum value (white)
				out_img( i, j ) = out_img.get_max_value();
			}
			else
			{
				// use computed value as new pixel intensity
				out_img( i, j ) = new_pixel;
			}
		}
	}
	out_img *= 1.15;  // change intensity values

	return( out_img );
}

/**
	Test draw functions on image.
	@param img Image to draw on
 */
void test_draw( Image& img )
{
	Image_Coordinate left( 0, 40 );
	Image_Coordinate right( 100, 16 );
	Image_Coordinate l2( 590, 40 );
	Image_Coordinate r2( 300, 180 );
	//img.draw( left );
	//img.draw( right );
	img.draw_line( left, right );
	img.draw_line( l2, r2 );

	Image_Region win( 50, 50, 55, 75 );
	img.draw2( win );

	Image_Coordinate center( 500, 601 );
	img.draw_circle( center, 50 );

	img.write( "draw.pgm" );
}

/**
	Test rotate function on image.
	@param img Image to rotate
 */
void test_rotate( Image& img )
{
	Image out_img_1 = img;
	Image out_img_2 = img;
	Image out_img_3 = img;
	//out_img.scale( 1.189, 1.189 );
	//out_img_1.scale( .841, .841 );
	//out_img_2.scale( 1, 1 );
	//out_img_3.scale( 2, 2 );

	// note that rotation is in radians
	Image_Coordinate coord( 0, 100);
	out_img_1.rotate( .4 );
	//out_img_2.rotate( .4, coord );
	out_img_3.rotate( 2 * 3.1415 / 3,//3.1415 / 4,
		Image_Coordinate( img.row_size() / 2, img.col_size() / 2 ) );

	// write the output image to the specified file
	out_img_1.write( "out1.pgm" );
	//out_img_2.write( "out2.pgm" );
	out_img_3.write( "out3.pgm" );

	//Image out_img_4 = img;
	//out_img_4 *= 10.3;
	//out_img_4.rotate( .079 );
	//out_img_4.rotate( .079, 1, 1 );//, img.col_size() / 2, img.row_size() / 2 );
	//out_img_4.rotate( -.0079 );
	//out_img_4.write( "out4.pgm" );
}

/**
	Test rotate function on image.
	@param img Image to rotate
 */
void test_rotate2( Image& img )
{
	// img.rotate( 3.1415 / 2, 0, 0 );
	// img.push_region( 20, 20, 80, 80 );
	// img.write( "out2.pgm" );
	img.rotate( 45 );
	// img.pop_region( );
	// img.rotate( 30, false );
}

/**
	Test statistics functions on image.
	@param img Image to stat
void test_stat( Image& img )
{
	img.write( "in.lsh" );
	Image out_img = img;
	Hist h( out_img );
	fprintf( stderr, "%u %lf %lf %lf %lf\n",
		h.get_total_count(),
		h.get_avg(),
		h.get_var(),
		h.get_dev(),
		h.get_med()
	);
	h.equalize( out_img );
	//h.expand( out_img );
	h.write_matrix( "hist.m" );
	out_img.write( "out.pgm" );

	Hist h2( out_img );
	h2.write_matrix( "hist2.m" );
}
 */

/**
	Test normalization function on image.
	@param img Image to test
 */
void test_normalize( Image& img )
{
	// copy every other pixel into the output image with a higher intensity
	for( unsigned i = 0; i != img.row_size(); i++ )
	{
		for( unsigned j = 0; j != img.col_size(); j++ )
		{
			// make even column values < 0
			if( j % 2 == 0 )
			{
				img(i, j) -= 300;
			}
		}
	}

	img.write( "unnorm.pgm", "a" );
	img.normalize();
	img.write( "norm.pgm", "a" );
}

/**
	Test Gaussian functions on image.
	@param img Image to test
 */
void test_gauss( Image& img )
{
	double sigma = 1.0;
	Vector vec = make_gauss_1D_filter( sigma );
	vec.write( "", 5 );

	fprintf( stderr, "\n" );

	Matrix mat = make_gauss_2D_filter( sigma );
	mat.write( "", 5 );

	fprintf( stderr, "\n" );

	Matrix log = make_log_filter( sigma );
	log.write( "", 5 );

	Matrix gx_filter = make_gauss_2D_dx_filter( 1.0 );
	Matrix gy_filter = gx_filter.transpose();
	gx_filter.write( "", 5 );
	gy_filter.write( "", 5 );
}

/**
	Test matrix functions.
 */
void test_matrix( )
{
	fprintf( stderr, "Solving Ax=b\n" );

	unsigned size = 2;
	Matrix A( size, size );
	for( unsigned i = 0; i != A.row_size(); ++i )
	{
		for( unsigned j = 0; j != A.col_size(); ++j )
		{
			if( i == j )
			{
				A(i,j) = 1;
			}
			else
			{
				A(i,j) = 0;
			}
			//A(i,j) = i + j;
		}
	}

	A = Matrix::random( size, size, 10 );
	/*
	A(0,0) = 8;
	A(0,1) = 3;
	A(0,2) = 3;
	A(1,0) = 9;
	A(1,1) = 3;
	A(1,2) = 1;
	A(2,0) = 4;
	A(2,1) = 0;
	A(2,2) = 6;
	 */

	fprintf( stderr, "Original A:\n" );
	A.write();
	fprintf( stderr, "\n" );

	Vector b( A.row_size() );
	for( unsigned i = 0; i != b.size(); ++i )
	{
		b(i) = i;
	}

	fprintf( stderr, "b:\n" );
	b.write();
	fprintf( stderr, "\n" );

	// save A and b 
	Matrix A_2 = A;
	Vector b_2 = b;

	// solve using original method
	Vector x( A.row_size() );
	Lin_Alg::solve_system( A, b, x );

	fprintf( stderr, "LU (new A):\n" );
	A.write( "", 5 );
	fprintf( stderr, "\n" );

	fprintf( stderr, "x (solution):\n" );
	x.write( "", 5 );
	fprintf( stderr, "\n" );

	fprintf( stderr, "A * x =\n" );
	b.write( "", 5 );
	fprintf( stderr, "\n" );

	// solve using new method
	Lin_Alg::solve_system( A_2, b_2 );

	fprintf( stderr, "LU (new A):\n" );
	A_2.write( "", 5 );
	fprintf( stderr, "\n" );

	fprintf( stderr, "x (solution):\n" );
	b.write( "", 5 );
	fprintf( stderr, "\n" );

	fprintf( stderr, "A * x =\n" );
	Matrix b_result = A * b;
	b_result.write( "", 5 );
	fprintf( stderr, "\n" );
}

/**
	Compare matrix versions for speed.
	Should have identical function in example/example.cpp
 */
void comp_matrix( )
{
	time_t start_time = time(NULL);

	unsigned size = 5000;
	//Matrix m = Matrix::random( size );
	Matrix m( size );

	for( unsigned i = 0; i != m.row_size(); ++i )
	{
		for( unsigned j = 0; j != m.col_size(); ++j )
		{
			m(i,j) = 100;
		}
	}
	// m *= m;
	//m.write();

	time_t end_time = time(NULL);
	fprintf( stderr, "Time: %u\n", unsigned(end_time - start_time) );
}
