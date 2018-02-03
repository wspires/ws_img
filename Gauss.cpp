/**
	@file   Gauss.cpp
	@author Wade Spires
	@date   2006/1/27
	@brief  Functions based on Gaussian.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Gauss.hpp"

using std::string;
using std::vector;

namespace ws_img
{

double DEFAULT_STD_SPREAD = 3;

/**
  	One-dimensional Gaussian function with mean mu:
		g(x) = (1 / sigma * sqrt(2 * pi) * e^(-(x - mu)^2 / 2*sigma^2 )

	@param[in] x Position
	@param[in] sigma Standard deviation
	@param[in] mu Mean
	@return value Value of Gaussian with std. dev. sigma at position x
 */
double
gauss( double x, double sigma, double mu )
{
	// 1 / (sigma * sqrt(2 * pi))
	// double factor = 1 / (sigma * 2.506628274631);
	double factor = 1;  // omit scaling factor

	double x_mean = x - mu;

	// e^( -(x - mu)^2 / (2 * sigma^2) )
	return( factor * exp( -(x_mean * x_mean) / (2 * sigma * sigma) ) );
}

/**
  	First derivative of one-dimensional Gaussian function:
		g'(x) = (-x / sigma^2) * g(x)

	@param[in] x Position
	@param[in] sigma Standard deviation
	@return value Value of Gaussian derivative with std. dev. sigma at position x
 */
double
gauss_dx( double x, double sigma )
{
	return( (-x / (sigma * sigma)) * gauss(x, sigma) );
}

/**
  	Second derivative of one-dimensional Gaussian function:
		g''(x) = (x^2 - 1) / sigma^2 * g(x)

	@param[in] x Position
	@param[in] sigma Standard deviation
	@return value Value of second Gaussian derivative with std. dev. sigma at
		position x
 */
double
gauss_dxdx( double x, double sigma )
{
	return( ((x * x) - 1) / (sigma * sigma) * gauss(x, sigma) );
}

/**
  	Two-dimensional Gaussian function with 0 mean:
		g(x,y) = (1 / 2 * pi * sigma^2) * e^(-(x^2 + y^2) / 2 * sigma^2 )

	@param[in] x Position
	@param[in] y Position
	@param[in] sigma Standard deviation
	@return value Value of Gaussian with std. dev. sigma at position (x,y)
 */
double
gauss_2D( double x, double y, double sigma )
{
	// double factor = 1 / (TWO_PI * sigma * sigma);
	double factor = 1;  // omit scaling factor

	return( factor * exp( -((x * x) + (y * y)) / (2 * sigma * sigma) ) );
}

/**
	First derivative with respect to x of Gaussian 2-D function.
	(-x / sigma^2) * e^( -(x^2 + y^2) / (2 * sigma^2) )
 */
double
gauss_2D_dx( double x, double y, double sigma )
{
	double sigma2 = sigma * sigma; // sigma^2
	return( (-x / sigma2) * exp( -((x * x) + (y * y)) / (2 * sigma2) ) );
}

/**
	Laplacian of Gaussian function.
		(-1 / pi * sigma^4 ) * (1 - (x^2 + y^2) / sigma^2)
			* e^( -(x^2 + y^2) / (2 * sigma^2) )
 */
double
laplace_gauss( double x, double y, double sigma )
{
	double sigma2     = sigma * sigma;      // sigma^2
	double two_sigma2 = 2 * sigma2;         // 2 * sigma^2
	double x2_plus_y2 = (x * x) + (y * y);  // x^2 + y^2

	if( sigma <= 0 )
	{
		err_quit( "laplace_gauss::Standard deviation sigma must be positive\n" );
	}

	// double scale_factor = -(1 / (PI * sigma2 * sigma2));
	// double scale_factor = 1;
	double scale_factor = -4;
	// double scale_factor = 255;
	// return( scale_factor *
			// (1 - (x2_plus_y2 / two_sigma2)) * exp( -x2_plus_y2 / two_sigma2 ) );
	return( scale_factor *
			((x2_plus_y2 - sigma2) / (sigma2 * sigma2)) * exp( -x2_plus_y2 / two_sigma2 ) );
}

/**
	Create 1D Gaussian filter using given sigma.

	@param[in] sigma Standard deviation of Gaussian--used to determine length of
		filter if no size is given (larger sigma implies longer filter)
	@param[in] std_spread Spread of standard deviation
 */
Vector
make_gauss_filter( double sigma, double std_spread )
{
	return( make_gauss_1D_filter( sigma, std_spread ) );
}

/**
	Create 1D Gaussian filter using given sigma.

	@param[in] sigma Standard deviation of Gaussian--used to determine length of
		filter if no size is given (larger sigma implies longer filter)
	@param[in] std_spread Spread of standard deviation
 */
Vector
make_gauss_1D_filter( double sigma, double std_spread )
{
	// set size to capture ~99.73% of the area under the Gaussian
	Vector::size_type size =
			static_cast<unsigned>( round( std_spread * sigma ) ) + 1;
	if( ws_tools::is_even( size ) )
	{
		++size;
	}
	Vector filter( size );

	// compute each value of the filter using 1D Gaussian (note that the Gaussian
	// is symmetric about its mean)
	Vector::size_type mid_point = filter.size() / 2;
	for( Vector::size_type i = 0; i != (mid_point + 1); ++i )
	{
		Vector::value_type g = gauss( i, sigma );
		filter( mid_point + i ) = g;
		filter( mid_point - i ) = g;
	}

	// normalize filter--we divide by the sum of the filter coefficients instead
	// of (1/sigma * sqrt(2*pi) that is associated with the Gaussian since we
	// are only sampling the Gaussian over a coarse integer grid, but we still
	// need the coefficients to sum to 1
	filter /= filter.sum();

	return( filter );
}

/**
	Create 1D Gaussian derivative filter using given sigma.
	@param[in] sigma Standard deviation of Gaussian--used to determine length of
		filter (larger sigma implies longer filter)
	@param[in] std_spread Spread of standard deviation
 */
Vector
make_gauss_dx_filter( double sigma, double std_spread )
{
	// set size to capture ~99.73% of the area under the Gaussian
	Vector::size_type size =
			static_cast<unsigned>( round( std_spread * sigma ) ) + 1;
	if( ws_tools::is_even( size ) )
	{
		++size;
	}
	Vector filter( size );

	// compute each value of the filter using 1D Gaussian (note that the Gaussian
	// is symmetric about its mean)
	Vector::size_type mid_point = filter.size() / 2;
	for( Vector::size_type i = 0; i != (mid_point + 1); ++i )
	{
		Vector::value_type g = gauss_dx( i, sigma );
		filter( mid_point + i ) = g;
		filter( mid_point - i ) = g;
	}

	// do not normalize since a derivative filter sums to 0

	return( filter );
}

/**
	Make first-derivative Gaussian kernel in the x direction (along the rows).

	Note: The kernel in the y direction (derivative along the columns) is given
	by filter.transpose().

	@param[in] sigma
	@param[in] std_spread Spread of standard deviation
	@retval filter
 */
Matrix
make_gauss_2D_dx_filter( double sigma, double std_spread )
{
	// set size to capture ~99.73% of the area under the Gaussian
	Matrix::size_type size =
			static_cast<unsigned>( round( std_spread * sigma ) ) + 1;
	if( ws_tools::is_even( size ) )
	{
		++size;
	}
	Matrix filter( size );

	Matrix::size_type mid_point = filter.row_size() / 2;
	for( Matrix::size_type i = 0; i != (mid_point + 1); ++i )
	{
		for( Matrix::size_type j = 0; j != (mid_point + 1); ++j )
		{
			double x = j;  // must convert to floating-point to take -x or -y
			double y = i;

			filter( mid_point - i, mid_point + j ) = gauss_2D_dx( x, y, sigma );

			// symmetric along columns
			filter( mid_point + i, mid_point + j ) =
					filter( mid_point - i, mid_point + j );

			// symmetric along rows after negation
			filter( mid_point - i, mid_point - j ) =
					-filter( mid_point - i, mid_point + j );

			filter( mid_point + i, mid_point - j ) =
					-filter( mid_point - i, mid_point + j );
		}
	}

	// do not normalize since a derivative filter sums to 0

	return( filter );
}

/**
	Create 1D Gaussian second derivative filter using given sigma.
	@param[in] sigma Standard deviation of Gaussian--used to determine length of
		filter (larger sigma implies longer filter)
	@param[in] std_spread Spread of standard deviation
	@retval filter
 */
Vector
make_gauss_dxdx_filter( double sigma, double std_spread )
{
	// set size to capture ~99.73% of the area under the Gaussian
	Vector::size_type size =
			static_cast<unsigned>( round( std_spread * sigma ) ) + 1;
	if( ws_tools::is_even( size ) )
	{
		++size;
	}
	Vector filter( size );

	// compute each value of the filter using 1D Gaussian (note that the Gaussian
	// is symmetric about its mean)
	Vector::size_type mid_point = filter.size() / 2;
	for( Vector::size_type i = 0; i != (mid_point + 1); ++i )
	{
		Vector::value_type g = gauss_dxdx( i, sigma );
		filter( mid_point + i ) = g;
		filter( mid_point - i ) = g;
	}

	// do not normalize since a derivative filter sums to 0

	return( filter );
}

/**
	Create 2D Gaussian filter using given sigma. Note that the Gaussian is a
	separable function; hence, filtering with a 2D Gaussian can be done more
	efficiently by filtering the rows and columns with 2 1D Gaussians.

	@param[in] sigma Standard deviation of Gaussian--used to determine length of
		filter if no size is given (larger sigma implies longer filter)
	@param[in] std_spread Spread of standard deviation
 */
Matrix
make_gauss_2D_filter( double sigma, double std_spread )
{
	// set size to capture ~99.73% of the area under the Gaussian
	Matrix::size_type size =
			static_cast<unsigned>( round( std_spread * sigma ) ) + 1;
	if( ws_tools::is_even( size ) )
	{
		++size;
	}
	Matrix filter( size );

	// compute each value of the filter using 1D Gaussian (note that the Gaussian
	// is symmetric about its mean)
	Matrix::size_type mid_point = filter.row_size() / 2;
	for( Matrix::size_type i = 0; i != (mid_point + 1); ++i )
	{
		for( Matrix::size_type j = 0; j != (mid_point + 1); ++j )
		{
			Matrix::value_type g = gauss_2D( i, j, sigma );
			filter( mid_point + i, mid_point + j ) = g;
			filter( mid_point + i, mid_point - j ) = g;
			filter( mid_point - i, mid_point + j ) = g;
			filter( mid_point - i, mid_point - j ) = g;
		}
	}

	// normalize filter
	filter /= filter.sum().sum();

	return( filter );
}

/**
	Create 2D Gaussian filter using given sigma. Note that the Gaussian is a
	separable function; hence, filtering with a 2D Gaussian can be done more
	efficiently by filtering the rows and columns with 2 1D Gaussians.

	Users should use the version that uses separable filters.

	TODO scaling problems

	@param[in] sigma Standard deviation of Gaussian--used to determine length of
		filter if no size is given (larger sigma implies longer filter)
	@param[in] std_spread Spread of standard deviation
 */
Matrix
make_log_filter( double sigma, double std_spread )
{
	// set size to capture ~99.73% of the area under the Gaussian
	Matrix::size_type size =
			static_cast<unsigned>( round( std_spread * sigma ) ) + 1;
	if( ws_tools::is_even( size ) )
	{
		++size;
	}
	Matrix filter( size );

	// compute each value of the filter using LoG function (note that the LoG
	// is symmetric about its center)
	Matrix::size_type mid_point = filter.row_size() / 2;
	for( Matrix::size_type i = 0; i != (mid_point + 1); ++i )
	{
		for( Matrix::size_type j = 0; j != (mid_point + 1); ++j )
		{
			Matrix::value_type log = laplace_gauss( i, j, sigma );
			filter( mid_point + i, mid_point + j ) = log;
			filter( mid_point + i, mid_point - j ) = log;
			filter( mid_point - i, mid_point + j ) = log;
			filter( mid_point - i, mid_point - j ) = log;
		}
	}

	// do not normalize since a derivative filter sums to 0
	fprintf( stderr, "%lf\n", laplace_gauss(0,1,1) );

	return( filter );
}

/**
	Create 2D Laplacian of Gaussian separable filters using the given sigma.

	The convolution is given by 4 1D convolutions:
		(g(x) * (g_yy(y) * I(x,y))) + (g_xx(x) * (g(y) * I(x,y)))
	In words, this is given by the following algorithm:
		1. Convolve image I with 2nd derivative of Gaussian mask g_yy(y) along
			each column;
		2. Convolve resultant image from 1. with Gaussian mask g(x) along each
			row--call this image I_y;
		3. Convolve original image I with Gaussian mask g(y) along each
			column; and
		4. Convolve resultant image from 3. with 2nd derivative of Gaussian mask
			g_xx(x) along each row--call this image I_x.
		5. Add I_x and I_y.

	@param[in] filter_1 Row filter
	@param[in] filter_2 Column filter
	@param[in] sigma Standard deviation of Gaussian--used to determine length of
		filter if no size is given (larger sigma implies longer filter)
	@param[in] std_spread Spread of standard deviation
 */
void
make_log_filter( Vector& filter_1, Vector& filter_2, double sigma,
		double std_spread )
{
	// set size to capture ~99.73% of the area under the Gaussian
	Matrix::size_type size =
			static_cast<unsigned>( round( std_spread * sigma ) ) + 1;
	if( ws_tools::is_even( size ) )
	{
		++size;
	}
	filter_1 = Vector( size );  // second derivative of Gaussian
	filter_2 = Vector( size );  // Gaussian

	double s2 = sigma * sigma;

	// compute LoG for each position
	int z = - static_cast<int>( size / 2 );
	for( Vector::size_type i = filter_1.vec_begin(); i != filter_1.vec_end();
			++i )
	{
		Vector::value_type t = z++;
		Vector::value_type a = (t * t) / s2;
		Vector::value_type b = exp (- a / 2);
		filter_1(i) = (1 - a) * b;
		filter_2(i) = b;
	}

	// Gaussian should technically sum to 1, but the outputs are too faint
	// filter_2 /= filter_2.sum();

	// derivative mask's coefficients should sum to zero
	Vector::value_type value = filter_1.sum() / filter_1.size();
	for( Vector::size_type i = filter_1.vec_begin(); i != filter_1.vec_end();
			++i )
	{
		filter_1(i) -= value;
	}
}

/**
	Apply 2D Laplacian of Gaussian separable filters with the given sigma to the
	input matrix.

	The convolution is given by 4 1D convolutions:
		(g(x) * (g_yy(y) * I(x,y))) + (g_xx(x) * (g(y) * I(x,y)))
	In words, this is given by the following algorithm:
		1. Convolve image I with 2nd derivative of Gaussian mask g_yy(y) along
			each column;
		2. Convolve resultant image from 1. with Gaussian mask g(x) along each
			row--call this image I_y;
		3. Convolve original image I with Gaussian mask g(y) along each
			column; and
		4. Convolve resultant image from 3. with 2nd derivative of Gaussian mask
			g_xx(x) along each row--call this image I_x.
		5. Add I_x and I_y.

	@param[in] A Input matrix
	@param[in] sigma Standard deviation of Gaussian--used to determine length of
		filter if no size is given (larger sigma implies longer filter)
	@param[in] std_spread Spread of standard deviation
 */
Matrix&
filter_log( Matrix& A, double sigma, double std_spread )
{
	// create LoG filters
	Vector f, g;
	make_log_filter( f, g, sigma, std_spread );

	// apply algorithm in the comments above
	A.convolve( f, g ) += Matrix( A ).convolve( g, f );
	return( A );
}

/**
	Apply 2D Laplacian of Gaussian separable filters with the given sigma to the
	input image.

	The convolution is given by 4 1D convolutions:
		(g(x) * (g_yy(y) * I(x,y))) + (g_xx(x) * (g(y) * I(x,y)))
	In words, this is given by the following algorithm:
		1. Convolve image I with 2nd derivative of Gaussian mask g_yy(y) along
			each column;
		2. Convolve resultant image from 1. with Gaussian mask g(x) along each
			row--call this image I_y;
		3. Convolve original image I with Gaussian mask g(y) along each
			column; and
		4. Convolve resultant image from 3. with 2nd derivative of Gaussian mask
			g_xx(x) along each row--call this image I_x.
		5. Add I_x and I_y.

	@param[in] A Input image
	@param[in] sigma Standard deviation of Gaussian--used to determine length of
		filter if no size is given (larger sigma implies longer filter)
	@param[in] std_spread Spread of standard deviation
	@param[in] erase_unfiltered_border Whether to set unfiltered border pixels to
		0 after applying the LoG filter (default is true)
 */
Image&
filter_log( Image& A, double sigma, double std_spread,
		bool erase_unfiltered_border )
{
	// create LoG filters
	Vector f, g;
	make_log_filter( f, g, sigma, std_spread );

	// apply algorithm in the comments above
	A.convolve( f, g ) += Image( A ).convolve( g, f );

	Vector::size_type mid = f.size() / 2;

	/*
		Set un-filtered border pixels to black
	 */
	if( erase_unfiltered_border )
	{
		Image::size_type M = A.sub_row_size();
		Image::size_type N = A.sub_col_size();

		// left border
		A.push_region( A.row_begin(), A.col_begin(), M, mid );
		for( Image::size_type i = A.row_begin(); i != A.row_end(); ++i )
		{
			for( Image::size_type j = A.col_begin(); j != A.col_end(); ++j )
			{
				A(i,j) = Image::Black_Pixel;
			}
		}
		A.pop_region();

		// top border
		A.push_region( A.row_begin(), A.col_begin(), mid, N );
		for( Image::size_type i = A.row_begin(); i != A.row_end(); ++i )
		{
			for( Image::size_type j = A.col_begin(); j != A.col_end(); ++j )
			{
				A(i,j) = Image::Black_Pixel;
			}
		}
		A.pop_region();

		// bottom border
		A.push_region( M - mid, A.col_begin(), mid, N );
		for( Image::size_type i = A.row_begin(); i != A.row_end(); ++i )
		{
			for( Image::size_type j = A.col_begin(); j != A.col_end(); ++j )
			{
				A(i,j) = Image::Black_Pixel;
			}
		}
		A.pop_region();

		// right border
		A.push_region( A.row_begin(), N - mid, M, mid );
		for( Image::size_type i = A.row_begin(); i != A.row_end(); ++i )
		{
			for( Image::size_type j = A.col_begin(); j != A.col_end(); ++j )
			{
				A(i,j) = Image::Black_Pixel;
			}
		}
		A.pop_region();
	}

	return( A );
}

} // namespace ws_img
