/**
	@file   Edge_Detector.cpp
	@author Wade Spires
	@date   2005/12/1
	@brief  Class Edge_Detector.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Edge_Detector.hpp"

using std::string;
using std::vector;

using namespace ws_img;

typedef Edge_Detector::pixel_type  pixel_type; 
typedef Edge_Detector::size_type   size_type; 

// edges are white and non-edges are black (used locally by functions but not by
// the user)
const pixel_type Edge_Detector::_Edge_Pixel     = Image::White_Pixel;
const pixel_type Edge_Detector::_Non_Edge_Pixel = Image::Black_Pixel;

// edges are white and non-edges are black (relative to a bilevel image so used
// by the user)
// const pixel_type Edge_Detector::Edge_Pixel     = 0;
// const pixel_type Edge_Detector::Non_Edge_Pixel = 1;
const pixel_type Edge_Detector::Edge_Pixel     = Image::White_Pixel;
const pixel_type Edge_Detector::Non_Edge_Pixel = Image::Black_Pixel;

const pixel_type Edge_Detector::DEFAULT_THRESHOLD = 150;

/**
	Detect edges in image using given threshold.
	@param[in] img_name Image name
	@param[in] threshold Threshold
	@retval edges Edge image 
 */
Image
Edge_Detector::detect( string& img_name, pixel_type threshold )
{
	Image edge_img( img_name );
	detect( edge_img, threshold );
	return( edge_img );
}

/**
	Detect edges in image using given threshold.
	@param[in,out] img Image to detect edges in--the image is overwritten with
		the edges
	@param[in] threshold Threshold
 */
void
Edge_Detector::detect( Image& img, pixel_type threshold )
{
	switch( _detector_type )
	{
		case Prewitt:
			prewitt( img, threshold );
			break;

		case Sobel:
			sobel( img, threshold );
			break;

		case Canny:
			canny( img, threshold );
			break;

		default:
			err_quit( "detect: Unsupported edge detector given\n" );
			break;
	}
}

/**
	Detect edges in image using the Prewitt edge detector.
	@param[in,out] img Image to find edges in
	@param[in] threshold Threshold
 */
void
Edge_Detector::prewitt( Image& img, pixel_type threshold )
{
	// construct Prewitt derivative filters
	Matrix x_filter(3);
	x_filter(0,0) =  1;
	x_filter(0,1) =  0;
	x_filter(0,2) = -1;
	x_filter(1,0) =  1;
	x_filter(1,1) =  0;
	x_filter(1,2) = -1;
	x_filter(2,0) =  1;
	x_filter(2,1) =  0;
	x_filter(2,2) = -1;

	Matrix y_filter(3);
	y_filter(0,0) =  1;
	y_filter(0,1) =  1;
	y_filter(0,2) =  1;
	y_filter(1,0) =  0;
	y_filter(1,1) =  0;
	y_filter(1,2) =  0;
	y_filter(2,0) = -1;
	y_filter(2,1) = -1;
	y_filter(2,2) = -1;

	// derivative images of gray-scale image
	Image Ix = img.to_gray();
	Image Iy = Ix;

	// filter image
	Ix.convolve( x_filter );
	Iy.convolve( y_filter );

	compute_edge_map( img, Ix._red_data, Iy._red_data, threshold );

	// convert to black-and-white (binary) image where the edges are white (0)
	// img.to_bilevel();
}

/**
	Detect edges in image using the Sobel edge detector.
	@param[in,out] img Image to find edges in
	@param[in] threshold Threshold
 */
void
Edge_Detector::sobel( Image& img, pixel_type threshold )
{
	// construct Sobel derivative filters
	Matrix x_filter(3);
	x_filter(0,0) =  1;
	x_filter(0,1) =  0;
	x_filter(0,2) = -1;
	x_filter(1,0) =  2;
	x_filter(1,1) =  0;
	x_filter(1,2) = -2;
	x_filter(2,0) =  1;
	x_filter(2,1) =  0;
	x_filter(2,2) = -1;

	Matrix y_filter(3);
	y_filter(0,0) =  1;
	y_filter(0,1) =  2;
	y_filter(0,2) =  1;
	y_filter(1,0) =  0;
	y_filter(1,1) =  0;
	y_filter(1,2) =  0;
	y_filter(2,0) = -1;
	y_filter(2,1) = -2;
	y_filter(2,2) = -1;

	// derivative images of gray-scale image
	Image Ix = img.to_gray();
	Image Iy = Ix;

	// filter image
	Ix.convolve( x_filter );
	Iy.convolve( y_filter );

	compute_edge_map( img, Ix._red_data, Iy._red_data, threshold );

	// convert to black-and-white (binary) image where the edges are white (0)
	// img.to_bilevel();
}

/**
  	Compute edge map.

	Algorithm:
	for each pixel position
		compute gradient magnitude using image derivatives Ix and Iy:
				sqrt( Ix^2 + Iy^2 )

		if gradient magnitude exceeds given threshold
			mark pixel position as an edge

	TODO change indices to be correct for different regions--currently assumes
	that Ix is the same size as all of *this

	@param[out] img Image
	@param[in] Ix Partial derivative of image in x direction
	@param[in] Iy Partial derivative of image in y direction
	@param[in] threshold Threshold to compute edges with
 */
void
Edge_Detector::compute_edge_map( Image& img, const Matrix& Ix, const Matrix& Iy,
	pixel_type threshold )
{
	// compute gradient magnitude
	for( size_type i = img.row_begin(); i != img.row_end(); ++i )
	{
		for( size_type j = img.col_begin(); j != img.col_end(); ++j )
		{
			// slower, more precise
			// pixel_type grad_mag = hypot( Ix(i,j), Iy(i,j) );

			// faster approximation
			pixel_type grad_mag = fabs( Ix(i,j) ) + fabs( Iy(i,j) );

			if( grad_mag > threshold )
			{
				img(i,j) = _Edge_Pixel;
			}
			else
			{
				img(i,j) = _Non_Edge_Pixel;
			}
		}
	}
}

/**
	Detect edges in image using the Canny edge detector.
	@param[in,out] img Image
	@param[in] low_threshold Low threshold
	@param[in] high_threshold High threshold
 */
void
Edge_Detector::canny( Image& img, pixel_type low_threshold,
	pixel_type high_threshold )
{
	// derivative images
	Image Ix, Iy;

	// gradient magnitude
	Matrix grad_mag;

	canny( img, Ix, Iy, grad_mag, low_threshold, high_threshold );
}

/**
	Detect edges in image using the Canny edge detector.
	Save intermediate representations, e.g., first derivative images.

	@param[in,out] img Image to find edges in
	@param[out] Ix First derivative in x direction
	@param[out] Iy First derivative in y direction
	@param[out] grad_mag Gradient magnitude of image
	@param[out] edge_map Edge map after non-maxima suppression
	@param[in] low_threshold Low threshold
	@param[in] high_threshold High threshold
 */
void
Edge_Detector::canny( Image& img, Image& Ix, Image& Iy, Matrix& grad_mag,
	pixel_type low_threshold, pixel_type high_threshold )
{
	Matrix edges;
	canny( img, Ix, Iy, grad_mag, edges, low_threshold, high_threshold );
}

// faster (hopefully) Canny implementation that tries to avoid creating
// temporary objects whenever possible
// #define PERFORMANCE

/**
	Detect edges in image using the Canny edge detector.
	Save intermediate representations, e.g., first derivative images.

	@param[in,out] img Image to find edges in
	@param[out] Ix First derivative in x direction
	@param[out] Iy First derivative in y direction
	@param[out] grad_mag Gradient magnitude of image
	@param[out] edge_map Edge map after non-maxima suppression
	@param[in] low_threshold Low threshold
	@param[in] high_threshold High threshold
 */
void
Edge_Detector::canny( Image& img, Image& Ix, Image& Iy, Matrix& grad_mag,
	Matrix& edge_map, pixel_type low_threshold, pixel_type high_threshold )
{
	// TODO Handle color properly: (1) compute each color channel separately and
	// add the 3 results, or (2) Use p. 337 of Dig. Img. Proc. to compute
	// gradient correctly

// faster Canny implementation (hopefully) that tries to avoid creating
// temporary objects whenever possible
#ifdef PERFORMANCE

	// construct Canny derivative filters (first derivative of Gaussian)
	static Matrix gx_filter = make_gauss_2D_dx_filter( 1, 2 );
	static Matrix gy_filter = make_gauss_2D_dx_filter( 1, 2 ).transpose();

	// derivative images of gray-scale image (we use img in place of Iy, so
	// users should not expect to have Iy computed with this version)
	Ix = Image( img.to_gray() ).convolve( gx_filter );
	img.convolve( gy_filter );

	// perform non-maxima suppression and hysteresis
	img._red_data = non_maxima_suppression( Ix._red_data, img._red_data,
			grad_mag );
	hysteresis( img._red_data, low_threshold, high_threshold );

	// don't handle color since we convert to gray-scale above

#else  // straight-forward Canny implementation

	if( img.sub_row_size() <= 3 || img.sub_col_size() <= 3 )
	{
		err_warn( "Image is too small to perform Canny edge detection\n" );
		return;
	}

	// derivative images of gray-scale image (don't convert img, though)
	Ix = Image( img ).to_gray();
	Iy = Ix;

	// construct Canny derivative filters (first derivative of Gaussian)
	static const Matrix gx_filter = make_gauss_2D_dx_filter( 1, 2 );
	static const Matrix gy_filter = make_gauss_2D_dx_filter( 1, 2 ).transpose();

	// apply filters
	Ix.convolve( gx_filter );
	Iy.convolve( gy_filter );

	// perform non-maxima suppression and hysteresis
	edge_map = non_maxima_suppression( Ix._red_data, Iy._red_data, grad_mag );

	// Image( edge_map ).write( "redges.pgm" );
	hysteresis( img, edge_map, low_threshold, high_threshold );

	// store edges in image
	// TODO should we assume that we are using the red channel?
	if( img.is_color() )
	{
		img._green_data = img._red_data;
		img._blue_data  = img._red_data;
	}

	// convert to black-and-white (binary) image where the edges are white (0)
	// img.to_bilevel();

#endif
}

/**
  	Perform non-maxima suppression.

	Algorithm:
		compute gradient magnitude:
			for each pixel position
				compute gradient magnitude using image derivatives Ix and Iy:
						sqrt( Ix^2 + Iy^2 )

		perform non-maxima suppression: 

	Implementation notes:
		For non-maxima suppression, we skip the first and last rows (and columns)
		since, for each pixel (i,j), we must access one of (i,j)'s neighbors:
			(i-1,j), (i,j-1), (i-1,j-1), ..., (i+1,j+1)
		which is undefined for border pixels.

	@param[in] Ix Partial derivative of image in x direction
	@param[in] Iy Partial derivative of image in y direction
 */
Matrix
Edge_Detector::non_maxima_suppression( const Matrix& Ix, const Matrix& Iy,
	Matrix& grad_mag )
{
	// compute gradient magnitude using image derivatives Ix and Iy:
	// sqrt( Ix^2 + Iy^2 )
	grad_mag.resize( Ix.row_size(), Ix.col_size(), false );
	for( size_type i = grad_mag.row_begin(); i != grad_mag.row_end(); ++i )
	{
		for( size_type j = grad_mag.col_begin(); j != grad_mag.col_end(); ++j )
		{
			// slower, more precise
			// grad_mag(i,j) = hypot( Ix(i,j), Iy(i,j) );

			// faster approximation
			grad_mag(i,j) = fabs( Ix(i,j) ) + fabs( Iy(i,j) );
		}
	}

	// initialize positions relative to current pixel in which to test for maxima
	size_type prev_pos_x = 0;
	size_type prev_pos_y = 0;
	size_type next_pos_x = 0;
	size_type next_pos_y = 0;

	// initialize all edge values to all grad. mag. values--some edges will be
	// pruned in the suppression stage next
	Matrix edge_map = grad_mag;

	// skip the first and last row/column since for each pixel (i,j) we access 
	// (i,j)'s neighbors
	grad_mag.set_region( grad_mag.row_begin() + 1, grad_mag.col_begin() + 1,
			grad_mag.sub_row_size() - 2, grad_mag.sub_col_size() - 2 );

	// perform non-maxima suppression
	for( size_type i = grad_mag.row_begin(); i != grad_mag.row_end(); ++i )
	{
		for( size_type j = grad_mag.col_begin(); j != grad_mag.col_end(); ++j )
		{
			// gradient direction at this point (avoid division by 0)
			// TODO determine if need tan or atan here: some confusion
			// double grad_dir = atan2( Iy(i,j), Ix(i,j) );
			// if -PI < theta < 0, convert angle so it is between 0 and 360 degrees
			// if( grad_dir < 0 )
			// {
				// grad_dir += 2 * PI;
			// }
			double grad_dir = Iy(i,j) / (Ix(i,j) != 0) ? Ix(i,j) : .00001;

			// use direction of gradient to determine positions to save results 
			update_position( grad_dir, i, j, &prev_pos_x, &prev_pos_y,
					&next_pos_x, &next_pos_y );

			// suppress non-maxima pixel by marking as non-edge
			if( grad_mag(i, j) <= grad_mag( prev_pos_y, prev_pos_x )
					|| grad_mag(i, j) <= grad_mag( next_pos_y, next_pos_x ) )
			{
				edge_map( i, j ) = _Non_Edge_Pixel;
			}
		}
	}

	grad_mag.reset_region();

	return( edge_map );
}

/**
  	Determine gradient direction
	Compute positions into magnitude matrix to check for non-maxima
	suppression by comparing tan( theta ) to different ranges of values.

	Note that the x-axis corresponds to the column position (index j) and the
	y-axis corresponds to the row position (index i).

	@param[in] tan_theta Tan
	@param[in] i Row coordinate of current pixel
	@param[in] j Column coordinate of current pixel
	@param[out] prev_pos_x Previous pixel coordinate relative to current pixel
	@param[out] prev_pos_y Previous pixel coordinate relative to current pixel
	@param[out] next_pos_x Next pixel coordinate relative to current pixel
	@param[out] next_pos_y Next pixel coordinate relative to current pixel
 */
void
Edge_Detector::update_position( double tan_theta,
		size_type i, size_type j,
		size_type* prev_pos_x, size_type* prev_pos_y,
		size_type* next_pos_x, size_type* next_pos_y )
{
	// ensure pointers are valid
	assert( prev_pos_x != NULL && prev_pos_y != NULL &&
		next_pos_x != NULL && next_pos_y != NULL );

	// compare tan of theta to ranges to determine gradient direction
	if( -.4142 < tan_theta && tan_theta <= .4142 )
	{
		assert( j > 0 );

		// horizontal direction (0 degrees)
		*prev_pos_y = i;
		*prev_pos_x = j - 1;
		*next_pos_y = i;
		*next_pos_x = j + 1;
	}
	else if( .4142 < tan_theta && tan_theta < 2.4142 )
	{
		assert( i > 0 && j > 0 );

		// diagonal direction (45 degrees)
		*prev_pos_y = i - 1;
		*prev_pos_x = j - 1;
		*next_pos_y = i + 1;
		*next_pos_x = j + 1;
	}
	else if( 2.4142 <= fabs( tan_theta ) )
	{
		assert( i > 0 );

		// vertical direction (90 degrees)
		*prev_pos_y = i - 1;
		*prev_pos_x = j;
		*next_pos_y = i + 1;
		*next_pos_x = j;
	}
	else if( -2.4142 < tan_theta && tan_theta <= -.4142 )
	{
		assert( i > 0 && j > 0 );

		// diagonal direction (135 degrees)
		*prev_pos_y = i - 1;
		*prev_pos_x = j + 1;
		*next_pos_y = i + 1;
		*next_pos_x = j - 1;
	}
}

/**
	Hysteresis.

	Algorithm:
		for each pixel
			if pixel > high threshold
				mark pixel as edge
			if pixel < low threshold
				mark pixel as non-edge
			if low threshold <= pixel <= high threshold
				if pixel is adjacent to an edge pixel
					mark pixel as edge
					recursively mark each in-between pixel adjacent to this pixel as
							an edge

		mark any remaining in-between pixels as non-edges

	@param[out] img Output image
	@param[in] edges Map of edges
	@param[in,out] edges Map of edges
	@param[in] low_threshold Low threshold
	@param[in] high_threshold High threshold
 */
void
Edge_Detector::hysteresis( Image& img, const Matrix& edges,
		pixel_type low_threshold, pixel_type high_threshold )
{
	// determine low and high thresholds
	// pixel_type high_thresh = threshold;
	// pixel_type low_thresh  = .75 * threshold;

	// threshold = edges.get_max();

	if( low_threshold > high_threshold )
	{
		low_threshold = high_threshold;
	}

	// skip the first and last row/column since for each pixel (i,j) we access 
	// (i,j)'s neighbors
	Image_Region reg( edges.row_begin() + 1, edges.col_begin() + 1,
			edges.sub_row_size() - 2, edges.sub_col_size() - 2 );

	// threshold
	for( size_type i = reg.row_begin(); i != reg.row_end(); ++i )
	{
		for( size_type j = reg.col_begin(); j != reg.col_end(); ++j )
		{
			if( edges(i,j) < low_threshold )
			{
				img(i,j) = _Non_Edge_Pixel;
			}
			else if( edges(i,j) > high_threshold )
			{
				img(i,j) = _Edge_Pixel;

				// hysteresis( edges, i, j, low_threshold, high_threshold );
			}
			else // low_threshold <= edges(i,j) <= high_threshold
			{
				// if any neighbor of (i,j) is an edge pixel
				if( img(i - 1, j - 1)     == _Edge_Pixel
					|| img(i - 1, j)       == _Edge_Pixel
					|| img(i - 1, j + 1)   == _Edge_Pixel
					|| img(i, j - 1)       == _Edge_Pixel
					|| edges(i, j + 1)     >  high_threshold
					|| edges(i + 1, j - 1) >  high_threshold
					|| edges(i + 1, j)     >  high_threshold
					|| edges(i + 1, j + 1) >  high_threshold
				)
				{
					// recursively update adjacent pixels
					hysteresis( img, edges, i, j, low_threshold, high_threshold );
				}
				else
				{
					// mark in-between pixel as non-edges
					img(i,j) = _Non_Edge_Pixel;
				}
			}
		}
	}
}

/**
	Recursive hysteresis procedure.

	TODO Possible infinite loop if threshold is too high. Figure out iterative
	way to do it

	@param[in,out] edges Map of edges
	@param[in] low_thresh Low threshold
	@param[in] high_thresh High threshold
 */
void
Edge_Detector::hysteresis( Image& img, const Matrix& edges,
	size_type i, size_type j,
	pixel_type low_thresh, pixel_type high_thresh )
{
	// pixel is between the low and high thresholds, so mark as an edge
	img(i,j) = _Edge_Pixel;

	// recursively try each direction until reaching a pixel that is either
	// an edge pixel, along the image border, or is not between the low/high
	// thresholds

	// go to left pixel
	const size_type left = j - 1;
	if( img(i,left) != _Edge_Pixel
			&& j != edges.col_begin()
			&& (low_thresh <= edges(i,left) && edges(i,left) <= high_thresh) )
	{
		hysteresis( img, edges, i, left, low_thresh, high_thresh );
	}

	// go to top pixel
	const size_type top = i - 1;
	if( img(top,j) != _Edge_Pixel
			&& i != edges.row_begin()
			&& (low_thresh <= edges(top,j) && edges(top,j) <= high_thresh) )
	{
		hysteresis( img, edges, top, j, low_thresh, high_thresh );
	}

	// go to right pixel
	const size_type right = j + 1;
	if( img(i,right) != _Edge_Pixel
			&& j != (edges.col_end() - 1)
			&& (low_thresh <= edges(i,right) && edges(i,right) <= high_thresh) )
	{
		hysteresis( img, edges, i, right, low_thresh, high_thresh );
	}

	// go to bottom pixel
	const size_type bottom = i + 1;
	if( img(bottom,j) != _Edge_Pixel
			&& i != (edges.row_end() - 1)
			&& (low_thresh <= edges(bottom,j) && edges(bottom,j) <= high_thresh) )
	{
		hysteresis( img, edges, bottom, j, low_thresh, high_thresh );
	}
}
