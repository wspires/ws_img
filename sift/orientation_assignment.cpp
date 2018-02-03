/**
	@file   orientation_assignment.cpp
	@author Wade Spires
	@date   2006/2/22
	@brief  Class orientation_assignment.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "orientation_assignment.hpp"

using std::string;
using std::vector;

// my namespace
using namespace ws_img;

typedef Image::size_type   size_type;
typedef Image::value_type  value_type;
typedef Image::pixel_type  pixel_type;
typedef Image_Coordinate   Img_Coord;

const ws_img::Matrix::size_type  NUM_BINS  = 36;
const ws_img::Matrix::value_type BIN_WIDTH = (2 * PI) / NUM_BINS;

// filter for all histogram filtering
ws_img::Matrix Gauss;

/**
	3. Orientation assignment: One or more orientations are assigned to each
	keypoint location based on local image gradient directions. All future
	operations are performed on image data that has been transformed relative to
	the assigned orientation, scale, and location for each feature, thereby
	providing invariance to these transformations.  

	@param[in,out] scale_space
	@param[in,out] keypoints
 */
Keypoint_List
orientation_assignment( const Scale_Space& scale_space,
	const Keypoint_List& keypoints )
{
	log_msg( "\nOrientation assignment\n" );

	Keypoint_List new_keypoints;

	const Matrix::size_type num_bins = NUM_BINS;

	// dimension of region to construct histogram
	const size_type reg_size = 16;

	// Gaussian filter
	Gauss = make_filter( reg_size );

	// examine each keypoint
	for( unsigned k = 0; k != keypoints.size(); ++k )
	{
		const size_type  y      = keypoints[k].img_coord.row;
		const size_type  x      = keypoints[k].img_coord.col;
		const size_type  scale  = keypoints[k].scale_pos;
		const size_type  octave = keypoints[k].octave_pos;
		const Image&     L      = scale_space[octave][scale].L; 

		// skip images in the last scale space image since it is empty
		if( scale == (scale_space[octave].size() - 1) )
		{
			continue;
		}

		// compute gradient magnitude and direction around given point
		Matrix mag, dir;
		grad_mag_dir( L, y, x, reg_size, mag, dir );

		if( mag.row_size() <= 1 && mag.row_size() <= 1 )
		{
			continue;
		}

// Image( mag ).normalize().write( "norm.pgm" );
// dir.write();
// Ix.write();
// mag.write();
// err_quit( "\n" );

		// const value_type sigma = scale_space[octave][scale].s;
// fprintf( stderr, "%u %u %u %u %lf\n", octave, scale, y, x, sigma );
// fprintf( stderr, "%u %u %u %u\n", octave, scale, y, x );

		// create histogram around point
		Vector hist = make_histogram( mag, dir, num_bins );

		// add keypoints to list using the peaks in the histogram
		assign_points_using_hist_peaks( hist, keypoints[k], new_keypoints,
				mag, dir );
// hist.write();
// err_quit( "\n" );
	}

	log_msg( "   Keypoints: %u, Old keypoints: %u\n",
			new_keypoints.size(), keypoints.size() );


	return( new_keypoints );
}

/**
	Compute gradient magnitude and direction of image.

	@param[in] img Image
	@param[out] mag Gradient magnitude
	@param[out] dir Gradient direction
 */
void
grad_mag_dir( const Image& img, size_type y, size_type x, size_type reg_size,
	Matrix& mag, Matrix& dir )
{
	// ensure a 17x17 region can fit centered at (x,y) (the additional row and
	// column is necessary since we are taking forward differences)
	size_type half_size = reg_size / 2;
	int row_begin     = y - half_size;
	int col_begin     = x - half_size;
	size_type row_end = y + half_size;
	size_type col_end = x + half_size;
	if( row_begin < 0 || col_begin < 0 || row_end >= img.sub_row_size()
			|| col_end >= img.sub_col_size() )
	{
		return;
	}
	Image_Region reg( row_begin, col_begin, reg_size, reg_size );

	// compute derivatives in 16x16 neighborhood around point (y,x)
	Matrix Ix( reg_size ), Iy( reg_size );
	for( size_type i = reg.row_begin(), ii = 0; i != reg.row_end(); ++i, ++ii )
	{
		for( size_type j = reg.col_begin(), jj = 0; j != reg.col_end(); ++j, ++jj )
		{
			Iy(ii,jj) = img(i + 1,j) - img(i,j);
			Ix(ii,jj) = img(i,j + 1) - img(i,j);
		}
	}

	// allocate results
	mag = dir = Matrix( Ix.sub_row_size(), Ix.sub_col_size() );

	// compute gradient magnitude
	for( size_type i = Ix.row_begin(); i != Ix.row_end(); ++i )
	{
		for( size_type j = Ix.col_begin(); j != Ix.col_end(); ++j )
		{
			// gradient magnitude: slower, more precise
			mag(i,j) = hypot( Ix(i,j), Iy(i,j) );

			// gradient magnitude: faster approximation
			// mag(i,j) = fabs( Ix(i,j) ) + fabs( Iy(i,j) );

			// gradient direction
			value_type theta = atan2( Iy(i,j), Ix(i,j) );

			// if -PI < theta < 0, convert angle to be between 0 and 2 * PI
			if( theta < 0 )
			{
				theta += 2 * PI;
			}

			// ensure angles are in the range [0, 2 * PI)
			if( theta >= 2 * PI )
			{
				theta = 0;
			}
			dir(i,j) = theta;
		}
	}
}

/**
	Compute gradient magnitude and direction.

	@param[in] grad_mag Gradient magnitude
	@param[in] grad_dir Gradient direction
	@param[in] num_bins Number of bins in histogram
 */
Vector
make_histogram( const Matrix& grad_mag, const Matrix& grad_dir,
	const Matrix::size_type num_bins )
{
	// allocate histogram with the given number of bins
	Vector hist( num_bins );

	// number of angles that are mapped to each bin
	double bin_width = (2 * PI) / num_bins;

// gauss.write();
// fprintf( stderr, "%u %u\n", gauss.row_size(), gauss.col_size() );
// err_quit( "\n" );

	// construct histogram of angles
	for( size_type i = grad_dir.row_begin(); i != grad_dir.row_end(); ++i )
	{
		for( size_type j = grad_dir.col_begin(); j != grad_dir.col_end(); ++j )
		{
			// value_type hist_count = grad_mag(i,j);
			value_type hist_count = grad_mag(i,j) * Gauss(i,j);
			value_type theta      = grad_dir(i,j);

			// increment histogram bin
			size_type bin = static_cast<size_type>( round( theta / bin_width ) );
			hist( bin ) += hist_count;
		}
	}

// hist.write();
// err_quit( "\n" );

	return( hist );
}

/**
	Add keypoints to list using the peaks in the histogram.
	@param[in] hist Histogram
	@param[in] old_keypoint Old keypoint
	@param[out] new_keypoints List of new keypoints
 */
void
assign_points_using_hist_peaks( const Vector& hist, const Keypoint& old_keypoint,
	Keypoint_List& new_keypoints, const Matrix& mag, const Matrix& dir )
{
	// find local peaks in histogram that are within 80% of the largest peak
	value_type max = hist.get_max();
	if( max == 0 )
	{
		return;
	}
	for( size_type i = hist.vec_begin(); i != hist.vec_end(); ++i )
	{
		// must be within 80% of the max
		if( hist(i) >= .8 * max )
		{
			// must be a local peak (i.e., greater than its adjacent neighbor)
			if( (i > hist.vec_begin() && hist(i) <= hist(i - 1))
					|| (i < (hist.vec_end() - 1) && hist(i) <= hist(i + 1)) )
			{
				continue;
			}

			// assign initial orientation
			double orientation = i * BIN_WIDTH;

			// fit parabola if we have 2 values to the left and right; otherwise,
			// just use the orientation of the bin directly
			if( i > hist.vec_begin() && i < (hist.vec_end() - 1) )
			{
				value_type center = fit_parabola( i - 1, i, i + 1,
						hist(i - 1), hist(i), hist(i + 1) );
				orientation = center * BIN_WIDTH;
			}

			// add keypoint with given orientation to new list
			Keypoint keypoint    = old_keypoint;
			keypoint.grad_mag    = mag;
			keypoint.grad_dir    = dir;
			keypoint.orientation = orientation;
			new_keypoints.push_back( keypoint );
		}
	}
}

/**
	Fit parabola to points.

	@param[in] x1
	@param[in] x2
	@param[in] x3
	@param[in] y1
	@param[in] y2
	@param[in] y3
	@retval x0 Location of the center parabola (position of max. value)
 */
double
fit_parabola( double x1, double x2, double x3, double y1, double y2, double y3 )
{
	// 1. solve Ax = b where A are positions under a parabola p(x)
	Matrix A(3, 3);
	A(0,0) = x1 * x1;
	A(0,1) = x1;
	A(0,2) = 1;
	A(1,0) = x2 * x2;
	A(1,1) = x2;
	A(1,2) = 1;
	A(2,0) = x3 * x3;
	A(2,1) = x3;
	A(2,2) = 1;

	// b = p(x) for all x above
	Vector b(3);
	b(0) = y1;
	b(1) = y2;
	b(2) = y3;

	// solve 3-by-3 (using inverse() defined in keypoint_localization.hpp)
	Matrix A_inv;
	if( !inverse(A, A_inv) )
	{
		// return the original center if the matrix is ill-conditioned
		// (not likely)
		return( x2 );
	}
	Vector x = A_inv * b;

	// 2. find center: x0 = -b / 2a since p'(x) = 2a * x + b
	double x0 = -x(1) / (2 * x(0));

	// 3. evaluate p(x0) = a x0^2 + b x0 + c
	// Note: we only need x0, not y0 = p(x0)
	// value_type p_x0 = (x(0) * x0 * x0) + (x(1) * x0) + x(2);

	return( x0 );
}

/**
	Create Gaussian filter with given size.

	@param[in] filter_size Size of Gaussian filter
	@retval gauss Gaussian filter
 */
Matrix
make_filter( size_type filter_size )
{
	Matrix gauss( filter_size );

	// version 1.
	// create a Gaussian mask whose size allows it to overlap the mag. matrix and
	// with the given std. dev. (by default, the returned filter excludes 0s, so
	// we must create the matrix gauss with the 0s included explicitly)
	Matrix g = make_gauss_2D_filter( 3 * Scale_Space::DEFAULT_SIGMA );
	// size_type gauss_center = g.row_size() / 2;
	for( size_type i = g.row_begin(); i != g.row_end(); ++i )
	{
		for( size_type j = g.col_begin(); j != g.col_end(); ++j )
		{
			// gauss( i + gauss_center, j + gauss_center ) = g(i,j);
			gauss(i, j) = g(i,j);
		}
	}

	// version 2.
	// create a Gaussian mask whose size allows it to overlap the mag. matrix,
	// we have sigma = 3 * 1.5 = 4.5 -> filter size = round(2 * 3 * 4.5) + 1
	// = 15 for a symmetric filter with a spread of 3 and of odd length
	// assert( grad_mag.row_size() == 15 && grad_mag.col_size() == 15 );
	// static const Matrix gauss = make_gauss_2D_filter( 3 * 1.5, 3 );

	return( gauss );
}
