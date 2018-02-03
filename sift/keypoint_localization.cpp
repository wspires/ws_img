/**
	@file   keypoint_localization.cpp
	@author Wade Spires
	@date   2006/2/12
	@brief  Functions for localizing keypoints.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "keypoint_localization.hpp"

using std::string;
using std::vector;

using namespace ws_img;

typedef Image::size_type   size_type;
typedef Image::value_type  value_type;
typedef Image::pixel_type  pixel_type;
typedef Image_Coordinate   Img_Coord;

/**
	2. Keypoint localization: At each candidate location, a detailed model is fit
	to determine location and scale. Keypoints are selected based on measures of
	their stability.  

	@param[in] keypoints Current set of keypoints
	@param[in] scale_space Scale space
	@param[in] threshold Threshold
	@retval new_keypoints New set of keypoints
 */
Keypoint_List
keypoint_localization( const Keypoint_List& keypoints,
		const Scale_Space& scale_space, double threshold, double r )
{
	log_msg( "\nLocalizing keypoints\n" );

	Keypoint_List new_keypoints;

	// examine each keypoint
	for( unsigned k = 0; k != keypoints.size(); ++k )
	{
		// get keypoint information
		const size_type y      = keypoints[k].img_coord.row;
		const size_type x      = keypoints[k].img_coord.col;
		const size_type scale  = keypoints[k].scale_pos;
		const size_type octave = keypoints[k].octave_pos;

		// verify that the keypoint position is valid
		assert( octave < scale_space.size() && scale < scale_space[octave].size()
				&& y < scale_space[octave][scale].D.row_size()
				&& x < scale_space[octave][scale].D.col_size()
		);

		// get DoG image that this keypoint lies in
		const Image& D = scale_space[ octave ][ scale ].D;

		// verify that a neighborhood exists around the point (y,x,sigma)
		if( y <= 0 || y >= (D.row_size() - 1)
				|| x <= 0 || x >= (D.col_size() - 1)
				|| scale <= 0 || scale >= (scale_space[octave].size() - 1) )
		{
			continue;
		}

		// get DoG images around this keypoint for computing derivatives in sigma
		// (B)ottom and (T)op
		const Image& B = scale_space[ octave ][ scale - 1 ].D;
		const Image& T = scale_space[ octave ][ scale + 1 ].D;

		assert( y > 0 && y < (D.row_size() - 1)
				&& x > 0 && x < (D.col_size() - 1) );

		// first order derivative--partials in y, x, and sigma directions
		Vector dx( 3 );
		dx(0) = D(y + 1, x) - D(y, x);
		dx(1) = D(y, x + 1) - D(y, x);
		dx(2) = T(y, x) - D(y, x);

		// Hessian matrix (second partials in y, x, and sigma directions) 
		Matrix H( 3, 3);

		// second partial derivative with respect to y^2, yx, ys
		H(0,0) = D(y + 1,x) - 2 * D(y,x) + D(y - 1,x);
		H(0,1) = (D(y + 1,x + 1) - D(y,x + 1)) - (D(y + 1,x) - D(y,x));
		H(0,2) = (T(y + 1,x) - T(y,x)) - (D(y + 1,x) - D(y,x));

		// second partial derivative with respect to xy, x^2, xs
		H(1,0) = (D(y + 1,x + 1) - D(y + 1,x)) - (D(y,x + 1) - D(y,x));
		H(1,1) = D(y,x + 1) - 2 * D(y,x) + D(y,x - 1);
		H(1,2) = (T(y,x + 1) - T(y,x)) - (D(y,x + 1) - D(y,x));

		// second partial derivative with respect to sy, sx, s^2
		H(2,0) = (T(y + 1,x) - D(y + 1,x)) - (T(y,x) - D(y,x));
		H(2,1) = (T(y,x + 1) - D(y,x + 1)) - (T(y,x) - D(y,x));
		H(2,2) = T(y,x) - 2 * D(y,x) + B(y,x);

		// solve the system H * x_hat = -dx for x_hat
		// we have 3 ways to solve this system

		// 1. solve system by finding the inverse directly; this is easy and fast
		// for a 3x3, but the matrix may be singular, so we the skip point if so
		Matrix H_inv;
		if( !inverse( H, H_inv ) )
		{
			continue;
		}
		Vector x_hat = H_inv * -dx;
		
		// 2. solve system using Gaussian elimination
		// Vector x_hat;
		// Lin_Alg::solve_system( H, -dx, x_hat );

		// 3. solve system using SVD; this method automatically applies
		// least-squares fitting if H is singular, but is much slower
		// Vector x_hat;
		// Lin_Alg::solve_system_SVD( H, -dx, x_hat );

		// threshold using the contrast: D(x_hat) = D + 1/2 * (dx * x_hat)
		double D_x_hat = fabs( D(y,x) + .5 * dx.dot_product(x_hat) );
		if( D_x_hat < threshold )
		{
			continue;  // skip point
		}

		// if the point moved too far
		double limit = 2;
		if( x_hat(0) > limit || x_hat(1) > limit || x_hat(2) > limit )
		{
			continue;  // skip point
		}

		// eliminate edge responses:
		// (1) find 2x2 Hessian H = [ Dxx Dxy; Dxy Dyy ]
		// (2) compute Tr(H) = Dxx + Dyy; Det(H) = DxxDyy - (Dxy)^2
		// (3) accept if Tr(H)^2 / Det(H) < (r + 1)^2 / r
		// Note: we have already computed these derivatives above
		value_type trace = H(0,0) + H(1,1);
		value_type det   = (H(0,0) * H(1,1)) - (H(0,1) * H(1,0));
		if( det <= 0 || (trace * trace) / det >= ((r + 1) * (r + 1)) / r )
		{
			continue;  // skip point
		}


		// add interpolated keypoint
		// compute new position (only changes if x_hat > .5)
		double new_row   = round( y + x_hat(0) );
		double new_col   = round( x + x_hat(1) );
		double new_scale = round( scale + x_hat(2) );

		// ensure that the new position is within the bounds of the scale space
		if( new_row < 0 )
		{
			new_row = 0;
		}
		else if( new_row >= D.sub_row_size() )
		{
			new_row = D.sub_row_size() - 1;
		}
		if( new_col < 0 )
		{
			new_col = 0;
		}
		else if( new_col >= D.sub_col_size() )
		{
			new_col = D.sub_col_size() - 1;
		}
		if( new_scale < 0 )
		{
			new_scale = 0;
		}
		else if( new_scale >= scale_space[octave].size() )
		{
			// don't use the last image since it has no DoG
			new_scale = scale_space[octave].size() - 2;
		}

		// save new keypoint
		Keypoint new_keypoint;
		new_keypoint.octave_pos    = octave;
		new_keypoint.scale_pos     = static_cast<size_type>( new_scale );
		new_keypoint.img_coord.row = static_cast<size_type>( new_row );
		new_keypoint.img_coord.col = static_cast<size_type>( new_col );
		new_keypoints.push_back( new_keypoint );
	}

	log_msg( "   Keypoints: %u, Non-keypoints: %u, Total: %u\n",
			new_keypoints.size(), keypoints.size() - new_keypoints.size(),
			keypoints.size()
		);

	return( new_keypoints );
}

/**
	Find inverse of 3-by-3 matrix A directly.
	@param[in] A Matrix to invert
	@param[out] A_inv Inverse of A
	@retval is_nonsingular Whether matrix is nonsingular (i.e., returns false if
		det(A) = 0; A_inv is invalid in this case)
 */
bool
inverse( const Matrix& A, Matrix& A_inv )
{
	if( A.sub_row_size() != 3 || A.sub_col_size() != 3 )
	{
		err_quit( "Matrix dimensions (%u x %u) must be 3-by-3",
				A.sub_row_size(), A.sub_col_size() );
	}

	// easier to deal with short names
	const Matrix::value_type a = A(0,0);
	const Matrix::value_type b = A(0,1);
	const Matrix::value_type c = A(0,2);
	const Matrix::value_type d = A(1,0);
	const Matrix::value_type e = A(1,1);
	const Matrix::value_type f = A(1,2);
	const Matrix::value_type g = A(2,0);
	const Matrix::value_type h = A(2,1);
	const Matrix::value_type i = A(2,2);

	// compute determinant
	double det = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
	if( det == 0 )
	{
		return( false );
	}

	// compute inverse
	A_inv = Matrix( 3, 3 );
	A_inv(0,0) = e * i - f * h;
	A_inv(0,1) = c * h - b * i;
	A_inv(0,2) = b * f - c * e;
	A_inv(1,0) = f * g - d * i;
	A_inv(1,1) = a * i - c * g;
	A_inv(1,2) = c * d - a * f;
	A_inv(2,0) = d * h - e * g;
	A_inv(2,1) = b * g - a * h;
	A_inv(2,2) = a * e - b * d;

	A_inv /= det;

	return( true );
}
