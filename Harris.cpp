/**
	@file   Harris.cpp
	@author Wade Spires
	@date   2006/4/8
	@brief  Class Harris.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Harris.hpp"

using std::string;
using std::vector;

using ws_img::Image;
using ws_img::Image_Region;
using ws_img::Matrix;
using ws_img::Vector;

typedef Image::size_type  size_type;
typedef Image::value_type value_type;

/**
  	Harris corner detector.

	TODO Return list of points instead of drawing on the image.

	@param[in,out] img Image to find corners on
	@param[in] reg_size Size of neighborhood
	@param[in] threshold Threshold on point to determine if its a corner
	@param[in] kappa Parameter
 */
void
harris( Image& img, size_type reg_size, value_type threshold,
		const double kappa )
{
	if( img.sub_row_size() <= reg_size || img.sub_col_size() <= reg_size )
	{
		return;
	}

	img.to_gray();

	// smooth image with Gaussian
	const Vector filter = ws_img::make_gauss_filter( 1.0 );
	img.convolve( filter, filter );

	// Matrix Ix( img.sub_row_size(), img.sub_col_size() );
	// Matrix Iy( img.sub_row_size(), img.sub_col_size() );
	Matrix IxIx( img.sub_row_size(), img.sub_col_size() );
	Matrix IxIy( img.sub_row_size(), img.sub_col_size() );
	Matrix IyIy( img.sub_row_size(), img.sub_col_size() );

	// compute derivative images
	Matrix::size_type ii = IxIx.row_begin();
	for( size_type i = img.row_begin(); i != (img.row_end() - 1); ++i, ++ii )
	{
		Matrix::size_type jj = IxIx.col_begin();
		for( size_type j = img.col_begin(); j != (img.col_end() - 1); ++j, ++jj )
		{
			// TODO should we compute 1st order or 2nd order derivative?

			// compute 1st order derivative in x and y directions
			Matrix::value_type Ix = img(i,j + 1) - img(i,j);
			Matrix::value_type Iy = img(i + 1,j) - img(i,j);
			IxIx(ii,jj) = Ix * Ix;
			IxIy(ii,jj) = Ix * Iy;
			IyIy(ii,jj) = Iy * Iy;
		}
	}

	// Ix.convolve( filter, filter );
	// Iy.convolve( filter, filter );

	// Image( Ix ).write( "Ix.pgm" );
	// Image( Iy ).write( "Iy.pgm" );

	if( ws_tools::is_even( reg_size ) )
	{
		++reg_size;
	}

	// limit the scanning region
	Image_Region scan_reg = IxIx.get_region();
	scan_reg.set_row_begin( reg_size / 2 );
	scan_reg.set_col_begin( reg_size / 2 );
	scan_reg.set_row_size( scan_reg.row_size() - reg_size - 1 );
	scan_reg.set_col_size( scan_reg.col_size() - reg_size - 1 );

	Image R_img( img.sub_row_size(), img.sub_col_size() );
	Image_Region box(0,0,5,5);

	Image_Region local_reg( reg_size, reg_size );
	for( size_type i = scan_reg.row_begin(); i != scan_reg.row_end(); ++i )
	{
		for( size_type j = scan_reg.col_begin(); j != scan_reg.col_end(); ++j )
		{
			Matrix::value_type sum_IxIx = 0;
			Matrix::value_type sum_IxIy = 0;
			Matrix::value_type sum_IyIy = 0;

			// look at local neighborhood
			local_reg.set_begin( i, j );
			for( size_type m = local_reg.row_begin();
					m != local_reg.row_end();
					++m )
			{
				for( size_type n = local_reg.col_begin();
						n != local_reg.col_end();
						++n )
				{
					sum_IxIx += IxIx(m,n);
					sum_IxIy += IxIy(m,n);
					sum_IyIy += IyIy(m,n);
				}
			}

			// compute determinant and trace of M = [ Ix^2 IxIxy; IxIy Iy^2 ]
			const value_type det   = (sum_IxIx * sum_IyIy) - (sum_IxIy * sum_IxIy);
			const value_type trace = (sum_IxIx + sum_IyIy);

			// compute corner measure
			const value_type R = det + (kappa * trace * trace);
			R_img(i,j) = R;
			
			// if( R > threshold )
			{
				// ++count;
				// box.set_row_begin(i);
				// box.set_col_begin(j);
				// img.draw2( box );
			}
		}
	}

	// img.to_color();
	R_img.normalize();
	// R_img.write( "R_img.pgm" );
	for( size_type i = scan_reg.row_begin(); i != scan_reg.row_end(); ++i )
	{
		for( size_type j = scan_reg.col_begin(); j != scan_reg.col_end(); ++j )
		{
			// R_img(i,j) is an interest point if it exceeds some threshold
			if( R_img(i,j) > threshold )
			{
				box.set_row_begin(i);
				box.set_col_begin(j);
				img.draw( box );

				// img.red(i,j)   = 255;
				// img.green(i,j) = 0;
				// img.blue(i,j)  = 0;
			}
		}
	}
	/*
	 */

	// img.write( "test.pgm" );
}
