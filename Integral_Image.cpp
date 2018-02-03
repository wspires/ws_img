/**
	@file   Integral_Image.cpp
	@author Wade Spires
	@date   2006/4/20
	@brief  Class Integral_Image.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Integral_Image.hpp"

using std::string;
using std::vector;

using namespace ws_img;

typedef Integral_Image::value_type value_type;

/**
	Construct integral image from given image.

	We store one additional row and column to simplify getting the area sum.
	@param[in] img Image to compute integral image of
 */
Integral_Image::Integral_Image( const Image& img )
: Image( img.sub_row_size() + 1, img.sub_col_size() + 1 )
{
	// initialize first row and column to 0
	for( size_type ii = row_begin(); ii != row_end(); ++ii )
	{
		at(ii,0) = 0;
	}
	for( size_type jj = col_begin(); jj != col_end(); ++jj )
	{
		at(0,jj) = 0;
	}

	// sum for all positions
	// i and j are input image positions; ii and jj are integral image positions
	for( size_type i = img.row_begin(), ii = (row_begin() + 1);
			i != img.row_end();
			++i, ++ii )
	{
		for( size_type j = img.col_begin(), jj = (col_begin() + 1);
				j != img.col_end();
				++j, ++jj )
		{
			at(ii,jj) = at(ii, jj - 1) + at(ii - 1, jj) + img(i,j)
				- at(ii - 1, jj - 1);
		}
	}
}
