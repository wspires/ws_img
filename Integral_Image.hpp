/**
	@file   Integral_Image.hpp
	@author Wade Spires
	@date   2006/4/20
	@brief  Class Integral_Image.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef INTEGRAL_IMAGE_HPP
#define INTEGRAL_IMAGE_HPP

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include <algorithm>
#include <string>
#include <vector>

// c headers
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// my headers
#include "Image.hpp"
#include "Image_Region.hpp"
#include "ws_tools.hpp"

namespace ws_img
{

/**
	@brief Integral_Image
 */
class Integral_Image : public Image
{

public:

	/**
		Default constructor.
	Integral_Image( )
	{ }
	 */

	/**
		Construct integral image from given image.
	 */
	Integral_Image( const Image& );

	/**
		Destructor does nothing since no member variables are dynamically
		allocated.
	 */
	virtual ~Integral_Image( )
	{ }

	/**
		Access sum of current region starting at row i and column j.
		The pixel sum of any upright rectangle r = (x, y, w, h) can be determined
		by 4 table lookups in the Summed Area Table (SAT):
			RecSum(r) = SAT(x - 1, y - 1)
				+ SAT(x + w - 1, y + h - 1)
				- SAT(x, y + h -1)
				- SAT(x + w - 1, y - 1)

		@param[in] i Starting row position
		@param[in] j Starting column position
		@param[in] h Height of region
		@param[in] w Width of region
		@retval area Area sum of region
	 */
	inline value_type get_area( size_type i, size_type j, size_type h,
			size_type w  ) const
	{
		// we store an additional, unused row and column at the start
		// of an integral image to simplify this computation, so we increment the
		// starting position by one
		// ++i;
		// ++j;

		// value_type area = at(i - 1, j - 1)
			// + at(i + h - 1, j + w - 1)
			// - at(i - 1, j + w - 1)
			// - at(i + h - 1, j - 1);

		// no need to subtract 1 since we store an additional row and column
		const value_type area = at(i, j)
			+ at(i + h, j + w)
			- at(i, j + w)
			- at(i + h, j);

		return( area );
	}

	/**
	  	Get area sum of region.

		@param[in] region Region to find area sum for
		@retval area Area sum of region
	 */
	inline value_type get_area( const Image_Region& region ) const
	{
		return(
			get_area( region.row_begin(), region.col_begin(), region.row_size(),
				region.col_size() )
		);
	}

	/*
		Override Image's size functions to compensate for the additional row and
		column that are stored for an Integral Image (and unused by users of the
		class).
	 */

	/**
		Get the number of rows in the image.
		@retval num_rows Number of rows in image
	 */
	inline size_type row_size( ) const
	{
		return( Image::row_size() - 1 );
	}

	/**
		Get the number of rows in the image.
		@retval height Height of image
	 */
	inline size_type height( ) const
	{ return( row_size() ); }

	/**
		Get the number of rows in the image.
		@retval num_rows Number of rows in image
	 */
	inline size_type get_num_rows( ) const
	{ return( row_size() ); }

	/**
		Get the number of columns in the image.
		@retval num_cols Number of columns in image
	 */
	inline size_type col_size( ) const
	{
		return( Image::col_size() - 1 );
	}

	/**
		Get the number of columns in the image.
		@retval width Width of image
	 */
	inline size_type width( ) const
	{ return( col_size() ); }

	/**
		Get the number of columns in the image.
		@retval num_cols Number of columns in image
	 */
	inline size_type get_num_cols( ) const
	{
		return( col_size() );
	}
};

} // namespace ws_img

#endif // INTEGRAL_IMAGE_HPP 
