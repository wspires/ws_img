/**
	@file   Image_Coordinate.hpp
	@author Wade Spires
	@date   2005/07/27
	@brief  Image coordinate.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef IMG_COORD_HPP
#define IMG_COORD_HPP

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <cstdio>

#include "ws_tools.hpp"
#include "Image.hpp"

namespace ws_img
{

/**
	@brief Row and column coordinate in an image
 */
struct Image_Coordinate
{
	unsigned  row;  //< Row coordinate
	unsigned  col;  //< Column coordinate

	Image_Coordinate( )
	: row(0), col(0)
	{ }

	Image_Coordinate( unsigned r, unsigned c )
	: row(r), col(c)
	{ }

	/**
		Compare two coordinates.
		@param[in] coord Coordinate to compare this coordinate to
	 */
	bool operator<( const Image_Coordinate& coord ) const
	{
		if( row < coord.row )
		{
			return( true );
		}
		else if( row == coord.row && col < coord.col )
		{
			return( true );
		}
		return( false );
	}
};

} // namespace ws_img

#endif // IMG_COORD_HPP
