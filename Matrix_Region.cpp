/**
	@file   Matrix_Region.cpp
	@author Wade Spires
	@date   2006/12/16
	@brief  Class Matrix_Region.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Matrix_Region.hpp"
#include "ws_tools.hpp"

using namespace ws_img;
using namespace ws_tools;

/**
	Construct region starting at position (0,0) that is has no rows or columns.
 */
Matrix_Region::Matrix_Region( )
: _row_begin(0), _col_begin(0), _row_size(0), _col_size(0),
	_row_end(0), _col_end(0)
{ }

/**
	Construct region with the specified number of rows and columns.
	@param[in] row_size Number of rows in the region
	@param[in] col_size Number of columns in the region
 */
Matrix_Region::Matrix_Region( size_type row_size, size_type col_size )
: _row_begin(0), _col_begin(0), _row_size(row_size), _col_size(col_size),
	_row_end(row_size), _col_end(col_size)
{ }

/**
	Construct region starting at given row and column positions and with the
	specified sizes.

	@param[in] row_begin Starting row position of the region
	@param[in] col_begin Starting col position of the region
	@param[in] row_size Number of rows in the region
	@param[in] col_size Number of columns in the region
 */
Matrix_Region::Matrix_Region( size_type row_begin, size_type col_begin,
	size_type row_size, size_type col_size )
: _row_begin(row_begin), _col_begin(col_begin),
	_row_size(row_size), _col_size(col_size),
	_row_end(row_begin + row_size), _col_end(col_begin + col_size)
{ }

/**
	Construct region from string--inverse of to_string().
	@param[in] region String representation of region
 */
Matrix_Region::Matrix_Region( const std::string& region )
: _row_begin(0), _col_begin(0), _row_size(0), _col_size(0),
	_row_end(0), _col_end(0)
{
	// format is [row_begin,col_begin]->[row_end,col_end]
	std::vector<std::string> words = split_string( region, "[,]->" );
	if( words.size() != 4 )
	{
		err_quit( "Invalid region: %s\n", region.c_str() );
	}

	_row_begin = size_type( string_to_double(words[0]) );
	_col_begin = size_type( string_to_double(words[1]) );
	_row_size  = size_type( string_to_double(words[2]) ) - _row_begin + 1;
	_col_size  = size_type( string_to_double(words[3]) ) - _col_begin + 1;

	_row_end = row_begin() + row_size();
	_col_end = col_begin() + col_size();
}
