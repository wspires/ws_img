/**
	@file   Matrix_Region.hpp
	@author Wade Spires
	@date   2005/10/28
	@brief  Matrix region.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef MATRIX_REGION_HPP
#define MATRIX_REGION_HPP

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include <cstdio>
#include <cmath>

namespace ws_img
{

/**
	@brief Region inside (possibly) larger matrix to use for operating on only a
	subsection of a matrix.
 */
class Matrix_Region
{

public:

	/// Type that can represent position and size of a matrix
	typedef unsigned size_type;

private:

	size_type  _row_begin;  //< Row position region starts in matrix region
	size_type  _col_begin;  //< Column position region starts in matrix region

	size_type  _row_size;   //< Number of rows in matrix region
	size_type  _col_size;   //< Number of columns in matrix region

	size_type  _row_end;    //< One past the last valid row in matrix region
	size_type  _col_end;    //< One past the last valid column in matrix region

public:

	/**
		Construct region starting at position (0,0) that is has no rows or columns.
	 */
	Matrix_Region( );

	/**
		Construct region with the specified number of rows and columns.
		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region
	 */
	Matrix_Region( size_type, size_type );

	/**
		Construct region starting at given row and column positions and with the
		specified sizes.

		@param[in] row_begin Starting row position of the region
		@param[in] col_begin Starting col position of the region
		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region
	 */
	Matrix_Region( size_type, size_type, size_type, size_type );

	/**
		Construct region from string--inverse of to_string().
		@param[in] region String representation of region
	 */
	Matrix_Region( const std::string& );

	/**
		Get starting row position of the region.
		@retval row_begin Starting row position of the region
	 */
	size_type row_begin( ) const;

	/**
		Get starting column position of the region.
		@retval col_begin Starting column position of the region
	 */
	size_type col_begin( ) const;

	/**
		Get number of rows in region.
		@retval row_size Number of rows in region
	 */
	size_type row_size( ) const;

	/**
		Get number of columns in region.
		@retval col_size Number of columns in region
	 */
	size_type col_size( ) const;

	/**
		Get ending row position of the region.
		@retval row_end Ending row position of the region
	 */
	size_type row_end( ) const;

	/**
		Get ending column position of the region.
		@retval col_end Ending column position of the region
	 */
	size_type col_end( ) const;

	/**
		Get left-most position of the region.
		@retval left Left-most position of the region
	 */
	size_type left( ) const;

	/**
		Get right-most valid position of the region.
		@retval right Right-most position of the region
	 */
	size_type right( ) const;

	/**
		Get top-most position of the region.
		@retval top Top-most position of the region
	 */
	size_type top( ) const;

	/**
		Get bottom-most valid position of the region.
		@retval bottom Bottom-most position of the region
	 */
	size_type bottom( ) const;

	/**
		Set starting row position of the region.
		@param[in] row_begin Starting row position of the region
		@retval this This region
	 */
	Matrix_Region& row_begin( size_type );

	/**
		Set starting row position of the region.
		@param[in] row_begin Starting row position of the region
		@retval this This region
	 */
	Matrix_Region& set_row_begin( size_type );

	/**
		Set starting column position of the region.
		@param[in] col_begin Starting column position of the region
		@retval this This region
	 */
	Matrix_Region& col_begin( size_type );

	/**
		Set starting column position of the region.
		@param[in] col_begin Starting column position of the region
		@retval this This region
	 */
	Matrix_Region& set_col_begin( size_type );

	/**
		Set number of rows in region.
		@param[in] row_size Number of rows in region
		@retval this This region
	 */
	Matrix_Region& row_size( size_type );

	/**
		Set number of rows in region.
		@param[in] row_size Number of rows in region
		@retval this This region
	 */
	Matrix_Region& set_row_size( size_type );

	/**
		Set number of columns in region.
		@param[in] col_size Number of columns in region
		@retval this This region
	 */
	Matrix_Region& col_size( size_type );

	/**
		Set number of columns in region.
		@param[in] col_size Number of columns in region
		@retval this This region
	 */
	Matrix_Region& set_col_size( size_type );

	/**
		Set the region's starting row and column positions.

		@param[in] row_begin Starting row position of the region
		@param[in] col_begin Starting col position of the region
		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region
		@retval this This region
	 */
	Matrix_Region& begin( size_type, size_type );

	/**
		Set the region's starting row and column positions.

		@param[in] row_begin Starting row position of the region
		@param[in] col_begin Starting col position of the region
		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region
		@retval this This region
	 */
	Matrix_Region& set_begin( size_type, size_type );

	/**
		Set the region's row and column sizes.

		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region
		@retval this This region
	 */
	Matrix_Region& size( size_type, size_type );

	/**
		Set the region's row and column sizes.

		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region
		@retval this This region
	 */
	Matrix_Region& set_size( size_type, size_type );

	/**
		Set the region's starting row and column positions and sizes.

		@param[in] row_begin Starting row position of the region
		@param[in] col_begin Starting col position of the region
		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region
		@retval this This region
	 */
	Matrix_Region& operator()( size_type, size_type, size_type, size_type );

	/**
		Set the region's starting row and column positions and sizes.

		@param[in] row_begin Starting row position of the region
		@param[in] col_begin Starting col position of the region
		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region
		@retval this This region
	 */
	Matrix_Region& set( size_type, size_type, size_type, size_type );

	/**
		Compare two regions for equality.

		@retval equal Whether the given region is the same as this region
	 */
	bool operator==( const Matrix_Region& ) const;

	/**
		Compare two regions for inequality.

		@retval equal Whether the given region is not the same as this region
	 */
	bool operator!=( const Matrix_Region& ) const;

	/**
		Get area size of region.
		@retval area Area of region
	 */
	size_type area( );

	/**
		Get perimeter size around region.
		@retval perimeter perimeter of region
	 */
	size_type perimeter( );

	/**
		Determine if this region lies inside another region.

		For example,
		----------
		|  ---   |
		|  | |   |
		|  ---   |
		|        |
		----------
		returns true.

		@param[in] region Matrix_Region
		@retval result Whether this region lies inside of the given region
    */
   bool inside( const Matrix_Region& ) const;

	/**
		Determine if this region overlaps with the given region such that this
		region is to the right of the given region.

		For example,
		  --------
		--|---   |
		| |  |   |
		| --------
		|    |
		------
		returns true.

		@param[in] region Matrix_Region
		@retval result Whether region overlaps with this region
    */
   bool overlap( const Matrix_Region& ) const;

	/**
		Convert region to a string representation.
		@retval region_string String representation of region
	 */
	std::string to_string( ) const;

	/**
		Perform rescaling on region, which is useful for algorithms that scan a
		matrix at multiple positions and scales.

		@param[in] scale_factor Amount matrix/region is scaled each time
		@param[in] scale_count Number of times scaling has been performed
	 */
	Matrix_Region& rescale( double, size_type );

	/**
		Perform rescaling on region. This is useful for detection algorithms that
		scan a matrix at multiple positions and scales.

		@param[in] scale_amount Amount the matrix/region has been scaled in total
	 */
	Matrix_Region& rescale( double );
};

/*
	Define inlined member functions.
 */

/**
	Get starting row position of the region.
	@retval row_begin Starting row position of the region
 */
inline Matrix_Region::size_type
Matrix_Region::row_begin( ) const
{
	return( _row_begin );
}

/**
	Get starting column position of the region.
	@retval col_begin Starting column position of the region
 */
inline Matrix_Region::size_type
Matrix_Region::col_begin( ) const
{
	return( _col_begin );
}

/**
	Get number of rows in region.
	@retval row_size Number of rows in region
 */
inline Matrix_Region::size_type
Matrix_Region::row_size( ) const
{
	return( _row_size );
}

/**
	Get number of columns in region.
	@retval col_size Number of columns in region
 */
inline Matrix_Region::size_type
Matrix_Region::col_size( ) const
{
	return( _col_size );
}

/**
	Get ending row position of the region.
	@retval row_end Ending row position of the region
 */
inline Matrix_Region::size_type
Matrix_Region::row_end( ) const
{
	return( _row_end );
}

/**
	Get ending column position of the region.
	@retval col_end Ending column position of the region
 */
inline Matrix_Region::size_type
Matrix_Region::col_end( ) const
{
	return( _col_end );
}

/**
	Get left-most position of the region.
	@retval left Left-most position of the region
 */
inline Matrix_Region::size_type
Matrix_Region::left( ) const
{
	return( col_begin() );
}

/**
	Get right-most valid position of the region.
	@retval right Right-most position of the region
 */
inline Matrix_Region::size_type
Matrix_Region::right( ) const
{
	return( col_end() - 1 );
}

/**
	Get top-most position of the region.
	@retval top Top-most position of the region
 */
inline Matrix_Region::size_type
Matrix_Region::top( ) const
{
	return( row_begin() );
}

/**
	Get bottom-most valid position of the region.
	@retval bottom Bottom-most position of the region
 */
inline Matrix_Region::size_type
Matrix_Region::bottom( ) const
{
	return( row_end() - 1 );
}

/**
	Set starting row position of the region.
	@param[in] row_begin Starting row position of the region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::row_begin( size_type new_row_begin )
{
	_row_begin = new_row_begin;
	_row_end = row_begin() + row_size();
	return( *this );
}

/**
	Set starting row position of the region.
	@param[in] row_begin Starting row position of the region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::set_row_begin( size_type new_row_begin )
{
	return( row_begin( new_row_begin ) );
}

/**
	Set starting column position of the region.
	@param[in] col_begin Starting column position of the region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::col_begin( size_type new_col_begin )
{
	_col_begin = new_col_begin;
	_col_end = col_begin() + col_size();
	return( *this );
}

/**
	Set starting column position of the region.
	@param[in] col_begin Starting column position of the region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::set_col_begin( size_type new_col_begin )
{
	return( col_begin( new_col_begin ) );
}

/**
	Set number of rows in region.
	@param[in] row_size Number of rows in region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::row_size( size_type new_row_size )
{
	_row_size = new_row_size;
	_row_end = row_begin() + row_size();
	return( *this );
}

/**
	Set number of rows in region.
	@param[in] row_size Number of rows in region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::set_row_size( size_type new_row_size )
{
	return( row_size( new_row_size ) );
}

/**
	Set number of columns in region.
	@param[in] col_size Number of columns in region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::col_size( size_type new_col_size )
{
	_col_size = new_col_size;
	_col_end = col_begin() + col_size();
	return( *this );
}

/**
	Set number of columns in region.
	@param[in] col_size Number of columns in region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::set_col_size( size_type new_col_size )
{
	return( col_size( new_col_size ) );
}

/**
	Set the region's starting row and column positions.

	@param[in] row_begin Starting row position of the region
	@param[in] col_begin Starting col position of the region
	@param[in] row_size Number of rows in the region
	@param[in] col_size Number of columns in the region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::begin( size_type new_row_begin, size_type new_col_begin )
{
	row_begin(new_row_begin);
	col_begin(new_col_begin);
	return( *this );
}

/**
	Set the region's starting row and column positions.

	@param[in] row_begin Starting row position of the region
	@param[in] col_begin Starting col position of the region
	@param[in] row_size Number of rows in the region
	@param[in] col_size Number of columns in the region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::set_begin( size_type new_row_begin, size_type new_col_begin )
{
	return( begin( new_row_begin, new_col_begin ) );
}

/**
	Set the region's row and column sizes.

	@param[in] row_size Number of rows in the region
	@param[in] col_size Number of columns in the region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::size( size_type new_row_size, size_type new_col_size )
{
	row_size(new_row_size);
	col_size(new_col_size);
	return( *this );
}

/**
	Set the region's row and column sizes.

	@param[in] row_size Number of rows in the region
	@param[in] col_size Number of columns in the region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::set_size( size_type new_row_size, size_type new_col_size )
{
	return( size( new_row_size, new_col_size ) );
}

/**
	Set the region's starting row and column positions and sizes.

	@param[in] row_begin Starting row position of the region
	@param[in] col_begin Starting col position of the region
	@param[in] row_size Number of rows in the region
	@param[in] col_size Number of columns in the region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::operator()( size_type new_row_begin, size_type new_col_begin,
	size_type new_row_size, size_type new_col_size )
{
	begin( new_row_begin, new_col_begin );
	size( new_row_size, new_col_size );
	return( *this );
}

/**
	Set the region's starting row and column positions and sizes.

	@param[in] row_begin Starting row position of the region
	@param[in] col_begin Starting col position of the region
	@param[in] row_size Number of rows in the region
	@param[in] col_size Number of columns in the region
	@retval this This region
 */
inline Matrix_Region&
Matrix_Region::set( size_type new_row_begin, size_type new_col_begin,
	size_type new_row_size, size_type new_col_size )
{
	return( operator()( new_row_begin, new_col_begin,
				new_row_size, new_col_size )
			);
}

/**
	Compare two regions for equality.

	@retval equal Whether the given region is the same as this region
 */
inline bool
Matrix_Region::operator==( const Matrix_Region& region ) const
{
	return(
			row_begin() == region.row_begin()
			&& col_begin() == region.col_begin()
			&& row_end() == region.row_end()
			&& col_end() == region.col_end()
	);
}

/**
	Compare two regions for inequality.

	@retval equal Whether the given region is not the same as this region
 */
inline bool
Matrix_Region::operator!=( const Matrix_Region& region ) const
{
	return( !(*this == region) );
}

/**
	Get area size of region.
	@retval area Area of region
 */
inline Matrix_Region::size_type
Matrix_Region::area( )
{
	return( row_size() * col_size() );
}

/**
	Get perimeter size around region.
	@retval perimeter perimeter of region
 */
inline Matrix_Region::size_type
Matrix_Region::perimeter( )
{
	return( (2 * row_size()) + (2 * col_size()) );
}

/**
	Determine if this region lies inside another region.

	For example,
	----------
	|  ---   |
	|  | |   |
	|  ---   |
	|        |
	----------
	returns true.

	Note: a region is inside of itself, so 'this->inside( *this )' returns true.

	@param[in] region Matrix_Region
	@retval result Whether this region lies inside of the given region
 */
inline bool
Matrix_Region::inside( const Matrix_Region& region ) const
{
	return(
			(row_begin() >= region.row_begin() && row_begin() < region.row_end())
			&& (col_begin() >= region.col_begin() && col_begin() < region.col_end())
			&& (row_end() >= region.row_begin() && row_end() <= region.row_end())
			&& (col_end() >= region.col_begin() && col_end() <= region.col_end())
	);
}

/**
	Determine if this region overlaps with the given region such that this
	region is to the right of the given region.

	For example,
	  --------
	--|---   |
	| |  |   |
	| --------
	|    |
	------
	returns true.

	@param[in] region Matrix_Region
	@retval result Whether region overlaps with this region
 */
inline bool
Matrix_Region::overlap( const Matrix_Region& region ) const
{
	// if regions span some similar columns, then they may overlap, so check
	// if the rows overlap next
	if( (col_begin() <= region.col_begin() && col_end() > region.col_begin() )
	|| (region.col_begin() < col_begin() && region.col_end() > col_begin()) )
	{
		// if this region's rows start above and extend below right
		// region's rows, then they overlap
		if( row_begin() <= region.row_begin()
			&& row_end() > region.row_begin() )
		{
			return( true );
		}

		// or if right region's rows start above and extend below this
		// region's rows, then they overlap
		else if( region.row_begin() < row_begin()
			&& region.row_end() > row_begin() )
		{
			return( true );
		}
	}

	return( false );
}

/**
	Convert region to a string representation.
	@retval region_string String representation of region
 */
inline std::string
Matrix_Region::to_string( ) const
{
	/*
	 TODO (we must determine the number of digits used to do this)
	// verify that we do not overflow the character buffer
	assert(
		ws_tools::count_digits( row_begin() )
			+ ws_tools::count_digits( col_begin() )
			+ ws_tools::count_digits( row_begin() + row_size() - 1 )
			+ ws_tools::count_digits( col_begin() + col_size() - 1 )
			+ std::string( "[,]->[,]" ).size() )
		< BUFSIZ
	);
	 */

	char buf[ BUFSIZ ];
	sprintf( buf, "[%u,%u]->[%u,%u]",
		row_begin(), col_begin(),
		row_begin() + row_size() - 1,
		col_begin() + col_size() - 1 );
	
	return( std::string(buf) );
}

/**
	Perform rescaling on region, which is useful for algorithms that scan a
	matrix at multiple positions and scales.

	@param[in] scale_factor Amount matrix/region is scaled each time
	@param[in] scale_count Number of times scaling has been performed
 */
inline Matrix_Region&
Matrix_Region::rescale( double scale_factor, size_type scale_count )
{
	double scale_amount = pow( scale_factor, scale_count );
	rescale( scale_amount );
	return( *this );
}

/**
	Perform rescaling on region. This is useful for detection algorithms that
	scan a matrix at multiple positions and scales.

	@param[in] scale_amount Amount the matrix/region has been scaled in total
 */
inline Matrix_Region&
Matrix_Region::rescale( double scale_amount )
{
	set_row_begin( size_type( row_begin() / scale_amount ) );
	set_col_begin( size_type( col_begin() / scale_amount ) );
	set_row_size(  size_type( row_size()  / scale_amount ) );
	set_col_size(  size_type( col_size()  / scale_amount ) );
	return( *this );
}

} // namespace ws_img

#endif // MATRIX_REGION_HPP
