/**
	@file   Matrix_Iterator.hpp
	@author Wade Spires
	@date   2005/10/28
	@brief  Definition of matrix iterator types.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef MATRIX_ITERATOR_HPP
#define MATRIX_ITERATOR_HPP

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include  <algorithm>
#include  <functional>
#include  <stdexcept>
#include  <iomanip>
#include  <iostream>
#include  <iterator>
#include  <string>
#include  <vector>

// c headers
#include  <cfloat>
#include  <climits>
#include  <cmath>

#include "Matrix_Region.hpp"
#include "ws_tools.hpp"

namespace ws_img
{

/**
	@brief Iterator to entire matrix.
 */
template <class _Traits>
class _matrix_iterator
{
	// use types defined in traits
	typedef typename _Traits::value_type        value_type;
	typedef typename _Traits::size_type         size_type;
	typedef typename _Traits::reference         reference;
	// typedef typename _Traits::iterator_category iterator_category;
	typedef typename _Traits::iterator          iterator;
	typedef typename _Traits::pointer           pointer;
	typedef typename _Traits::difference_type   difference_type;

	iterator  _iter;       //< Iterator column in matrix
	size_type _pos;        //< Column position in matrix
	size_type _col_size;   //< Number of columns in matrix region
	size_type _row_count;  //< Number of rows traversed

	/**
		Amount to offset iterator by when at end of column in the region:
			_offset = (total # columns in matrix - last column of region)
				+ number of columns before region
			= number of columns after region + number of columns before region
	 */
	difference_type _offset;

public:

	/**
		Default constructor constructs an invalid iterator.
	 */
	_matrix_iterator( )
	: _iter(0), _pos(0), _col_size(0), _row_count(0), _offset(0)
	{ }

	/**
		Construct iterator at given column.
		@param[in] iter Iterator to position in matrix region
		@param[in] col_begin First column of region
		@param[in] col_end One passed the end of the column of region
		@param[in] full_col_size Number of columns in entire matrix
			(not the region)
		@param[in] offset Number of columns in matrix
	 */
	_matrix_iterator( const iterator iter,
		const size_type col_begin = 0,
		const size_type col_end = 0,
		const size_type full_col_size = 0 )
	: _iter(iter), _pos(0), _col_size(col_end - col_begin),
		_row_count(1),
		_offset((full_col_size - col_end) + col_begin)
	{ }

	/**
		Increment iterator position by given offset.
		TODO Remove for loop by figuring out how to compute the next position
		arithmetically. See implementation notes for operator++ as to why this is
		non-trivial (i.e., not simply '_iter += offset * _offset;').

		@param[in] offset Amount to offset iterator by
		@retval iterator Iterator to column
	 */
	inline _matrix_iterator& operator+=( const difference_type offset )
	{
		for( difference_type i = 0; i != offset; ++i )
		{
			++(*this);
		}
		return( *this );
	}

	/**
		Increment iterator position by given offset.
		@param[in] offset Amount to offset iterator by
		@retval col_iterator Iterator to column
	 */
	inline _matrix_iterator operator+( const difference_type offset ) const
	{
		return( (_matrix_iterator( *this ) += offset) );
	}

	/**
		Prefix operator to increment iterator position to next position in the
		matrix region.

		Implementation note:
		The implementation is somewhat complicated as we have two cases to
		handle for determining the next iterator position in the region:
		(i) The next iterator position is in the current row of the region.
			Just increment the position by one to get to the next column.
		(ii) The next iterator position is at the start of the next row in
			the region. Increment the iterator by _offset where _offset is the
			number of columns before and after the region (this is set in the
			constructor).
		Although this implies that we should have a conditional statement to
		check for each possibility, we can rearrange the data to handle these
		cases using arithmetic only, which should be more efficient. The version
		using a conditional statement is also given below (commented out) for
		reference and in case this method is either incorrect or less efficient.

		Using pointer arithmetic may be less efficient than using conditional
		statements on architectures that support branch predication (this is
		different from branch _prediction_), such as Intel's IA-64 platform.

		@param[in] offset Amount to offset iterator by
		@retval col_iterator Iterator to column
	 */
	inline _matrix_iterator operator++( )
	{
		// determine if at a new row in region and increment row count if so
		++_pos;
		bool is_new_row = !bool( (_row_count * _col_size) - _pos);
		_row_count += difference_type( is_new_row );  // cast to integral type

		// amount to increment iterator by: either the offset amount if at a
		// new row or just 1 (next column position) if still in current row of
		// the region 
		difference_type offset = (difference_type( is_new_row ) * _offset) + 1;

		_iter += offset;

		// fprintf( stderr, "%u %u => offset %u\n",
			// ((_row_count + 1) * _col_size),
			// ((_row_count + 1) * _col_size) - _pos,
			// offset );

		// short conditional version: increment iterator (must reset _pos?!)
		// _iter += (++_pos == _col_size) ? 1 : (_offset + 1);

		// long conditional version: increment iterator
		// ++_pos;
		// if( _pos == _col_size )
		// {
			// _pos = 0;
			// _iter += _offset;
		// }
		// ++_iter;

		return( *this );
	}

	/**
		Postfix operator to increment iterator position by 1 to go to next row
		of the same column.
		@param[in] unused Unused variable to satisfy operator++ syntax requirement
		@retval col_iterator Iterator to column
	 */
	inline _matrix_iterator operator++( int unused )
	{
		_matrix_iterator iter_copy = *this;
		++(*this);
		return( iter_copy );
	}

	/**
		Increment iterator position by given offset.
		@param[in] offset Amount to offset iterator by
		@retval iterator Iterator to column
	 */
	inline _matrix_iterator& operator-=( const difference_type offset )
	{
		return( (*this += -offset) );
	}

	/**
		Increment iterator position by given offset.
		@param[in] offset Amount to offset iterator by
		@retval col_iterator Iterator to column
	 */
	inline _matrix_iterator operator-( const difference_type offset ) const
	{
		return( (_matrix_iterator( *this ) -= offset) );
	}

	/**
		Prefix operator to increment iterator position by 1 to go to next row of
		the same column.
		@param[in] offset Amount to offset iterator by
		@retval col_iterator Iterator to column
	 */
	inline _matrix_iterator operator--( )
	{
		return( (*this -= 1) );
	}

	/**
		Postfix operator to increment iterator position by 1 to go to next row
		of the same column.
		@param[in] offset Amount to offset iterator by
		@retval col_iterator Iterator to column
	 */
	inline _matrix_iterator operator--( int unused )
	{
		_matrix_iterator iter_copy = *this;
		--_iter;
		return( iter_copy );
	}

	/**
		Dereference iterator to access value pointed to by iterator.
		@retval reference Reference to column element
	 */
	inline reference operator*( ) const
	{
		return( *_iter );
	}

	/**
		Access pointer at current iterator position.
		@retval pointer Pointer to column element
	 */
	inline pointer operator->( ) const
	{
		return( _iter );
	}

	/**
		Compare iterators.
		@param[in] rhs Iterator on right-handside of this iterator
		@retval less_than Whether this iterator position is less than the
			right-handside iterator position
	 */
	inline bool operator<( const _matrix_iterator& rhs ) const
	{
		assert( _pos == rhs._pos );	
		return( _iter < rhs._iter );
	}

	/**
		Compare iterators.
		@param[in] rhs Iterator on right-handside of this iterator
		@retval equal Whether this iterator position is equal to the
			right-handside iterator position
	 */
	inline bool operator!=( const _matrix_iterator& rhs ) const
	{
		assert( _pos == rhs._pos );	
		return( _iter != rhs._iter );
	}
};

/**
	@brief Iterator to row or column in matrix.

	In the comments, below we will refer to this class being an iterator for
	columns only for clarity since the same idea applies to rows directly.
 */
template <class _Traits>
class _row_col_iterator
{
	// use types defined in traits
	typedef typename _Traits::size_type         size_type;
	typedef typename _Traits::reference         reference;
	// typedef typename _Traits::iterator_category iterator_category;
	typedef typename _Traits::iterator          iterator;
	typedef typename _Traits::pointer           pointer;
	typedef typename _Traits::difference_type   difference_type;

	iterator  _iter;    //< Iterator column in matrix
	size_type _pos;     //< Column position in matrix

	/// Amount to offset iterator by with each increment
	difference_type _offset;

public:

	/**
		Default construct constructs an invalid iterator.
	 */
	_row_col_iterator( )
	: _iter(0), _pos(0), _offset(0)
	{ }

	/**
		Construct iterator at given column.
		@param[in] iter Iterator to column in matrix
		@param[in] pos Column position of iterator
		@param[in] offset Number of columns in matrix
	 */
	_row_col_iterator( const iterator iter, const difference_type pos = 0,
		const difference_type offset = 0 )
	: _iter(iter), _pos(pos), _offset(offset)
	{ }

	/**
		Increment iterator position by given offset.
		@param[in] offset Amount to offset iterator by
		@retval iterator Iterator to column
	 */
	inline _row_col_iterator& operator+=( const difference_type offset )
	{
		_iter += offset * _offset;
		return( *this );
	}

	/**
		Increment iterator position by given offset.
		@param[in] offset Amount to offset iterator by
		@retval col_iterator Iterator to column
	 */
	inline _row_col_iterator operator+( const difference_type offset ) const
	{
		return( (_row_col_iterator( *this ) += offset) );
	}

	/**
		Prefix operator to increment iterator position by 1 to go to next row of
		the same column.
		@param[in] offset Amount to offset iterator by
		@retval col_iterator Iterator to column
	 */
	inline _row_col_iterator operator++( )
	{
		return( (*this += 1) );
	}

	/**
		Postfix operator to increment iterator position by 1 to go to next row
		of the same column.
		@param[in] offset Amount to offset iterator by
		@retval col_iterator Iterator to column
	 */
	inline _row_col_iterator operator++( int unused )
	{
		_row_col_iterator iter_copy = *this;
		++(*this);
		return( iter_copy );
	}

	/**
		Increment iterator position by given offset.
		@param[in] offset Amount to offset iterator by
		@retval iterator Iterator to column
	 */
	inline _row_col_iterator& operator-=( const difference_type offset )
	{
		return( (*this += -offset) );
	}

	/**
		Increment iterator position by given offset.
		@param[in] offset Amount to offset iterator by
		@retval col_iterator Iterator to column
	 */
	inline _row_col_iterator operator-( const difference_type offset ) const
	{
		return( (_row_col_iterator( *this ) -= offset) );
	}

	/**
		Prefix operator to increment iterator position by 1 to go to next row of
		the same column.
		@param[in] offset Amount to offset iterator by
		@retval col_iterator Iterator to column
	 */
	inline _row_col_iterator operator--( )
	{
		return( (*this -= 1) );
	}

	/**
		Postfix operator to increment iterator position by 1 to go to next row
		of the same column.
		@param[in] offset Amount to offset iterator by
		@retval col_iterator Iterator to column
	 */
	inline _row_col_iterator operator--( int unused )
	{
		_row_col_iterator iter_copy = *this;
		--(*this);
		return( iter_copy );
	}

	/**
		Subtract two iterators to yield their distance apart.
		The result should be s.t. '(iter += (rhs - iter)) == rhs is be true
		@param[in] rhs Right-hand side of this operator
		@retval col_iterator Iterator to column
	 */
	inline difference_type operator-( const _row_col_iterator rhs ) const
	{
		// bad
		if( _offset != rhs._offset )
		{
			return(0);
		}

		assert( _offset != 0 );
		return( (_iter - rhs._iter) / _offset );
	}

	/**
		Dereference iterator to access value pointed to by iterator.
		@retval reference Reference to column element
	 */
	inline reference operator*( ) const
	{
		return( *_iter );
	}

	/**
		Access pointer at current iterator position.
		@retval pointer Pointer to column element
	 */
	inline pointer operator->( ) const
	{
		return( _iter );
	}

	/**
		Compare iterators.
		@param[in] rhs Iterator on right-handside of this iterator
		@retval less_than Whether this iterator position is less than the
			right-handside iterator position
	 */
	inline bool operator<( const _row_col_iterator& rhs ) const
	{
		assert( _pos == rhs._pos );	
		return( _iter < rhs._iter );
	}

	/**
		Compare iterators.
		@param[in] rhs Iterator on right-handside of this iterator
		@retval equal Whether this iterator position is equal to the
			right-handside iterator position
	 */
	inline bool operator!=( const _row_col_iterator& rhs ) const
	{
		assert( _pos == rhs._pos );	
		return( _iter != rhs._iter );
	}
};

} // namespace ws_img

#endif // MATRIX_ITERATOR_HPP
