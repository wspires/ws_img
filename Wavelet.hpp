/**
	@file   Wavelet.hpp
	@author Wade Spires
	@date   2005/11/21
	@brief  Class Wavelet.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef WAVELET_HPP
#define WAVELET_HPP

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

// GSL headers
#include <gsl/gsl_block.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_wavelet2d.h>

// my headers
#include "Matrix.hpp"
#include "Vector.hpp"
#include "Wavelet_Coefficients.hpp"
#include "Wavelet_Subband.hpp"
#include "ws_tools.hpp"

#ifndef NULL
#define NULL 0
#endif

namespace ws_img
{

/**
	@brief Wavelet
 */
class Wavelet
{

public:

	typedef Matrix::value_type        value_type;
	typedef Matrix::size_type         size_type;

	enum Decomp_Type { Standard, Pyramid };
	enum Wavelet_Type {
		Haar,
		Daub_4, Daub_6, Daub_8, Daub_10,
		Daub_12, Daub_14, Daub_16, Daub_18, Daub_20,
		Bspline_103, Bspline_105,
		Bspline_202, Bspline_204, Bspline_206, Bspline_208,
		Bspline_301, Bspline_303, Bspline_305, Bspline_307, Bspline_309
	};

	/**
		Default constructor.
	 */
   Wavelet( )
	: _wavelet(NULL), _wavelet_type(Haar), _workspace(NULL)
	{
		_wavelet = get_wavelet( _wavelet_type );
		set_wavelet( _wavelet_type );
	}

	/**
		Construct wavelet using given type.
		@param[in] wavelet_type Type of wavelet transform to perform
	 */
   Wavelet( const Wavelet_Type wavelet_type )
	: _wavelet(NULL), _wavelet_type(wavelet_type), _workspace(NULL)
	{
		_wavelet = get_wavelet( _wavelet_type );
		set_wavelet( _wavelet_type );
	}

	/**
		Construct wavelet using given name.
		@param[in] wavelet_name name of wavelet transform to perform
	 */
   Wavelet( const std::string& wavelet_name )
	: _wavelet(NULL), _workspace(NULL)
	{
		_wavelet_type = get_wavelet_type( wavelet_name );
		_wavelet      = get_wavelet( _wavelet_type );
		set_wavelet( _wavelet_type );
	}

	/**
		Copy constructor.
	 */
   Wavelet( const Wavelet& wavelet )
	: _wavelet(NULL), _wavelet_type(wavelet._wavelet_type), _workspace(NULL)
	{
		if( wavelet._wavelet != NULL )
		{
			_wavelet = get_wavelet( wavelet._wavelet_type );
		}
	}

	/**
		Assignment operator.
		@param[in] wavelet Wavelet to copy
	  */
	Wavelet& operator=( const Wavelet& wavelet )
	{
		if( &wavelet != this )
		{
			if( wavelet._wavelet != NULL )
			{
				_wavelet_type = wavelet._wavelet_type;
				_wavelet      = get_wavelet( _wavelet_type );
			}
		}

		return( *this );
	}


	/**
		Destructor does nothing since no member variables are dynamically
		allocated.
	 */
   ~Wavelet( )
	{
		if( _wavelet != NULL )
		{
			gsl_wavelet_free( _wavelet );
		}
		if( _workspace != NULL )
		{
			gsl_wavelet_workspace_free( _workspace );
		}
	}

	void forward_transform( Matrix&, Decomp_Type = Standard ) const;
	void inverse_transform( Matrix&, Decomp_Type = Standard ) const;

	Wavelet_Subbands get_subbands( Matrix&, unsigned ) const;

	static Wavelet_Type get_wavelet_type( const std::string& );
	static std::string  get_wavelet_type( const Wavelet_Type& );

	void print_filters( ) const;

protected:

	void forward_standard( Matrix& ) const;
	void forward_pyramid( Matrix& ) const;

	void inverse_standard( Matrix& ) const;
	void inverse_pyramid( Matrix& ) const;
	
	gsl_matrix*  matrix_to_gsl( const Matrix& ) const;
	void         gsl_to_matrix( gsl_matrix*, Matrix& ) const;

	gsl_wavelet* get_wavelet( const Wavelet_Type ) const;

	void set_wavelet( const Wavelet_Type wavelet_type, bool = false );

	/**
		Get number of levels the signal can be decomposed.
		@param[in] num_levels Number of requested levels
		@param[in] size Size/length of signal to transform
		@retval num_levels Number of levels that signal can be transformed
	  */
	inline unsigned
	get_num_levels( unsigned num_levels, unsigned size ) const
	{
		// compute how many times the signal can be split in half and take the
		// smaller between this value and the number of requested levels
		int n = (size == 0) ? 0 : int( floor( log2(size) ) );
		return( (n < 0) ? 0 : std::min<int>( n, num_levels ) );
	}

private:

	/// GSL wavelet
	gsl_wavelet*  _wavelet;

	/// Type of wavelet
	Wavelet_Type  _wavelet_type;

	/**
		TODO
		Workspace for performing tranformation--we can possibly save this space
		to avoid reallocation
	 */
	gsl_wavelet_workspace* _workspace;

	std::string name;
	const double* _h1;
	const double* _g1;
	const double* _h2;
	const double* _g2;
	size_t _nc;
	size_t _offset;

	void wavelet2d_nstransform( Matrix&, unsigned ) const;

	// void wavelet2d_nstransform( const my_wavelet_type*, double*, size_t, size_t,
			// size_t, Vector& );
	int binary_logn( const unsigned ) const;
	// void dwt_step( const my_wavelet_type*, double*, size_t, size_t, Vector& );
	// void inverse_dwt_step( const my_wavelet_type*, double*, size_t, size_t,
			// Vector& );
};

} // namespace ws_img

#endif // WAVELET_HPP 
