/**
	@file   Lift_Scheme.hpp
	@author Wade Spires
	@date   2005/10/16
	@brief  Class that applies the lifting scheme to one- and two-dimensional
	signals.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef LIFT_SCHEME_HPP
#define LIFT_SCHEME_HPP

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include <algorithm>
#include <vector>

// c headers
#include <cstdio>
#include <cstdlib>
#include <cmath>

// my headers
#include "Matrix.hpp"
#include "Vector.hpp"
#include "Linear_Algebra.hpp"
#include "Wavelet_Subband.hpp"
#include "ws_tools.hpp"

namespace ws_img
{

/**
	@brief Lifting scheme.
 */
// TODO make template
// template <unsigned N, unsigned N_tilde>
class Lift_Scheme
{

public:

	typedef Matrix::size_type               size_type;
	typedef Matrix::row_col_iterator        iterator;
	typedef Matrix::const_row_col_iterator  const_iterator;

	/// Default number of transform levels to use
	static const unsigned DEFAULT_LEVELS = 100000;

protected:

	/**
		Get number of coefficients in current iteration.
		@param[in] size Size/length of signal
		@param[in] step Current step size
		@retval num_coef Number of coefficients
	 */
	inline unsigned
	get_num_coef( unsigned size, unsigned step )
	{
		assert( step != 0 );
		return( (size + step - 1) / step );
	}

	/**
		Get number of levels the signal can be decomposed.
		@param[in] num_levels Number of requested levels
		@param[in] size Size/length of signal to transform
		@retval num_levels Number of levels that signal can be transformed 
	 */
	static inline unsigned
	get_num_levels( unsigned num_levels, unsigned size )
	{
		// compute how many times the signal can be split in half and take the
		// smaller between this value and the number of requested levels
		int n = (size == 0) ? 0 : int( floor( log2(size) ) );
		return( (n < 0) ? 0 : std::min<int>( n, num_levels ) );
	}

	/**
		Get number of levels the signal can be decomposed.
		calculate number of iterations num_levels, which is the maximum value
		s.t.
			2^n * (N - 1) <= L < 2^(n + 1) * (N - 1) + 1
		where L = max( width, height )
		and N = max( # dual vanish. moments, # real vanish. moments )
		=> num_levels = floor( log2( (L - 1) / (N - 1) ) )

		@param[in] num_levels Number of requested levels
		@param[in] size Size/length of signal to transform
		@retval num_levels Number of levels that signal can be transformed 
	 */
	static inline unsigned
	get_num_levels( unsigned num_levels, unsigned size, unsigned N,
		unsigned N_tilde )
	{
		assert( int(N) > 0 && int(N_tilde) > 0 );

		unsigned max_N = std::max<unsigned>(N, N_tilde) - 1;
		return( get_num_levels( num_levels, (size - 1) / max_N ) );
	}

	/**
		Get number of levels the signal can be decomposed for Haar.
		The only difference is that ceil() is used instead of floor().
		@param[in] num_levels Number of requested levels
		@param[in] size Size/length of signal to transform
		@retval num_levels Number of levels that signal can be transformed 
	 */
	inline unsigned
	get_num_levels_Haar( unsigned num_levels, unsigned size )
	{
		// compute how many times the signal can be split in half and take the
		// smaller between this value and the number of requested levels
		int n = (size == 0) ? 0 : int( ceil( log2(size) ) );
		return( (n < 0) ? 0 : std::min<int>( n, num_levels ) );
	}

	/**
		@brief Lifting filter.
	 */
	class Lift_Filter : public Matrix
	{
		unsigned n_;

	public:

		/**
			Filter to handle boundary cases.
			@param[in] num_moments Number of moments
		 */
		Lift_Filter( unsigned num_moments = 2 )
		: Matrix( (num_moments / 2) + 1, num_moments ), n_(num_moments)
		{
			Vector x_vec( n_ );  // position for samples to interpolate

			// verify that N is valid
			if( ws_tools::is_odd(n_) || n_ == 0 )
			{
				Vector empty_vec; 
				//return( empty_vec );
			}

			// sample positions are symmetric around 0, except start at .5 offsets
			x_vec(0) = .5 * (1 - int(n_));  // = (-N / 2) + .5
			for( unsigned j = 1; j != n_; ++j )
			{
				x_vec(j) = x_vec(j - 1) + 1.0;
			}

			// fill in filter matrix in opposite order
			unsigned row_pos = real_row_size() - 1;
			unsigned col_pos = n_ - 1;

			for( unsigned i = 0; i != real_row_size(); ++i )
			{
				col_pos = n_ - 1;  // reset column position

				for( unsigned j = 0; j != n_; ++j )
				{
					// values to interpolate: set location of current coefficient
					// to 1 (all others are 0 by default)
					Vector y_vec( n_ );
					y_vec(j) = 1;

					// fit interpolating polynomial using Neville's algorithm
					// and evaluate it at row i
					at(row_pos, col_pos) = neville( x_vec, y_vec, double(i) );

					--col_pos;
				}
				--row_pos;
			}
		}

		/**
			Destructor does nothing since no member variables are dynamically
			allocated.
		 */
		~Lift_Filter( )
		{ }

		void write( ) const
		{
			fprintf( stderr, "Filter coefficients with N = %u\n", col_size() );
			for( unsigned i = 0; i != real_row_size(); ++i )
			{
				fprintf( stderr, "%d to the left and %d to the right\n",
					i, col_size() - i );

				for( unsigned j = 0; j != col_size(); ++j )
				{
					fprintf( stderr, "%lf ", at(i, j) );
				}
				fprintf( stderr, "\n" );
			}
		}

	protected:

		/**
			Given n points of the form (x,y), use the Neville algorithm to find
			the value at val of the polynomial of degree (n - 1) that interpolates
			the points (x,y).

			@param[in] x Sample position
			@param[in] y Sample value
			@param[in] val Value to evaluate polynomial at
		 */
		double
		neville( const Vector& x, const Vector& y, double val )
		{
			Vector vec( n_ );

			for( unsigned i = 0; i != n_; ++i )
			{
				vec(i) = y(i);

				for( int j = (i - 1); j >= 0; --j )
				{
					double x_diff = x(i) - x(j);

					if( round(x_diff) == 0 )
					{
						err_quit( "%lf (~%d) is too close to zero\n",
								x_diff, round(x_diff) );
					}

					vec(j) = vec(j + 1) + (vec(j + 1) - vec(j))
						* (val - x(i)) / x_diff;
				}
			}

			return( vec(0) );
		}
	};

	/**
		@brief Moments in x and y directions.
	 */
	struct Moments
	{
		// moments in x and y directions
		Matrix X;  //< Moments in x direction
		Matrix Y;  //< Moments in y direction

		typedef Matrix::size_type size_type;

		/**
			Default constructor.
		 */
		Moments( )
		{ }

		/**
			Initialize values of the integral moment matrices. The moments are
			equal to k^i where k is the coefficient index and i is the moment (0,
			1, 2, ...).  In case the dimensions of X and Y are different, one
			matrix for each direction is allocated.

			@param[in] width Signal's width
			@param[in] height Signal's height
			@param[in] N_tilde Number of dual vanishing moments
		 */
		Moments( size_type width, size_type height, size_type N_tilde )
		{
			set( width, height, N_tilde );
		}

		/**
			Initialize values of the integral moment matrices. The moments are
			equal to k^i where k is the coefficient index and i is the moment (0,
			1, 2, ...).  In case the dimensions of X and Y are different, one
			matrix for each direction is allocated.

			@param[in] width Signal's width
			@param[in] height Signal's height
			@param[in] N_tilde Number of dual vanishing moments
		 */
		void
		set( size_type width, size_type height, size_type N_tilde )
		{
			X.resize(width, N_tilde);
			Y.resize(height, N_tilde);

			if( ws_tools::is_odd(N_tilde) )
			{
				// err_quit( "Constructing moments: N_tilde (%u) cannot be odd\n",
					// N_tilde );
			}

			// initialize moments in X direction
			for( size_type i = 0; i != X.real_row_size(); ++i )
			{
				for( size_type j = 0; j != X.real_col_size(); ++j )
				{
					if( i == 0 && j == 0 )
					{
						// pow(0,0) is undefined, so set to 1
						X(i, j) = 1.0;
						// X(i, j) = 0.0;
					}
					else // moment = i^j
					{
						X(i, j) = pow( i, j );
						// X(i, j) = pow( j, i );  // maybe try this
					}
				}
			}

			// initialize moments in Y direction
			for( size_type i = 0; i != Y.real_row_size(); ++i )
			{
				for( size_type j = 0; j != Y.real_col_size(); ++j )
				{
					if( i == 0 && j == 0 )
					{
						// pow(0,0) is undefined, so set to 1
						Y(i, j) = 1.0;
						// Y(i, j) = 0.0;
					}
					else // moment = i^j
					{
						Y(i, j) = pow( i, j );
						// Y(i, j) = pow( j, i );
					}
				}
			}
		}

		/**
			Destructor does nothing since no member variables are dynamically
			allocated.
		 */
		~Moments( )
		{ }
	};

	/**
		@brief Lifting coefficients.
	 */
	struct Lifts
	{
		Matrix X;  //< Lifting coefficients in x direction
		Matrix Y;  //< Lifting coefficients in y direction

		/**
			Default constructor.
		 */
		Lifts( )
		{ }

		/**
			Allocate and initialize lifting matrices.
			@param[in] width Signal's width
			@param[in] height Signal's height
			@param[in] N_tilde Number of dual vanishing moments
		 */
		Lifts( unsigned width, unsigned height, unsigned N, unsigned N_tilde,
			unsigned num_levels )
		{
			set( width, height, N, N_tilde, num_levels );
		}

		/**
			Allocate and initialize lifting matrices.
			@param[in] width Signal's width
			@param[in] height Signal's height
			@param[in] N_tilde Number of dual vanishing moments
		 */
		void
		set( unsigned width, unsigned height, unsigned N, unsigned N_tilde,
			unsigned num_levels )
		{
			unsigned row_levels = get_num_levels( num_levels, width, N, N_tilde );
			unsigned d_num_X = width / 2;
			if( d_num_X > 0 && row_levels > 0 )
			{
				X.resize( row_levels, d_num_X * N_tilde );
			}

			// number of differences in y direction
			unsigned col_levels = get_num_levels( num_levels, height, N, N_tilde);
			unsigned d_num_Y = height / 2;
			if( d_num_Y > 0 && col_levels > 0 )
			{
				Y.resize( col_levels, d_num_Y * N_tilde );
			}
		}

		/**
			Delete allocated lifting matrices.
		 */
		virtual ~Lifts( )
		{ }
	};

public:

	/**
		Default constructor.	
	 */
	Lift_Scheme( )
	: n_(2), n_tilde_(2), _filter(2),
		_current_row_size(0), _current_col_size(0)
	{ }

	/**
		Construct lifting parameters using given values for each parameter.
		@param[in] n Number of real vanishing moments
		@param[in] n_tilde Number of dual vanishing moments
	 */
	Lift_Scheme( unsigned n, unsigned N_tilde )
	: n_(n), n_tilde_(N_tilde), _filter(n),
		_current_row_size(0), _current_col_size(0)
	{ }

	/**
		Destructor does nothing since no member variables are dynamically
		allocated.
	 */
	~Lift_Scheme( )
	{ }

	// lifting functions
	void lift( Vector&, unsigned = DEFAULT_LEVELS );
	void lift( Matrix&, unsigned = DEFAULT_LEVELS );

	void lift_Haar( Vector&, unsigned = DEFAULT_LEVELS );
	void lift_Haar( Matrix&, unsigned = DEFAULT_LEVELS );

	// inverse lifting functions
	void inverse_lift( Vector&, unsigned = DEFAULT_LEVELS );
	void inverse_lift( Matrix&, unsigned = DEFAULT_LEVELS );

	void inverse_lift_Haar( Vector&, unsigned = DEFAULT_LEVELS );
	void inverse_lift_Haar( Matrix&, unsigned = DEFAULT_LEVELS );

	void lift_int( Vector&, unsigned = DEFAULT_LEVELS );
	void inverse_lift_int( Vector&, unsigned = DEFAULT_LEVELS );

	// functions for conversion into the standard octave format
	void octave_form( Matrix&, unsigned = DEFAULT_LEVELS ) const;
	void inverse_octave_form( Matrix&, unsigned = DEFAULT_LEVELS ) const;

	Wavelet_Subbands get_subbands( Matrix&, unsigned = DEFAULT_LEVELS ) const;

	/**
		Get N.
		@retval N
	 */
	inline unsigned get_N( ) const
	{
		return( n_ );
	}

	/**
		Get N tilde.
		@retval N tilde
	 */
	inline unsigned get_N_tilde( ) const
	{
		return( n_tilde_ );
	}

protected:

	void lift( Vector&, Lifts&, Moments&, unsigned );
	void lift( Matrix&, Lifts&, Moments&, unsigned );
	void lift_Haar_1D( iterator, unsigned, unsigned );

	void inverse_lift( Vector&, Lifts&, Moments&, unsigned );
	void inverse_lift( Matrix&, Lifts&, Moments&, unsigned );
	void inverse_lift_Haar_1D( iterator, unsigned, unsigned );

	// auxiliary functions
	void calc_moments_lifts( unsigned, unsigned, unsigned, Moments&, Lifts& );

	void calc_moment( Matrix&, unsigned, unsigned );

	// void calc_lift( Vector&, Matrix&, unsigned, unsigned );
	void calc_lift( iterator, Matrix&, unsigned, unsigned );

	void Predict( iterator, unsigned, size_type );
	void Predict_row( Matrix&, const size_type, const unsigned,
			const std::vector<unsigned>& );
	void Predict_col( Matrix&, const size_type, const unsigned,
			const std::vector<unsigned>& );
	void Update( iterator, unsigned, const_iterator, size_type );
	void Update_row( Matrix&, const size_type, const unsigned,
			const Matrix&, const size_type, const std::vector<unsigned>& );
	void Update_col( Matrix&, const size_type, const unsigned,
			const Matrix&, const size_type, const std::vector<unsigned>& );

	std::vector<unsigned> get_coef_cases( unsigned, unsigned, unsigned );
	void invert_filter_lift( Lifts& );

	void lift_int( Vector&, Lifts&, Moments&, unsigned );
	void inverse_lift_int( Vector&, Lifts&, Moments&, unsigned );
	void Predict_int( Vector&, unsigned );
	void Update_int( Vector&, unsigned, const_iterator,
			size_type );

	void octave_form_1D( iterator, size_type ) const;
	void inverse_octave_form_1D( iterator, size_type ) const;

	void print_filter_moments_lifts( const Lift_Filter&, const Moments&,
			const Lifts& );

private:

	unsigned     n_;        //< Number of real vanishing moments
	unsigned     n_tilde_;  //< Number of dual vanishing moments

	Lift_Filter  _filter;   //< Filter coefficients

	Lifts        _lifts;    //< Lifting coefficients
	Moments      _moments;  //< Moments

	size_type    _current_row_size;
	size_type    _current_col_size;
};

} // namespace ws_img

#endif // LIFT_SCHEME_HPP 
