/**
	@file   Wavelet_Subband.hpp
	@author Wade Spires
	@date   2005/11/26
	@brief  Class Wavelet_Subband.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef WAVELET_SUBBAND_HPP
#define WAVELET_SUBBAND_HPP

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include <algorithm>
#include <map>
#include <string>
#include <vector>

// c headers
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// my headers
#include "Matrix.hpp"
#include "ws_tools.hpp"

namespace ws_img
{

/**
	@brief Wavelet_Level_Band Decomposition level and band position for a
	wavelet transform.
 */
struct Wavelet_Level_Band
{
	/// Specific subband type--used for finding particular band in Subbands map
	enum Wavelet_Band {

		// 1D wavelet bands: low-pass and high-pass
		L, H,

		// 2D wavelet bands: low-pass, low-pass; low-pass, high-pass; etc.
		LL, LH, HL, HH,

		// bands for horizontal and vertical differences for LL band
		LL_H, LL_V,

		// bands for overcomplete transformation: even and odd phases
		LL_EE, LH_EE, HL_EE, HH_EE,
		LL_EO, LH_EO, HL_EO, HH_EO,
		LL_OE, LH_OE, HL_OE, HH_OE,
		LL_OO, LH_OO, HL_OO, HH_OO
	};

	/// Decomposition level: level 1 is from first decomposition, ...
	unsigned  _level;

	/// Subband (LL, LH, etc.) at a given level
	enum Wavelet_Band _band;

	/**
		Construct level 1 LL band by default.
	 */
	Wavelet_Level_Band( )
	: _level(1), _band(LL)
	{ }

	/**
		Construct level 1 LL band by default.
	 */
	Wavelet_Level_Band( const unsigned level, const Wavelet_Band& band )
	: _level( level ), _band( band )
	{ }

	/**
		Destructor does nothing since no member variables are dynamically
		allocated.
	 */
	~Wavelet_Level_Band( )
	{ }

	/**
		Compare two Wavelet_Level_Band structures for inequality--first compare
		level, then compare band type.

		@param rhs Wavelet_Level_Band on right-hand side of < operator.
	 */
	inline bool operator<( const Wavelet_Level_Band& rhs ) const
	{
		if( _level < rhs._level )
		{
			return( true );
		}
		else if( _level == rhs._level )
		{
			if( _band < rhs._band )
			{
				return( true );
			}
		}

		return( false );
	}

	/**
		Convert level-band to a string representation.
		@retval Converted string
	 */
	inline std::string to_string( ) const
	{
		std::string str =
			ws_tools::int_to_string(_level) + "_" + to_string( _band );
		return( str );
	}

	inline Wavelet_Band get_band( ) const
	{
		return( _band );
	}

	inline unsigned get_level( ) const
	{
		return( _level );
	}

	/**
		Convert string to its band representation.
		@param str String to convert
		@retval band Band represented by str
	 */
	Wavelet_Band to_band( const std::string& str )
	{
		Wavelet_Band band = LL;

		// TODO Verify this is correct: copied from old Wavelet code
		if     ( str == "LL" )   band = LL;
		else if( str == "LH" )   band = LH;
		else if( str == "HL" )   band = HL;
		else if( str == "HH" )   band = HH;
		else if( str == "LLH" || str == "LL_H" )  band = LL_H;
		else if( str == "LLV" || str == "LL_V" )  band = LL_V;
		else
		{
			err_quit( "Unable to convert string to band\n" );
		}

		return( band );
	}

protected:

	/**
		Convert string to a level-band representation.
		@param str String to convert
		@retval band Wavelet_Band represented by str
	 */
	std::string to_string( const Wavelet_Band& band ) const
	{
		std::string str;

		switch( band )
		{
			case L:     str = "L";   break;
			case H:     str = "H";   break;

			case LL:    str = "LL";  break;
			case LH:    str = "LH";  break;
			case HL:    str = "HL";  break;
			case HH:    str = "HH";  break;

			// we omit the '_' since other code relies on not seeing a '_'
			case LL_V:  str = "LLV"; break;
			case LL_H:  str = "LLH"; break;
			default:    str = "";    break;
		}
		return( str );
	}
};

/**
	@brief Wavelet_Subband from wavelet decomposition.
	Can be used as an image but also contains information about what
	orientation and level the image came from in a wavelet decomposition.
 */
class Wavelet_Subband : public ws_img::Matrix
{
	Wavelet_Level_Band _level_band;

public:

	/**
		Create default subband using default image.
	 */
	Wavelet_Subband( )
	: ws_img::Matrix()
	{ }

	/**
		Create subband using given image.
		@param matrix Matrix to set subband to
	 */
	Wavelet_Subband( const ws_img::Matrix& matrix )
	: ws_img::Matrix(matrix)
	{ }

	/**
		Create subband using given image and band-level.
		@param matrix Matrix to set subband to
		@param b Band-level to set this subband's band-level to
	 */
	Wavelet_Subband( const ws_img::Matrix& matrix, const unsigned level,
		const Wavelet_Level_Band::Wavelet_Band& band )
	: ws_img::Matrix(matrix), _level_band(level, band)
	{ }

	/**
		Create subband using given image and band-level.
		@param matrix Matrix to set subband to
		@param b Band-level to set this subband's band-level to
	 */
	Wavelet_Subband( const ws_img::Matrix& matrix, const Wavelet_Level_Band& b )
	: ws_img::Matrix(matrix), _level_band(b)
	{ }

	/**
		Destructor does nothing since no member variables are dynamically
		allocated.
	 */
	~Wavelet_Subband( )
	{ }

	/**
		Set wavelet band-level that subband came from.
		@param b Band-level to set this subband's band-level to
	 */
	inline void set_level_band( const Wavelet_Level_Band& level_band )
	{
		_level_band = level_band;
	}

	/**
		Set wavelet band-level that subband came from.
		@param L Level to set this subband's band-level to
		@param B Band to set this subband's band-level to
	 */
	inline void set_level_band( const unsigned level,
		const Wavelet_Level_Band::Wavelet_Band& band )
	{
		_level_band = Wavelet_Level_Band( level, band );
	}

	/**
		Get this subband's band-level.
		@retval _level_band Band-level for this subband
	 */
	inline Wavelet_Level_Band get_level_band( ) const
	{
		return( _level_band );
	}

	// void octave_form( Matrix&, unsigned = DEFAULT_LEVELS );
	// void inverse_octave_form( Matrix&, unsigned = DEFAULT_LEVELS );

protected:
	// void octave_form_1D( iterator, size_type );
	// void inverse_octave_form_1D( iterator, size_type );

	// void print_filter_moments_lifts( const Lift_Filter&, const Moments&,
			// const Lifts& );
};

/**
	@brief Collection of subband images indexed by its decomposition level and
	subband (i.e., a Wavelet_Level_Band structure).
 */
typedef std::map< Wavelet_Level_Band, Wavelet_Subband > Wavelet_Subbands;

} // namespace ws_img

#endif // WAVELET_SUBBAND_HPP 
