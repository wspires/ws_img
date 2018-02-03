/**
	@file   Edge_Detector.hpp
	@author Wade Spires
	@date   2005/12/1
	@brief  Class Edge_Detector.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef EDGE_DETECTOR_HPP
#define EDGE_DETECTOR_HPP

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
#include "Gauss.hpp"
#include "Matrix.hpp"
#include "ws_tools.hpp"

namespace ws_img
{

/**
	@brief Edge_Detector
 */
class Edge_Detector
{
	
public:

	typedef Image::pixel_type  pixel_type;
	typedef Image::size_type   size_type;

	/// Types of edge detectors
	enum Detector_Type { Prewitt, Sobel, Canny };

	static const pixel_type Edge_Pixel;      //< Value used for edge pixels
	static const pixel_type Non_Edge_Pixel;  //< Value used for non-edge pixels

	static const pixel_type DEFAULT_THRESHOLD;

	/**
		Default constructor.
	 */
   Edge_Detector( )
	: _detector_type( Sobel )
	{ }

	/**
		Construct edge detector using given type.
	 */
   Edge_Detector( const Detector_Type detector_type )
	: _detector_type( detector_type )
	{ }

	/**
		Destructor does nothing since no member variables are dynamically
		allocated.
	 */
   ~Edge_Detector( )
	{ }

	Image detect( std::string&, pixel_type = DEFAULT_THRESHOLD );
	void  detect( Image&, pixel_type = DEFAULT_THRESHOLD );

	static void prewitt( Image&, pixel_type = DEFAULT_THRESHOLD );
	static void sobel( Image&, pixel_type = DEFAULT_THRESHOLD );
	static void canny( Image&, pixel_type = 40, pixel_type = 70 );

	static void canny( Image&, Image&, Image&, Matrix&, 
			pixel_type = 40, pixel_type = 70 );

	static void canny( Image&, Image&, Image&, Matrix&, Matrix&,
			pixel_type = 40, pixel_type = 70 );

	// prewitt and sobel function
	static void compute_edge_map( Image&, const Matrix&, const Matrix&,
			pixel_type );

	// canny functions
	static Matrix non_maxima_suppression( const Matrix&, const Matrix&,
			Matrix& );
	static void hysteresis( Image&, const Matrix&, pixel_type, pixel_type );

protected:

	static void hysteresis( Image&, const Matrix&, size_type, size_type,
			pixel_type, pixel_type );
	static void update_position( double, size_type, size_type,
		size_type*, size_type*, size_type*, size_type* );

private:

	// edges are on and non-edges are off (used locally by functions but not by
	// the user); we have to define two different sets of values since when we
	// convert to bilevel, the values are flipped
	static const pixel_type _Edge_Pixel;
	static const pixel_type _Non_Edge_Pixel;

	/// Type of edge detector
	Detector_Type _detector_type;
};

} // namespace ws_img

#endif // EDGE_DETECTOR_HPP 
