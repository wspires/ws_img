/**
	@file   Keypoint.hpp
	@author Wade Spires
	@date   2006/1/30
	@brief  Class Keypoint.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef _KEYPOINT_HPP
#define _KEYPOINT_HPP

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include <algorithm>
#include <string>
#include <vector>

// c headers
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// my headers
#include "ws_tools.hpp"

// Image headers
#include "Image.hpp"
#include "Image_Coordinate.hpp"
#include "Vector.hpp"

#ifndef PI
#define PI M_PI
#endif

/**
	@brief Keypoint
 */
class Keypoint
{

typedef ws_img::Image::size_type  size_type;
typedef ws_img::Image_Coordinate  Img_Coord;
	
public:

	/**
		Default constructor.
	 */
   Keypoint( )
	{ }

	/**
		Destructor does nothing since no member variables are dynamically
		allocated.
	 */
   ~Keypoint( )
	{ }

	size_type  scale_pos;    //< Scale point was found in
	size_type  octave_pos;   //< Octave point was found in

	Img_Coord       img_coord;    //< Spatial location
	double          scale;        //< Image scale
	double          orientation;  //< Orientation
	ws_img::Vector  descriptors;  //< Set of descriptors

	ws_img::Matrix  grad_mag;  //< Gradient magnitude
	ws_img::Matrix  grad_dir;  //< Gradient direction

	std::string to_string( ) const;
};

/**
	@brief Keypoint list
 */
class Keypoint_List : public std::vector<Keypoint>
{

public:

	/**
		Default constructor.
	 */
   Keypoint_List( )
	{ }

	/**
		Destructor does nothing since no member variables are dynamically
		allocated.
	 */
   ~Keypoint_List( )
	{ }

	void write( const std::string& );
};

#endif // _KEYPOINT_HPP 
