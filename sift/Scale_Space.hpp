/**
	@file   Scale_Space.hpp
	@author Wade Spires
	@date   2006/1/28
	@brief  Class Scale_Space.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef _SCALE_SPACE_HPP
#define _SCALE_SPACE_HPP

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

// tools headers
#include "ws_tools.hpp"

// Image headers
#include "Gauss.hpp"
#include "Image.hpp"

/**
	@brief Scale_Space_Image
 */
struct Scale_Space_Image
{
	ws_img::Image  L;  //< Gaussian smoothed image
	ws_img::Image  D;  //< Difference of Gaussian image
	double         s;  //< Sigma used at this scale
};

typedef std::vector<Scale_Space_Image> Octave;

/**
	@brief Scale_Space
 */
class Scale_Space : public std::vector<Octave>
{
	
public:

	static const double    INITIAL_SIGMA;         //< Initial std. deviation
	static const double    DEFAULT_SIGMA;         //< Default std. deviation
	static const unsigned  DEFAULT_NUM_SCALES;    //< Default scales per octave
	static const unsigned  DEFAULT_NUM_OCTAVES;   //< Default number of octaves

	/**
		Default constructor.
	 */
   Scale_Space( )
	: _sigma( DEFAULT_SIGMA )
	{ }

	/**
		Construct scale-space using given image and parameters.
	 */
   Scale_Space( const ws_img::Image&, double = DEFAULT_SIGMA,
			unsigned = DEFAULT_NUM_OCTAVES,
			unsigned = DEFAULT_NUM_SCALES );

	/**
		Destructor does nothing since no member variables are dynamically
		allocated.
	 */
   ~Scale_Space( )
	{ }

	void write( const std::string& );

protected:
	unsigned get_num_octaves( const ws_img::Image&, unsigned, unsigned );
	void scale_img( ws_img::Image&, double, double );

private:
	double _sigma;  //< Initial sigma
};

#endif // _SCALE_SPACE_HPP 
