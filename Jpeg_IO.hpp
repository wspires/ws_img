/**
	@file   Jpeg_IO.hpp
	@author Wade Spires
	@date   2005/12/8
	@brief  Class Jpeg_IO reads and writes JPEG images.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef JPEG_IO_HPP
#define JPEG_IO_HPP

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include <string>

// my headers
#include "Image.hpp"

namespace ws_img
{

class Image;  // forward declaration

/**
	@brief Jpeg_IO
 */
class Jpeg_IO
{
	
public:

	/**
		Default constructor.
	 */
   Jpeg_IO( )
	{ }

	/**
		Destructor does nothing since no member variables are dynamically
		allocated.
	 */
   ~Jpeg_IO( )
	{ }

	static void write( const Image&, const std::string&, int );
	static int read( Image&, const std::string& );
};

} // namespace ws_img

#endif // JPEG_IO_HPP 
