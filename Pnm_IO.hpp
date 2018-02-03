/**
	@file   Pnm_IO.hpp
	@author Wade Spires
	@date   2006/12/16
	@brief  Class Pnm_IO.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef PNM_IO_HPP
#define PNM_IO_HPP

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include <string>

namespace ws_img
{

class Image;

/**
	@brief Pnm_IO
 */
class Pnm_IO
{

private:
	
	/**
		Default constructor.
	 */
   Pnm_IO( );

	/**
		Destructor does nothing since no member variables are dynamically
		allocated.
	 */
   virtual ~Pnm_IO( );

public:

	/**
		Type of image data: gray-scale ascii, gray-scale binary, color ascii, etc.
		This is used internally for reading images.
	 */
	enum Image_Type { PBM_ASCII, PBM_BINARY, PGM_ASCII, PGM_BINARY,
		PPM_ASCII, PPM_BINARY };

	/**
	  Determine if file is a pnm image or not by its name. Can be passed to
	  ws_tools::dir_traverse() as a file name filter.
	  @param[in] file_name Name of file to check
	 */
	static bool is_pnm( const std::string& );

	/**
		Read image file.
		@param[in] image_file_name Name of image file to read image from
	 */
	static void read( Image&, const std::string& );

	/**
		Write contents of raster image to file.
		@param[in] image_file_name Name of image file to write image to
		@param[in] mode File mode to write as: binary or ascii (binary is default)
	 */
	static void write( const Image&, const std::string& = "",
		const std::string = "b" );

protected:  // internal functions for read()

	/**
		Read header information from image file.
		@param[in] image_file_fp image file pointer
	 */
	static void read_header( FILE*, unsigned*, unsigned*, unsigned*,
			Image_Type* );

	/**
      Remove comments from string. Comments begin with the '#' character.
      @param[in,out] str String to remove comments from
    */
   static void remove_comments( std::string& );

	/**
      Remove white space from string. White space consists of spaces, tabs, and
      newlines.
      @param[in,out] str String to remove white space from
    */
   static void remove_white_space( std::string& );

	/**
		Read magic ID at start of PGM to determine if data is ascii or binary.
		@param[in] line Line from 
		@param[in] pos_1 Position in line of non-whitespace character
		@param[in] pos_2 Position after pos_1 containing next whitespace
		(or npos)
	 */
	static Image_Type read_magic_id( const std::string&,
		const std::string::size_type, const std::string::size_type );

	/**
		Read unsigned integer from header of PGM (row, column, or max number).
		@param[in] line Line from 
		@param[in] pos_1 Position in line of non-whitespace character
		@param[in] pos_2 Position after pos_1 containing next whitespace
		(or npos)
	 */
	static unsigned read_integer( const std::string&,
		const std::string::size_type, const std::string::size_type );

	/**
		Determine type of image.
		@param[in] type Type of image (internal format)
	 */
	static void set_image_type( Image&, const Image_Type& );

	/**
		Read in ascii image data from image file.
		@param[in] image_file_fp File pointer to image file currently reading
	 */
	static void read_ascii_data( Image&, FILE* );

	/**
		Read in binary image data from image file.

		All image data in the file should be a sequence of bytes with no
		intervening white-space. The byte sequence is interpreted differently
		depending on the image type:
		- For bilevel images, pixels are packed into bytes where each bit
		represents a single pixel value (0 -> white, 1 ->black). The last byte in
		a row will be padded at the end if the number of columns is not divisible
		by 8.
		- For gray-scale images, each byte represents a single pixel value in the
		range [0, 255].
		- For color images, pixels are stored as a 3-byte triplet (r,g,b)
		representing the colors red, green, and blue with each byte for a single
		color channel a value in the range [0, 255].

		TODO Handle both big- and little-endian if necessary.

		@param[in] image_file_fp File pointer to image file currently reading
	 */
	static void read_binary_data( Image&, FILE* );

	/**
		Determines if image is an binary image.
		@param[in] img_type Image type (internal format)
		@retval is_binary Whether image is binary--true or false
	 */
	static bool is_binary( const Image_Type& );

	/**
		Determines if image is an ascii image.
		@param[in] img_type Image type (internal format)
		@retval is_ascii Whether image is ascii--true or false
	 */
	static bool is_ascii( const Image_Type& );

protected:  // internal functions for write()

	/**
		Write header information to image file.
		@param[in] image_file_fp File pointer to image file currently writing
	 */
	static void write_binary_header( const Image&, FILE* );

	/**
		Write header information to image file.
		@param[in] image_file_fp File pointer to image file currently writing
	 */
	static void write_ascii_header( const Image&, FILE* );

	/**
		Get time stamp to record time that image was written. 
		@retval time_stamp Current time in string format
	 */
	static std::string get_time_stamp( );

	/**
		Write image data to image file in binary form.
		@param[in] image_file_fp File pointer to image file currently writing
	 */
	static void write_binary_data( const Image&, FILE* );

	/**
		Write image data to image file in ascii form.
		Note: Only 70 characters are written per line, which adds to the code's
		complexity (all '*_width' variables are for this purpose); however, this
		conforms to the formal specification for an ascii image.

		@param[in] image_file_fp File pointer to image file currently writing
	 */
	static void write_ascii_data( const Image&, FILE* );
};

} // namespace ws_img

#endif // PNM_IO_HPP
