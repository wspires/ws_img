/**
   @file   Image.hpp
   @author Wade Spires
   @date   2005/08/02
   @brief  Image class for reading, writing, and processing grayscale PGM and
		color PPM images.

	Functions for inputting and outputting an image are defined in this header
	file. Functions for processing images are defined in the source file
	Image.cpp.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef IMAGE_HPP
#define IMAGE_HPP

#define NDEBUG  // comment this line out to enable assert() debugging calls
#include <cassert>

// c++ headers
#include <algorithm>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

// c headers
#include <cctype>
#include <cerrno>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>

// whether to use boost library's smart pointer
// #define HAVE_BOOST
#ifdef HAVE_BOOST
	#include  <boost/shared_ptr.hpp>
#endif // HAVE_BOOST

// local headers
#include "ws_tools.hpp"
#include "Image_Coordinate.hpp"
#include "Image_Region.hpp"
#include "Jpeg_IO.hpp"
#include "Lift_Scheme.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"
#include "Wavelet.hpp"
#include "Wavelet_Subband.hpp"

namespace ws_img
{

#ifdef HAVE_BOOST
class Image;
// typedef boost::shared_ptr<Image> Img_Ptr;
#endif // HAVE_BOOST

/**
	@brief Class representing PBM/PGM/PPM images.
 */
class Image
{
	// forward declarations
	friend class Jpeg_IO;
	friend class Edge_Detector;

public:  // Data type definitions

	typedef Matrix                   Gray_Image;
	typedef Matrix::value_type       pixel_type;
	typedef pixel_type               value_type;
	typedef Matrix::size_type        size_type;
	typedef Matrix::reference        reference;
	typedef Matrix::const_reference  const_reference;
	typedef Matrix::pointer          pointer;
	typedef Matrix::const_pointer    const_pointer;
	typedef Matrix::difference_type  difference_type;

	typedef Matrix::iterator            iterator;
	typedef Matrix::const_iterator      const_iterator;
	typedef Matrix::row_iterator        row_iterator;
	typedef Matrix::const_row_iterator  const_row_iterator;
	typedef Matrix::col_iterator        col_iterator;
	typedef Matrix::const_col_iterator  const_col_iterator;

	typedef unsigned char  uchar;           //< Abbreviation for unsigned char

	static const pixel_type Black_Pixel;    //< Black pixel value
	static const pixel_type White_Pixel;    //< White pixel value

	static const pixel_type INVALID_PIXEL;  //< Invalid pixel value
	static const pixel_type DEFAULT_MAX;    //< Default maximum pixel value

	/// Type of image: black-white, gray, and color
	enum Image_Type { Bilevel, Gray, Color };

	/// Type of color channel: red, green, and blue
	enum Color_Type { Red, Green, Blue };

private:  // Member variables

	/**
		Type of image data: gray-scale ascii, gray-scale binary, color ascii, etc.
		This is used internally for reading images.
	 */
	enum _Image_Type { PBM_ASCII, PBM_BINARY, PGM_ASCII, PGM_BINARY,
		PPM_ASCII, PPM_BINARY };

	/* Declare member variables */

	// matrices that store the actual image data
	Gray_Image   _red_data;    //< Red color data (also used for gray-scale)
	Gray_Image   _green_data;  //< Green color data
	Gray_Image   _blue_data;   //< Blue color data

	/**
		Points to currently active color channel, which determines the value
		returned by at().
	 */
	Gray_Image*  _current_data;

	/// Maximum possible pixel value in image
	pixel_type   _max_value;

	/// Type of image--gray or color
	Image_Type   _img_type;

public:  // Constructors and destructors

	/**
		Default constructor that constructs an image with a single black pixel.
	 */
	Image( )
	: _red_data( 1, 1 ), _current_data( &_red_data ), _max_value(DEFAULT_MAX),
	_img_type( Gray )
	{ }

	/**
		Copy an image using the given Image object.

		Note: we must supply a copy constructor since, otherwise, only a shallow
		copy of _current_data would be created that pointed to the input image's
		data.

		@param[in] img Image object to copy
	 */
	Image( const Image& img )
	: _red_data( img._red_data ), _green_data( img._green_data ),
	_blue_data( img._blue_data ),
	_max_value( img._max_value ), _img_type( img._img_type )
	{
		// determine which color channel to point to (default is red)
		_current_data = &_red_data;
		if( img._current_data == &img._green_data )
		{
			_current_data = &_green_data;
		}
		else if( img._current_data == &img._blue_data )
		{
			_current_data = &_blue_data;
		}
	}

	/**
		Copy an image using the given Image object.

		Note: we must supply operator= since, otherwise, only a shallow
		copy of _current_data would be created that pointed to the input image's
		data.

		@param[in] img Image object to copy
		@retval this This image
	 */
	Image& operator=( const Image& img )
	{
		if( &img == this )
		{
			return( *this );
		}

		// copy data
		_red_data   = img._red_data;
		_green_data = img._green_data;
		_blue_data  = img._blue_data;
		_max_value  = img._max_value;
		_img_type   = img._img_type;

		// determine which color channel to point to (default is red)
		_current_data = &_red_data;
		if( img._current_data == &img._green_data )
		{
			_current_data = &_green_data;
		}
		else if( img._current_data == &img._blue_data )
		{
			_current_data = &_blue_data;
		}

		return( *this );
	}

	/**
		Construct an image by reading the specified file.
		@param[in] image_file_name Name of image file
	 */
	Image( const std::string& image_file_name )
	: _current_data( &_red_data ), _max_value(DEFAULT_MAX), _img_type( Gray )
	{
		read( image_file_name );  // defer work to read function
	}

	/**
		Construct an image using the specified raw image data.
		Image's size assumed to be greater than 0 (otherwise, segfault).
		@param[in] gray_img_data Gray-scale image data
	 */
	Image( const Gray_Image& gray_img_data )
	: _red_data( gray_img_data ), _current_data( &_red_data ),
	_max_value(DEFAULT_MAX), _img_type( Gray )
	{ }

	/**
		Assign the image data to the specified raw image data.

		Note: Supplying the operator= is not only more efficient (prevents
		construction of temporaries) but also correct since, otherwise, only a
		shallow copy of _current_data would be created that pointed to the input
		image's data.

		@param[in] gray_img_data Gray-scale image data
		@retval this This image
	 */
	Image& operator=( const Gray_Image& gray_img_data )
	{
		_red_data     = gray_img_data;
		_max_value    = DEFAULT_MAX;
		_img_type     = Gray;
		_current_data = &_red_data;

		// delete unused channels
		_green_data = Gray_Image();
		_blue_data  = Gray_Image();

		return( *this );
	}

	/**
		Assign the image data to the specified value.

		@param[in] val New value to assign to each pixel
		@retval this This image
	 */
	Image& operator=( value_type val )
	{
		*_current_data = val;
		return( *this );
	}

	/**
		Construct an empty image with the specified number of rows and columns.
		@param[in] num_rows Number of rows in new image
		@param[in] num_cols Number of columns in new image
	Image( const size_type num_rows, const size_type num_cols,
		const Image_Type& type = Gray )
	: _max_value(DEFAULT_MAX), _img_type(type)
	{
		allocate_image( num_rows, num_cols );
	}
	 */

	/**
		Construct an empty image with the specified number of rows and columns.

		@param[in] num_rows Number of rows in new image
		@param[in] num_cols Number of columns in new image
	 */
	explicit Image( const size_type num_rows, const size_type num_cols )
	: _red_data( num_rows, num_cols ),
		_current_data( &_red_data ),
		_max_value(DEFAULT_MAX), _img_type(Gray)
	{ }

	/**
		Construct an empty image with the specified number of rows and columns.

		Constructs matrix of data in initialization list since it is twice as
		fast as calling allocate_image().  We can either (1) go ahead and
		construct green/blue images in the initialization as well or (2)
		separately construct them if the image is indeed color. (1) is slower if
		the image is just gray-scale and faster if the image is color, while (2)
		is faster if the image is gray-scale and slower if the image is
		color. Since we have a constructor that takes in the rows and columns only
		with the type assumed to be gray-level, we assume a color image is being
		created; otherwise, a user could use the other version for greater
		efficiency.  

		@param[in] num_rows Number of rows in new image
		@param[in] num_cols Number of columns in new image
	 */
	explicit Image( const size_type num_rows, const size_type num_cols,
		const Image_Type& type )
	: _red_data( num_rows, num_cols ),
		_green_data( num_rows, num_cols ),
		_blue_data( num_rows, num_cols ),
		_current_data( &_red_data ),
		_max_value(DEFAULT_MAX), _img_type(type)
	{
		// create other color channels
		if( is_color() )
		{
			// _green_data = Gray_Image( num_rows, num_cols );
			// _blue_data  = Gray_Image( num_rows, num_cols );
		}
	}

	/**
		Construct an image with the specified number of rows and columns with
		each value initialized to the given value.

		Constructs matrix of data in initialization list since it is twice as
		fast as calling allocate_image().  We can either (1) go ahead and
		construct green/blue images in the initialization as well or (2)
		separately construct them if the image is indeed color. (1) is slower if
		the image is just gray-scale and faster if the image is color, while (2)
		is faster if the image is gray-scale and slower if the image is
		color.

		@param[in] num_rows Number of rows in new image
		@param[in] num_cols Number of columns in new image
	 */
	explicit Image( const size_type num_rows, const size_type num_cols,
		const Image_Type& type, const pixel_type value )
	: _red_data( num_rows, num_cols, value ),
		_green_data( num_rows, num_cols, value ),
		_blue_data( num_rows, num_cols, value ),
		_current_data( &_red_data ),
		_max_value(DEFAULT_MAX), _img_type(type)
	{
		// create other color channels
		if( is_color() )
		{
			// _green_data = Gray_Image( num_rows, num_cols, value );
			// _blue_data  = Gray_Image( num_rows, num_cols, value );
		}
	}

	/**
		Construct a single color image by combining 3 separate RGB color channels.
		@param[in] red_img Image representing red values
		@param[in] green_img Image representing green values
		@param[in] blue_img Image representing blue values
	 */
	Image( const Image& red_img, const Image& green_img, const Image& blue_img )
	: _max_value(red_img.get_max_value()), _img_type(Color)
	{
		// verify that image dimensions agree
		assert( red_img.sub_row_size() == green_img.sub_row_size() );
		assert( red_img.sub_row_size() == blue_img.sub_row_size() );
		assert( red_img.sub_col_size() == green_img.sub_col_size() );
		assert( red_img.sub_col_size() == blue_img.sub_col_size() );

		// Note: we don't have to check the image type since we can just use the
		// current channel of each image

 		allocate_image( red_img.sub_row_size(), red_img.sub_col_size() );

		// copy each channel to image
		for( size_type i = row_begin(); i != red_img.row_end(); ++i )
		{
			for( size_type j = col_begin(); j != red_img.col_end(); ++j )
			{
				red(i,j) = red_img(i,j);
			}
		}
		for( size_type i = row_begin(); i != green_img.row_end(); ++i )
		{
			for( size_type j = col_begin(); j != green_img.col_end(); ++j )
			{
				green(i,j) = green_img(i,j);
			}
		}
		for( size_type i = row_begin(); i != blue_img.row_end(); ++i )
		{
			for( size_type j = col_begin(); j != blue_img.col_end(); ++j )
			{
				blue(i,j) = blue_img(i,j);
			}
		}
	}

	/**
		Default destructor does nothing since no member variables are
		dynamically allocated.
	 */
	virtual ~Image( )
	{ }

public:  // Public accessor functions

	/**
		Return number of rows in image, ignoring the current image region.
		@retval num_rows Number of rows
	 */
	size_type real_row_size( ) const;

	/**
		Get the number of rows in the image.
		@retval num_rows Number of rows in image
	 */
	size_type row_size( ) const;

	/**
		Get the number of rows in the image.
		@retval height Height of image
	 */
	size_type height( ) const;

	/**
		Get the number of rows in the image.
		@retval num_rows Number of rows in image
	 */
	size_type get_num_rows( ) const;
	
	/**
		Return number of columns in image, ignoring the current image region.
		@retval num_cols Number of columns
	 */
	size_type real_col_size( ) const;

	/**
		Get the number of columns in the image.
		@retval num_cols Number of columns in image
	 */
	size_type col_size( ) const;

	/**
		Get the number of columns in the image.
		@retval width Width of image
	 */
	size_type width( ) const;

	/**
		Get the number of columns in the image.
		@retval num_cols Number of columns in image
	 */
	size_type get_num_cols( ) const;

	/**
		Return number of pixels in image, ignoring the current image region.

		@retval num_pixels Number of pixels
	 */
	size_type real_size( ) const;

	/**
		Return number of pixels in image.

		The size actually refers to the current submatrix region.

		@retval num_pixels Number of pixels
	 */
	size_type size( ) const;

	/**
		Get the maximum pixel value that can be written by the image.

		Note that larger values can be stored in the image, which is useful for
		many image processing algorithms. However, when the image is written to
		file by a call to write(), any value larger than get_max_value() will be
		written as get_max_value().

		@retval _max_value Maximum pixel value image can hold
	 */
	pixel_type get_max_value( ) const;

	/**
		Get the maximum pixel value that can be written by the image.

		Note that larger values can be stored in the image, which is useful for
		many image processing algorithms. However, when the image is written to
		file by a call to write(), any value larger than get_max_value() will be
		written as get_max_value().

		@retval _max_value Maximum pixel value image can hold
	 */
	pixel_type max_value( ) const;

	/**
		Get the maximum pixel value that can be written by the image.

		Note that larger values can be stored in the image, which is useful for
		many image processing algorithms. However, when the image is written to
		file by a call to write(), any value larger than get_max_value() will be
		written as get_max_value().

		@retval _max_value Maximum pixel value image can hold
	 */
	pixel_type max_value( pixel_type );

	/**
		Get image type.
		@retval img_type Image type
	 */
	Image_Type get_type( ) const;

	/**
		Get image type.
		@retval img_type Image type
	 */
	Image_Type type( ) const;

	/**
		Change type of image.
		@param[in] new_type Type of image
	 */
	void type( const Image_Type& );

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	reference at( size_type, size_type );

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	const_reference at( size_type, size_type ) const;

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	reference operator()( size_type, size_type );

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	const_reference operator()( size_type, size_type ) const;

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	reference operator()( const Image_Coordinate& );

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	const_reference operator()( const Image_Coordinate& ) const;

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	reference gray( size_type, size_type );

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	const_reference gray( size_type, size_type ) const;

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	reference gray( const Image_Coordinate& );

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	const_reference gray( const Image_Coordinate& ) const;

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	reference red( size_type, size_type );

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	const_reference red( size_type, size_type ) const;

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	reference red( const Image_Coordinate& );

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	const_reference red( const Image_Coordinate& ) const;

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	reference green( size_type, size_type );

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	const_reference green( size_type, size_type ) const;

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	reference green( const Image_Coordinate& );

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	const_reference green( const Image_Coordinate& ) const;

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	reference blue( size_type, size_type );

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	const_reference blue( size_type, size_type ) const;

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	reference blue( const Image_Coordinate& );

	/**
		Access row i and column j in image.
		@retval ref Reference to row i and column j in image
	 */
	const_reference blue( const Image_Coordinate& ) const;

	/**
		Get pointer to internal data.
		@retval _matrix Internal data
	 */
	pointer data( );

	/**
		Get read-only pointer to internal data.
		@retval _matrix Internal data
	 */
	const_pointer data( ) const;

	/**
		Convert image to an array (allocated with 'new').
		The returned array must be dellocated by a call to 'delete[]'.
		@retval _matrix Internal data
	 */
	double* to_array( ) const;

	/**
		Determine if file is an image or not by its name. Can be passed to
		ws_tools::dir_traverse() as a file name filter.
		@param[in] file_name Name of file to check
	 */
	static bool is_image( const std::string& );

	/**
		Determine if file is a PNM image or not by its name. Can be passed to
		ws_tools::dir_traverse() as a file name filter.
		@param[in] file_name Name of file to check
	 */
	static bool is_pnm( const std::string& );

	/**
		Determine if file is a JPEG image or not by its name. Can be passed to
		ws_tools::dir_traverse() as a file name filter.
		@param[in] file_name Name of file to check
	 */
	static bool is_jpg( const std::string& );

	/**
		Determine if file is a MATLAB matrix or not by its name. Can be passed to
		ws_tools::dir_traverse() as a file name filter.
		@param[in] file_name Name of file to check
	 */
	static bool is_matrix( const std::string& );

	/**
		Determine if file is an video or not by its name. Can be passed to
		ws_tools::dir_traverse() as a file name filter.
		@param[in] file_name Name of file to check
	 */
	static bool is_video( const std::string& );

	/**
		Determines if image is a bilevel (black and white) image.
		@retval is_bilevel Whether image is bilevel--true or false
	 */
	bool is_bilevel( ) const;

	/**
		Determines if image is a gray-scale image.
		@retval is_gray Whether image is gray-scale--true or false
	 */
	bool is_gray( ) const;

	/**
		Determines if image is a color image.
		@retval is_color Whether image is in color--true or false
	 */
	bool is_color( ) const;

	Image& sort( );
	Image& reverse_sort( );

	/**
		Read image from file.
		@param[in] image_file_name Name of image file to read image from
	 */
	void read( const std::string& image_file_name )
	{
		// read PNM image
		if( is_pnm( image_file_name ) )
		{
			read_pnm( image_file_name );
		}

		// read JPEG image
		else if( is_jpg( image_file_name ) )
		{
			read_jpg( image_file_name );
		}
		else if( is_matrix( image_file_name ) )
		{
			// read matrix into the current channel; hence, separate matrices can
			// be read in that represent each color
			_current_data->read( image_file_name );
		}
		else // default to reading image if no file extension exists
		{
			// read file
			// read_pnm( image_file_name );

			// don't read file
			err_quit( "Unable to open image file '%s'\n", 
				image_file_name.c_str() );
		}
	}

	/**
		Write image to file.
		@param[in] image_file_name Name of image file to write image to
		@param[in] mode File mode to write as: binary or ascii (binary is default)
	 */
	void write( const std::string& image_file_name = "",
		const std::string mode = "b" ) const
	{
		// write as pnm image
		if( is_pnm( image_file_name ) )
		{
			write_pnm( image_file_name, mode );
		}

		// write as JPEG image
		else if( is_jpg( image_file_name ) )
		{
			write_jpg( image_file_name );
		}

		// write as matrix
		else if( is_matrix( image_file_name ) )
		{
			if( is_gray() )
			{
				_red_data.write( image_file_name );
				// _current_data->write( image_file_name );
			}

			// write separate matrices for each color channel
			else if( is_color() )
			{
				_red_data.write( "r_" + image_file_name );
				_green_data.write( "g_" + image_file_name );
				_blue_data.write( "b_" + image_file_name );
			}
		}
		else // default to writing as image if file type is undetermined
		{
			write_pnm( image_file_name, mode );
		}
	}

	/**
		Convert pixel value to an unsigned char that is suitable for printing.
		If value is negative, return 0; if value is larger than the maximum
		allowable pixel value, then return the maximum pixel value.
		Otherwise, round the value to the nearest integer/character.

		Note: this is inlined since it is called by most file writers.

		@param pixel Pixel to convert to unsigned character
		@param min_pixel Minimum value pixel may be
		@param max_pixel Maximum value pixel may be
		@retval pixel Converted pixel
	 */
	uchar pixel_to_uchar( pixel_type, pixel_type, pixel_type ) const;

protected:  // Internal member functions

	/**
		Allocate space for image.
		@param[in] num_rows Number of rows in new image
		@param[in] num_cols Number of columns in new image
	 */
	void allocate_image( size_type num_rows, size_type num_cols  )
	{
		// _red_data = Gray_Image( num_rows, num_cols );
		// if( is_color() )
		// {
			// _green_data = Gray_Image( num_rows, num_cols );
			// _blue_data  = Gray_Image( num_rows, num_cols );
		// }

		_red_data.resize( num_rows, num_cols, false );
		if( is_color() )
		{
			_green_data.resize( num_rows, num_cols, false );
			_blue_data.resize( num_rows, num_cols, false );
		}

		// use the red channel by default
		_current_data = &_red_data;
	}

	/**
		Write contents of raster image to PNM file.
		@param[in] image_file_name Name of image file to write image to
		@param[in] mode File mode to write as: binary or ascii (binary is default)
	 */
	void write_pnm( const std::string& = "", const std::string = "b" ) const;

	/**
		Read contents of PNM file into raster image.
		@param[in] image_file_name Name of image file to read image from
	 */
	void read_pnm( const std::string& );

	/**
		Write contents of raster image to JPEG file.
		@param[in] image_file_name Name of image file to write image to
	 */
	void write_jpg( const std::string& image_file_name = "" ) const;

	/**
		Read contents of JPEG file into raster image.
		@param[in] image_file_name Name of image file to read image from
	 */
	void read_jpg( const std::string& image_file_name = "" );

public:  // Miscellaneous public functions

	// public functions defined in Image.cpp
	Image& to_bilevel( );
	Image& to_gray( );
	Image& to_color( );

	void set_channel( const Color_Type& );

	void filter_gauss( );
	Image& convolve_row( const Vector&, Image& );
	Image& convolve_row( const Vector& );
	Image& convolve_col( const Vector&, Image& );
	Image& convolve_col( const Vector& );
	Image& convolve( const Vector&, const Vector& );
	Image& convolve( const Matrix& );

	Image transpose( );

	Image& crop( );
	Image& crop( const Image_Region& );

	Image& resize( size_type, size_type, bool = true );
	Image& scale( double );
	Image& scale( double, double );

	Image& nointerp_scale( double, double );

	Image& rotate( const double, const bool = false );
	Image& rotate( const double, const size_type, const size_type,
			const bool = false );

	Image& rotate( double, const Image_Coordinate&, bool = true );
	void rotate( const Image_Coordinate&, const Image_Coordinate&,
		double, double, int*, int* );

	Image& transform( double );

	// binary operators for two images
	Image& operator+=( const Image& );
	Image  operator+( const Image& ) const;
	Image& operator-=( const Image& );
	Image  operator-( const Image& ) const;
	Image& operator*=( const Image& );
	Image  operator*( const Image& ) const;

	// binary operators for image and value
	Image& operator+=( const value_type );
	Image  operator+( const value_type ) const;
	friend Image operator+( const value_type, Image& );

	Image& operator-=( const value_type );
	Image  operator-( const value_type ) const;
	friend Image operator-( const value_type, Image& );

	Image& operator*=( const value_type );
	Image  operator*( const value_type ) const;
	friend Image operator*( const value_type, Image& );

	Image& operator/=( const value_type );
	Image  operator/( const value_type ) const;

	/**
		Cast to matrix.
		@retval ref Constant reference to the current color channel
	 */
	operator const Matrix&() const
	{
		return( *_current_data );
	}

	Image& normalize( const pixel_type = DEFAULT_MAX );
	void stat( double*, double*, double* ) const;

	value_type get_min( ) const;
	value_type get_max( ) const;
	void get_min_max( value_type*, value_type* ) const;

	void sub_image( const Image_Region&, Image& );

	// drawing functions
	void draw( const Image_Coordinate&,
		const Image::pixel_type& = White_Pixel );
	void draw_line( const Image_Coordinate&,
		const Image_Coordinate&, const Image::pixel_type& = White_Pixel );

	void draw( const Image_Region&, const Image::pixel_type& = White_Pixel );
	void draw2( const Image_Region& );

	void draw_circle( const Image_Coordinate&, const double,
		const Image::pixel_type& = White_Pixel );

	// Fourier transform functions
	void fft( Image&, Image& ) const;
	void ifft( Image&, Image& );

	// lifting functions
	void lift( unsigned, unsigned, unsigned = 10000 );
	void inverse_lift( unsigned, unsigned, unsigned = 10000 );

	void lift( Lift_Scheme&, unsigned = 10000 );
	void inverse_lift( Lift_Scheme&, unsigned = 10000 );

	void octave_form( const Lift_Scheme&, unsigned = 10000 );
	void inverse_octave_form( const Lift_Scheme&, unsigned = 10000 );

	Wavelet_Subbands get_subbands( const Lift_Scheme&, unsigned = 10000 );

	// wavelet functions
	typedef Wavelet::Wavelet_Type Wavelet_Type;
	typedef Wavelet::Decomp_Type  Decomp_Type;

	void wavelet_transform( const Wavelet&, Decomp_Type = Wavelet::Standard );
	void inverse_wavelet_transform( const Wavelet&,
			Decomp_Type = Wavelet::Standard );

	void wavelet_transform( Wavelet_Type = Wavelet::Haar,
			Decomp_Type = Wavelet::Standard );
	void inverse_wavelet_transform( Wavelet_Type = Wavelet::Haar,
			Decomp_Type = Wavelet::Standard );

	Wavelet_Subbands get_subbands( const Wavelet&, unsigned = 10000 );

public: // Integral image functions

	Image to_integral( );
	Image& to_integral( Image& );

	/**
		Access sum of current region starting at row i and column j.
		The pixel sum of any upright rectangle r = (x, y, w, h) can be determined
		by 4 table lookups in the Summed Area Table (SAT):
			RecSum(r) = SAT(x - 1, y - 1)
				+ SAT(x + w - 1, y + h - 1)
				- SAT(x, y + h -1)
				- SAT(x + w - 1, y - 1)

		The function to_integral() is assumed to have been called prior to
		calling this function.  

		@param[in] i Starting row position
		@param[in] j Starting column position
		@param[in] h Height of region
		@param[in] w Width of region
		@retval area Area sum of region
	 */
	inline value_type get_area( size_type i, size_type j,
			size_type h, size_type w  ) const
	{
		// we store an additional, unused row and column at the start
		// of an integral image to simplify this computation, so we increment the
		// starting position by one
		// ++i;
		// ++j;

		// value_type area = at(i - 1, j - 1)
			// + at(i + h - 1, j + w - 1)
			// - at(i - 1, j + w - 1)
			// - at(i + h - 1, j - 1);

		// no need to subtract 1 since we store an additional row and column
		const value_type area = at(i, j)
			+ at(i + h, j + w)
			- at(i, j + w)
			- at(i + h, j);

		return( area );
	}

	/**
	  	Get area sum of region.

		The function to_integral() is assumed to have been called prior to
		calling this function.  

		@param[in] region Region to find area sum for
		@retval area Area sum of region
	 */
	inline value_type get_area( const Image_Region& region ) const
	{
		return(
			get_area( region.row_begin(), region.col_begin(), region.row_size(),
				region.col_size() )
		);
	}

public:  // Iterator functions

	/**
		Get iterator to start of image.
		@retval iter Iterator to row
	 */
	inline iterator begin( )
	{
		return( _current_data->begin() );
	}

	/**
		Get iterator to start of image.
		@retval iter Iterator to row
	 */
	inline const_iterator const_begin( ) const
	{
		return( _current_data->const_begin() );
	}

	/**
		Get iterator to end of image.
		@retval iter Iterator to row
	 */
	inline iterator end( )
	{
		return( _current_data->end() );
	}

	/**
		Get iterator to end of image.
		@retval iter Iterator to row
	 */
	inline const_iterator const_end( ) const
	{
		return( _current_data->const_end() );
	}

	void clear( )
	{
		_red_data.clear();
		if( is_color() )
		{
			_green_data.clear();
			_blue_data.clear();
		}
	}

public:  // Region functions

	/**
		Get starting row position of the matrix.
		
		The position actually refers to the current submatrix region, which is the
		first row of the matrix (0) unless explicitly changed.

		@retval row_size Number of rows in region
	 */
	inline size_type row_begin( ) const
	{
		return( _current_data->row_begin() );
	}

	/**
		Get starting column position of the matrix.

		The position actually refers to the current submatrix region, which is the
		first column of the matrix (0) unless explicitly changed.

		@retval col_size Number of columns in region
	 */
	inline size_type col_begin( ) const
	{
		return( _current_data->col_begin() );
	}

	/**
		Return last row position in region.
		@retval row_end 
	 */
	inline size_type row_end( ) const
	{ 
		return( _current_data->row_end() );
	}

	/**
		Return last column position in region.
		@retval num_rows Number of rows
	 */
	inline size_type col_end( ) const
	{ 
		return( _current_data->col_end() );
	}

	/**
		Return number of rows in submatrix marked by region.

		The size actually refers to the current submatrix region, which is the
		entire matrix unless explicitly changed.

		@retval num_rows Number of rows
	 */
	inline size_type sub_row_size( ) const
	{ 
		return( _current_data->sub_row_size() );	
	}

	/**
		Return number of columns in submatrix marked by region.

		The size actually refers to the current submatrix region, which is the
		entire matrix unless explicitly changed.

		@retval num_cols Number of columns
	 */
	inline size_type sub_col_size( ) const
	{ 
		return( _current_data->sub_col_size() );	
	}

	/**
		Get current image region.
		@retval region Current image region
	 */
	const Image_Region& get_region( ) const
	{
		return( _current_data->get_region() );
	}

	/**
		Get current image region.
		@retval region Current image region
	 */
	const Image_Region& region( ) const
	{
		return( _current_data->region() );
	}

	/**
		Set image region to the given region.

		@param[in] row_begin Starting row position of the region
		@param[in] col_begin Starting col position of the region
		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region

		@retval this This image
	 */
	const Image& set_region( size_type row_begin, size_type col_begin,
		size_type row_size, size_type col_size )
	{
		return( region( row_begin, col_begin, row_size, col_size ) );
	}

	/**
		Set image region to the given region.

		@param[in] row_begin Starting row position of the region
		@param[in] col_begin Starting col position of the region
		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region

		@retval this This image
	 */
	const Image& region( size_type row_begin, size_type col_begin,
		size_type row_size, size_type col_size )
	{
		const Image_Region new_region( row_begin, col_begin, row_size, col_size );
		return( region( new_region ) );
	}

	/**
		Set image region to the given region.
		@param[in] new_region New image region

		@retval this This image
	 */
	const Image& set_region( const Image_Region& new_region )
	{
		return( region( new_region ) );
	}

	/**
		Set image region to the given region.
		@param[in] new_region New image region

		@retval this This image
	 */
	const Image& region( const Image_Region& new_region )
	{
		_red_data.set_region( new_region );
		if( is_color() )
		{
			_green_data.set_region( new_region );
			_blue_data.set_region( new_region );
		}
		return( *this );
	}

	/**
		Save current region.
		This is useful if only the current region needs to be saved before
		changing to multiple regions that do not need saving.
	 */
	void save_region( )
	{
		_red_data.save_region();
		if( is_color() )
		{
			_green_data.save_region();
			_blue_data.save_region();
		}
	}

	/**
		Save current region and set new region.

		@param[in] row_begin Starting row position of the region
		@param[in] col_begin Starting col position of the region
		@param[in] row_size Number of rows in the region
		@param[in] col_size Number of columns in the region
	 */
	void push_region( size_type row_begin, size_type col_begin,
		size_type row_size, size_type col_size )
	{
		Image_Region region( row_begin, col_begin, row_size, col_size );
		push_region( region );
	}

	/**
		Save current region and set new region.
		@param[in] region New image region
	 */
	void push_region( const Image_Region& region )
	{
		_red_data.push_region( region );
		if( is_color() )
		{
			_green_data.push_region( region );
			_blue_data.push_region( region );
		}
	}

	/**
		Restore last saved region.
	 */
	void pop_region( )
	{
		_red_data.pop_region();
		if( is_color() )
		{
			_green_data.pop_region();
			_blue_data.pop_region();
		}
	}

	/**
		Restore last saved region.
	 */
	void restore_region( )
	{
		_red_data.restore_region();
		if( is_color() )
		{
			_green_data.restore_region();
			_blue_data.restore_region();
		}
	}

	/**
		Reset image region to the entire image.
	 */
	void reset_region( )
	{
		_red_data.reset_region();
		if( is_color() )
		{
			_green_data.reset_region();
			_blue_data.reset_region();
		}
	}

};

/*
	Define inlined member functions.
 */

/**
	Return number of rows in image, ignoring the current image region.

	@retval num_rows Number of rows
 */
inline Image::size_type Image::real_row_size( ) const
{ return( _current_data->real_row_size() ); }

/**
	Get the number of rows in the image.
	@retval num_rows Number of rows in image
 */
inline Image::size_type Image::row_size( ) const
{ return( _current_data->row_size() ); }

/**
	Get the number of rows in the image.
	@retval height Height of image
 */
inline Image::size_type Image::height( ) const
{ return( row_size() ); }

/**
	Get the number of rows in the image.
	@retval num_rows Number of rows in image
 */
inline Image::size_type Image::get_num_rows( ) const
{ return( row_size() ); }

/**
	Return number of columns in image, ignoring the current image region.

	@retval num_cols Number of columns
 */
inline Image::size_type Image::real_col_size( ) const
{ return( _current_data->real_col_size() ); }

/**
	Get the number of columns in the image.
	@retval num_cols Number of columns in image
 */
inline Image::size_type Image::col_size( ) const
{
	assert( row_size() >= 1 );
	return( _current_data->col_size() );
}

/**
	Get the number of columns in the image.
	@retval width Width of image
 */
inline Image::size_type Image::width( ) const
{ return( col_size() ); }

/**
	Get the number of columns in the image.
	@retval num_cols Number of columns in image
 */
inline Image::size_type Image::get_num_cols( ) const
{
	return( col_size() );
}

/**
	Return number of pixels in image, ignoring the current image region.

	@retval num_pixels Number of pixels
 */
inline Image::size_type
Image::real_size( ) const
{
	return( _current_data->real_size() );
}

/**
	Return number of pixels in image.

	The size actually refers to the current submatrix region.

	@retval num_pixels Number of pixels
 */
inline Image::size_type
Image::size( ) const
{
	return( _current_data->size() );
}

/**
	Get the maximum pixel value that can be written by the image.

	Note that larger values can be stored in the image, which is useful for
	many image processing algorithms. However, when the image is written to
	file by a call to write(), any value larger than get_max_value() will be
	written as get_max_value().

	@retval _max_value Maximum pixel value image can hold
 */
inline Image::pixel_type Image::get_max_value( ) const
{
	return( max_value() );
}

/**
	Get the maximum pixel value that can be written by the image.

	Note that larger values can be stored in the image, which is useful for
	many image processing algorithms. However, when the image is written to
	file by a call to write(), any value larger than get_max_value() will be
	written as get_max_value().

	@retval _max_value Maximum pixel value image can hold
 */
inline Image::pixel_type Image::max_value( ) const
{
	return( _max_value );
}

/**
	Get the maximum pixel value that can be written by the image.

	Note that larger values can be stored in the image, which is useful for
	many image processing algorithms. However, when the image is written to
	file by a call to write(), any value larger than get_max_value() will be
	written as get_max_value().

	@retval _max_value Maximum pixel value image can hold
 */
inline Image::pixel_type Image::max_value( pixel_type new_max_value )
{
	_max_value = new_max_value;
	return( max_value() );
}

/**
	Get image type.
	@retval img_type Image type
 */
inline Image::Image_Type Image::get_type( ) const
{ return( _img_type ); }

/**
	Get image type.
	@retval img_type Image type
 */
inline Image::Image_Type Image::type( ) const
{
	return( _img_type );
}

/**
	Change type of image.
	@param[in] new_type Type of image
 */
inline void Image::type( const Image_Type& new_type )
{
	_img_type = new_type;
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::reference Image::at( size_type i, size_type j )
{
	assert( _current_data == &_red_data || _current_data == &_green_data
			|| _current_data == &_blue_data );

	return( _current_data->at(i, j) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::const_reference Image::at( size_type i, size_type j ) const
{
	return( _current_data->at(i, j) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::reference Image::operator()( size_type i, size_type j )
{
	return( at(i, j) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::const_reference Image::operator()(
		size_type i, size_type j ) const
{
	return( at(i, j) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::reference Image::operator()( const Image_Coordinate& img_coord )
{
	return( at( img_coord.row, img_coord.col ) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::const_reference Image::operator()(
		const Image_Coordinate& img_coord ) const
{
	return( at( img_coord.row, img_coord.col ) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::reference Image::gray( size_type i, size_type j )
{
	return( at(i, j) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::const_reference Image::gray( size_type i, size_type j ) const
{
	return( at(i, j) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::reference Image::gray( const Image_Coordinate& img_coord )
{
	return( at( img_coord.row, img_coord.col ) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::const_reference Image::gray(
		const Image_Coordinate& img_coord ) const
{
	return( at( img_coord.row, img_coord.col ) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::reference Image::red( size_type i, size_type j )
{
	return( _red_data(i, j) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::const_reference Image::red( size_type i, size_type j ) const
{
	return( _red_data(i, j) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::reference Image::red( const Image_Coordinate& img_coord )
{
	return( _red_data( img_coord.row, img_coord.col ) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::const_reference Image::red(
		const Image_Coordinate& img_coord ) const
{
	return( _red_data( img_coord.row, img_coord.col ) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::reference Image::green( size_type i, size_type j )
{
	return( _green_data(i, j) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::const_reference Image::green( size_type i, size_type j ) const
{
	return( _green_data(i, j) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::reference Image::green( const Image_Coordinate& img_coord )
{
	return( _green_data( img_coord.row, img_coord.col ) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::const_reference Image::green(
		const Image_Coordinate& img_coord ) const
{
	return( _green_data( img_coord.row, img_coord.col ) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::reference Image::blue( size_type i, size_type j )
{
	return( _blue_data(i, j) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::const_reference Image::blue( size_type i, size_type j ) const
{
	return( _blue_data(i, j) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::reference Image::blue( const Image_Coordinate& img_coord )
{
	return( _blue_data( img_coord.row, img_coord.col ) );
}

/**
	Access row i and column j in image.
	@retval ref Reference to row i and column j in image
 */
inline Image::const_reference Image::blue(
		const Image_Coordinate& img_coord ) const
{
	return( _blue_data( img_coord.row, img_coord.col ) );
}

/**
	Get pointer to internal data.
	@retval _matrix Internal data
 */
inline Image::pointer Image::data( )
{
	return( _current_data->data() );
}

/**
	Get read-only pointer to internal data.
	@retval _matrix Internal data
 */
inline Image::const_pointer Image::data( ) const
{
	return( _current_data->data() );
}

/**
	Convert image to an array (allocated with 'new').
	The returned array must be dellocated by a call to 'delete[]'.
	@retval _matrix Internal data
 */
inline double* Image::to_array( ) const
{
	return( _current_data->to_array() );
}

/**
	Convert pixel value to an unsigned char that is suitable for printing.
	If value is negative, return 0; if value is larger than the maximum
	allowable pixel value, then return the maximum pixel value.
	Otherwise, round the value to the nearest integer/character.

	Note: this is inlined since it is called by most image file writers for each
	pixel.  

	@param pixel Pixel to convert to unsigned character
	@param min_pixel Minimum value pixel may be
	@param max_pixel Maximum value pixel may be
	@retval pixel Converted pixel
 */
inline Image::uchar
Image::pixel_to_uchar( pixel_type pixel,
		pixel_type min_pixel, pixel_type max_pixel ) const
{
	// if value < 0, write 0; if value is too large, write the max
	pixel = std::max<pixel_type>( pixel, min_pixel );
	pixel = std::min<pixel_type>( pixel, max_pixel );

	// round pixel to nearest integer before casting to unsigned char
	return( (uchar) round(pixel) );
}

/**
	Determine if file is an image or not by its name. Can be passed to
	ws_tools::dir_traverse() as a file name filter.
	@param[in] file_name Name of file to check
 */
inline bool
Image::is_image( const std::string& file_name )
{
	std::string ext = ws_tools::get_ext_name( file_name );
	if( is_pnm( file_name ) || is_jpg( file_name ) )
	{
		return( true );
	}
	return( false );
}

/**
	Determine if file is a PNM image or not by its name. Can be passed to
	ws_tools::dir_traverse() as a file name filter.
	@param[in] file_name Name of file to check
 */
inline bool
Image::is_pnm( const std::string& file_name )
{
	// TODO 'incomplete type' compilere error
	// return( Pnm_IO::is_pnm( file_name ) );

	std::string ext = ws_tools::get_ext_name( file_name );
	ws_tools::to_lower( ext );
	if( ext == "pgm" || ext == "ppm" || ext == "pbm" )
	{
		return( true );
	}
	return( false );
}

/**
	Determine if file is a JPEG image or not by its name. Can be passed to
	ws_tools::dir_traverse() as a file name filter.
	@param[in] file_name Name of file to check
 */
inline bool
Image::is_jpg( const std::string& file_name )
{
	std::string ext = ws_tools::get_ext_name( file_name );
	ws_tools::to_lower( ext );
	if( ext == "jpg" || ext == "jpeg" )
	{
		return( true );
	}
	return( false );
}

/**
	Determine if file is a MATLAB matrix or not by its name. Can be passed to
	ws_tools::dir_traverse() as a file name filter.
	@param[in] file_name Name of file to check
 */
inline bool
Image::is_matrix( const std::string& file_name )
{
	return( Matrix::is_matrix(file_name) );
}

/**
	Determine if file is an video or not by its name. Can be passed to
	ws_tools::dir_traverse() as a file name filter.
	@param[in] file_name Name of file to check
 */
inline bool
Image::is_video( const std::string& file_name )
{
	std::string ext = ws_tools::get_ext_name( file_name );
	ws_tools::to_lower( ext );
	if( ext == "mpg" || ext == "mpeg" || ext == "mp4" )
	{
		return( true );
	}
	else if( ext == "avi" || ext == "mov" || ext == "wmv" )
	{
		return( true );
	}
	return( false );
}

/**
	Determines if image is a bilevel (black and white) image.
	@retval is_bilevel Whether image is bilevel--true or false
 */
inline bool Image::is_bilevel( ) const
{
	return( type() == Bilevel );
}

/**
	Determines if image is a gray-scale image.
	@retval is_gray Whether image is gray-scale--true or false
 */
inline bool Image::is_gray( ) const
{
	return( type() == Gray );
}

/**
	Determines if image is a color image.
	@retval is_color Whether image is in color--true or false
 */
inline bool Image::is_color( ) const
{
	return( type() == Color );
}

inline Image&
Image::sort( )
{
	// _current_data->sort();

	_red_data.sort();
	if( is_color() )
	{
		_green_data.sort();
		_blue_data.sort();
	}
	return( *this );
}

inline Image&
Image::reverse_sort( )
{
	// _current_data->reverse_sort();

	_red_data.reverse_sort();
	if( is_color() )
	{
		_green_data.reverse_sort();
		_blue_data.reverse_sort();
	}

	return( *this );
}

/**
	Directory for reading image files.
 */
class Img_Dir : public ws_tools::Dir
{

public:

	Img_Dir( const std::string& dir = ".", bool recurse = true )
	{
		// always read, never write, to the directory
		open( dir, Dir::Read, recurse );
	}

	virtual ~Img_Dir( )
	{ }

protected:

	virtual bool filter( const std::string& file_name )
	{
		return( Image::is_image( file_name ) );
	}

};

} // namespace ws_img

#endif // IMAGE_HPP
