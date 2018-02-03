/**
	@file   Pnm_IO.cpp
	@author Wade Spires
	@date   2006/12/16
	@brief  Class Pnm_IO.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Pnm_IO.hpp"

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

// local headers
#include "Image.hpp"

using namespace ws_img;

using std::string;
using std::vector;

typedef Image::size_type   size_type;
typedef Image::value_type  value_type;
typedef Image::pixel_type  pixel_type;
typedef Image::uchar       uchar;

/**
	Default constructor.
 */
Pnm_IO::Pnm_IO( )
{ }

/**
	Destructor does nothing since no member variables are dynamically
	allocated.
 */
Pnm_IO::~Pnm_IO( )
{ }

/**
	Determine if file is a pnm image or not by its name. Can be passed to
	ws_tools::dir_traverse() as a file name filter.
	@param[in] file_name Name of file to check
 */
bool
Pnm_IO::is_pnm( const string& file_name )
{
	return( Image::is_pnm( file_name ) );

	/*
	string ext = ws_tools::get_ext_name( file_name );
	ws_tools::to_lower( ext );
	if( ext == "pgm" || ext == "ppm" || ext == "pbm" )
	{
		return( true );
	}
	return( false );
	 */
}

/**
	Read image file.
	@param[in] image_file_name Name of image file to read image from
 */
void
Pnm_IO::read( Image& img, const string& image_file_name )
{
	FILE*        image_file_fp;
	size_type    num_rows  = 0;
	size_type    num_cols  = 0;
	unsigned     max_value = 0;
	Image_Type   img_type  = PBM_ASCII;

	// open in binary mode
	if( (image_file_fp = fopen( image_file_name.c_str(), "rb" )) == NULL )
	{
		err_quit( "Unable to open image file '%s'\n",
				image_file_name.c_str() );
	}

	// read header information; allocate images; read in data

	read_header( image_file_fp, &max_value, &num_rows, &num_cols, &img_type );

	set_image_type( img, img_type );
	img.max_value( static_cast<Image::pixel_type>( max_value ) );

	// allocate an image with the correct size (no interpolation when resizing)
	img.resize( num_rows, num_cols, false );

	if( is_ascii(img_type) )
	{
		read_ascii_data( img, image_file_fp );
	}
	else if( is_binary(img_type) )
	{
		read_binary_data( img, image_file_fp );
	}

	if( fclose( image_file_fp ) != 0 )
	{
		err_warn( "Problems closing file '%s'\n", image_file_name.c_str() );
	}
}

/**
	Write contents of raster image to file.
	@param[in] image_file_name Name of image file to write image to
	@param[in] mode File mode to write as: binary or ascii (binary is default)
 */
void
Pnm_IO::write( const Image& img, const string& image_file_name,
	const string mode )
{
	FILE* image_file_fp = stdout;  // use stdout by default
	
	if( image_file_name != "" )
	{
		if( (image_file_fp = fopen( image_file_name.c_str(), "wb" )) == NULL )
		{
			err_quit( "Unable to open image file '%s'\n", 
					image_file_name.c_str() );
		}
	}

	if( mode == "b" )
	{
		write_binary_header( img, image_file_fp );
		write_binary_data( img, image_file_fp );
	}
	else if( mode == "a" )
	{
		write_ascii_header( img, image_file_fp );
		write_ascii_data( img, image_file_fp );
	}

	if( fflush( image_file_fp ) != 0 )
	{
		err_quit( "Unable to write to image file '%s'\n", 
				image_file_name.c_str() );
	}

	if( image_file_name != "" )
	{
		if( fclose( image_file_fp ) != 0 )
		{
			err_warn( "Problems closing file '%s'\n", image_file_name.c_str() );
		}
	}
}

/**
	Read header information from image file.
	@param[in] image_file_fp image file pointer
 */
void
Pnm_IO::read_header( FILE* image_file_fp, unsigned* max_value,
	size_type* num_rows, size_type* num_cols, Image_Type* img_type )
{
	// temporary place to store line read in from file
	char tmp_str[ BUFSIZ ];
	const char* white_space = " \t\r\n";

	// flags to determine which fields have been read in
	bool have_read_magic_id  = false;
	bool have_read_rows      = false;
	bool have_read_cols      = false;
	bool have_read_max_value = false;

	do // read lines until all parts of header have been read (or error)
	{
		// read in line
		fgets( tmp_str, BUFSIZ, image_file_fp );

		// convert to C++ string; remove comments from line
		string line = tmp_str;
		remove_comments( line );

		// skip empty lines
		if( line.empty() )
		{
			continue;
		}

		// find first word in line
		string::size_type pos_1 = line.find_first_not_of( white_space );
		string::size_type pos_2 = line.find_first_of( white_space, pos_1 );

		// skip if line is all whitespace
		if( pos_1 == string::npos ) 
		{
			continue;
		}

		// try to read P1, P2, ...
		if( !have_read_magic_id )
		{
			assert( img_type != 0 );
			*img_type = read_magic_id( line, pos_1, pos_2 );
			have_read_magic_id = true;

			// PBM images do not store a max value, so we prevent looping back
			// and attempting to read the image data after reading the number
			// of rows
			if( *img_type == PBM_ASCII || *img_type == PBM_BINARY )
			{
				*max_value = 1;
				have_read_max_value = true;
			}

			// try to get next word or try next line
			pos_1 = line.find_first_not_of( white_space, pos_2 );
			pos_2 = line.find_first_of( white_space, pos_1 );
			if( pos_1 == string::npos ) continue;
		}

		if( !have_read_cols )
		{
			assert( num_cols != 0 );
			*num_cols = read_integer( line, pos_1, pos_2 );
			have_read_cols = true;

			// try to get next word or try next line
			pos_1 = line.find_first_not_of( white_space, pos_2 );
			pos_2 = line.find_first_of( white_space, pos_1 );

			if( pos_1 == string::npos ) continue;
		}

		if( !have_read_rows )
		{
			assert( num_rows != 0 );
			*num_rows = read_integer( line, pos_1, pos_2 );
			have_read_rows = true;

			// try to get next word or try next line
			pos_1 = line.find_first_not_of( white_space, pos_2 );
			pos_2 = line.find_first_of( white_space, pos_1 );

			if( pos_1 == string::npos ) continue;
		}

		if( !have_read_max_value )
		{
			assert( *img_type != PBM_ASCII && *img_type != PBM_BINARY );

			*max_value = read_integer( line, pos_1, pos_2 );
			have_read_max_value = true;
		}
	}
	while( !( have_read_magic_id && have_read_cols && have_read_rows 
		&& have_read_max_value ) );
}

/**
	Remove comments from string. Comments begin with the '#' character.
	@param[in,out] str String to remove comments from
 */
void
Pnm_IO::remove_comments( string& str )
{
	// find start of comment
	string::size_type comment_pos = str.find_first_of( '#' );

	// nothing to do if no comments
	if( comment_pos == string::npos )
	{
		return;
	}

	// erase comments
	str.erase( comment_pos );
}

/**
	Remove white space from string. White space consists of spaces, tabs, and
	newlines.
	@param[in,out] str String to remove white space from
 */
void
Pnm_IO::remove_white_space( string& str )
{
	// points one past the end of the unremoved elements after remove()
	string::iterator new_end;

	// remove spaces
	new_end = std::remove( str.begin(), str.end(), ' ' );
	str.assign( str.begin(), new_end );

	// remove tabs
	new_end = std::remove( str.begin(), str.end(), '\t' );
	str.assign( str.begin(), new_end );

	// remove newlines
	new_end = std::remove( str.begin(), str.end(), '\n' );
	str.assign( str.begin(), new_end );
}

/**
	Read magic ID at start of PGM to determine if data is ascii or binary.
	@param[in] line Line from 
	@param[in] pos_1 Position in line of non-whitespace character
	@param[in] pos_2 Position after pos_1 containing next whitespace
	(or npos)
 */
Pnm_IO::Image_Type
Pnm_IO::read_magic_id( const string& line,
	const string::size_type pos_1,
	const string::size_type pos_2 )
{
	if( ( pos_2 == string::npos && line.length() - pos_1 != 2 )
		|| ( pos_2 != string::npos && ( pos_2 - pos_1 ) != 2 ) )
	{
		err_quit( "Magic code 'P[1-6]' does not exist\n" );
	}

	if( line.at( pos_1 ) != 'P')
	{
		err_quit( "Magic code 'P[1-6]' does not exist\n" );
	}

	// determine type of data
	if( line.at( pos_1 + 1 ) == '1' ) 
	{
		return( PBM_ASCII );
	}
	else if( line.at( pos_1 + 1 ) == '2' ) 
	{
		return( PGM_ASCII );
	}
	else if( line.at( pos_1 + 1 ) == '3' )
	{
		return( PPM_ASCII );
	}
	else if( line.at( pos_1 + 1 ) == '4' ) 
	{
		return( PBM_BINARY );
	}
	else if( line.at( pos_1 + 1 ) == '5' )
	{
		return( PGM_BINARY );
	}
	else if( line.at( pos_1 + 1 ) == '6' )
	{
		return( PPM_BINARY );
	}
	else
	{
		err_quit( "Magic code 'P[1-6]' does not exist\n" );
	}

	// return something to quell compiler warnings
	return( PGM_BINARY );
}

/**
	Read unsigned integer from header of PGM (row, column, or max number).
	@param[in] line Line from 
	@param[in] pos_1 Position in line of non-whitespace character
	@param[in] pos_2 Position after pos_1 containing next whitespace
	(or npos)
 */
unsigned
Pnm_IO::read_integer( const string& line,
	const string::size_type pos_1,
	const string::size_type pos_2 )
{
	string tmp_str;

	// assign tmp_str to word from line
	if( ( pos_2 == string::npos ) )
	{
		tmp_str.assign( line, pos_1, line.length() - pos_1 );
	}
	else
	{
		tmp_str.assign( line, pos_1, pos_2 - pos_1 );
	}

	// convert string to integer
	int parsed_int = atoi( tmp_str.c_str() ); 

	// verify it's a valid value
	if( parsed_int <= 0 )
	{
		err_quit( "Field %d in header is invalid\n", parsed_int );
	}

	return( parsed_int );
}

/**
	Determine type of image.
	@param[in] type Type of image (internal format)
 */
void
Pnm_IO::set_image_type( Image& img, const Image_Type& type )
{
	switch( type )
	{
		case PBM_ASCII:
		case PBM_BINARY:
			img.type( Image::Bilevel );
			break;

		case PGM_ASCII:
		case PGM_BINARY:
			img.type( Image::Gray );
			break;

		case PPM_ASCII:
		case PPM_BINARY:
			img.type( Image::Color );
			break;

		default:  // should not get here
			err_quit( "Pnm_IO::Unknown image type\n" );
	}
}

/**
	Read in ascii image data from image file.
	@param[in] image_file_fp File pointer to image file currently reading
 */
void
Pnm_IO::read_ascii_data( Image& img, FILE* image_file_fp )
{
	// determine image type
	if( img.is_bilevel() )
	{
		// read each pixel from file and store in raster
		for( size_type i = img.row_begin(); i != img.row_end(); ++i )
		{
			for( size_type j = img.col_begin(); j != img.col_end(); ++j )
			{
				int pixel_value;
				fscanf( image_file_fp, "%d", &pixel_value );

				// 0 -> white, 1 -> black
				// note: we only check the right-most bit to determine the pixel
				pixel_value = (pixel_value & 01) ? 1 : 0;

				img(i,j) = (pixel_type) pixel_value;
			}
		}
	}
	else if( img.is_gray() )
	{
		// read each pixel from file and store in raster
		for( size_type i = img.row_begin(); i != img.row_end(); ++i )
		{
			for( size_type j = img.col_begin(); j != img.col_end(); ++j )
			{
				int pixel_value;
				fscanf( image_file_fp, "%d", &pixel_value );
				img(i,j) = (pixel_type) pixel_value;
			}
		}
	}
	else if( img.is_color() )
	{
		// read each pixel from file and store in raster
		for( size_type i = img.row_begin(); i != img.row_end(); ++i )
		{
			for( size_type j = img.col_begin(); j != img.col_end(); ++j )
			{
				int red_value;
				int green_value;
				int blue_value;

				fscanf( image_file_fp, "%d", &red_value );
				fscanf( image_file_fp, "%d", &green_value );
				fscanf( image_file_fp, "%d", &blue_value );

				img.red(i, j)   = (pixel_type) red_value;
				img.green(i, j) = (pixel_type) green_value;
				img.blue(i, j)  = (pixel_type) blue_value;
			}
		}
	}
}

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
void
Pnm_IO::read_binary_data( Image& img, FILE* image_file_fp )
{
	// determine image type
	if( img.is_bilevel() )
	{
		const unsigned bits_per_byte = 8;

		// number of bits in the last byte of a row, which may not equal 
		// the number of bits per byte if the last byte was padded
		const unsigned num_bits_in_last_byte =
				(img.sub_col_size() % bits_per_byte);

		// number of full bytes in one row to read from the file (not counting
		// the last byte if it has padding bits)
		const unsigned num_full_bytes = (img.sub_col_size() / bits_per_byte);

		// read each row
		for( size_type i = img.row_begin(); i != img.row_end(); ++i )
		{
			size_type j = img.col_begin();

			// read each byte except for the last one
			for( unsigned b = 0; b != num_full_bytes; ++b )
			{
				uchar read_in_byte = fgetc( image_file_fp );

				// handle each bit in byte (from most significant bit to least)
				for( unsigned current_bit_pos = bits_per_byte;
						current_bit_pos != 0; --current_bit_pos )
				{
					// examine current bit to determine pixel value
					if( (read_in_byte >> (current_bit_pos - 1)) & 01 )
					{
						img(i,j) = 1;
					}
					else
					{
						img(i,j) = 0;
					}

					++j;  // move to next column in image
				}
			}

			// handle last byte if some bits bits are used for padding
			if( num_bits_in_last_byte != 0 )
			{
				// handle each bit in byte up to the amount of padding
				uchar read_in_byte = fgetc( image_file_fp );
				for( unsigned current_bit_pos = bits_per_byte;
						current_bit_pos != (bits_per_byte - num_bits_in_last_byte);
						--current_bit_pos )
				{
					// examine current bit to determine pixel value
					if( (read_in_byte >> (current_bit_pos - 1)) & 01 )
					{
						img(i,j) = 1;
					}
					else
					{
						img(i,j) = 0;
					}

					++j;  // move to next column in image
				}
			}
		}
	}
	else if( img.is_gray() )
	{
		// read each pixel from file and store in raster
		for( size_type i = img.row_begin(); i != img.row_end(); ++i )
		{
			for( size_type j = img.col_begin(); j != img.col_end(); ++j )
			{
				const uchar pixel_value = fgetc( image_file_fp );
				img(i,j) = (pixel_type) pixel_value;
			}
		}
	}
	else if( img.is_color() )
	{
		// read each pixel from file and store in raster
		for( size_type i = img.row_begin(); i != img.row_end(); ++i )
		{
			for( size_type j = img.col_begin(); j != img.col_end(); ++j )
			{
				// read each color channel
				const uchar red_value   = fgetc( image_file_fp );
				const uchar green_value = fgetc( image_file_fp );
				const uchar blue_value  = fgetc( image_file_fp );

				img.red(i,j)   = (pixel_type) red_value;
				img.green(i,j) = (pixel_type) green_value;
				img.blue(i,j)  = (pixel_type) blue_value;
			}
		}
	}
}

/**
	Determines if image is an binary image.
	@param[in] img_type Image type (internal format)
	@retval is_binary Whether image is binary--true or false
 */
bool
Pnm_IO::is_binary( const Image_Type& img_type )
{
	return( img_type == PBM_BINARY || img_type == PGM_BINARY
		|| img_type == PPM_BINARY );
}

/**
	Determines if image is an ascii image.
	@param[in] img_type Image type (internal format)
	@retval is_ascii Whether image is ascii--true or false
 */
bool
Pnm_IO::is_ascii( const Image_Type& img_type )
{
	return( img_type == PBM_ASCII || img_type == PGM_ASCII
		|| img_type == PPM_ASCII );
}

/**
	Write header information to image file.
	@param[in] image_file_fp File pointer to image file currently writing
 */
void
Pnm_IO::write_binary_header( const Image& img, FILE* image_file_fp )
{
	// embed time stamp into comment line
	string time_stamp = get_time_stamp();

	if( img.is_bilevel() )
	{
		fprintf( image_file_fp, "P4\n" );
		fprintf( image_file_fp, "# Created by Wade Spires on %s",
			time_stamp.c_str() );
		fprintf( image_file_fp, "%d %d\n",
				img.sub_col_size(), img.sub_row_size() );
		// bilevel images do not have a maximum value in the file
	}
	else if( img.is_gray() )
	{
		fprintf( image_file_fp, "P5\n" );
		fprintf( image_file_fp, "# Created by Wade Spires on %s",
			time_stamp.c_str() );
		fprintf( image_file_fp, "%d %d\n",
				img.sub_col_size(), img.sub_row_size() );
		fprintf( image_file_fp, "%u\n", unsigned( img.max_value() ) );
	}
	else if( img.is_color() )
	{
		fprintf( image_file_fp, "P6\n" );
		fprintf( image_file_fp, "# Created by Wade Spires on %s",
			time_stamp.c_str() );
		fprintf( image_file_fp, "%d %d\n",
				img.sub_col_size(), img.sub_row_size() ); 
		fprintf( image_file_fp, "%u\n", unsigned( img.max_value() ) );
	}
}

/**
	Write header information to image file.
	@param[in] image_file_fp File pointer to image file currently writing
 */
void
Pnm_IO::write_ascii_header( const Image& img, FILE* image_file_fp )
{
	string time_stamp = get_time_stamp();

	if( img.is_bilevel() )
	{
		fprintf( image_file_fp, "P1\n" );
		fprintf( image_file_fp, "# Created by Wade Spires on %s",
			time_stamp.c_str() );
		fprintf( image_file_fp, "%d %d\n",
				img.sub_col_size(), img.sub_row_size() );
		// bilevel images do not have a maximum value in the file
	}
	else if( img.is_gray() )
	{
		fprintf( image_file_fp, "P2\n" );
		fprintf( image_file_fp, "# Created by Wade Spires on %s",
			time_stamp.c_str() );
		fprintf( image_file_fp, "%d %d\n",
				img.sub_col_size(), img.sub_row_size() );
		fprintf( image_file_fp, "%u\n", unsigned( img.max_value() ) );
	}
	else if( img.is_color() )
	{
		fprintf( image_file_fp, "P3\n" );
		fprintf( image_file_fp, "# Created by Wade Spires on %s",
			time_stamp.c_str() );
		fprintf( image_file_fp, "%d %d\n",
				img.sub_col_size(), img.sub_row_size() ); 
		fprintf( image_file_fp, "%u\n", unsigned( img.max_value() ) );
	}
}

/**
	Get time stamp to record time that image was written. 
	@retval time_stamp Current time in string format
 */
string
Pnm_IO::get_time_stamp( )
{
	// get time stamp
	time_t current_time = time( NULL );
	char time_str[ BUFSIZ ];

	assert( BUFSIZ >= 26 );

	// store time in given string--ctime_r() is thread-safe, unlike ctime()
	// TODO Is ctime_r available for Windows? Perhaps we need a #ifdef?
	ctime_r( &current_time, time_str );

	// if ctime_r() is not available, but ctime() is (note: not thread-safe)
	// strncpy( time_str, ctime( &current_time ), BUFSIZ );

	// if ctime() is not available, we must write a newline
	// time_str[0] = '\n';
	// time_str[1] = '\0';

	return( string( time_str ) );
}

/**
	Write image data to image file in binary form.
	@param[in] image_file_fp File pointer to image file currently writing
 */
void
Pnm_IO::write_binary_data( const Image& img, FILE* image_file_fp )
{
	const pixel_type min_pixel = 0;
	const pixel_type max_pixel = img.max_value();

	if( img.is_bilevel() )
	{
		const unsigned bits_per_byte = 8;

		for( size_type i = img.row_begin(); i != img.row_end(); i++ )
		{
			// initialize start of next row: fill_byte is byte to write to image
			// file where each bit represents one pixel; current_bit_pos is
			// the current bit in fill_byte being set
			uchar fill_byte = 0;
			unsigned current_bit_pos = 1 << (bits_per_byte - 1);

			for( size_type j = img.col_begin(); j != img.col_end(); j++ )
			{
				// if byte is completely filled, write it to file, clear the
				// byte value, and reset the bit position (note: this will not
				// be called for the last byte, so we must write after the loop)
				if( current_bit_pos == 0 )
				{
					fwrite( &fill_byte, sizeof( uchar ), 1, image_file_fp );
					fill_byte = 0;
					current_bit_pos = 1 << (bits_per_byte - 1);
				}

				uchar pixel = img.pixel_to_uchar( img(i,j),
						min_pixel, max_pixel );
				if( pixel != 0 )
				{
					// insert a '1' bit at the current position in the byte
					fill_byte |= current_bit_pos;
				}
				current_bit_pos >>= 1;  // move to next bit
			}

			// write last byte in the row
			// note: last byte is automatically padded with 0s if the width is
			// not divisible by 8
			fwrite( &fill_byte, sizeof( uchar ), 1, image_file_fp );
		}
	}
	else if( img.is_gray() )
	{
		for( size_type i = img.row_begin(); i != img.row_end(); i++ )
		{
			for( size_type j = img.col_begin(); j != img.col_end(); j++ )
			{
				// write pixel to file
				const uchar pixel = img.pixel_to_uchar( img(i,j),
						min_pixel, max_pixel );
				fwrite( &pixel, sizeof( uchar ), 1, image_file_fp );
			}
		}
	}
	else if( img.is_color() )
	{
		for( size_type i = img.row_begin(); i != img.row_end(); i++ )
		{
			for( size_type j = img.col_begin(); j != img.col_end(); j++ )
			{
				const uchar red_pixel   = img.pixel_to_uchar( img.red(i,j),
						min_pixel, max_pixel );
				const uchar green_pixel = img.pixel_to_uchar( img.green(i,j),
						min_pixel, max_pixel );
				const uchar blue_pixel  = img.pixel_to_uchar( img.blue(i,j),
						min_pixel, max_pixel );

				// write pixels to file
				fwrite( &red_pixel, sizeof( uchar ), 1, image_file_fp );
				fwrite( &green_pixel, sizeof( uchar ), 1, image_file_fp );
				fwrite( &blue_pixel, sizeof( uchar ), 1, image_file_fp );
			}
		}
	}
}

/**
	Write image data to image file in ascii form.
	Note: Only 70 characters are written per line, which adds to the code's
	complexity (all '*_width' variables are for this purpose); however, this
	conforms to the formal specification for an ascii image.

	@param[in] image_file_fp File pointer to image file currently writing
 */
void
Pnm_IO::write_ascii_data( const Image& img, FILE* image_file_fp )
{
	const pixel_type min_pixel = 0;
	const pixel_type max_pixel = img.max_value();

	const size_type max_line_width = 70; // max. characters allowed on line
	size_type line_width  = 0;           // current width of line

	if( img.is_bilevel() || img.is_gray() )
	{
		for( size_type i = img.row_begin(); i != img.row_end(); i++ )
		{
			for( size_type j = img.col_begin(); j != img.col_end(); j++ )
			{
				// convert pixel to writeable type
				uchar pixel = img.pixel_to_uchar( img(i,j),
						min_pixel, max_pixel );

				// determine how much space the values will need on the line
				size_type pixel_width = 0;
				if( pixel >= 100 )     pixel_width = 4;
				else if( pixel >= 10 ) pixel_width = 3;
				else                   pixel_width = 2;

				// determine if need a new line
				line_width += pixel_width;
				if( line_width >= max_line_width )
				{
					fprintf( image_file_fp, "\n" );
					line_width = pixel_width;
				}

				// write pixel to file
				fprintf( image_file_fp, "%u ", pixel );
			}
		}
	}
	else if( img.is_color() )
	{
		for( size_type i = img.row_begin(); i != img.row_end(); i++ )
		{
			for( size_type j = img.col_begin(); j != img.col_end(); j++ )
			{
				// convert each pixel to writeable type
				uchar red_pixel   = img.pixel_to_uchar( img.red(i,j),
						min_pixel, max_pixel );
				uchar green_pixel = img.pixel_to_uchar( img.green(i,j),
						min_pixel, max_pixel );
				uchar blue_pixel  = img.pixel_to_uchar( img.blue(i,j),
						min_pixel, max_pixel );

				// determine how much space the values will need on the line
				size_type red_width = 0;
				if( red_pixel >= 100 )     red_width = 4;
				else if( red_pixel >= 10 ) red_width = 3;
				else                       red_width = 2;

				size_type green_width = 0;
				if( green_pixel >= 100 )     green_width = 4;
				else if( green_pixel >= 10 ) green_width = 3;
				else                         green_width = 2;

				size_type blue_width = 0;
				if( blue_pixel >= 100 )     blue_width = 4;
				else if( blue_pixel >= 10 ) blue_width = 3;
				else                        blue_width = 2;

				size_type total_width = red_width + green_width + blue_width;

				// determine if need a new line
				line_width += total_width;
				if( line_width >= max_line_width )
				{
					fprintf( image_file_fp, "\n" );
					line_width = total_width;
				}

				// write pixels to file
				fprintf( image_file_fp, "%u ", red_pixel );
				fprintf( image_file_fp, "%u ", green_pixel );
				fprintf( image_file_fp, "%u ", blue_pixel );
			}
		}
	}
}
