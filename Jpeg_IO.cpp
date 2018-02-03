/**
	@file   Jpeg_IO.cpp
	@author Wade Spires
	@date   2005/12/8
	@brief  Class Jpeg_IO reads and writes JPEG images.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Jpeg_IO.hpp"

// c++ headers
#include <algorithm>
#include <vector>

// c headers
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

// libjpeg headers
#if HAVE_LIBJPEG
#include <jpeglib.h>
#include <setjmp.h>
#endif

// my headers
#include "ws_tools.hpp"

using std::string;
using std::vector;

using namespace ws_img;

/**
	Write image to file in the JPEG format.
	This code closely follows example.c included with libjpeg found at
	http://www.ijg.org/

	@param[in] img Image to write
	@param[in] file_name Name of file to write
	@param[in] quality Image quality in the range [0,100] where higher values
		imply better quality but less compression
 */
void
Jpeg_IO::write( const Image& img, const std::string& file_name, int quality )
{
#if HAVE_LIBJPEG
	const Image::pixel_type min_pixel = 0;
	const Image::pixel_type max_pixel = img.max_value();

	// copy of image for libjpeg to use
	JSAMPLE* image_buffer = NULL;

	// set number of color components per pixel and assign the color space
	if( img.is_bilevel() || img.is_gray() )
	{
		// allocate single channel
		image_buffer = (JSAMPLE*) malloc(
				sizeof(JSAMPLE) * img.sub_row_size() * img.sub_col_size() );

		// copy image into image buffer
		unsigned k = 0;
		for( Image::size_type i = img.row_begin(); i != img.row_end(); ++i )
		{
			for( Image::size_type j = img.col_begin(); j != img.col_end(); ++j )
			{
				image_buffer[k++] = img.pixel_to_uchar( img(i,j),
						min_pixel, max_pixel );
			}
		}
	}
	else // img.is_color() == true
	{
		// allocate 3 color channels
		image_buffer = (JSAMPLE*) malloc(
				3 * sizeof(JSAMPLE) * img.sub_row_size() * img.sub_col_size() );

		// copy image into image buffer
		unsigned k = 0;
		for( Image::size_type i = img.row_begin(); i != img.row_end(); ++i )
		{
			for( Image::size_type j = img.col_begin(); j != img.col_end(); ++j )
			{
				image_buffer[k++] = img.pixel_to_uchar( img.red(i,j),
						min_pixel, max_pixel );
				image_buffer[k++] = img.pixel_to_uchar( img.green(i,j),
						min_pixel, max_pixel );
				image_buffer[k++] = img.pixel_to_uchar( img.blue(i,j),
						min_pixel, max_pixel );
			}
		}
	}

	// This struct contains the JPEG compression parameters and pointers to
	// working space (which is allocated as needed by the JPEG library).
	// It is possible to have several such structures, representing multiple
	// compression/decompression processes, in existence at once.  We refer
	// to any one struct (and its associated working data) as a "JPEG object".
	struct jpeg_compress_struct cinfo;

	// This struct represents a JPEG error handler.  It is declared separately
	// because applications often want to supply a specialized error handler
	// (see the second half of this file for an example).  But here we just
	// take the easy way out and use the standard error handler, which will
	// print a message on stderr and call exit() if compression fails.
	// Note that this struct must live as long as the main JPEG parameter
	// struct, to avoid dangling-pointer problems.
	struct jpeg_error_mgr jerr;

	FILE*     out_fp;          // target file
	JSAMPROW  row_pointer[1];  // pointer to JSAMPLE row[s]
	int       row_stride;      // physical row width in image buffer

	/* Step 1: allocate and initialize JPEG compression object */

	// We have to set up the error handler first, in case the initialization
	// step fails.  (Unlikely, but it could happen if you are out of memory.)
	// This routine fills in the contents of struct jerr, and returns jerr's
	// address which we place into the link field in cinfo.
	cinfo.err = jpeg_std_error( &jerr );

	// Now we can initialize the JPEG compression object.
	jpeg_create_compress( &cinfo );

	/* Step 2: specify data destination (eg, a file) */
	/* Note: steps 2 and 3 can be done in either order. */

	// Here we use the library-supplied code to send compressed data to a
	// stdio stream.  You can also write your own code to do something else.
	// VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
	// requires it in order to write binary files.
	if( (out_fp = fopen( file_name.c_str(), "wb")) == NULL )
	{
		err_quit( "Unable to open file '%s'\n", file_name.c_str() );
	}
	jpeg_stdio_dest( &cinfo, out_fp );

	/* Step 3: set parameters for compression */

	// First we supply a description of the input image.
	// Four fields of the cinfo struct must be filled in:
	cinfo.image_height = img.sub_row_size();  // image width and height
	cinfo.image_width  = img.sub_col_size();

	// set number of color components per pixel and assign the color space
	if( img.is_bilevel() || img.is_gray() )
	{
		cinfo.input_components = 1;
		cinfo.in_color_space   = JCS_GRAYSCALE;
		row_stride = img.sub_col_size();  // JSAMPLEs per row in image_buffer
	}
	else // img.is_color() == true
	{
		cinfo.input_components = 3;
		cinfo.in_color_space   = JCS_RGB;
		row_stride = 3 * img.sub_col_size();  // JSAMPLEs per row in image_buffer
	}

	// Now use the library's routine to set default compression parameters.
	// (You must set at least cinfo.in_color_space before calling this,
	// since the defaults depend on the source color space.)
	jpeg_set_defaults( &cinfo );

	// Now you can set any non-default parameters you wish to.
	// Here we just illustrate the use of quality (quantization table) scaling:
	jpeg_set_quality( &cinfo, quality, TRUE /* limit to baseline-JPEG values*/ );

	/* Step 4: Start compressor */

	// TRUE ensures that we will write a complete interchange-JPEG file.
	// Pass TRUE unless you are very sure of what you're doing.
	jpeg_start_compress( &cinfo, TRUE );

	/* Step 5: while (scan lines remain to be written) */
	/*           jpeg_write_scanlines(...); */


	// Here we use the library's state variable cinfo.next_scanline as the
	// loop counter, so that we don't have to keep track ourselves.
	// To keep things simple, we pass one scanline per call; you can pass
	// more if you wish, though.
	while( cinfo.next_scanline < cinfo.image_height )
	{
		// jpeg_write_scanlines expects an array of pointers to scanlines.
		// Here the array is only one element long, but you could pass
		// more than one scanline at a time if that's more convenient.
		row_pointer[0] = &image_buffer[ cinfo.next_scanline * row_stride ];
		(void) jpeg_write_scanlines( &cinfo, row_pointer, 1 );
	}

	/* Step 6: Finish compression */

	jpeg_finish_compress( &cinfo );

	// After finish_compress, we can close the output file.
	fclose( out_fp );

	/* Step 7: release JPEG compression object */

	// This is an important step since it will release a good deal of memory.
	jpeg_destroy_compress( &cinfo );

	free( image_buffer );

#else
	err_warn( "Unable to read image '%s': libjpeg (http://www.ijg.org/)"
			" is not installed.\n", file_name.c_str() );
	return;
#endif
}

#if HAVE_LIBJPEG

/*
	Here's the extended error handler struct:
 */
struct my_error_mgr
{
	struct jpeg_error_mgr  pub;            // "public" fields
	jmp_buf                setjmp_buffer;  // for return to caller
};

typedef struct my_error_mgr* my_error_ptr;

/*
	Here's the routine that will replace the standard error_exit method:
 */
METHODDEF(void)
my_error_exit( j_common_ptr cinfo )
{
	// cinfo->err really points to a my_error_mgr struct, so coerce pointer
	my_error_ptr myerr = (my_error_ptr) cinfo->err;

	// Always display the message.
	// We could postpone this until after returning, if we chose.
	(*cinfo->err->output_message)( cinfo );

	// Return control to the setjmp point
	longjmp( myerr->setjmp_buffer, 1 );
}

#endif

/**
	Read image from file in the JPEG format.
	This code closely follows example.c included with libjpeg found at
	http://www.ijg.org/

	@param[in] img Image to write
	@param[in] file_name Name of file to write
 */
int
Jpeg_IO::read( Image& img, const std::string& file_name )
{
#if HAVE_LIBJPEG

	// This struct contains the JPEG decompression parameters and pointers to
	// working space (which is allocated as needed by the JPEG library).
	struct jpeg_decompress_struct cinfo;

	// We use our private extension JPEG error handler.
	// Note that this struct must live as long as the main JPEG parameter
	// struct, to avoid dangling-pointer problems.
	struct my_error_mgr jerr;

	FILE*       in_fp;        // input file
	JSAMPARRAY  buffer;       // output row buffer
	int         row_stride;   // physical row width in output buffer

	// In this example, we open the input file before doing anything else,
	// so that the setjmp() error recovery below can assume the file is open.
	// VERY IMPORTANT: use "b" option to fopen() if you are on a machine that
	// requires it in order to read binary files.
	if( (in_fp = fopen(file_name.c_str(), "rb")) == NULL )
	{
		err_quit( "Unable to open file '%s'\n", file_name.c_str() );
	}

	/* Step 1: allocate and initialize JPEG decompression object */

	// We set up the normal JPEG error routines, then override error_exit.
	cinfo.err = jpeg_std_error(&jerr.pub);
	jerr.pub.error_exit = my_error_exit;

	// Establish the setjmp return context for my_error_exit to use.
	if( setjmp( jerr.setjmp_buffer ) )
	{
		// If we get here, the JPEG code has signaled an error.
		// We need to clean up the JPEG object, close the input file, and return.
		jpeg_destroy_decompress( &cinfo );
		fclose( in_fp );
		return( 0 );
	}
	// Now we can initialize the JPEG decompression object.
	jpeg_create_decompress( &cinfo );

	/* Step 2: specify data source (eg, a file) */

	jpeg_stdio_src( &cinfo, in_fp );

	/* Step 3: read file parameters with jpeg_read_header() */

	// We can ignore the return value from jpeg_read_header since
	// (a) suspension is not possible with the stdio data source, and
	// (b) we passed TRUE to reject a tables-only JPEG file as
	// an error.
	// See libjpeg.doc for more info.
	(void) jpeg_read_header( &cinfo, TRUE );

	/* Step 4: set parameters for decompression */

	// In this example, we don't need to change any of the defaults set by
	// jpeg_read_header(), so we do nothing here.

	/* Step 5: Start decompressor */

	(void) jpeg_start_decompress( &cinfo );
	// We can ignore the return value since suspension is not possible
	// with the stdio data source.

	// We may need to do some setup of our own at this point before reading
	// the data.  After jpeg_start_decompress() we have the correct scaled
	// output image dimensions available, as well as the output colormap
	// if we asked for color quantization.
	// In this example, we need to make an output work buffer of the right size.
	// JSAMPLEs per row in output buffer
	row_stride = cinfo.output_width * cinfo.output_components;
	// Make a one-row-high sample array that will go away when done with image
	buffer = (*cinfo.mem->alloc_sarray)
			( (j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1 );

	/* Step 6: while (scan lines remain to be read) */
	/*           jpeg_read_scanlines(...); */

	// allocate image (gray-scale or color)
	if( cinfo.output_components == 1 )
	{
		img._img_type = Image::Gray;
	}
	else if( cinfo.output_components == 3 )
	{
		img._img_type = Image::Color;
	}
	else
	{
		err_quit( "Jpeg_IO::read:Unknown image type\n" );
	}
	img._max_value = Image::DEFAULT_MAX;
	img.allocate_image( cinfo.output_height, cinfo.output_width );

	Image::size_type i = 0;
	if( img.is_bilevel() || img.is_gray() )
	{
		// Here we use the library's state variable cinfo.output_scanline as the
		// loop counter, so that we don't have to keep track ourselves.
		while( cinfo.output_scanline < cinfo.output_height )
		{
			// jpeg_read_scanlines expects an array of pointers to scanlines.
			// Here the array is only one element long, but you could ask for
			// more than one scanline at a time if that's more convenient.
			(void) jpeg_read_scanlines( &cinfo, buffer, 1 );

			// copy buffer into row of image
			unsigned buf_j = 0;
			for( Image::size_type j = img.col_begin(); j != img.col_end(); ++j )
			{
				img(i,j) = buffer[0][buf_j++];
			}
			++i;
		}
	}
	else // img.is_color() == true
	{
		// Here we use the library's state variable cinfo.output_scanline as the
		// loop counter, so that we don't have to keep track ourselves.
		while( cinfo.output_scanline < cinfo.output_height )
		{
			// jpeg_read_scanlines expects an array of pointers to scanlines.
			// Here the array is only one element long, but you could ask for
			// more than one scanline at a time if that's more convenient.
			(void) jpeg_read_scanlines( &cinfo, buffer, 1 );

			// copy buffer into row of image
			unsigned buf_j = 0;
			for( Image::size_type j = img.col_begin(); j != img.col_end(); ++j )
			{
				img.red(i,j)   = buffer[0][buf_j++];
				img.green(i,j) = buffer[0][buf_j++];
				img.blue(i,j)  = buffer[0][buf_j++];
			}
			++i;
		}
	}

	/* Step 7: Finish decompression */

	(void) jpeg_finish_decompress( &cinfo );
	// We can ignore the return value since suspension is not possible
	// with the stdio data source.

	/* Step 8: Release JPEG decompression object */

	// This is an important step since it will release a good deal of memory.
	jpeg_destroy_decompress( &cinfo );

	// After finish_decompress, we can close the input file.
	// Here we postpone it until after no more JPEG errors are possible,
	// so as to simplify the setjmp error logic above.  (Actually, I don't
	// think that jpeg_destroy can do an error exit, but why assume anything...)
	fclose( in_fp );

	// At this point you may want to check to see whether any corrupt-data
	// warnings occurred (test whether jerr.pub.num_warnings is nonzero).

	// successfully finished
	return( 1 );

#else
	err_warn( "Unable to read image '%s': libjpeg (http://www.ijg.org/)"
			" is not installed.\n", file_name.c_str() );
	return( 0 );
#endif
}
