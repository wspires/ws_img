/**
	@file   Video_IO.cpp
	@author Wade Spires
	@date   2006/3/24
	@brief  Class Video_IO.

	This code is based on the FFmpeg tutorial written by Martin Böhme at
	http://www.inb.uni-luebeck.de/~boehme/using_libavcodec.html
	with changes to make object-oriented and interface with my image library.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Video_IO.hpp"

using std::string;
using std::vector;

using namespace ws_img;

typedef Image::size_type  size_type;
typedef Image::value_type value_type;

/**
	Default constructor.
 */
Video_IO::Video_IO( )
: _format_ctx(0), _codec_ctx(0), _codec(0), _frame(0), _converted_frame(0),
_buffer(0), _frame_num(0)
{ }

/**
	Construct decoder using the given video file.
	@param[in] file_name Video file to decode
	@param[in] image_type Image type--either Gray or Color for gray-scale or
		color images
 */
Video_IO::Video_IO( const std::string& file_name, Image::Image_Type image_type )
: _format_ctx(0), _codec_ctx(0), _codec(0), _frame(0), _converted_frame(0),
_buffer(0), _frame_num(0)
{
	// register all formats and codecs
	av_register_all();

	// open video file
	if( av_open_input_file( &_format_ctx, file_name.c_str(), NULL, 0, NULL)
			!= 0 )
	{
		err_quit( "Unable to open video file '%s'\n", file_name.c_str() );
	}

	// retrieve stream information
	if( av_find_stream_info( _format_ctx ) < 0 )
	{
		err_quit( "Unable to find stream information for video file '%s'\n",
				file_name.c_str() );
	}

	// find the first video stream
	_video_stream = -1;
	for( int i = 0; i < _format_ctx->nb_streams; ++i )
	{
		if( _format_ctx->streams[i]->codec->codec_type == CODEC_TYPE_VIDEO )
		{
			_video_stream = i;
			break;
		}
	}
	if( _video_stream == -1 )  // if no video stream was found
	{
		err_quit( "Unable to find video stream for video file '%s'\n",
				file_name.c_str() );
	}

	// get a pointer to the codec context for the video stream
	_codec_ctx = _format_ctx->streams[_video_stream]->codec;

	// find the decoder for the video stream
	if( (_codec = avcodec_find_decoder( _codec_ctx->codec_id )) == NULL )
	{
		err_quit( "Unable to find codec for video file '%s'\n",
				file_name.c_str() );
	}

	// open codec
	if( avcodec_open(_codec_ctx, _codec) < 0 )
	{
		err_quit( "Unable to open codec for video file '%s'\n",
				file_name.c_str() );
	}

	// allocate video frame
	_frame = avcodec_alloc_frame();

	// allocate an AVFrame structure
	if( (_converted_frame = avcodec_alloc_frame()) == NULL )
	{
		err_quit( "Unable to allocate frame structure for video file '%s'\n",
				file_name.c_str() );
	}

	if( image_type == Image::Gray )
	{
		_pixel_format = PIX_FMT_GRAY8;
	}
	else if( image_type == Image::Color )
	{
		_pixel_format = PIX_FMT_RGB24;
	}

	// determine required buffer size and allocate buffer
	unsigned num_bytes = avpicture_get_size( _pixel_format, _codec_ctx->width,
			_codec_ctx->height );
	_buffer = new uint8_t[num_bytes];

	// assign appropriate parts of buffer to image planes in _converted_frame
	avpicture_fill( (AVPicture*) _converted_frame, _buffer, _pixel_format,
			_codec_ctx->width, _codec_ctx->height );
}

/**
	Destructor closes files and deallocates buffers and codecs.
 */
Video_IO::~Video_IO( )
{
	cleanup();
}

/**
	Get next image frame.
	@param[out] img Next image frame from video
	@retval finished Whether any more frames are left to extract
 */
bool
Video_IO::next_frame( Image& img )
{
	AVPacket packet;
	while( av_read_frame( _format_ctx, &packet ) >= 0 )
	{
		// if this is a packet from the video stream
		if( packet.stream_index == _video_stream )
		{
			// decode video frame
			int is_frame_finished;
			avcodec_decode_video( _codec_ctx, _frame, &is_frame_finished,
					packet.data, packet.size );

			// if we got a complete video frame
			if( is_frame_finished )
			{
				++_frame_num;

				// convert the image from its native format to RGB (or gray-scale)
				img_convert( (AVPicture*) _converted_frame, _pixel_format, 
						(AVPicture*) _frame, _codec_ctx->pix_fmt, _codec_ctx->width, 
						_codec_ctx->height );

				// create image from frame
				if( _pixel_format == PIX_FMT_GRAY8 )
				{
					make_gray_frame( _converted_frame, _codec_ctx->height,
							_codec_ctx->width, img );
				}
				else if( _pixel_format == PIX_FMT_RGB24 )
				{
					make_color_frame( _converted_frame, _codec_ctx->height,
							_codec_ctx->width, img );
				}
				return( true );
			}
		}

		// free the packet that was allocated by av_read_frame
		av_free_packet( &packet );
	}

	// if last frame was processed, stop
	cleanup();
	return( false );
}

/**
	Create a gray-scale image from the given frame.
	@param[in] frame
	@param[in] height
	@param[in] width
	@param[out] img
 */
void
Video_IO::make_gray_frame( const AVFrame* frame, int height, int width,
		Image& img ) const
{
	size_type row_size = static_cast<size_type>( height );
	size_type col_size = static_cast<size_type>( width );

	// ensure the frame has the correct size and is in gray-scale
	if( img.sub_row_size() != row_size || img.sub_col_size() != col_size )
	{
		img.resize( row_size, col_size, false );
	}
	img.to_gray();

	// copy image data from the buffer into the image
	unsigned k = 0;
	for( size_type i = img.row_begin(); i != img.row_end(); ++i )
	{
		for( size_type j = img.col_begin(); j != img.col_end(); ++j )
		{
			img(i,j) = static_cast<value_type>( frame->data[0][ k++ ] );
		}
	}
}

/**
	Create a color image from the given frame.
	@param[in] frame
	@param[in] height
	@param[in] width
	@param[out] img
 */
void
Video_IO::make_color_frame( const AVFrame* frame, int height, int width,
		Image& img ) const
{
	size_type row_size = static_cast<size_type>( height );
	size_type col_size = static_cast<size_type>( width );

	// ensure the frame has the correct size and is in color
	if( img.sub_row_size() != row_size || img.sub_col_size() != col_size )
	{
		img.resize( row_size, col_size, false );
	}
	img.to_color();

	// copy image data from the buffer into the image
	unsigned k = 0;
	for( size_type i = img.row_begin(); i != img.row_end(); ++i )
	{
		for( size_type j = img.col_begin(); j != img.col_end(); ++j )
		{
			img.red(i,j)   = static_cast<value_type>( frame->data[0][ k++ ] );
			img.green(i,j) = static_cast<value_type>( frame->data[0][ k++ ] );
			img.blue(i,j)  = static_cast<value_type>( frame->data[0][ k++ ] );
		}
	}
}

/**
	Close files and deallocate buffers and codecs.
 */
void
Video_IO::cleanup( )
{
	// free the RGB image
	if( _buffer != 0 )
	{
		delete [] _buffer;
	}
	_buffer = 0;
	if( _converted_frame != 0 )
	{
		av_free(_converted_frame);
	}
	_converted_frame = 0;

	// free the YUV frame
	if( _frame != 0 )
	{
		av_free(_frame);
	}
	_frame = 0;

	// close the codec
	if( _codec_ctx != 0 )
	{
		avcodec_close(_codec_ctx);
	}
	_codec_ctx = 0;

	// close the video file
	if( _format_ctx != 0 )
	{
		av_close_input_file(_format_ctx);
	}
	_format_ctx = 0;
}
