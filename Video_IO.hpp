/**
	@file   Video_IO.hpp
	@author Wade Spires
	@date   2006/3/24
	@brief  Class Video_IO.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef VIDEO_IO_HPP
#define VIDEO_IO_HPP

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
#include "Image.hpp"

// libav headers
#include <ffmpeg/avcodec.h>
#include <ffmpeg/avformat.h>

// local headers

namespace ws_img
{

/**
	@brief Video_IO
 */
class Video_IO
{
	
public:

   Video_IO( const std::string&,
			ws_img::Image::Image_Type = ws_img::Image::Gray );

private: // the user must open a video (don't allow a default or copying)
   Video_IO( );
   Video_IO( const Video_IO& );
   Video_IO& operator=( const Video_IO& );

public:
   virtual ~Video_IO( );

	bool next_frame( ws_img::Image& );

	/**
		Get number of frames decoded so far.
		@retval num_frames Number of frames decoded so far
	 */
	inline unsigned get_frame_num( ) const
	{ return( _frame_num ); }

protected:

	void make_color_frame( const AVFrame*, int, int, ws_img::Image& ) const;
	void make_gray_frame( const AVFrame*, int, int, ws_img::Image& ) const;
	void cleanup( );

private:

	AVFormatContext*  _format_ctx;       //< Format context
	int               _video_stream;     //< Video stream
	AVCodecContext*   _codec_ctx;        //< Codec context
	AVCodec*          _codec;            //< Codec
	AVFrame*          _frame;            //< Frame
	AVFrame*          _converted_frame;  //< Frame converted from native to RGB
	uint8_t*          _buffer;           //< Buffer
	PixelFormat       _pixel_format;     //< Pixel format
	unsigned          _frame_num;        //< Current number of frames

};

} // namespace ws_img

#endif // VIDEO_IO_HPP 
