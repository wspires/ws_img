/**
	@file   Keypoint.cpp
	@author Wade Spires
	@date   2006/1/30
	@brief  Class Keypoint.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Keypoint.hpp"

using std::string;
using std::vector;

using namespace ws_img;

/**
	Convert keypoint into a string representation.
	@retval str Keypoint converted to a string
 */
std::string
Keypoint::to_string( ) const
{
	return( ws_tools::format_string( "(%u,%u,%u,%u)",
			octave_pos, scale_pos, img_coord.row, img_coord.col ) );
}

/**
	Write keypoint descriptors to file.

	@param[in] keypoints List of keypoints
	@param[in] file_name File to write points to
 */
void
Keypoint_List::write( const string& file_name )
{
	FILE* file_fp;
	if( (file_fp = fopen( file_name.c_str(), "w" )) == NULL )
	{
		err_quit( "Unable to open file '%s'\n", file_name.c_str() );
	}

	// get number of keypoint descriptors
	unsigned num_descriptors = 0;
	if( size() > 0 )
	{
		num_descriptors = at(0).descriptors.size();
	}

	// write number of points and number of descriptors per point
	fprintf( file_fp, "%u %u\n", (unsigned) size(), num_descriptors );

	// write each point
	for( unsigned i = 0; i != size(); ++i )
	{
		const Keypoint& key = at(i);

		if( key.descriptors.size() != num_descriptors )
		{
			err_quit( "Number of descriptors do not match\n" );
		}

		// shift angle so that it is in the range [-PI, PI] like Lowe's angles
		double orientation = key.orientation;
		if( orientation > PI )
		{
			orientation -= 2 * PI;
		}

		fprintf( file_fp, "%.2lf %.2lf %.2lf %.3lf\n",
				static_cast<double>( key.img_coord.row ),
				static_cast<double>( key.img_coord.col ),
				key.scale, orientation );

		// write each descriptor
		unsigned k = 0;
		for( Vector::size_type j = 0; j != key.descriptors.size(); ++j )
		{
			// fprintf( file_fp, " %lf", key.descriptors(j) );
			fprintf( file_fp, " %1.lf", round(key.descriptors(j)) );
			++k;

			// ensure that only 20 descriptors are on each line
			if( k == 20 )
			{
				fprintf( file_fp, "\n" );
				k = 0;
			}
		}
		fprintf( file_fp, "\n" );
	}

	fclose( file_fp );
}
