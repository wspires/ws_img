/**
	@file   K_Means.hpp
	@author Wade Spires
	@date   2006/10/20
	@brief  Class K_Means.

	Copyright (c) 2006 Wade Spires. All rights reserved.
 */

#ifndef K_MEANS_HPP
#define K_MEANS_HPP

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

// image headers
#include "Image.hpp"
#include "Matrix.hpp"

namespace ws_img
{

/**
	@brief K_Means
 */
class K_Means
{
	
public:

	/**
		Default constructor.
   K_Means( )
	{ }
	 */

	/**
		Destructor does nothing since no member variables are dynamically
		allocated.
   virtual ~K_Means( )
	{ }
	 */

	static ws_img::Matrix
	k_means( const unsigned, const ws_img::Matrix& );

protected:

	static void
	select_initial_means( const ws_img::Matrix&, ws_img::Matrix& );

	static void
	classify_samples( const ws_img::Matrix&, const ws_img::Matrix&,
				std::vector<ws_img::Matrix::size_type>& );

	static void
	select_new_means( const ws_img::Matrix&,
				const std::vector<ws_img::Matrix::size_type>&, ws_img::Matrix& );

	static double
	compute_difference( const ws_img::Matrix&, const ws_img::Matrix& );

	static void
	matrix_to_image( const ws_img::Matrix&, const ws_img::Matrix::size_type,
				const ws_img::Matrix::size_type, const std::string& );
};

} // namespace ws_img

#endif // K_MEANS_HPP
