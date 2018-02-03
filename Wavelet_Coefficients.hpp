/**
	@file   Wavelet_Coefficients.hpp
	@author Wade Spires
	@date   2006/2/1
	@brief  Coefficients for popular wavelets.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#ifndef WAVELET_COEFFICIENTS_HPP
#define WAVELET_COEFFICIENTS_HPP

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

// local headers

namespace ws_img
{

// sqrt( 1 / 2 ): used in some coefficients
static const double SQRT_1_OVER_2 = .70710678118654752440;

/**
	@brief Coefficients for Haar wavelets
 */
struct Haar_Coefficients
{
	static const double ch_2[2];
	static const double cg_2[2];
};

/**
	@brief Coefficients for Daubechies' wavelets
 */
struct Daub_Coefficients
{
	static const double h_4[4];
	static const double g_4[4];
	static const double h_6[6];
	static const double g_6[6];
	static const double h_8[8];
	static const double g_8[8];
	static const double h_10[10];
	static const double g_10[10];
	static const double h_12[12];
	static const double g_12[12];
	static const double h_14[14];
	static const double g_14[14];
	static const double h_16[16];
	static const double g_16[16];
	static const double h_18[18];
	static const double g_18[18];
	static const double h_20[20];
	static const double g_20[20];
};

/**
	@brief Coefficients for B-Spline wavelets
 */
struct Bspline_Coefficients
{
	static const double h1_103[6];
	static const double g2_103[6];
	static const double h1_105[10];
	static const double g2_105[10];
	static const double g1_1[10];
	static const double h2_1[10];
	static const double h1_202[6];
	static const double g2_202[6];
	static const double h1_204[10];
	static const double g2_204[10];
	static const double h1_206[14];
	static const double g2_206[14];
	static const double h1_208[18];
	static const double g2_208[18];
	static const double h2_2[18];
	static const double g1_2[18];
	static const double h1_301[4];
	static const double g2_301[4];
	static const double h1_303[8];
	static const double g2_303[8];
	static const double h1_305[12];
	static const double g2_305[12];
	static const double h1_307[16];
	static const double g2_307[16];
	static const double h1_309[20];
	static const double g2_309[20];
	static const double h2_3[20];
	static const double g1_3[20];
};

} // namespace ws_img

#endif // WAVELET_COEFFICIENTS_HPP 
