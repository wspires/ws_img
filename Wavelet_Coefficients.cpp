/**
	@file   Wavelet_Coefficients.cpp
	@author Wade Spires
	@date   2006/2/1
	@brief  Coefficients for popular wavelets.

	Copyright 2006 Wade Spires.
	Distributed under the GNU Lesser General Public License, Version 2.1.
	(See accompanying file LICENSE.txt or copy at
	http://www.gnu.org/licenses/lgpl.txt)
 */

#include "Wavelet_Coefficients.hpp"

using namespace ws_img;

/*
	Coefficients for Haar wavelets.
 */
const double Haar_Coefficients::ch_2[2] =
{
	SQRT_1_OVER_2,
	SQRT_1_OVER_2
};

const double Haar_Coefficients::cg_2[2] =
{
	SQRT_1_OVER_2,
	-(SQRT_1_OVER_2)
};

/*
	Coefficients for Daubechies wavelets of extremal phase are from
	I. Daubechies, "Orthonormal Bases of Compactly Supported Wavelets",
	Communications on Pure and Applied Mathematics, 41 (1988) 909--996
	(table 1).
	Additional digits have been obtained using the Mathematica package
	Daubechies.m by Tong Chen & Meng Xu available at
	http://www.cwp.mines.edu/wavelets/.
 */
const double Daub_Coefficients::h_4[4] =
{
	0.48296291314453414337487159986,
	0.83651630373780790557529378092,
	0.22414386804201338102597276224,
	-0.12940952255126038117444941881
};

const double Daub_Coefficients::g_4[4] =
{
	-0.12940952255126038117444941881,
	-0.22414386804201338102597276224,
	0.83651630373780790557529378092,
	-0.48296291314453414337487159986
};

const double Daub_Coefficients::h_6[6] =
{
	0.33267055295008261599851158914,
	0.80689150931109257649449360409,
	0.45987750211849157009515194215,
	-0.13501102001025458869638990670,
	-0.08544127388202666169281916918,
	0.03522629188570953660274066472
};

const double Daub_Coefficients::g_6[6] =
{
	0.03522629188570953660274066472,
	0.08544127388202666169281916918,
	-0.13501102001025458869638990670,
	-0.45987750211849157009515194215,
	0.80689150931109257649449360409,
	-0.33267055295008261599851158914
};

const double Daub_Coefficients::h_8[8] =
{
	0.23037781330889650086329118304,
	0.71484657055291564708992195527,
	0.63088076792985890788171633830,
	-0.02798376941685985421141374718,
	-0.18703481171909308407957067279,
	0.03084138183556076362721936253,
	0.03288301166688519973540751355,
	-0.01059740178506903210488320852
};

const double Daub_Coefficients::g_8[8] =
{
	-0.01059740178506903210488320852,
	-0.03288301166688519973540751355,
	0.03084138183556076362721936253,
	0.18703481171909308407957067279,
	-0.02798376941685985421141374718,
	-0.63088076792985890788171633830,
	0.71484657055291564708992195527,
	-0.23037781330889650086329118304
};

const double Daub_Coefficients::h_10[10] =
{
	0.16010239797419291448072374802,
	0.60382926979718967054011930653,
	0.72430852843777292772807124410,
	0.13842814590132073150539714634,
	-0.24229488706638203186257137947,
	-0.03224486958463837464847975506,
	0.07757149384004571352313048939,
	-0.00624149021279827427419051911,
	-0.01258075199908199946850973993,
	0.00333572528547377127799818342
};

const double Daub_Coefficients::g_10[10] =
{
	0.00333572528547377127799818342,
	0.01258075199908199946850973993,
	-0.00624149021279827427419051911,
	-0.07757149384004571352313048939,
	-0.03224486958463837464847975506,
	0.24229488706638203186257137947,
	0.13842814590132073150539714634,
	-0.72430852843777292772807124410,
	0.60382926979718967054011930653,
	-0.16010239797419291448072374802
};

const double Daub_Coefficients::h_12[12] =
{
	0.11154074335010946362132391724,
	0.49462389039845308567720417688,
	0.75113390802109535067893449844,
	0.31525035170919762908598965481,
	-0.22626469396543982007631450066,
	-0.12976686756726193556228960588,
	0.09750160558732304910234355254,
	0.02752286553030572862554083950,
	-0.03158203931748602956507908070,
	0.00055384220116149613925191840,
	0.00477725751094551063963597525,
	-0.00107730108530847956485262161
};

const double Daub_Coefficients::g_12[12] =
{
	-0.00107730108530847956485262161,
	-0.00477725751094551063963597525,
	0.00055384220116149613925191840,
	0.03158203931748602956507908070,
	0.02752286553030572862554083950,
	-0.09750160558732304910234355254,
	-0.12976686756726193556228960588,
	0.22626469396543982007631450066,
	0.31525035170919762908598965481,
	-0.75113390802109535067893449844,
	0.49462389039845308567720417688,
	-0.11154074335010946362132391724
};

const double Daub_Coefficients::h_14[14] =
{
	0.07785205408500917901996352196,
	0.39653931948191730653900039094,
	0.72913209084623511991694307034,
	0.46978228740519312247159116097,
	-0.14390600392856497540506836221,
	-0.22403618499387498263814042023,
	0.07130921926683026475087657050,
	0.08061260915108307191292248036,
	-0.03802993693501441357959206160,
	-0.01657454163066688065410767489,
	0.01255099855609984061298988603,
	0.00042957797292136652113212912,
	-0.00180164070404749091526826291,
	0.00035371379997452024844629584
};

const double Daub_Coefficients::g_14[14] =
{
	0.00035371379997452024844629584,
	0.00180164070404749091526826291,
	0.00042957797292136652113212912,
	-0.01255099855609984061298988603,
	-0.01657454163066688065410767489,
	0.03802993693501441357959206160,
	0.08061260915108307191292248036,
	-0.07130921926683026475087657050,
	-0.22403618499387498263814042023,
	0.14390600392856497540506836221,
	0.46978228740519312247159116097,
	-0.72913209084623511991694307034,
	0.39653931948191730653900039094,
	-0.07785205408500917901996352196
};

const double Daub_Coefficients::h_16[16] =
{
	0.05441584224310400995500940520,
	0.31287159091429997065916237551,
	0.67563073629728980680780076705,
	0.58535468365420671277126552005,
	-0.01582910525634930566738054788,
	-0.28401554296154692651620313237,
	0.00047248457391328277036059001,
	0.12874742662047845885702928751,
	-0.01736930100180754616961614887,
	-0.04408825393079475150676372324,
	0.01398102791739828164872293057,
	0.00874609404740577671638274325,
	-0.00487035299345157431042218156,
	-0.00039174037337694704629808036,
	0.00067544940645056936636954757,
	-0.00011747678412476953373062823
};

const double Daub_Coefficients::g_16[16] =
{
	-0.00011747678412476953373062823,
	-0.00067544940645056936636954757,
	-0.00039174037337694704629808036,
	0.00487035299345157431042218156,
	0.00874609404740577671638274325,
	-0.01398102791739828164872293057,
	-0.04408825393079475150676372324,
	0.01736930100180754616961614887,
	0.12874742662047845885702928751,
	-0.00047248457391328277036059001,
	-0.28401554296154692651620313237,
	0.01582910525634930566738054788,
	0.58535468365420671277126552005,
	-0.67563073629728980680780076705,
	0.31287159091429997065916237551,
	-0.05441584224310400995500940520
};

const double Daub_Coefficients::h_18[18] =
{
	0.03807794736387834658869765888,
	0.24383467461259035373204158165,
	0.60482312369011111190307686743,
	0.65728807805130053807821263905,
	0.13319738582500757619095494590,
	-0.29327378327917490880640319524,
	-0.09684078322297646051350813354,
	0.14854074933810638013507271751,
	0.03072568147933337921231740072,
	-0.06763282906132997367564227483,
	0.00025094711483145195758718975,
	0.02236166212367909720537378270,
	-0.00472320475775139727792570785,
	-0.00428150368246342983449679500,
	0.00184764688305622647661912949,
	0.00023038576352319596720521639,
	-0.00025196318894271013697498868,
	0.00003934732031627159948068988
};

const double Daub_Coefficients::g_18[18] =
{
	0.00003934732031627159948068988,
	0.00025196318894271013697498868,
	0.00023038576352319596720521639,
	-0.00184764688305622647661912949,
	-0.00428150368246342983449679500,
	0.00472320475775139727792570785,
	0.02236166212367909720537378270,
	-0.00025094711483145195758718975,
	-0.06763282906132997367564227483,
	-0.03072568147933337921231740072,
	0.14854074933810638013507271751,
	0.09684078322297646051350813354,
	-0.29327378327917490880640319524,
	-0.13319738582500757619095494590,
	0.65728807805130053807821263905,
	-0.60482312369011111190307686743,
	0.24383467461259035373204158165,
	-0.03807794736387834658869765888
};

const double Daub_Coefficients::h_20[20] =
{
	0.02667005790055555358661744877,
	0.18817680007769148902089297368,
	0.52720118893172558648174482796,
	0.68845903945360356574187178255,
	0.28117234366057746074872699845,
	-0.24984642432731537941610189792,
	-0.19594627437737704350429925432,
	0.12736934033579326008267723320,
	0.09305736460357235116035228984,
	-0.07139414716639708714533609308,
	-0.02945753682187581285828323760,
	0.03321267405934100173976365318,
	0.00360655356695616965542329142,
	-0.01073317548333057504431811411,
	0.00139535174705290116578931845,
	0.00199240529518505611715874224,
	-0.00068585669495971162656137098,
	-0.00011646685512928545095148097,
	0.00009358867032006959133405013,
	-0.00001326420289452124481243668
};

const double Daub_Coefficients::g_20[20] =
{
	-0.00001326420289452124481243668,
	-0.00009358867032006959133405013,
	-0.00011646685512928545095148097,
	0.00068585669495971162656137098,
	0.00199240529518505611715874224,
	-0.00139535174705290116578931845,
	-0.01073317548333057504431811411,
	-0.00360655356695616965542329142,
	0.03321267405934100173976365318,
	0.02945753682187581285828323760,
	-0.07139414716639708714533609308,
	-0.09305736460357235116035228984,
	0.12736934033579326008267723320,
	0.19594627437737704350429925432,
	-0.24984642432731537941610189792,
	-0.28117234366057746074872699845,
	0.68845903945360356574187178255,
	-0.52720118893172558648174482796,
	0.18817680007769148902089297368,
	-0.02667005790055555358661744877
};

/*
	Coefficients are from A. Cohen, I. Daubechies, and J.-C. Feauveau;
	"Biorthogonal Bases of Compactly Supported Wavelets", Communications
	on Pure and Applied Mathematics, 45 (1992) 485--560 (table 6.1).

	Note the following errors in table 1:

	N = 2, N~ = 4, m0~
	the second term in z^-1 (45/64 z^-1) should be left out.

	N = 3, N~ = 7, m0~
	the term 336z^-3 should read 363z^-3.
 */
const double Bspline_Coefficients::h1_103[6] =
{
	-0.0883883476483184405501055452631,
	0.0883883476483184405501055452631,
	SQRT_1_OVER_2,
	SQRT_1_OVER_2,
	0.0883883476483184405501055452631,
	-0.0883883476483184405501055452631
};

const double Bspline_Coefficients::g2_103[6] =
{
	-0.0883883476483184405501055452631,
	-0.0883883476483184405501055452631,
	SQRT_1_OVER_2,
	-(SQRT_1_OVER_2),
	0.0883883476483184405501055452631,
	0.0883883476483184405501055452631
};

const double Bspline_Coefficients::h1_105[10] =
{
	0.0165728151840597076031447897368,
	-0.0165728151840597076031447897368,
	-0.1215339780164378557563951247368,
	0.1215339780164378557563951247368,
	SQRT_1_OVER_2,
	SQRT_1_OVER_2,
	0.1215339780164378557563951247368,
	-0.1215339780164378557563951247368,
	-0.0165728151840597076031447897368,
	0.0165728151840597076031447897368
};

const double Bspline_Coefficients::g2_105[10] =
{
	0.0165728151840597076031447897368,
	0.0165728151840597076031447897368,
	-0.1215339780164378557563951247368,
	-0.1215339780164378557563951247368,
	SQRT_1_OVER_2,
	-(SQRT_1_OVER_2),
	0.1215339780164378557563951247368,
	0.1215339780164378557563951247368,
	-0.0165728151840597076031447897368,
	-0.0165728151840597076031447897368
};

const double Bspline_Coefficients::g1_1[10] =
{
	0.0, 0.0, 0.0, 0.0,
	SQRT_1_OVER_2,
	-(SQRT_1_OVER_2),
	0.0, 0.0, 0.0, 0.0
};

const double Bspline_Coefficients::h2_1[10] =
{
	0.0, 0.0, 0.0, 0.0,
	SQRT_1_OVER_2,
	SQRT_1_OVER_2,
	0.0, 0.0, 0.0, 0.0
};

const double Bspline_Coefficients::h1_202[6] =
{
	-0.1767766952966368811002110905262,
	0.3535533905932737622004221810524,
	1.0606601717798212866012665431573,
	0.3535533905932737622004221810524,
	-0.1767766952966368811002110905262,
	0.0
};

const double Bspline_Coefficients::g2_202[6] =
{
	0.0,
	-0.1767766952966368811002110905262,
	-0.3535533905932737622004221810524,
	1.0606601717798212866012665431573,
	-0.3535533905932737622004221810524,
	-0.1767766952966368811002110905262
};

const double Bspline_Coefficients::h1_204[10] =
{
	0.0331456303681194152062895794737,
	-0.0662912607362388304125791589473,
	-0.1767766952966368811002110905262,
	0.4198446513295125926130013399998,
	0.9943689110435824561886873842099,
	0.4198446513295125926130013399998,
	-0.1767766952966368811002110905262,
	-0.0662912607362388304125791589473,
	0.0331456303681194152062895794737,
	0.0
};

const double Bspline_Coefficients::g2_204[10] =
{
	0.0,
	0.0331456303681194152062895794737,
	0.0662912607362388304125791589473,
	-0.1767766952966368811002110905262,
	-0.4198446513295125926130013399998,
	0.9943689110435824561886873842099,
	-0.4198446513295125926130013399998,
	-0.1767766952966368811002110905262,
	0.0662912607362388304125791589473,
	0.0331456303681194152062895794737
};

const double Bspline_Coefficients::h1_206[14] =
{
	-0.0069053396600248781679769957237,
	0.0138106793200497563359539914474,
	0.0469563096881691715422435709210,
	-0.1077232986963880994204411332894,
	-0.1698713556366120029322340948025,
	0.4474660099696121052849093228945,
	0.9667475524034829435167794013152,
	0.4474660099696121052849093228945,
	-0.1698713556366120029322340948025,
	-0.1077232986963880994204411332894,
	0.0469563096881691715422435709210,
	0.0138106793200497563359539914474,
	-0.0069053396600248781679769957237,
	0.0
};

const double Bspline_Coefficients::g2_206[14] =
{
	0.0,
	-0.0069053396600248781679769957237,
	-0.0138106793200497563359539914474,
	0.0469563096881691715422435709210,
	0.1077232986963880994204411332894,
	-0.1698713556366120029322340948025,
	-0.4474660099696121052849093228945,
	0.9667475524034829435167794013152,
	-0.4474660099696121052849093228945,
	-0.1698713556366120029322340948025,
	0.1077232986963880994204411332894,
	0.0469563096881691715422435709210,
	-0.0138106793200497563359539914474,
	-0.0069053396600248781679769957237,
};

const double Bspline_Coefficients::h1_208[18] =
{
	0.0015105430506304420992449678146,
	-0.0030210861012608841984899356291,
	-0.0129475118625466465649568669819,
	0.0289161098263541773284036695929,
	0.0529984818906909399392234421792,
	-0.1349130736077360572068505539514,
	-0.1638291834340902345352542235443,
	0.4625714404759165262773590010400,
	0.9516421218971785225243297231697,
	0.4625714404759165262773590010400,
	-0.1638291834340902345352542235443,
	-0.1349130736077360572068505539514,
	0.0529984818906909399392234421792,
	0.0289161098263541773284036695929,
	-0.0129475118625466465649568669819,
	-0.0030210861012608841984899356291,
	0.0015105430506304420992449678146,
	0.0
};

const double Bspline_Coefficients::g2_208[18] =
{
	0.0,
	0.0015105430506304420992449678146,
	0.0030210861012608841984899356291,
	-0.0129475118625466465649568669819,
	-0.0289161098263541773284036695929,
	0.0529984818906909399392234421792,
	0.1349130736077360572068505539514,
	-0.1638291834340902345352542235443,
	-0.4625714404759165262773590010400,
	0.9516421218971785225243297231697,
	-0.4625714404759165262773590010400,
	-0.1638291834340902345352542235443,
	0.1349130736077360572068505539514,
	0.0529984818906909399392234421792,
	-0.0289161098263541773284036695929,
	-0.0129475118625466465649568669819,
	0.0030210861012608841984899356291,
	0.0015105430506304420992449678146,
};

const double Bspline_Coefficients::h2_2[18] =
{
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.3535533905932737622004221810524,
	0.7071067811865475244008443621048,
	0.3535533905932737622004221810524,
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

const double Bspline_Coefficients::g1_2[18] =
{
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	-0.3535533905932737622004221810524,
	0.7071067811865475244008443621048,
	-0.3535533905932737622004221810524,
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

const double Bspline_Coefficients::h1_301[4] =
{
	-0.3535533905932737622004221810524,
	1.0606601717798212866012665431573,
	1.0606601717798212866012665431573,
	-0.3535533905932737622004221810524
};

const double Bspline_Coefficients::g2_301[4] =
{
	0.3535533905932737622004221810524,
	1.0606601717798212866012665431573,
	-1.0606601717798212866012665431573,
	-0.3535533905932737622004221810524
};

const double Bspline_Coefficients::h1_303[8] =
{
	0.0662912607362388304125791589473,
	-0.1988737822087164912377374768420,
	-0.1546796083845572709626847042104,
	0.9943689110435824561886873842099,
	0.9943689110435824561886873842099,
	-0.1546796083845572709626847042104,
	-0.1988737822087164912377374768420,
	0.0662912607362388304125791589473
};

const double Bspline_Coefficients::g2_303[8] =
{
	-0.0662912607362388304125791589473,
	-0.1988737822087164912377374768420,
	0.1546796083845572709626847042104,
	0.9943689110435824561886873842099,
	-0.9943689110435824561886873842099,
	-0.1546796083845572709626847042104,
	0.1988737822087164912377374768420,
	0.0662912607362388304125791589473
};

const double Bspline_Coefficients::h1_305[12] =
{
	-0.0138106793200497563359539914474,
	0.0414320379601492690078619743421,
	0.0524805814161890740766251675000,
	-0.2679271788089652729175074340788,
	-0.0718155324642587329469607555263,
	0.9667475524034829435167794013152,
	0.9667475524034829435167794013152,
	-0.0718155324642587329469607555263,
	-0.2679271788089652729175074340788,
	0.0524805814161890740766251675000,
	0.0414320379601492690078619743421,
	-0.0138106793200497563359539914474
};

const double Bspline_Coefficients::g2_305[12] =
{
	0.0138106793200497563359539914474,
	0.0414320379601492690078619743421,
	-0.0524805814161890740766251675000,
	-0.2679271788089652729175074340788,
	0.0718155324642587329469607555263,
	0.9667475524034829435167794013152,
	-0.9667475524034829435167794013152,
	-0.0718155324642587329469607555263,
	0.2679271788089652729175074340788,
	0.0524805814161890740766251675000,
	-0.0414320379601492690078619743421,
	-0.0138106793200497563359539914474
};

const double Bspline_Coefficients::h1_307[16] =
{
	0.0030210861012608841984899356291,
	-0.0090632583037826525954698068873,
	-0.0168317654213106405344439270765,
	0.0746639850740189951912512662623,
	0.0313329787073628846871956180962,
	-0.3011591259228349991008967259990,
	-0.0264992409453454699696117210896,
	0.9516421218971785225243297231697,
	0.9516421218971785225243297231697,
	-0.0264992409453454699696117210896,
	-0.3011591259228349991008967259990,
	0.0313329787073628846871956180962,
	0.0746639850740189951912512662623,
	-0.0168317654213106405344439270765,
	-0.0090632583037826525954698068873,
	0.0030210861012608841984899356291
};

const double Bspline_Coefficients::g2_307[16] =
{
	-0.0030210861012608841984899356291,
	-0.0090632583037826525954698068873,
	0.0168317654213106405344439270765,
	0.0746639850740189951912512662623,
	-0.0313329787073628846871956180962,
	-0.3011591259228349991008967259990,
	0.0264992409453454699696117210896,
	0.9516421218971785225243297231697,
	-0.9516421218971785225243297231697,
	-0.0264992409453454699696117210896,
	0.3011591259228349991008967259990,
	0.0313329787073628846871956180962,
	-0.0746639850740189951912512662623,
	-0.0168317654213106405344439270765,
	0.0090632583037826525954698068873,
	0.0030210861012608841984899356291
};

const double Bspline_Coefficients::h1_309[20] =
{
	-0.0006797443727836989446602355165,
	0.0020392331183510968339807065496,
	0.0050603192196119810324706421788,
	-0.0206189126411055346546938106687,
	-0.0141127879301758447558029850103,
	0.0991347824942321571990197448581,
	0.0123001362694193142367090236328,
	-0.3201919683607785695513833204624,
	0.0020500227115698857061181706055,
	0.9421257006782067372990864259380,
	0.9421257006782067372990864259380,
	0.0020500227115698857061181706055,
	-0.3201919683607785695513833204624,
	0.0123001362694193142367090236328,
	0.0991347824942321571990197448581,
	-0.0141127879301758447558029850103,
	-0.0206189126411055346546938106687,
	0.0050603192196119810324706421788,
	0.0020392331183510968339807065496,
	-0.0006797443727836989446602355165
};

const double Bspline_Coefficients::g2_309[20] =
{
	0.0006797443727836989446602355165,
	0.0020392331183510968339807065496,
	-0.0050603192196119810324706421788,
	-0.0206189126411055346546938106687,
	0.0141127879301758447558029850103,
	0.0991347824942321571990197448581,
	-0.0123001362694193142367090236328,
	-0.3201919683607785695513833204624,
	-0.0020500227115698857061181706055,
	0.9421257006782067372990864259380,
	-0.9421257006782067372990864259380,
	0.0020500227115698857061181706055,
	0.3201919683607785695513833204624,
	0.0123001362694193142367090236328,
	-0.0991347824942321571990197448581,
	-0.0141127879301758447558029850103,
	0.0206189126411055346546938106687,
	0.0050603192196119810324706421788,
	-0.0020392331183510968339807065496,
	-0.0006797443727836989446602355165
};

const double Bspline_Coefficients::h2_3[20] =
{
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	0.1767766952966368811002110905262,
	0.5303300858899106433006332715786,
	0.5303300858899106433006332715786,
	0.1767766952966368811002110905262,
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};

const double Bspline_Coefficients::g1_3[20] =
{
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
	-0.1767766952966368811002110905262,
	0.5303300858899106433006332715786,
	-0.5303300858899106433006332715786,
	0.1767766952966368811002110905262,
	0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
};