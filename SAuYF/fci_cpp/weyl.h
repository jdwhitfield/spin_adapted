/********************************************************************************
 * 
 * weyl.h, weyl.cpp
 *
 * These files are for generating the weyl tableaux systematically.
 *
 * JDWhitfield
 * Dartmouth 2017
 *
 ********************************************************************************/
//start of header guard
#ifndef WEYL
#define WEYL


int
num_weyl(int M, int N, int ms);

int
get_next_tableau(const int M, const std::vector<int> frame, std::vector<int>& tableau);

int 
get_init_tableau(const int M, const std::vector<int> frame, std::vector<int>&tableau);

#endif
