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
nchoosek(int , int);

int
num_weyl(int M, int N, int ms);

int
get_next_tableau(const int M, const std::vector<int> frame, std::vector<int>& tableau);

int 
get_init_tableau(const int M, const std::vector<int> frame, std::vector<int>&tableau);

int
tableau_pos(int row, int col, std::vector<int> frame_rows);

bool
row_col_to_pos(int& pos, const std::vector<int> frame, const int row, const int col);

void
pos_to_row_col(const int pos,const std::vector<int> frame, int& row, int& col);

void
print_tableau(const std::vector<int> tableau, const std::vector<int> frame);

int
tableau_pos(int row, int col, std::vector<int> frame_rows);

#endif
