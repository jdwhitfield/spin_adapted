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

//
/// fac[k] = k!
static constexpr std::array<int64_t,21> fac = {{1L, 1L, 2L, 6L, 24L, 120L, 720L, 5040L, 40320L, 362880L, 3628800L, 39916800L,
                                                    479001600L, 6227020800L, 87178291200L, 1307674368000L, 20922789888000L,
                                                    355687428096000L, 6402373705728000L, 121645100408832000L,
                                                    2432902008176640000L}};
 



#endif
