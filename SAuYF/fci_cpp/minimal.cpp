#include<stdio.h>
#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<assert.h>
#include<iomanip>
#include<random>
#include<chrono>
#include<algorithm>
#include"parser.h"

int 
main()
{

    const char* basis_fname="basis_funcs";

    libint2::BasisSet shells = parse_basisfile(basis_fname);

    for(auto s: shells)
	    std::cout << s << "\n";

    return 0;
}
