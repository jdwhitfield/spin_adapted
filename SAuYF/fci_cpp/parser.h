/******************************************************************************
 * 
 * parser.cpp
 *
 * This file has functions to parse inputs from various file formats for use in
 * other parts of the code base e.g. ci_matrix.cpp
 *
 * JDWhitfield 2017
 *
 *****************************************************************************/

//
// //start of header guard
#ifndef PARSER_LIB
#define PARSER_LIB

// Eigen matrix algebra library
#include <Eigen/Dense>
// Libint Gaussian integrals library
#include <libint2.hpp>


// any source file that includes this header will see this, and only this, var
// but it must be defined in one and only one of the files that's being complied 
extern const char *PeriodicTable[18];


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Matrix;  // import dense, dynamically sized Matrix type from Eigen;
                 // this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
                 // to meet the layout of the integrals returned by the Libint integral library
                 

//parse function declarations
libint2::BasisSet parse_basisfile(const char *);
std::vector<libint2::Atom>  	    parse_nucfile(const char* nuc_fname, int &);

//data for getting basis functions
const int NMAX=100;
const int ZMAX=18;
const char sto3g[]="!  STO-3G  EMSL  Basis Set Exchange Library   7/2/15 6:25 AM\n! Elements                             References\n! --------                             ----------\n!  H - Ne: W.J. Hehre, R.F. Stewart and J.A. Pople, J. Chem. Phys. 2657 (1969).\n! Na - Ar: W.J. Hehre, R. Ditchfield, R.F. Stewart, J.A. Pople,\n!          J. Chem. Phys.  2769 (1970).\n! K,Ca - : W.J. Pietro, B.A. Levy, W.J. Hehre and R.F. Stewart,\n! Ga - Kr: J. Am. Chem. Soc. 19, 2225 (1980).\n! Sc - Zn: W.J. Pietro and W.J. Hehre, J. Comp. Chem. 4, 241 (1983) + Gaussian.\n!  Y - Cd: W.J. Pietro and W.J. Hehre, J. Comp. Chem. 4, 241 (1983). + Gaussian\n!   \n\n\n****\nH     0 \nS   3   1.00\n      3.42525091             0.15432897       \n      0.62391373             0.53532814       \n      0.16885540             0.44463454       \n****\nHe     0 \nS   3   1.00\n      6.36242139             0.15432897       \n      1.15892300             0.53532814       \n      0.31364979             0.44463454       \n****\nLi     0 \nS   3   1.00\n     16.1195750              0.15432897       \n      2.9362007              0.53532814       \n      0.7946505              0.44463454       \nSP   3   1.00\n      0.6362897             -0.09996723             0.15591627       \n      0.1478601              0.39951283             0.60768372       \n      0.0480887              0.70011547             0.39195739       \n****\nBe     0 \nS   3   1.00\n     30.1678710              0.15432897       \n      5.4951153              0.53532814       \n      1.4871927              0.44463454       \nSP   3   1.00\n      1.3148331             -0.09996723             0.15591627       \n      0.3055389              0.39951283             0.60768372       \n      0.0993707              0.70011547             0.39195739       \n****\nB     0 \nS   3   1.00\n     48.7911130              0.15432897       \n      8.8873622              0.53532814       \n      2.4052670              0.44463454       \nSP   3   1.00\n      2.2369561             -0.09996723             0.15591627       \n      0.5198205              0.39951283             0.60768372       \n      0.1690618              0.70011547             0.39195739       \n****\nC     0 \nS   3   1.00\n     71.6168370              0.15432897       \n     13.0450960              0.53532814       \n      3.5305122              0.44463454       \nSP   3   1.00\n      2.9412494             -0.09996723             0.15591627       \n      0.6834831              0.39951283             0.60768372       \n      0.2222899              0.70011547             0.39195739       \n****\nN     0 \nS   3   1.00\n     99.1061690              0.15432897       \n     18.0523120              0.53532814       \n      4.8856602              0.44463454       \nSP   3   1.00\n      3.7804559             -0.09996723             0.15591627       \n      0.8784966              0.39951283             0.60768372       \n      0.2857144              0.70011547             0.39195739       \n****\nO     0 \nS   3   1.00\n    130.7093200              0.15432897       \n     23.8088610              0.53532814       \n      6.4436083              0.44463454       \nSP   3   1.00\n      5.0331513             -0.09996723             0.15591627       \n      1.1695961              0.39951283             0.60768372       \n      0.3803890              0.70011547             0.39195739       \n****\nF     0 \nS   3   1.00\n    166.6791300              0.15432897       \n     30.3608120              0.53532814       \n      8.2168207              0.44463454       \nSP   3   1.00\n      6.4648032             -0.09996723             0.15591627       \n      1.5022812              0.39951283             0.60768372       \n      0.4885885              0.70011547             0.39195739       \n****\nNe     0 \nS   3   1.00\n    207.0156100              0.15432897       \n     37.7081510              0.53532814       \n     10.2052970              0.44463454       \nSP   3   1.00\n      8.2463151             -0.09996723             0.15591627       \n      1.9162662              0.39951283             0.60768372       \n      0.6232293              0.70011547             0.39195739       \n****\nNa     0 \nS   3   1.00\n    250.7724300              0.1543289673     \n     45.6785110              0.5353281423     \n     12.3623880              0.4446345422     \nSP   3   1.00\n     12.0401930             -0.09996722919          0.1559162750     \n      2.7978819              0.39951282610          0.6076837186     \n      0.9099580              0.70011546890          0.3919573931     \nSP   3   1.00\n      1.4787406             -0.2196203690           0.01058760429    \n      0.4125649              0.2255954336           0.59516700530    \n      0.1614751              0.9003984260           0.46200101200    \n****\nMg     0 \nS   3   1.00\n    299.2374000              0.1543289673     \n     54.5064700              0.5353281423     \n     14.7515800              0.4446345422     \nSP   3   1.00\n     15.1218200             -0.09996722919          0.1559162750     \n      3.5139870              0.39951282610          0.6076837186     \n      1.1428570              0.70011546890          0.3919573931     \nSP   3   1.00\n      1.3954480             -0.2196203690           0.01058760429    \n      0.3893260              0.2255954336           0.59516700530    \n      0.1523800              0.9003984260           0.46200101200    \n****\nAl     0 \nS   3   1.00\n    351.4214767              0.1543289673     \n     64.01186067             0.5353281423     \n     17.32410761             0.4446345422     \nSP   3   1.00\n     18.89939621            -0.09996722919          0.1559162750     \n      4.391813233            0.39951282610          0.6076837186     \n      1.428353970            0.70011546890          0.3919573931     \nSP   3   1.00\n      1.3954482930          -0.2196203690           0.01058760429    \n      0.3893265318           0.2255954336           0.59516700530    \n      0.1523797659           0.9003984260           0.46200101200    \n****\nSi     0 \nS   3   1.00\n    407.7975514              0.1543289673     \n     74.28083305             0.5353281423     \n     20.10329229             0.4446345422     \nSP   3   1.00\n     23.19365606            -0.09996722919          0.1559162750     \n      5.389706871            0.39951282610          0.6076837186     \n      1.752899952            0.70011546890          0.3919573931     \nSP   3   1.00\n      1.4787406220          -0.2196203690           0.01058760429    \n      0.4125648801           0.2255954336           0.59516700530    \n      0.1614750979           0.9003984260           0.46200101200    \n****\nP     0 \nS   3   1.00\n    468.3656378              0.1543289673     \n     85.31338559             0.5353281423     \n     23.08913156             0.4446345422     \nSP   3   1.00\n     28.03263958            -0.09996722919          0.1559162750     \n      6.514182577            0.39951282610          0.6076837186     \n      2.118614352            0.70011546890          0.3919573931     \nSP   3   1.00\n      1.7431032310          -0.2196203690           0.01058760429    \n      0.4863213771           0.2255954336           0.59516700530    \n      0.1903428909           0.9003984260           0.46200101200    \n****\nS     0 \nS   3   1.00\n    533.1257359              0.1543289673     \n     97.10951830             0.5353281423     \n     26.28162542             0.4446345422     \nSP   3   1.00\n     33.32975173            -0.09996722919          0.1559162750     \n      7.745117521            0.39951282610          0.6076837186     \n      2.518952599            0.70011546890          0.3919573931     \nSP   3   1.00\n      2.0291942740          -0.2196203690           0.01058760429    \n      0.5661400518           0.2255954336           0.59516700530    \n      0.2215833792           0.9003984260           0.46200101200    \n****\nCl     0 \nS   3   1.00\n    601.3456136              0.1543289673     \n    109.5358542              0.5353281423     \n     29.64467686             0.4446345422     \nSP   3   1.00\n     38.96041889            -0.09996722919          0.1559162750     \n      9.053563477            0.39951282610          0.6076837186     \n      2.944499834            0.70011546890          0.3919573931     \nSP   3   1.00\n      2.1293864950          -0.2196203690           0.01058760429    \n      0.5940934274           0.2255954336           0.59516700530    \n      0.2325241410           0.9003984260           0.46200101200    \n****\nAr     0 \nS   3   1.00\n    674.4465184              0.1543289673     \n    122.8512753              0.5353281423     \n     33.24834945             0.4446345422     \nSP   3   1.00\n     45.16424392            -0.09996722919          0.1559162750     \n     10.49519900             0.39951282610          0.6076837186     \n      3.413364448            0.70011546890          0.3919573931     \nSP   3   1.00\n      2.6213665180          -0.2196203690           0.01058760429    \n      0.7313546050           0.2255954336           0.59516700530    \n      0.2862472356           0.9003984260           0.46200101200    \n****\nK     0 \nS   3   1.00\n    771.5103681              0.1543289673     \n    140.5315766              0.5353281423     \n     38.03332899             0.4446345422     \nSP   3   1.00\n     52.40203979            -0.0999672292           0.1559162750     \n     12.17710710             0.3995128261           0.6076837186     \n      3.960373165            0.7001154689           0.3919573931     \nSP   3   1.00\n      3.651583985           -0.2196203690           0.0105876043     \n      1.018782663            0.2255954336           0.5951670053     \n      0.3987446295           0.9003984260           0.4620010120     \nSP   3   1.00\n      0.5039822505          -0.3088441215          -0.1215468600     \n      0.1860011465           0.0196064117           0.5715227604     \n      0.08214006743          1.1310344420           0.5498949471     \n****\nCa     0 \nS   3   1.00\n    854.0324951              0.1543289673     \n    155.5630851              0.5353281423     \n     42.10144179             0.4446345422     \nSP   3   1.00\n     59.56029944            -0.0999672292           0.1559162750     \n     13.84053270             0.3995128261           0.6076837186     \n      4.501370797            0.7001154689           0.3919573931     \nSP   3   1.00\n      4.374706256           -0.2196203690           0.0105876043     \n      1.220531941            0.2255954336           0.5951670053     \n      0.477707930            0.9003984260           0.4620010120     \nSP   3   1.00\n      0.4558489757          -0.3088441215          -0.1215468600     \n      0.1682369410           0.0196064117           0.5715227604     \n      0.0742952070           1.1310344420           0.5498949471     \n****";

const char basis321g[]="!  3-21G  EMSL  Basis Set Exchange Library   10/13/15 11:18 AM\n! Elements                             References\n! --------                             ----------\n!  H - Ne: J.S. Binkley, J.A. Pople, W.J. Hehre, J. Am. Chem. Soc 102 939 (1980)\n! Na - Ar: M.S. Gordon, J.S. Binkley, J.A. Pople, W.J. Pietro and W.J. Hehre, \n!          J. Am. Chem. Soc. 104, 2797 (1983).\n!  K - Ca: K.D. Dobbs, W.J. Hehre, J. Comput. Chem. 7, 359 (1986). \n! Ga - Kr: K.D. Dobbs, W.J. Hehre, J. Comput. Chem. 7, 359 (1986).\n! Sc - Zn: K.D. Dobbs, W.J. Hehre, J. Comput. Chem. 8, 861 (1987). \n!  Y - Cd: K.D. Dobbs, W.J. Hehre, J. Comput. Chem. 8, 880 (1987). \n! Cs     : A 3-21G quality set derived from the Huzinage MIDI basis sets.\n!          E.D. Glendening and D. Feller, J. Phys. Chem. 99, 3060 (1995)\n!   \n\n\n****\nH     0 \nS   2   1.00\n      5.4471780              0.1562850        \n      0.8245470              0.9046910        \nS   1   1.00\n      0.1831920              1.0000000        \n****\nHe     0 \nS   2   1.00\n     13.6267000              0.1752300        \n      1.9993500              0.8934830        \nS   1   1.00\n      0.3829930              1.0000000        \n****\nLi     0 \nS   3   1.00\n     36.8382000              0.0696686        \n      5.4817200              0.3813460        \n      1.1132700              0.6817020        \nSP   2   1.00\n      0.5402050             -0.2631270              0.1615460        \n      0.1022550              1.1433900              0.9156630        \nSP   1   1.00\n      0.0285650              1.0000000              1.0000000        \n****\nBe     0 \nS   3   1.00\n     71.8876000              0.0644263        \n     10.7289000              0.3660960        \n      2.2220500              0.6959340        \nSP   2   1.00\n      1.2954800             -0.4210640              0.2051320        \n      0.2688810              1.2240700              0.8825280        \nSP   1   1.00\n      0.0773500              1.0000000              1.0000000        \n****\nB     0 \nS   3   1.00\n    116.4340000              0.0629605        \n     17.4314000              0.3633040        \n      3.6801600              0.6972550        \nSP   2   1.00\n      2.2818700             -0.3686620              0.2311520        \n      0.4652480              1.1994400              0.8667640        \nSP   1   1.00\n      0.1243280              1.0000000              1.0000000        \n****\nC     0 \nS   3   1.00\n    172.2560000              0.0617669        \n     25.9109000              0.3587940        \n      5.5333500              0.7007130        \nSP   2   1.00\n      3.6649800             -0.3958970              0.2364600        \n      0.7705450              1.2158400              0.8606190        \nSP   1   1.00\n      0.1958570              1.0000000              1.0000000        \n****\nN     0 \nS   3   1.00\n    242.7660000              0.0598657        \n     36.4851000              0.3529550        \n      7.8144900              0.7065130        \nSP   2   1.00\n      5.4252200             -0.4133010              0.2379720        \n      1.1491500              1.2244200              0.8589530        \nSP   1   1.00\n      0.2832050              1.0000000              1.0000000        \n****\nO     0 \nS   3   1.00\n    322.0370000              0.0592394        \n     48.4308000              0.3515000        \n     10.4206000              0.7076580        \nSP   2   1.00\n      7.4029400             -0.4044530              0.2445860        \n      1.5762000              1.2215600              0.8539550        \nSP   1   1.00\n      0.3736840              1.0000000              1.0000000        \n****\nF     0 \nS   3   1.00\n    413.8010000              0.0585483        \n     62.2446000              0.3493080        \n     13.4340000              0.7096320        \nSP   2   1.00\n      9.7775900             -0.4073270              0.2466800        \n      2.0861700              1.2231400              0.8523210        \nSP   1   1.00\n      0.4823830              1.0000000              1.0000000        \n****\nNe     0 \nS   3   1.00\n    515.7240000              0.0581430        \n     77.6538000              0.3479510        \n     16.8136000              0.7107140        \nSP   2   1.00\n     12.4830000             -0.4099220              0.2474600        \n      2.6645100              1.2243100              0.8517430        \nSP   1   1.00\n      0.6062500              1.0000000              1.0000000        \n****\nNa     0 \nS   3   1.00\n    547.6130000              0.0674911        \n     82.0678000              0.3935050        \n     17.6917000              0.6656050        \nSP   3   1.00\n     17.5407000             -0.1119370              0.1282330        \n      3.7939800              0.2546540              0.4715330        \n      0.9064410              0.8444170              0.6042730        \nSP   2   1.00\n      0.5018240             -0.2196600              0.0090665        \n      0.0609458              1.0891200              0.9972020        \nSP   1   1.00\n      0.0244349              1.0000000              1.0000000        \n****\nMg     0 \nS   3   1.00\n    652.8410000              0.0675982        \n     98.3805000              0.3917780        \n     21.2996000              0.6666610        \nSP   3   1.00\n     23.3727000             -0.1102460              0.1210140        \n      5.1995300              0.1841190              0.4628100        \n      1.3150800              0.8963990              0.6069070        \nSP   2   1.00\n      0.6113490             -0.3611010              0.0242633        \n      0.1418410              1.2150500              0.9866730        \nSP   1   1.00\n      0.0464011              1.0000000              1.0000000        \n****\nAl     0 \nS   3   1.00\n    775.7370000              0.0668347        \n    116.9520000              0.3890610        \n     25.3326000              0.6694680        \nSP   3   1.00\n     29.4796000             -0.1079020              0.1175740        \n      6.6331400              0.1462450              0.4611740        \n      1.7267500              0.9237300              0.6055350        \nSP   2   1.00\n      0.9461600             -0.3203270              0.0519383        \n      0.2025060              1.1841200              0.9726600        \nSP   1   1.00\n      0.0639088              1.0000000              1.0000000        \n****\nSi     0 \nS   3   1.00\n    910.6550000              0.0660823        \n    137.3360000              0.3862290        \n     29.7601000              0.6723800        \nSP   3   1.00\n     36.6716000             -0.1045110              0.1133550        \n      8.3172900              0.1074100              0.4575780        \n      2.2164500              0.9514460              0.6074270        \nSP   2   1.00\n      1.0791300             -0.3761080              0.0671030        \n      0.3024220              1.2516500              0.9568830        \nSP   1   1.00\n      0.0933392              1.0000000              1.0000000        \n****\nP     0 \nS   3   1.00\n   1054.9000000              0.0655410        \n    159.1950000              0.3840360        \n     34.5304000              0.6745410        \nSP   3   1.00\n     44.2866000             -0.1021300              0.1108510        \n     10.1019000              0.0815920              0.4564950        \n      2.7399700              0.9697880              0.6069360        \nSP   2   1.00\n      1.2186500             -0.3714950              0.0915820        \n      0.3955460              1.2709900              0.9349240        \nSP   1   1.00\n      0.1228110              1.0000000              1.0000000        \n****\nS     0 \nS   3   1.00\n   1210.6200000              0.0650070        \n    182.7470000              0.3820400        \n     39.6673000              0.6765450        \nSP   3   1.00\n     52.2236000             -0.1003100              0.1096460        \n     11.9629000              0.0650880              0.4576490        \n      3.2891100              0.9814550              0.6042610        \nSP   2   1.00\n      1.2238400             -0.2860890              0.1647770        \n      0.4573030              1.2280600              0.8708550        \nSP   1   1.00\n      0.1422690              1.0000000              1.0000000        \n****\nCl     0 \nS   3   1.00\n   1376.4000000              0.0645827        \n    207.8570000              0.3803630        \n     45.1554000              0.6781900        \nSP   3   1.00\n     60.8014000             -0.0987639              0.1085980        \n     13.9765000              0.0511338              0.4586820        \n      3.8871000              0.9913370              0.6019620        \nSP   2   1.00\n      1.3529900             -0.2224010              0.2192160        \n      0.5269550              1.1825200              0.8223210        \nSP   1   1.00\n      0.1667140              1.0000000              1.0000000        \n****\nAr     0 \nS   3   1.00\n   1553.7100000              0.0641707        \n    234.6780000              0.3787970        \n     51.0121000              0.6797520        \nSP   3   1.00\n     70.0453000             -0.0974661              0.1076190        \n     16.1473000              0.0390569              0.4595760        \n      4.5349200              0.9999160              0.6000410        \nSP   2   1.00\n      1.5420900             -0.1768660              0.2556870        \n      0.6072670              1.1469000              0.7898420        \nSP   1   1.00\n      0.1953730              1.0000000              1.0000000        \n****\n\n";

#endif

