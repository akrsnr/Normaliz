#include <iostream>
#include <algorithm>
#include <vector>
#include <libnormaliz/cone.h>
#include <iomanip>
#include "Eigen/Dense"
#include "rref.hpp"

#include <flint/fmpzxx.h>
#include <unordered_set>
#include "fmpz.h"
#include "fmpz_mat.h"

using namespace Eigen;
using std::vector;
using namespace libnormaliz;

typedef long long Integer;

const int INITIAL_VALUE = 999999;

void printComponents(const vector< vector<Integer> >& v, const string &type) {
    std::string delim;
    int index = 0;
    std::cout << type;
    for ( auto const& i : v ) {
        delim = "";
        std::cout << "\n";
        std::cout << "[";
        for ( auto const& j : i ) {
            std::cout << delim << j;
            delim = "\t";
        }
        std::cout << "],   ";
    }
    std::cout << "\n\n";
}


// Return a vector of basis elements of "self" as seperate matrices.
template<typename T>
vector<vector <T> > eigenTOvector(const MatrixXf &A) {
    vector<MatrixXf> transposedMatrices;
    long cols = A.cols();
    long rows = A.rows();
    vector<vector <T> > v {static_cast<unsigned long>(rows), std::vector<T>(cols, INITIAL_VALUE) };

    for (size_t i = 0; i < rows; ++i) {
        transposedMatrices.emplace_back(const_cast<MatrixXf&>(A).row(i).transpose());
    }

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            v[i][j] = (transposedMatrices.at(i).coeff(j, 0));
        }
    }
    //printComponents(v, "Ray");
    std::cout << std::endl;
    return v;

}


vector< vector<Integer > > hermiteNormalForm(const MatrixXf& A) {
    //std::cout << "original matrix \n" << A;
    MatrixXf B = A;
    B.transposeInPlace();
    const vector< vector< Integer> >& transposedMatrix = eigenTOvector<Integer>(B);
    unsigned long row = transposedMatrix.size();
    unsigned long col = transposedMatrix.at(0).size();
    //std::cerr << "row  " << row << "  col  " << col << '\n';


    fmpz_mat_t M {};
    fmpz_mat_t M_null {};
    fmpz_mat_t M_null_transpose {};
    fmpz_mat_t M_hermit {};


    fmpz_mat_init(M, row,col);
    fmpz_mat_init(M_null, col,col);
    fmpz_mat_init(M_null_transpose, col,col);
    fmpz_mat_init(M_hermit, col,col);

    /* Fill FLINT matrix */
    for (size_t i = 0; i < row; i++) {
        for (size_t j = 0; j < col; j++) {
            //std::cout << transposedMatrix.at(i).at(j) << "  ";
            fmpz_set_si(fmpz_mat_entry(M, i, j), transposedMatrix.at(i).at(j));
        }
        //std::cout << "\n";
    }


    /* Calculation hermite normal form of given matrix */
    long dimension = fmpz_mat_nullspace(M_null, M);
    fmpz_mat_transpose(M_null_transpose, M_null);
    fmpz_mat_hnf(M_hermit, M_null_transpose);

    /*
     flint_printf(" M = \n");
    fmpz_mat_print_pretty(M);

    std::cout << "dimension = " << dimension << std::endl;
    flint_printf(" \nNull space = \n\n");
    fmpz_mat_print_pretty(M_null);

    flint_printf(" \nNull space transpose = \n\n");
    fmpz_mat_print_pretty(M_null_transpose);

    flint_printf(" \nHermit Form = \n\n");
    fmpz_mat_print_pretty(M_hermit);
     */

    /* Up to here, we found hermite normal form of given matrix */


    /* In HNF, sometimes the library doesn't give simplified result. For example,
     * [2 0 2 0]
     * [0 2 0 2]
     * something similar matrices can be exist.
     * I write a code snippet overcoming the problem simply.
     * I put all elements in a 1D vector and 2D vector.
     * I convert negative numbers into positive numbers in 1D vector.
     */
    vector<Integer> duplicateFind {};
    vector< vector<Integer> > hermit {static_cast<unsigned long>(dimension), vector<Integer>(col, INITIAL_VALUE) };
    for (int i = 0; i < dimension; i++) {
        for (int j = 0; j < col; j++) {
            fmpz *val = fmpz_mat_entry(M_hermit, i, j);
            //std::cout << *val << "  ";
            hermit[i][j] = *val;
            if (*val < 0)
                duplicateFind.push_back(*val * -1);
            else
                duplicateFind.push_back(*val);
        }
        //std::cout << "\n";
    }

    /* There cannot be same elements in a set
     * To get rid of sorting performance penalty, I put 1D elements in unordered_set
     * Then again convert to vector
     */
    std::unordered_set<int> s;
    for (int i : duplicateFind)
        s.insert(i);
    duplicateFind.assign( s.begin(), s.end() );

    /*
    for (auto const& v : duplicateFind) {
        std::cout << v << "  ";
    }
    */

    /*  If the vector now contains only two elements and smaller one is 0
     * I divide all values by the max value
     */
    if (duplicateFind.size() == 2) {
        long max;
        long min;
        if (duplicateFind[0] < duplicateFind[1]) {
            max = duplicateFind[1];
            min = duplicateFind[0];
        } else {
            max = duplicateFind[0];
            min = duplicateFind[1];
        }
        if (max == 0 || min == 0)
            for (int i = 0; i < dimension; i++) {
                for (int j = 0; j < col; j++) {
                    hermit[i][j] /= max;
                }
            }
    }

    /* Free memory */
    fmpz_mat_clear(M_hermit);
    fmpz_mat_clear(M_null_transpose);
    fmpz_mat_clear(M_null);
    fmpz_mat_clear(M);

    /* return hermit 2d vector, echelon basis matrix */
    return hermit;
}


Cone<Integer> createU(const MatrixXf &A) {
    long columnSize = A.cols();
    long rowSize = A.rows();

    const vector<vector <Integer> >& matrixAInequalities = hermiteNormalForm(A);

    std::cout << "Hermit Normal Form\n";
    for (auto const& v : matrixAInequalities) {
        for (auto const i : v) {
            std::cout << i << "  ";
        }
        std::cout << '\n';
    }


    std::cout << "\nIdentity Matrix " << rowSize << " x " << rowSize << std::endl;
    // Row-sized identity matrix to get positive (orthant)
    MatrixXf identity = MatrixXf::Identity(rowSize, rowSize);
    std::cout << identity << std::endl;


    const vector<vector <Integer> >& identityMatrixRays = eigenTOvector<Integer>(identity);

    Type::InputType type = Type::cone;
    Cone<Integer> coneMatrixA = Cone<Integer>(type, matrixAInequalities);
    Cone<Integer> coneIdentity = Cone<Integer>(type, identityMatrixRays);

    /* OMITTED FOR NOW */
    const vector< vector<Integer> >& supportHyperPlanesA = coneMatrixA.getSupportHyperplanes();
    const vector< vector<Integer> >& supportHyperPlanesIdentity = coneIdentity.getSupportHyperplanes();
    /* OMITTED FOR NOW */

    vector< vector<Integer> > ineqsResultingCone;
    vector< vector<Integer> > equationsResultingCone;

    /* Pairs of A(kerneled), argument, Matrix
     *
       Type::inequalities
       Type::equations
       Type::congruences
     */
    map< InputType , vector< vector<Integer> > > pairsA =
        coneMatrixA.getConstraints();

    /* Pairs of Identity Matrix
     *
       Type::inequalities
       Type::equations
       Type::congruences
     */
    map< InputType , vector< vector<Integer> > > pairsIdentity =
        coneIdentity.getConstraints();

    /*    -------   NULL-SPACED MATRIX EQUATIONS and INEQUALITIES    -------   */

    /* Equations Matrix A(kerneled) */
    vector< vector<Integer> > equationsA = pairsA[Type::equations];
    std::cout << "Equations A(kerneled) vector size: " << equationsA.size();
    for (auto& v : equationsA) {
        std::reverse(std::begin(v), std::end(v));
    }
    printComponents(equationsA, "");

    /* Inequalities Matrix A(kerneled) */
    const vector< vector<Integer> >& ineqA = pairsA[Type::inequalities];
    std::cout << "Inequalities A(kerneled) vector size: " << ineqA.size();
    printComponents(ineqA, "");

    /*    -------   IDENTITY MATRIX EQUATIONS and INEQUALITIES    -------   */

    /* Equations identity */
    const vector< vector<Integer> >& equationsIdentity = pairsIdentity[Type::equations];
    std::cout << "Equations Identity vector size: " << equationsIdentity.size();
    printComponents(equationsIdentity, "");


    /* Inequalities identity */
    const vector< vector<Integer> >& ineqIdentity = pairsIdentity[Type::inequalities];
    std::cout << "Inequalities Identity vector size: " << ineqIdentity.size();
    printComponents(ineqIdentity, "");

    /*    -------   COMBINE ALL OBTAINED EQUATIONS and INEQUALITIES    -------   */

    /* Gathering "Equations" */
    equationsResultingCone.reserve(equationsA.size() + equationsIdentity.size());
    equationsResultingCone.insert( equationsResultingCone.end(), equationsA.begin(), equationsA.end() );
    equationsResultingCone.insert( equationsResultingCone.end(), equationsIdentity.begin(), equationsIdentity.end() );
    //printComponents(equationsResultingCone, "All equations");


    /* Gathering "Inequalities" */
    ineqsResultingCone.reserve(ineqA.size() + ineqIdentity.size());
    ineqsResultingCone.insert( ineqsResultingCone.end(), ineqA.begin(), ineqA.end() );
    ineqsResultingCone.insert( ineqsResultingCone.end(), ineqIdentity.begin(), ineqIdentity.end() );
    //printComponents(ineqsResultingCone, "All inequalities");


    /*                                   INTERSECTION                                       */
    /*    -------   CREATE A NEW CONE FROM COMBINED EQUATIONS and INEQUALITIES    -------   */

    Cone<Integer> resultingCone = Cone<Integer>(Type::equations, equationsResultingCone,
                                                Type::inequalities, ineqsResultingCone);

    resultingCone.compute(ConeProperty::Generators);
    const vector<vector <Integer> >& gens = resultingCone.getGenerators();

    return resultingCone;

}


int main() {


    //A << 0, 0, 1,   0, 1, 0,   1, 0, 0,   -1, 0, 0,    0, 0, -1,   0, -1, 0;
    //A << 0, 0, 1,   0, 1, 0,   1, 0, 0,   -1, 0, 0,    0, 0, -1,   0, -1, 0;


    /*
//pdf
    MatrixXf A{10, 3};
    A <<
      1, 0, 1 ,
        1, 0, 0 ,
        0, 1, 1 ,
        0, 1, 0 ,
        0, 0, 1 ,
        -1, 0, 0 ,
        0, 0, -1 ,
        0, -1, 1 ,
        0, -1, 0 ,
        -1, 0, 1;
*/

/*
    // CUBE
    MatrixXf A{6,3};
    A << 0, 0, 1,   0, 1, 0,   1, 0, 0,   -1, 0, 0,    0, 0, -1,   0, -1, 0;
*/

    /*
    //cross_polytope(2)
    MatrixXf A{4,2};
    A << -1, -1, -1, 1, 1, 1, 1, -1;
    */



      //cross_polytope(3)
    MatrixXf A{8,3};
    A << 	-1, 1, 1,
            -1, -1, 1,
            -1, -1, -1,
            -1, 1, -1,
             1, 1, -1,
             1, 1, 1,
             1, -1, 1,
             1, -1, -1;



     //cuboctahedron()
/*
    MatrixXf A{14,3};
    A << 0, 1, 0,
        1, 0, 0,
        -1, -1, 1,
        -1, 0, 0,
        -1, 1, 1,
        0, -1, 0,
        -1, 1, -1,
        0, 0, 1,
        0, 0, -1,
        1, -1, -1,
        -1, -1, -1,
        1, -1, 1,
        1, 1, -1,
        1, 1, 1;
*/


    Cone<Integer> resultingCone = createU(A);

    std::cout << "\n\n\n\n ~~ ~~ From now on main()-definitions are being run ~~ ~~" << std::endl;

    map< InputType , vector< vector<Integer> > > TypesResultingCone =
        resultingCone.getConstraints();

    printComponents(TypesResultingCone[Type::equations], "Equations of Intersected resulting Cone");
    printComponents(TypesResultingCone[Type::inequalities], "Inequalities of Intersected resulting Cone");






    //std::cout<< "\n\nCONE TEST\n";

    /*
     vector<vector <Integer> > data =  { {1,1,1}, {2,2,2} };
    Type::InputType type = Type::cone;
    Cone<Integer> MyCone = Cone<Integer>(type, data);
     */

    /* std::cout<< "CONE TEST END\n";

     MatrixXf b{6, 1};
     b << 1, 1, 1, 1, 1, 1;
     std::cout << "Here is the initial matrix d:\n" << b << std::endl;

     */

    /*
    FullPivLU<MatrixXf> luB(A);
    MatrixXf B_null_space = luB.kernel();
    B_null_space.transposeInPlace();
    std::cout << "Null space Transposed_B:\n" << B_null_space << std::endl;
     */

    //eigenTOvector(A);



    return 0;
}