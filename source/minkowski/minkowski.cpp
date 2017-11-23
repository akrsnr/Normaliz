#include <iostream>
#include <vector>
#include <libnormaliz/cone.h>
#include <iomanip>
#include "Eigen/Dense"
#include "rref.hpp"

using namespace Eigen;
using std::vector;
using namespace libnormaliz;

typedef long long Integer;


void printComponents(const vector< vector<Integer> >& v, string type) {
    std::string delim;
    int index = 0;
    for ( auto const& i : v ) {
        delim = "";
        std::cout << type << " Component-" << ++index << std::endl;
        std::cout << "(";
        for ( auto const& j : i ) {
            std::cout << delim << j;
            delim = ", ";
        }
        std::cout << ")" << std::endl;
    }
    std::cout << std::endl;
}


// Return a vector of basis elements of "self" as seperate matrices.
template<typename T>
vector<vector <T> > findGens(const MatrixXf& A) {
    vector<MatrixXf> transposedMatrices;
    long cols = A.cols();
    long rows = A.rows();
    vector<vector <T> > v {static_cast<unsigned long>(rows), std::vector<T>(cols, 999999) };


    /*
     v.push_back(A);
    std::cout << A.cols() << std::endl;
    std::cout << A.col(1).transpose() << std::endl;
     */

    //std::cout << "Col num:" << A.cols() << std::endl;
    for (size_t i = 0; i < rows; ++i) {
        transposedMatrices.emplace_back(const_cast<MatrixXf&>(A).row(i).transpose());
        //std::cout << "Ray " << i + 1 << ": \n" << transposedMatrices.at(i) << std::endl;
    }


    for (size_t i = 0; i < rows; ++i) {
        //std::cout << "Ray " << i + 1 << ": \n";
        for (size_t j = 0; j < cols; ++j) {
            //std::cout << transposedMatrices.at(i).coeff(j, 0) << std::endl;
            v[i][j] = (transposedMatrices.at(i).coeff(j, 0));
        }
    }

    /*
      for (size_t i = 0; i < rows; ++i) {
        std::cout << "Ray " << i + 1 << ": \n";
        for (size_t j = 0; j < cols; ++j) {
            std::cout << v.at(i).at(j) << std::endl;
        }
    }
     */

    //printComponents(v, "Ray");

    std::cout << std::endl;

    return v;

}

/* Here, "in case" */
void removeDuplicatePairs(std::vector<std::vector<Integer> >& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}


Cone<Integer> createU(MatrixXf A) {
    //MatrixXf B = A;
    long columnSize = A.cols();
    long rowSize = A.rows();

    std::cout << "Here is the initial matrix m:\n" << A << std::endl;

    std::cout << "The matrix m is of size "
              << columnSize << "x" << rowSize << std::endl;
    std::cout << "It has " << A.size() << " coefficients" << std::endl;

    A.transposeInPlace();
    std::cout << "and after being transposed:\n" << A << std::endl;


    /*

    MatrixXf B{test.size(), test.at(0).size()};

    for ( const auto &row : test )
    {
        for ( const auto &s : row ) {
            std::cout << s << '\t';
        }
        std::cout << std::endl;
    }

    for (int i = 0; i < test.size(); ++i) {
        for (int j = 0; j < test.at(i).size(); ++j) {
            B.row(i).coeffRef(j) = test.at(i).at(j);
        }
    }

    */






    FullPivLU<MatrixXf> lu(A);
    MatrixXf A_null_space = lu.kernel();
    std::cout << "Null space:\n" << A_null_space << std::endl;
    A_null_space.transposeInPlace();
    std::cout << "Null space Transposed_A but different ouput from ZZ's:\n" << A_null_space;

    std::cout << "soner test\n" << std::endl;

    vector<vector <double> > test = findGens<double>(A_null_space);
    to_reduced_row_echelon_form(test);
    vector<vector <Integer> > intVec {test.size(), vector<Integer>(test.at(0).size(), 999999) };


    for (int i = 0; i < test.size(); ++i) {
        std::cerr << "size = " << i << std::endl;
        for (int j = 0; j < test.at(i).size(); ++j) {
            std::cerr << "j = " << j << std::endl;
            intVec[i][j] = static_cast<Integer>(test.at(i).at(j));
        }
    }


    for ( const auto &row : intVec ) {
        for ( const auto &s : row ) {
            std::cout << s << '\t';
        }
        std::cout << std::endl;
    }

    std::cout << "soner test end\n" << std::endl;

    std::cout << std::endl << std::endl;
    std::cout << std::endl << std::endl;


    std::cout << "\n\nRays:" << std::endl;
    //vector<MatrixXf> argumentMatrixRays =
    //vector<vector <Integer> > argumentMatrixRays = findGens<Integer>(A_null_space);
    vector<vector <Integer> > argumentMatrixRays = intVec;
    /* *** We have null spaced matrix from now on *** */


    std::cout << "\nIdentity Matrix " << rowSize << " x " << rowSize << std::endl;
    // Row-sized identity matrix to get positive (orthant)
    MatrixXf identity = MatrixXf::Identity(rowSize, rowSize);
    std::cout << identity << std::endl;


    std::cout << "\nIdentity Matrix's Rays " << rowSize << " x " << rowSize << std::endl;
    //vector<MatrixXf> identityMatrixRays =
    vector<vector <Integer> > identityMatrixRays = findGens<Integer>(identity);


    Type::InputType type = Type::cone;
    Cone<Integer> coneMatrixA = Cone<Integer>(type, argumentMatrixRays);
    Cone<Integer> coneIdentity = Cone<Integer>(type, identityMatrixRays);

    /* OMITTED FOR NOW */
    const vector< vector<Integer> >& supportHyperPlanesA = coneMatrixA.getSupportHyperplanes();
    const vector< vector<Integer> >& supportHyperPlanesIdentity = coneIdentity.getSupportHyperplanes();
    /* OMITTED FOR NOW */

    vector< vector<Integer> > ineqsResultingCone;
    vector< vector<Integer> > equationsResultingCone;

    /* Pairs A(kerneled)
     *
       Type::inequalities
       Type::equations
       Type::congruences
     */
    map< InputType , vector< vector<Integer> > > pairsA =
        coneMatrixA.getConstraints();

    /* Pairs Identity
     *
       Type::inequalities
       Type::equations
       Type::congruences
     */
    map< InputType , vector< vector<Integer> > > pairsIdentity =
        coneIdentity.getConstraints();


    /* Equations Matrix A(kerneled) */
    const vector< vector<Integer> >& equationsA = pairsA[Type::equations];
    std::cout << "Equations A(kerneled) vector size: " << equationsA.size() << std::endl;
    printComponents(equationsA, "Equations A(kerneled)");

    /* Equations identity */
    const vector< vector<Integer> >& equationsIdentity = pairsIdentity[Type::equations];
    std::cout << "Equations Identity vector size: " << equationsIdentity.size() << std::endl;
    printComponents(equationsIdentity, "Equations Identity");

    /* Gathering "Equations" */
    equationsResultingCone.reserve(equationsA.size() + equationsIdentity.size());
    equationsResultingCone.insert( equationsResultingCone.end(), equationsA.begin(), equationsA.end() );
    equationsResultingCone.insert( equationsResultingCone.end(), equationsIdentity.begin(), equationsIdentity.end() );

    /* Inequalities Matrix A(kerneled) */
    const vector< vector<Integer> >& ineqA = pairsA[Type::inequalities];
    std::cout << "Inequalities A(kerneled) vector size: " << ineqA.size() << std::endl;
    printComponents(ineqA, "Inequalities A(kerneled)");


    /* Inequalities identity */
    const vector< vector<Integer> >& ineqIdentity = pairsIdentity[Type::inequalities];
    std::cout << "Inequalities Identity vector size: " << ineqIdentity.size() << std::endl;
    printComponents(ineqIdentity, "Inequalities A(kerneled)");

    /* Gathering "Inequalities" */
    ineqsResultingCone.reserve(ineqA.size() + ineqIdentity.size());
    ineqsResultingCone.insert( ineqsResultingCone.end(), ineqA.begin(), ineqA.end() );
    ineqsResultingCone.insert( ineqsResultingCone.end(), ineqIdentity.begin(), ineqIdentity.end() );



    Cone<Integer> resultingCone = Cone<Integer>(Type::equations, equationsResultingCone,
                                                Type::inequalities, ineqsResultingCone);


    std::cout << "Number of Excluded Faces: " << resultingCone.getNrDeg1Elements();

    /*
    map< InputType , vector< vector<Integer> > > TypesResultingCone =
                                                        resultingCone.getConstraints();

    printComponents(TypesResultingCone[Type::inequalities], "Inequalities of Resulting Cone");
    printComponents(TypesResultingCone[Type::equations], "Equations of Resulting Cone");
     */

    return resultingCone;

}


int main() {

    MatrixXf A{10, 3};
    //A << 0, 0, 1,   0, 1, 0,   1, 0, 0,   -1, 0, 0,    0, 0, -1,   0, -1, 0;
    //A << 0, 0, 1,   0, 1, 0,   1, 0, 0,   -1, 0, 0,    0, 0, -1,   0, -1, 0;



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



    /*
    MatrixXf A{3, 3};
    A << 1,2,-3,
        2,-1,4,
        4,3,-2;
*/

/*
    A.transposeInPlace();
    std::cout << "and after being transposed:\n" << A << std::endl;
    FullPivLU<MatrixXf> lu(A);
    MatrixXf A_null_space = lu.kernel();
    std::cout << "Null space:\n" << A_null_space << std::endl;
    A_null_space.transposeInPlace();
    std::cout << "Null space Transposed_A but different ouput from ZZ's:\n" << A_null_space;
*/
    /*HouseholderQR<MatrixXf> qr(A);
    cout << "soner \n" << qr.matrixQR().transpose();*/

    /*HouseholderQR<MatrixXf> lux(A);
    MatrixXf A_null_space = lux*/





     Cone<Integer> resultingCone = createU(A);

    std::cout << "\n ~~ ~~ From now on main()-definitions are being run ~~ ~~" << std::endl;

    map< InputType , vector< vector<Integer> > > TypesResultingCone =
        resultingCone.getConstraints();

    printComponents(TypesResultingCone[Type::inequalities], "Inequalities of Resulting Cone");
    printComponents(TypesResultingCone[Type::equations], "Equations of Resulting Cone");




    std::cout<< "\n\nCONE TEST\n";

    /*
     vector<vector <Integer> > data =  { {1,1,1}, {2,2,2} };
    Type::InputType type = Type::cone;
    Cone<Integer> MyCone = Cone<Integer>(type, data);
     */

    std::cout<< "CONE TEST END\n";

    MatrixXf b{6, 1};
    b << 1, 1, 1, 1, 1, 1;
    std::cout << "Here is the initial matrix d:\n" << b << std::endl;

    /*
    FullPivLU<MatrixXf> luB(A);
    MatrixXf B_null_space = luB.kernel();
    B_null_space.transposeInPlace();
    std::cout << "Null space Transposed_B:\n" << B_null_space << std::endl;
     */

    //findGens(A);



    return 0;
}