#include <iostream>
#include <vector>
#include <libnormaliz/cone.h>
#include "Eigen/Dense"

#include <string>

#include <cassert>




using namespace Eigen;
using std::vector;
using namespace libnormaliz;

typedef long long Integer;


// Return a vector of basis elements of "self" as seperate matrices.
vector<vector <Integer> > findGens(const MatrixXf& A) {
    vector<MatrixXf> transposedMatrices;
    long cols = A.cols();
    long rows = A.rows();
    vector<vector <Integer> > v {static_cast<unsigned long>(rows), std::vector<Integer>(cols, 999999) };


    /*
     v.push_back(A);
    std::cout << A.cols() << std::endl;
    std::cout << A.col(1).transpose() << std::endl;
     */

    //std::cout << "Col num:" << A.cols() << std::endl;
    for (size_t i = 0; i < rows; ++i) {
        transposedMatrices.emplace_back(A.row(i).transpose());
        //std::cout << "Ray " << i + 1 << ": \n" << transposedMatrices.at(i) << std::endl;
    }


    /*for (auto const &i : transposedMatrices) {
        std::cout << "soner\n" << i << std::endl;
    }*/


    for (size_t i = 0; i < rows; ++i) {
        //std::cout << "Ray " << i + 1 << ": \n";
        for (size_t j = 0; j < cols; ++j) {
            //std::cout << transposedMatrices.at(i).coeff(j, 0) << std::endl;
            v[i][j] = static_cast<long long>(transposedMatrices.at(i).coeff(j, 0));
        }
    }

    for (size_t i = 0; i < rows; ++i) {
        std::cout << "Ray " << i + 1 << ": \n";
        for (size_t j = 0; j < cols; ++j) {
            std::cout << v.at(i).at(j) << std::endl;
        }
    }

    std::cout << std::endl;

    return v;

}

void removeDuplicatePairs(std::vector<std::vector<Integer> >& v) {
    std::sort(v.begin(), v.end());
    v.erase(std::unique(v.begin(), v.end()), v.end());
}

void printComponents(const vector< vector<Integer> >& v, string type) {
    int index = 0;
    for ( auto const& i : v ) {
        std::cout << type << " " << ++index << std::endl;
        for ( auto const& j : i ) {
            std::cout << j << std::endl;
        }
    }
    std::cout << std::endl;
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

    FullPivLU<MatrixXf> lu(A);
    MatrixXf A_null_space = lu.kernel();
    std::cout << "Null space:\n" << A_null_space << std::endl;
    A_null_space.transposeInPlace();
    std::cout << "Null space Transposed_A but different ouput from ZZ's:\n" << A_null_space;


    std::cout << "\n\nRays:" << std::endl;
    //vector<MatrixXf> argumentMatrixRays =
    vector<vector <Integer> > argumentMatrixRays = findGens(A_null_space);


    std::cout << "\nIdentity Matrix " << rowSize << " x " << rowSize << std::endl;
    // Row-sized identity matrix to get positive (orthant)
    MatrixXf identity = MatrixXf::Identity(rowSize, rowSize);
    std::cout << identity << std::endl;


    std::cout << "\nIdentity Matrix's Rays " << rowSize << " x " << rowSize << std::endl;
    //vector<MatrixXf> identityMatrixRays =
    vector<vector <Integer> > identityMatrixRays = findGens(identity);


    Type::InputType type = Type::cone;
    Cone<Integer> coneMatrixA = Cone<Integer>(type, argumentMatrixRays);
    Cone<Integer> coneIdentity = Cone<Integer>(type, identityMatrixRays);

    vector< vector<Integer> > ineqsResultingCone;
    vector< vector<Integer> > equationsResultingCone;


    const vector< vector<Integer> >& supportHyperPlanesA = coneMatrixA.getSupportHyperplanes();
    const vector< vector<Integer> >& supportHyperPlanesIdentity = coneIdentity.getSupportHyperplanes();

    //std::cout << supportHyperPlanesA << std::endl;
/*
    int index = 0;
    for ( auto const& i : supportHyperPlanesA ) {
        std::cout << "SupportHyperPlane Matrix A " << ++index << std::endl;
        for ( auto const& j : i ) {
            std::cout << j << std::endl;
        }
    }
    std::cout << std::endl;

    index = 0;
    for ( auto const& i : supportHyperPlanesIdentity ) {
        std::cout << "SupportHyperPlane Identity " << ++index << std::endl;
        for ( auto const& j : i ) {
            std::cout << j << std::endl;
        }
    }
    std::cout << std::endl;
*/

    /* Equations A */
    map< InputType , vector< vector<Integer> > > equationsA =
        coneMatrixA.getConstraints();

    const vector< vector<Integer> >& temp = equationsA[Type::equations];
    std::cout << "Equations A vector size: " << temp.size() << std::endl;


    int index = 0;
    for ( auto const& i : temp ) {
        // Gathering all matrices, here by "Equations" of the given matrix A
        equationsResultingCone.push_back(i);
        std::cout << "Equations A " << ++index << std::endl;
        for ( auto const& j : i ) {
            std::cout << j << std::endl;
        }
    }
    std::cout << std::endl;


    /* Inequalities A */
    const vector< vector<Integer> >& ineqA = equationsA[Type::inequalities];
    std::cout << "Inequalities A vector size: " << ineqA.size() << std::endl;


    index = 0;
    for ( auto const& i : ineqA ) {
        // Gathering all matrices, here by "Inequalities" of the given matrix A
        ineqsResultingCone.push_back(i);
        std::cout << "Inequalities A " << ++index << std::endl;
        for ( auto const& j : i ) {
            std::cout << j << std::endl;
        }
    }
    std::cout << std::endl;


    /* Equations identity */
    map< InputType , vector< vector<Integer> > > equationsIdentity =
        coneIdentity.getConstraints();

    const vector< vector<Integer> >& temp2 = equationsIdentity[Type::equations];
    std::cout << "Identity equality vector size: " << temp2.size() << std::endl;
    index = 0;
    for ( auto const& i : temp2 ) {
        // Gathering all matrices, here by Equations of the Identity Matrix
        equationsResultingCone.push_back(i);
        std::cout << "Equations Identity " << ++index << std::endl;
        for ( auto const& j : i ) {
            std::cout << j << std::endl;
        }
    }
    std::cout << std::endl;

    /* Inequalities identity */
    const vector< vector<Integer> >& ineqIdentity = equationsIdentity[Type::inequalities];
    std::cout << "Identity inequalities vector size: " << ineqIdentity.size() << std::endl;
    index = 0;
    for ( auto const& i : ineqIdentity ) {
        // Gathering all matrices, here by "Inequalities" of the Identity Matrix
        ineqsResultingCone.push_back(i);
        std::cout << "Inequalities Identity " << ++index << std::endl;
        for ( auto const& j : i ) {
            std::cout << j << std::endl;
        }
    }
    std::cout << std::endl;



    Cone<Integer> resultingCone = Cone<Integer>(Type::equations, equationsResultingCone,
                                                Type::inequalities, ineqsResultingCone);

    map< InputType , vector< vector<Integer> > > TypesResultingCone =
        resultingCone.getConstraints();

    printComponents(TypesResultingCone[Type::inequalities], "Inequalities of Resulting Cone");
    printComponents(TypesResultingCone[Type::equations], "Equations of Resulting Cone");


    return resultingCone;

}





int main() {
    MatrixXf A{6, 3};
    A << 0, 0, 1, 0, 1, 0, 1, 0, 0, -1, 0, 0, 0, 0, -1, 0, -1, 0;

    createU(A);



    std::cout << "\n ~~ ~~ From now on main()-definitions are being run ~~ ~~" << std::endl;

    std::cout<< "CONE TEST\n";

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