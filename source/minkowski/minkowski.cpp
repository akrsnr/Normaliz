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


const int INITIAL_VALUE = 999999;


void printComponents(const vector< vector<Integer> >& v, const string &type, bool onlyVector = false) {
    std::string delim;
    int index = 0;
    std::cout << type;
    for ( auto const& i : v ) {
        delim = "";
        /*if (onlyVector) {
            //std::cout << type << " Component-" << ++index << std::endl;
            //onlyVector = false;
        }*/
        std::cout << "\n";
        std::cout << "(";
        for ( auto const& j : i ) {
            std::cout << delim << j;
            delim = ", ";
        }
        std::cout << "),   ";
        /*if (onlyVector) {
            std::cout << std::endl;
        }*/
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


vector<vector <Integer> > fromRREFtoVectorInteger(const vector<vector <double> >& v ) {
    vector<vector <Integer> > intVec {v.size(), vector<Integer>(v.at(0).size(), INITIAL_VALUE) };

    for (int i = 0; i < v.size(); ++i) {
        for (int j = 0; j < v.at(i).size(); ++j) {
            intVec[i][j] = static_cast<Integer >( v.at(i).at(j) );
        }
    }
    return  intVec;
}


Cone<Integer> createU(MatrixXf A) {
    long columnSize = A.cols();
    long rowSize = A.rows();

    /* Irregular-ordered right kernel calculation */
    A.transposeInPlace();
    FullPivLU<MatrixXf> lu(A);
    MatrixXf A_null_space = lu.kernel();
    A_null_space.transposeInPlace();

    /* Row reduce echelon form to make it regular */
    vector<vector <double> > rref = eigenTOvector<double>(A_null_space);
    to_reduced_row_echelon_form(rref);
    vector<vector <Integer> > intVecA = fromRREFtoVectorInteger(rref);

    std::cout << "RREF form of kernel\n";
    for ( const auto &row : intVecA ) {
        for ( const auto &s : row ) {
            std::cout << s << '\t';
        }
        std::cout << std::endl;
    }


    const vector<vector <Integer> >& matrixAInequalities = intVecA;
    /* *** We have null spaced and rref'ed matrix from now on *** */


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


    /* Equations Matrix A(kerneled) */
    const vector< vector<Integer> >& equationsA = pairsA[Type::equations];
    //std::cout << "Equations A(kerneled) vector size: " << equationsA.size() << std::endl;
    //printComponents(equationsA, "Equations A(kerneled)");

    /* Equations identity */
    const vector< vector<Integer> >& equationsIdentity = pairsIdentity[Type::equations];
    //std::cout << "Equations Identity vector size: " << equationsIdentity.size() << std::endl;
    //printComponents(equationsIdentity, "Equations Identity");

    /* Gathering "Equations" */
    equationsResultingCone.reserve(equationsA.size() + equationsIdentity.size());
    equationsResultingCone.insert( equationsResultingCone.end(), equationsA.begin(), equationsA.end() );
    equationsResultingCone.insert( equationsResultingCone.end(), equationsIdentity.begin(), equationsIdentity.end() );

    /* Inequalities Matrix A(kerneled) */
    const vector< vector<Integer> >& ineqA = pairsA[Type::inequalities];
    //std::cout << "Inequalities A(kerneled) vector size: " << ineqA.size() << std::endl;
    //printComponents(ineqA, "-test-Inequalities A(kerneled)", true);


    /* Inequalities identity */
    const vector< vector<Integer> >& ineqIdentity = pairsIdentity[Type::inequalities];
    //std::cout << "Inequalities Identity vector size: " << ineqIdentity.size() << std::endl;
    //printComponents(ineqIdentity, "Inequalities A(kerneled)");

    /* Gathering "Inequalities" */
    ineqsResultingCone.reserve(ineqA.size() + ineqIdentity.size());
    ineqsResultingCone.insert( ineqsResultingCone.end(), ineqA.begin(), ineqA.end() );
    ineqsResultingCone.insert( ineqsResultingCone.end(), ineqIdentity.begin(), ineqIdentity.end() );


    Cone<Integer> resultingCone = Cone<Integer>(Type::equations, equationsResultingCone,
                                                Type::inequalities, ineqsResultingCone);

    /* AFTER TUESDAY, 28th MEETING */

    /*
     * Here, I find generators of intersected(resulting) cone.
     */
    resultingCone.compute(ConeProperty::Generators);
    const vector<vector <Integer> >& gens = resultingCone.getGenerators();

    /*
     * ************* Documentation says **************************
     * Unless you are interested a single result, we recommend
     * to use compute since the data asked for can then be computed in a single run.
     * ***********************************************************
     * Here, I use the generators of intersected(resulting) cone to create a cone anew.
     * Frankly, the generators are used to recreate a cone as inequalities
     */
    Cone<Integer> reCreatebyGenerators = Cone<Integer>(Type::inequalities, gens);
    //reCreatebyGenerators.compute(ConeProperty::Generators);
    const vector<vector <Integer> >& gens2 = reCreatebyGenerators.getGenerators();
    printComponents(gens2, "Recreated cone's Generators", true);



    return resultingCone;

}


int main() {


    //A << 0, 0, 1,   0, 1, 0,   1, 0, 0,   -1, 0, 0,    0, 0, -1,   0, -1, 0;
    //A << 0, 0, 1,   0, 1, 0,   1, 0, 0,   -1, 0, 0,    0, 0, -1,   0, -1, 0;


//pdf
    MatrixXf A{11, 3};
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
        -1, 0, 1,
        -1,0,-1;



    /*
    MatrixXf A{3, 3};
    A << 1,2,-3,
        2,-1,4,
        4,3,-2;
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

    /*MatrixXf A{8,3};
    A << 	-1, 1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, 1, -1, 1, 1,
        1, 1, -1, 1, 1, -1, -1;*/


    Cone<Integer> resultingCone = createU(A);

    std::cout << "\n\n\n\n ~~ ~~ From now on main()-definitions are being run ~~ ~~" << std::endl;

    map< InputType , vector< vector<Integer> > > TypesResultingCone =
        resultingCone.getConstraints();

    printComponents(TypesResultingCone[Type::inequalities], "Inequalities of Intersected resulting Cone", true);
    printComponents(TypesResultingCone[Type::equations], "Equations of Intersected resulting Cone", true);








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