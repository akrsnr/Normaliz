/*
 * Normaliz
 * Copyright (C) 2007-2014  Winfried Bruns, Bogdan Ichim, Christof Soeger
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * As an exception, when this program is distributed through (i) the App Store
 * by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or (iii) Google Play
 * by Google Inc., then that store may impose any digital rights management,
 * device limits and/or redistribution restrictions that are required by its
 * terms of service.
 */

//---------------------------------------------------------------------------
#ifndef MATRIX_HPP
#define MATRIX_HPP
//---------------------------------------------------------------------------


#include <vector>
#include <list>
#include <iostream>
#include <string>

#include <libQnormaliz/libQnormaliz.h>
#include <libQnormaliz/Qinteger.h>
#include <libQnormaliz/Qconvert.h>

//---------------------------------------------------------------------------

namespace libQnormaliz {
using std::list;
using std::vector;
using std::string;

template<typename Number> class Matrix {

    template<typename> friend class Matrix;
    template<typename> friend class Lineare_Transformation;
    template<typename> friend class Sublattice_Representation;
    
    // public:

    size_t nr;
    size_t nc;
    vector< vector<Number> > elem;

//---------------------------------------------------------------------------
//              Private routines, used in the public routines
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//                      Rows and columns exchange
//---------------------------------------------------------------------------

    void exchange_rows(const size_t& row1, const size_t& row2);      //row1 is exchanged with row2
    void exchange_columns(const size_t& col1, const size_t& col2); // col1 is exchanged with col2

//---------------------------------------------------------------------------
//              Row and column reduction
//---------------------------------------------------------------------------
    // return value false undicates failure because of overflow
    // for all the routines below
    
    // reduction via integer division and elemntary transformations
    bool reduce_row(size_t corner);      //reduction by the corner-th row
    bool reduce_row (size_t row, size_t col); // corner at position (row,col)
            
    // replaces two columns by linear combinations of them
    bool linear_comb_columns(const size_t& col,const size_t& j,
            const Number& u,const Number& w,const Number& v,const Number& z);
    
    // column reduction with gcd method
    bool gcd_reduce_column (size_t corner, Matrix<Number>& Right);
    
//---------------------------------------------------------------------------
//                      Work horses
//---------------------------------------------------------------------------

    // takes |product of the diagonal elem|
    Number compute_vol(bool& success);  
        
    // Solve system with coefficients and right hand side in one matrix, using elementary transformations
    // solution replaces right hand side
    bool solve_destructive_inner(bool ZZinvertible, Number& denom);

    // asembles the matrix of the system (left side the submatrix of mother given by key
    // right side from column vectors pointed to by RS
    // both in a single matrix    
    void solve_system_submatrix_outer(const Matrix<Number>& mother, const vector<key_t>& key, const vector<vector<Number>* >& RS,
         Number& denom, bool ZZ_invertible, bool transpose, size_t red_col, size_t sign_col,
         bool compute_denom=true, bool make_sol_prime=false);
                    
    size_t row_echelon_inner_elem(bool& success); // does the work and checks for overflows
    // size_t row_echelon_inner_bareiss(bool& success, Number& det);
    // NOTE: Bareiss cannot be used if z-invertible transformations are needed
    
    size_t row_echelon(bool& success); // transforms this into row echelon form and returns rank
    size_t row_echelon(bool& success, Number& det); // computes also |det|
    size_t row_echelon(bool& success, bool do_compute_vol, Number& det); // chooses elem (or bareiss in former time)
    
    // reduces the rows a matrix in row echelon form upwards, from left to right
    bool reduce_rows_upwards();
    size_t row_echelon_reduce(bool& success); // combines row_echelon and reduce_rows_upwards
    
    // computes rank and index simultaneously, returns rank
    Number full_rank_index(bool& success);
    
    vector<key_t> max_rank_submatrix_lex_inner(bool& success) const;
    
    // A version of invert that circumvents protection and leaves it to the calling routine
    Matrix invert_unprotected(Number& denom, bool& sucess) const;
    
    bool SmithNormalForm_inner(size_t& rk, Matrix<Number>& Right);
    

//---------------------------------------------------------------------------
//                      Pivots for rows/columns operations
//---------------------------------------------------------------------------

    vector<long> pivot(size_t corner); //Find the position of an element x with
    //0<abs(x)<=abs(y) for all y!=0 in the right-lower submatrix of this
    //described by an int corner

    long pivot_column(size_t col);  //Find the position of an element x with
    //0<abs(x)<=abs(y) for all y!=0 in the lower half of the column of this
    //described by an int col
    
    long pivot_column(size_t row,size_t col); //in column col starting from row
    
//---------------------------------------------------------------------------
//                     Helpers for linear systems
//---------------------------------------------------------------------------

    Matrix bundle_matrices(const Matrix<Number>& Right_side)const;
    Matrix extract_solution() const;
    vector<vector<Number>* > row_pointers();
    void customize_solution(size_t dim, Number& denom, size_t red_col, 
                     size_t sign_col, bool make_sol_prime);
                    
public:

size_t row_echelon_inner_bareiss(bool& success, Number& det);

    vector<vector<Number>* > submatrix_pointers(const vector<key_t>& key);     
  
//---------------------------------------------------------------------------

//---------------------------------------------------------------------------
//                      Construction and destruction
//---------------------------------------------------------------------------

    Matrix();
    Matrix(size_t dim);                           //constructor of unit matrix
    Matrix(size_t row, size_t col);                 //main constructor, all entries 0
    Matrix(size_t row, size_t col, Number value); //constructor, all entries set to value
    Matrix(const vector< vector<Number> >& elem); //constuctor, elem=elem
    Matrix(const list< vector<Number> >& elems);

//---------------------------------------------------------------------------
//                             Data access
//---------------------------------------------------------------------------

    void write_column(size_t col, const vector<Number>& data); //write a column
    void print(const string& name, const string& suffix) const;         //  writes matrix into name.suffix
    void print_append(const string& name,const string& suffix) const;  // the same, but appends matrix
    void print(std::ostream& out) const;          // writes matrix to the stream
    void pretty_print(std::ostream& out, bool with_row_nr=false) const;  // writes matrix in a nice format to the stream                   // read a row
    size_t nr_of_rows() const;                       // returns nr
    size_t nr_of_columns() const;                   // returns nc
    void set_nr_of_columns(size_t c);
    /* generates a pseudo random matrix for tests, entries form 0 to mod-1 */
    void random(int mod=3);

    void set_zero(); // sets all entries to 0

    /* returns a submatrix with rows corresponding to indices given by
     * the entries of rows, Numbering from 0 to n-1 ! */
    Matrix submatrix(const vector<key_t>& rows) const;
    Matrix submatrix(const vector<int>& rows) const;
    Matrix submatrix(const vector<bool>& rows) const;

    void swap (Matrix<Number>& x);

	// returns the permutation created by sorting the rows with a grading function
    // or by 1-norm if computed is false
    vector<key_t> perm_sort_by_degree(const vector<key_t>& key, const vector<Number>& grading, bool computed) const;
    vector<key_t> perm_by_weights(const Matrix<Number>& Weights, vector<bool> absolute);
    
    void select_submatrix(const Matrix<Number>& mother, const vector<key_t>& rows);
    void select_submatrix_trans(const Matrix<Number>& mother, const vector<key_t>& rows);

    Matrix& remove_zero_rows(); // remove zero rows, modifies this

    // resizes the matrix to the given number of rows/columns
    // if the size shrinks it will keep all its allocated memory
    // useful when the size varies
    void resize(size_t nr_rows);
    void resize(size_t nr_rows, size_t nr_cols);
    void resize_columns(size_t nr_cols);
    void Shrink_nr_rows(size_t new_nr_rows);

    vector<Number> diagonal() const;     //returns the diagonale of this
                                  //this should be a quadratic matrix
    size_t maximal_decimal_length() const;    //return the maximal number of decimals
                                      //needed to write an entry
                                      
    vector<size_t> maximal_decimal_length_columnwise() const; // the same per column

    void append(const Matrix& M); // appends the rows of M to this
    void append(const vector<vector<Number> >& M); // the same, but for another type of matrix
    void append(const vector<Number>& v); // append the row v to this
    void append_column(const vector<Number>& v); // append the column v to this
    void remove_row(const vector<Number>& row); // removes all appearances of this row, not very efficient!
    void remove_duplicate_and_zero_rows();

    inline const Number& get_elem(size_t row, size_t col) const {
        return elem[row][col];
    }
    inline const vector< vector<Number> >& get_elements() const {
        return elem;
    }
    inline vector<Number> const& operator[] (size_t row) const {
        return elem[row];
    }
    inline vector<Number>& operator[] (size_t row) { 
        return elem[row];
    }
    void set_nc(size_t col){
        nc=col;
    }
    void set_nr(size_t rows){
        nc=rows;
    }

//---------------------------------------------------------------------------
//                  Basic matrices operations
//---------------------------------------------------------------------------

    Matrix add(const Matrix& A) const;                       // returns this+A
    Matrix multiplication(const Matrix& A) const;          // returns this*A
    Matrix multiplication(const Matrix& A, long m) const;// returns this*A (mod m)
    Matrix<Number> multiplication_cut(const Matrix<Number>& A, const size_t& c) const; // returns 
    // this*(first c columns of A)
    bool equal(const Matrix& A) const;             // returns this==A
    bool equal(const Matrix& A, long m) const;     // returns this==A (mod m)
    Matrix transpose() const;                     // returns the transpose of this
    
    bool is_diagonal() const;

//---------------------------------------------------------------------------
//                          Scalar operations
//---------------------------------------------------------------------------

    void scalar_multiplication(const Number& scalar);  //this=this*scalar
    void scalar_division(const Number& scalar);
    //this=this div scalar, all the elem of this must be divisible with the scalar
    void reduction_modulo(const Number& modulo);     //this=this mod scalar
    Number matrix_gcd() const; //returns the gcd of all elem
    void simplify_rows();         //each row of this is reduced by its gcd, 
                                          // vector of gcds returned
    void make_cols_prime(size_t from_col, size_t to_col);   
             // the columns of this in the specified range are reduced by their gcd

    Matrix multiply_rows(const vector<Number>& m) const;  //returns matrix were row i is multiplied by m[i]

//---------------------------------------------------------------------------
//                          Vector operations
//---------------------------------------------------------------------------

   void MxV(vector<Number>& result, const vector<Number>& v) const;//result = this*V
   vector<Number> MxV(const vector<Number>& v) const;//returns this*V
   vector<Number> VxM(const vector<Number>& v) const;//returns V*this
   vector<Number> VxM_div(const vector<Number>& v, const Number& divisor,bool& success) const; // additionally divides by divisor

//---------------------------------------------------------------------------
//                          Matrix operations
//           --- these are more complicated algorithms ---
//---------------------------------------------------------------------------

// Normal forms

// converts this to row echelon form over ZZ and returns rank, GMP protected, uses only elementary transformations over ZZ
    size_t row_echelon();

    // public version of row_echelon_reduce, GMP protected, uses only elementary transformations over ZZ
    size_t row_echelon_reduce();

    // transforms matrix into lower triangular form via column transformations
    // assumes that rk is the rank and that the matrix is zero after the first rk rows
    // Right = Right*(column transformation of this call)
    bool column_trigonalize(size_t rk, Matrix<Number>& Right);
    
    // combines row_echelon_reduce and column_trigonalize
    // returns column transformation matrix
    Matrix<Number> row_column_trigonalize(size_t& rk, bool& success);
    
// rank and determinant

    size_t rank() const; //returns rank
    Number full_rank_index() const; // returns index of full rank sublattice
    size_t rank_submatrix(const vector<key_t>& key) const; //returns rank of submarix defined by key
    
    // returns rank of submatrix of mother. "this" is used as work space    
    size_t rank_submatrix(const Matrix<Number>& mother, const vector<key_t>& key);
 
    // vol stands for |det|
    Number vol() const;
    Number vol_submatrix(const vector<key_t>& key) const;
    Number vol_submatrix(const Matrix<Number>& mother, const vector<key_t>& key);
    
// find linearly indepenpendent submatrix of maximal rank

    vector<key_t>  max_rank_submatrix_lex() const; //returns a vector with entries
    //the indices of the first rows in lexicographic order of this forming
    //a submatrix of maximal rank.
    
// Solution of linear systems with square matrix
  
    // In the following routines, denom is the absolute value of the determinant of the
    // left side matrix.
    // If the diagonal is asked for, ZZ-invertible transformations are used.
    // Otherwise ther is no restriction on the used algorithm
    
    //The diagonal of left hand side after transformation into an upper triangular matrix
    //is saved in diagonal, denom is |determinant|.
    
    // System with "this" as left side
    Matrix solve(const Matrix& Right_side, Number& denom) const;
    Matrix solve(const Matrix& Right_side, vector< Number >& diagonal, Number& denom) const;
    // solve the system this*Solution=denom*Right_side. 

    // system is defined by submatrix of mother given by key (left side) and column vectors pointed to by RS (right side)
    // NOTE: this is used as the matrix for the work     
    void solve_system_submatrix(const Matrix& mother, const vector<key_t>& key, const vector<vector<Number>* >& RS,
         vector< Number >& diagonal, Number& denom, size_t red_col, size_t sign_col);
    void solve_system_submatrix(const Matrix& mother, const vector<key_t>& key, const vector<vector<Number>* >& RS,
         Number& denom, size_t red_col, size_t sign_col, bool compute_denom=true, bool make_sol_prime=false);
    // the left side gets transposed
    void solve_system_submatrix_trans(const Matrix& mother, const vector<key_t>& key, const vector<vector<Number>* >& RS,
         Number& denom, size_t red_col, size_t sign_col);
        
                    
// For non-square matrices
                    
    // The next two solve routines do not require the matrix to be square.
    // However, we want rank = number of columns, ensuring unique solvability
    
    vector<Number> solve_rectangular(const vector<Number>& v, Number& denom) const;
    // computes solution vector for right side v, solution over the rationals
    // matrix needs not be quadratic, but must have rank = number of columns
    // with denominator denom. 
    // gcd of denom and solution is extracted !!!!!
    
    vector<Number> solve_ZZ(const vector<Number>& v) const;
    // computes solution vector for right side v
    // insists on integrality of the solution
                    
// homogenous linear systems

    Matrix<Number> kernel () const;
    // computes a ZZ-basis of the solutions of (*this)x=0
    // the basis is formed by the ROWS of the returned matrix
                    
// inverse matrix
                    
    //this*Solution=denom*I. "this" should be a quadratic matrix with nonzero determinant. 
    Matrix invert(Number& denom) const;
    
    void invert_submatrix(const vector<key_t>& key, Number& denom, Matrix<Number>& Inv, 
                bool compute_denom=true, bool make_sol_prime=false) const;
                    
// find linear form that is constant on the rows 

    vector<Number> find_linear_form () const;
    // Tries to find a linear form which gives the same value an all rows of this
    // this should be a m x n matrix (m>=n) of maxinal rank
    // returns an empty vector if there does not exist such a linear form
  
    vector<Number> find_linear_form_low_dim () const;
    //same as find_linear_form but also works with not maximal rank
    //uses a linear transformation to get a full rank matrix
    
// normal forms
        
    Matrix AlmostHermite(size_t& rk);
    // Converts "this" into lower trigonal column Hermite normal form, returns column 
    // transformation matrix
    // Almost: elements left of diagonal are not reduced mod diagonal 
    
    // Computes Smith normal form and returns column transformation matrix
    Matrix SmithNormalForm(size_t& rk);
    
//for simplicial subcones

    // computes support hyperplanes and volume
    void simplex_data(const vector<key_t>& key, Matrix<Number>& Supp, Number& vol, bool compute_vol) const; 
    
// Sorting of rows
    
    Matrix& sort_by_weights(const Matrix<Number>& Weights, vector<bool> absolute);
    Matrix& sort_lex();
    void order_rows_by_perm(const vector<key_t>& perm);
    
// solve homogeneous congruences
    
    Matrix<Number> solve_congruences(bool& zero_modulus) const;
    
// saturate sublattice
    
    void saturate();

// find the indices of the rows in which the linear form L takes its max and min values
    
    vector<key_t> max_and_min(const vector<Number>& L, const vector<Number>& norm) const;
    
// try to sort the rows in such a way that the extreme points of the polytope spanned by the rows come first
    
    size_t extreme_points_first(const vector<Number> norm=vector<Number>(0));

// find an inner point in the cone spanned by the rows of the matrix
    
    vector<Number> find_inner_point();
};
//class end *****************************************************************

//---------------------------------------------------------------------------
//                  Utilities
//---------------------------------------------------------------------------

template<typename Number> class order_helper {
    
public:
    
    vector<Number> weight;
    key_t index;
    vector<Number>* v;
};

template<typename Number>
vector<vector<Number> > to_matrix(const vector<Number>& v){
    
    vector<vector<Number> > mat(1);
    mat[0]=v;
    return mat;    
}

template<typename Number>
Matrix<Number>  readMatrix(const string project);

//---------------------------------------------------------------------------
//                  Conversion between integer types
//---------------------------------------------------------------------------

template<typename ToType, typename FromType>
void convert(Matrix<ToType>& to_mat, const Matrix<FromType>& from_mat);

template<typename Number>
void mat_to_mpz(const Matrix<Number>& mat, Matrix<mpz_class>& mpz_mat);

template<typename Number>
void mat_to_Int(const Matrix<mpz_class>& mpz_mat, Matrix<Number>& mat);

template<typename Number>
void mpz_submatrix(Matrix<mpz_class>& sub, const Matrix<Number>& mother, const vector<key_t>& selection);

template<typename Number>
void mpz_submatrix_trans(Matrix<mpz_class>& sub, const Matrix<Number>& mother, const vector<key_t>& selection);

} // namespace

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------
