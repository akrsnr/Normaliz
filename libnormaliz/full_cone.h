/*
 * Normaliz 2.5
 * Copyright (C) 2007-2010  Winfried Bruns, Bogdan Ichim, Christof Soeger
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
 */

#ifndef FULL_CONE_H
#define FULL_CONE_H

#include <set>
#include <list>
#include <bitset>

#include "libnormaliz.h"
#include "integer.h"
#include "matrix.h"
#include "simplex.h"
#include "cone_dual_mode.h"

namespace libnormaliz {

template<typename Integer> class Cone;
template<typename Integer> class Full_Cone;

//generates a lifted cone with the lower part simplicial, needed for computing the triangulation by lifting
template<typename Integer>
void lift(Full_Cone<Integer>& Lifted, Matrix<Integer> Extreme_Generators);


template<typename Integer>
class Full_Cone {
	int dim;
	int nr_gen;
	int hyp_size;
	
	bool pointed;
	bool ht1_generated;
	bool ht1_extreme_rays;
	bool ht1_hilbert_basis;
	bool integrally_closed;
	bitset<ConeProperty::EnumSize> is_Computed;
	vector<Integer> Linear_Form;
	Integer multiplicity;
	Matrix<Integer> Generators;
	vector<bool> Extreme_Rays;
	list<vector<Integer> > Support_Hyperplanes;
	list< Simplex<Integer> > Triangulation;
	list<vector<Integer> > Hilbert_Basis;
	list<vector<Integer> > Homogeneous_Elements;
	vector<Integer> H_Vector;
	vector<Integer> Hilbert_Polynomial;

	friend void lift<Integer>(Full_Cone<Integer>&, Matrix<Integer>);
	friend class Cone<Integer>;

/* ---------------------------------------------------------------------------
 *				Private routines, used in the public routines
 * ---------------------------------------------------------------------------
 */
	void add_hyperplane(const int & size, const vector<Integer> & positive_gen, const vector<Integer> & negative_gen);
	void transform_values(const int & size, const vector<int> & test_key);
	void add_simplex(const int & new_generator, const int & size, const vector<int> & col, const vector<int> & col_inv);

	/* adds a new element to the Hilbert basis */
	void reduce_and_insert(const vector<Integer> & new_element);
	/* adds a new element to the Hilbert basis
	 * faster as above, provided the scalar products are precomputed, which needs more memory */
	void reduce_and_insert_speed(const vector<Integer> & new_element);

	/* Returns true if new_element is reducible versus the elements in Ired */
	bool is_reducible(list<vector<Integer> > & Ired, const vector<Integer> & new_element);
	bool is_reducible(list<vector<Integer> *> & Ired, const vector<Integer> & new_element);

	/* reduce Red versus Ired */
	void reduce(list<vector<Integer> > & Ired, list<vector<Integer> > & Red, const int & size);

	/* adds a matrix with new elements to the Hilbert basis */
	void reduce_and_insert(const Matrix<Integer> & New_Elements);
	/* adds a list with new elements to the Hilbert basis */
	void reduce_and_insert(const list<vector<Integer> > & New_Elements);


	/* to be used with a shelling in order to add to each simples the maximal new face */
	void find_new_face();
	
	/* */
	Simplex<Integer> find_start_simplex() const;

	/* compute triangulations of the not compressed, not simplicial pieces and add them to Triangulation*/
	void process_non_compressed(list<vector<int> > & non_compressed);
	/* reduce the Candidates against itself and stores the remaining elements in Hilbert_Basis */
	void global_reduction(set<vector<Integer> > & Candidates);
	/* computes a degree function, s.t. every generator has value >0 */
	vector<Integer> compute_degree_function() const;

	void compute_support_hyperplanes(const bool do_partial_triang = false);
	void compute_support_hyperplanes_triangulation();
	void support_hyperplanes_partial_triangulation();
	void compute_support_hyperplanes_pyramid(const bool do_triang = false);
	void support_hyperplane_common();
	void compute_extreme_rays();
	void compute_hilbert_basis();
	void compute_ht1_elements();
	void compute_hilbert_polynomial();
	void compute_hilbert_basis_polynomial();

	void check_pointed();
	void check_ht1_generated();
	void check_ht1_extreme_rays();
	void check_ht1_hilbert_basis();
	void check_integrally_closed();

	void compute_multiplicity();
	bool low_part_simplicial();
	void line_shelling();
	void triangulation_lift();
	/* computes the e vector using the h vector */
	vector<Integer> compute_e_vector();
	/* computes the Hilbert polynomial using the h-vector */
	void compute_polynomial();

	/* support hyperplanes computation for a dynamic lifting
	 * adjusts the lifting if necessary, used in dual algorithm */
	void support_hyperplanes_dynamic();


    /* constructor for recursively generated subcones
     * int i is a dummy parameter to distinguish it from the standard constructor */
    Full_Cone(Matrix<Integer> M, int i);

public:
	Full_Cone();
	Full_Cone(Matrix<Integer> M);            //main constructor
	Full_Cone(const Cone_Dual_Mode<Integer> &C);
	Full_Cone(const Full_Cone<Integer> & C); //copy constructor
	~Full_Cone();                   //destructor

/*---------------------------------------------------------------------------
 *						Data access
 *---------------------------------------------------------------------------
 */
	void print() const;                //to be modified, just for tests
	int getDimension() const;        //returns dimension
	int getNrGenerators() const;    //returns the number of generators
	bool isPointed() const;
	bool isHt1ExtremeRays() const;     //returns homogeneous
	bool isHt1HilbertBasis() const;
	bool isIntegrallyClosed() const;
	vector<Integer> getLinearForm() const; //returns the linear form
	Integer getMultiplicity() const; //returns multiplicity
	const Matrix<Integer>& getGenerators() const;
	vector<bool> getExtremeRays() const;
	Matrix<Integer> getSupportHyperplanes() const;
	Matrix<Integer> getTriangulation() const;
	/* read the triangulation and the volume of each simplex,
	 * the volume is saved on the last column
	 * the vectors corresponding to the generators of each simplex are sorted */
	Matrix<Integer> getTriangulationVolume() const;
	Matrix<Integer> getHilbertBasis() const;
	Matrix<Integer> getHt1Elements() const;
	vector<Integer> getHVector() const;
	vector<Integer> getHilbertPolynomial() const;
	
	bool isComputed(ConeProperty::Enum prop) const; 


/*---------------------------------------------------------------------------
 *				Computation Methods
 *---------------------------------------------------------------------------
 */
	void support_hyperplanes();
	void support_hyperplanes_pyramid();
	void support_hyperplanes_triangulation();
	void support_hyperplanes_triangulation_pyramid();
	void triangulation_hilbert_basis();
	void hilbert_basis();
	void hilbert_polynomial();
	void hilbert_basis_polynomial();
	void ht1_elements();

	/* computes the multiplicity of the ideal in case of a Rees algebra
	 * (not the same as the multiplicity of the semigroup) */
	Integer primary_multiplicity() const;

	/* completes the computation when a Cone_Dual_Mode is given */
	void dual_mode();

	/* checks if the cone is compressed, support hyperplanes must be computed */
	bool check_compressed();

	void error_msg(string s) const;
};
//class end *****************************************************************
//---------------------------------------------------------------------------

}

//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------------