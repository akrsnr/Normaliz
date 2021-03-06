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

#ifndef CONE_PROPERTY_H_
#define CONE_PROPERTY_H_

#include <bitset>
#include <ostream>

namespace libQnormaliz {

/* An enumeration of things, that can be computed for a cone.
 * The namespace prevents interfering with other names.
 * Remember to change also the string conversion if you change this enum!
 */
namespace ConeProperty {
    enum Enum {
        //
        // goals that can be computed (or are defined by input data)
        //
        // matrix valued
        Generators,
        ExtremeRays,
        VerticesOfPolyhedron,
        SupportHyperplanes,
        HilbertBasis,
        ModuleGenerators,
        Deg1Elements,
        ModuleGeneratorsOverOriginalMonoid,
        Sublattice, 
        ExcludedFaces,
        OriginalMonoidGenerators,
        MaximalSubspace,
        Equations, // new
        Congruences, // new
        //vector valued
        Grading,
        Dehomogenization,
        WitnessNotIntegrallyClosed,
        // Cardinalities
        TriangulationSize,
        // Number valued,        
        TriangulationDetSum,
        ReesPrimaryMultiplicity,
        GradingDenom, // new
        UnitGroupIndex, // new
        InternalIndex, // new
        ExternalIndex, // new
        // rational valued
        Multiplicity,
        // dimensions
        RecessionRank,
        AffineDim,
        ModuleRank,
        Rank, // new
        EmbeddingDim, // new      
        // boolean valued 
        IsPointed,
        IsDeg1ExtremeRays,
        IsDeg1HilbertBasis,
        IsIntegrallyClosed,
        IsReesPrimary,
        IsInhomogeneous, // new        
        // complex structures
        Triangulation,
        HilbertSeries,
        InclusionExclusionData,
        StanleyDec,        
        ClassGroup,        
        NumberHull,
        ConeDecomposition,
        HilbertQuasiPolynomial,
                //
        // integer type for computations
        //
        BigInt,
        //
        // algorithmic variants
        //
        DefaultMode,
        Approximate,
        BottomDecomposition,
        NoBottomDec,       
        DualMode,
        PrimalMode, //new
        Symmetrize, // new
        NoSymmetrization, // new
        KeepOrder,
        HSOP,
        //
        // checking properties of already computed data
        // (cannot be used as a computation goal)
        //
        IsTriangulationNested,  //new
        IsTriangulationPartial,  //new
        
        EnumSize // this has to be the last entry, to get the number of entries in the enum
    }; // remember to change also the string conversion function if you change this enum
}

class ConeProperties {
public:
    /* Constructors */
    ConeProperties();
    ConeProperties(ConeProperty::Enum);
    ConeProperties(ConeProperty::Enum, ConeProperty::Enum);
    ConeProperties(ConeProperty::Enum, ConeProperty::Enum, ConeProperty::Enum);
    ConeProperties(const std::bitset<ConeProperty::EnumSize>&);

    /* set properties */
    ConeProperties& set(ConeProperty::Enum, bool value=true);
    ConeProperties& set(const std::string s, bool value=true);
    ConeProperties& set(ConeProperty::Enum, ConeProperty::Enum);
    ConeProperties& set(ConeProperty::Enum, ConeProperty::Enum, ConeProperty::Enum);
    ConeProperties& set(const ConeProperties&);

    /* reset (=unset) properties */
    ConeProperties& reset(ConeProperty::Enum Property);
    ConeProperties& reset(const ConeProperties&);
    ConeProperties& reset_compute_options();

    /* test which/how many properties are set */
    bool test(ConeProperty::Enum Property) const;
    bool any() const;
    bool none() const;
    size_t count () const;

    /* return the restriction of this to the goals / options */
    ConeProperties goals();
    ConeProperties options();

    /* the following methods are used internally */
    void set_preconditions();    // activate properties which are needed implicitly
    void prepare_compute_options(bool inhomogeneous);
    void check_sanity(bool inhomogeneous);
    void check_Q_permissible();

    /* print it in a nice way */
    friend std::ostream& operator<<(std::ostream&, const ConeProperties&);


private:
    std::bitset<ConeProperty::EnumSize> CPs;

};

// conversion to/from strings
bool isConeProperty(ConeProperty::Enum& cp, const std::string& s);
ConeProperty::Enum toConeProperty(const std::string&);
const std::string& toString(ConeProperty::Enum);
std::ostream& operator<<(std::ostream&, const ConeProperties&);

}

#endif /* CONE_PROPERTY_H_ */
