(* ::Package:: *)

(* copy this module to a directory in $Path.  Then invoke with <<GA20` *)
BeginPackage[ "GA20`" ]

(* Must reference any global symbol (or some of them) before Unprotecting it, since it may not have
   been loaded:

   http://mathematica.stackexchange.com/a/137007/10
 *)
{D, TraditionalForm, DisplayForm, StandardForm, Format};

Unprotect[
   Bivector,
   BivectorSelection,
   DisplayForm,
   Format,
   GradeSelection,
   Scalar,
   ScalarProduct,
   ScalarSelection,
   ScalarValue,
   StandardForm,
   TeXForm ,
   TraditionalForm,
   Vector,
   VectorSelection,
   e
]

ClearAll[
   Bivector,
   BivectorSelection,
   GAdisplay,
   GradeSelection,
   Scalar,
   ScalarProduct,
   ScalarSelection,
   ScalarValue,
   Vector,
   VectorSelection,
   antisymmetric,
   binaryOperator,
   bivectorQ,
   bladeQ,
   bold,
   directProduct,
   displayMapping,
   esub,
   grade,
   oneTeXForm,
   pauliGradeSelect,
   pauliGradeSelect0,
   pauliGradeSelect1,
   pauliGradeSelect2,
   pmagnitude,
   scalarQ,
   signedSymmetric,
   symmetric,
   vectorQ
]

GA20::usage = "GA20: An implementation of Euclidean (CL(2,0)) Geometric Algebra.

Pauli matrices are used to represent the algebraic elements.  This provides an efficient and compact representation
of the entire algebraic space.

Internally, a multivector is represented by a pair (grade, pauli-representation).  The grade portion will be
obliterated when adding objects that have different grade, or multiplying vectors or bivectors.  When
it is available, certain operations can be optimized.  Comparison ignores the cached grade if it exists.

Elements of the algebra can be constructed with one of

   Scalar[ v ]
   Vector[ v, n ]
   Bivector[ v ]

Example:

   m = Scalar[ Sin[ x ] ] + Vector[ Log[ z ], 3 ]
   m // StandardForm

> e[ 3 ] Log[ z ] + Sin[ x ]

A few operators are provided:
   ==         Compare two multivectors, ignoring the cached grade if any.
   m1 + m2
   m1 - m2
   - m
   st * vb    Scalars can multiply vectors and bivectors in any order
   vb1 ** vb1 Vectors and bivectors when multiplied have to use the NonCommutativeMultiply operator, but any grade object may also.
   m1 . m2    Dot product.  The functional form Dot[ m1, m2 ] may also be used.
   m1 ^ m2   Wedgeproduct.  Enter with m1 [ Esc ]^[ Esc ] m2.  The functional form Wedge[ m1, m2 ]
   <m>        Scalar selection.  Enter with [ Esc ]<[ Esc ] m [ Esc ]>[ Esc ].  The functional form ScalarValue[ m ] may also be used.  This returns the numeric (or expression) value of the scalar grade of the multivector, and not a grade[ ] object.
   <m1,m2>    Scalar product.  Enter with [ Esc ]<[ Esc ] m1,m2 [ Esc ]>[ Esc ].  The functional form ScalarProduct[ m1, m2 ] may also be used.  This returns the numeric (or expression) value of the scalar product of the multivectors, and not a grade[ ] object.

   Functions provided:

   - GradeSelection
   - ScalarSelection
   - VectorSelection
   - BivectorSelection
   - ScalarValue, < m >
   - ScalarProduct, < m1, m2 >

The following built-in methods are overridden:

   - TraditionalForm
   - DisplayForm
   - StandardForm
   - TeXForm

Internal functions:

   - scalarQ
   - vectorQ
   - bivectorQ
   - bladeQ
   - gradeAnyQ
   - notGradeQ

TODO:

1) How to get better formatted output by default without using one of TraditionalForm, DisplayForm, StandardForm ?

2) Can a package have options (i.e. to define the name of the e[ ] operator used in StandardForm that represents a basis vector).

3) proper packaging stuff:  private for internals.
" ;

grade::usage = "grade.  (internal) An upvalue type that represents a CL(2,0) algebraic element as a pair {grade, v}, where v is a sum of products of Pauli matrices.  These matrices may be scaled by arbitrary numeric or symbolic factors." ;
Scalar::usage = "Scalar[ v ] constructs a scalar grade quantity with value v." ;
Scalar[ v_ ] := grade[ 0, v IdentityMatrix[ 2 ] ] ;
Vector::usage = "Vector[ v, n ], where n = {1,2} constructs a vector grade quantity with value v in direction n." ;
Vector[ v_, k_Integer /; k == 1 ] := grade[ 1, v PauliMatrix[ 1 ] ] ;
Vector[ v_, k_Integer /; k == 2 ] := grade[ 1, v PauliMatrix[ 3 ] ] ;
Bivector::usage = "Bivector[ v ], constructs a bivector grade quantity with value v in the plane e1,e2." ;
Bivector[ v_ ] := grade[ 2, v PauliMatrix[ 1 ] . PauliMatrix[ 3 ] ] ;

(*Begin[ "`Private`" ]*)
gradeQ::usage = "gradeQ[ m, n ] tests if the multivector m is of grade n.  n = -1 is used internally to represent values of more than one grade." ;
gradeQ[ m_grade, n_Integer ] := ((m // First) == n)
scalarQ::usage = "scalarQ[ m ] tests if the multivector m is of grade 0 (scalar)" ;
scalarQ[ m_grade ] := gradeQ[ m, 0 ]
vectorQ::usage = "vectorQ[ m ] tests if the multivector m is of grade 1 (vector)" ;
vectorQ[ m_grade ] := gradeQ[ m, 1 ]
bivectorQ::usage = "bivectorQ[ m ] tests if the multivector m is of grade 2 (bivector)" ;
bivectorQ[ m_grade ] := gradeQ[ m, 2 ]
bladeQ::usage = "bladeQ[ m ] tests if the multivector is of a single grade." ;
bladeQ[ m_grade ] := ((m // First) >= 0)
gradeAnyQ::usage = "gradeAnyQ[ ].  predicate pattern match for grade[ _ ]" ;
gradeAnyQ[ m_grade ] := True
gradeAnyQ[ _ ] := False
notGradeQ::usage = "notGradeQ[ ].  predicate pattern match for !grade[ ]" ;
notGradeQ[ v_ ] := Not[ gradeAnyQ[ v ] ]

directProduct[ t_, v1_, v2_ ] := grade[ t, (v1 // Last).(v2 // Last) ] ;
signedSymmetric[ t_, v1_, v2_, s_ ] :=
  Module[ {a = (v1 // Last), b = (v2 // Last)},
   grade[ t, (a.b + s b.a)/2 ] ] ;
symmetric[ t_, v1_, v2_ ] := signedSymmetric[ t, v1, v2, 1 ] ;
antisymmetric[ t_, v1_, v2_ ] := signedSymmetric[ t, v1, v2, -1 ] ;

(*These operator on just the Pauli matrix portions x of \
pauliGradeSelect[ ,x ]*)
pauliGradeSelect[ {{a_,  _}, { _, d_}}, 0 ] := ( (1/2) (a + d) IdentityMatrix[ 2 ] )
pauliGradeSelect[ {{a_, b_}, {c_, d_}}, 1 ] := ( (1/2) (b + c) PauliMatrix[ 1 ] + (1/2)(a - d) PauliMatrix[ 3 ] )
pauliGradeSelect[ {{a_, b_}, {c_, d_}}, 2 ] := ( (1/2) (c - b) PauliMatrix[ 1 ].PauliMatrix[ 3 ] )


pauliGradeSelect0 := pauliGradeSelect[ #, 0 ] & ;
pauliGradeSelect1 := pauliGradeSelect[ #, 1 ] & ;
pauliGradeSelect2 := pauliGradeSelect[ #, 2 ] & ;
(*End[ "`Private`" ]*)

GradeSelection::usage = "GradeSelection[ m, k ] selects the grade k elements from the multivector m.  The selected result is represented internally as a grade[ ] type (so scalar selection is not just a number)." ;
GradeSelection[ m_?scalarQ, 0 ] := m ;
GradeSelection[ m_?vectorQ, 1 ] := m ;
GradeSelection[ m_?bivectorQ, 2 ] := m ;
GradeSelection[ m_, k_Integer /; k >= 0 && k <= 2 ] := grade[ k, pauliGradeSelect[ m // Last, k ] ] ;
ScalarSelection::usage = "ScalarSelection[ m ] selects the grade 0 (scalar) elements from the multivector m.  The selected result is represented internally as a grade[ ] type (not just a number or an expression)." ;
ScalarSelection := GradeSelection[ #, 0 ] & ;
VectorSelection::usage = "VectorSelection[ m ] selects the grade 1 (vector) elements from the multivector m.  The selected result is represented internally as a grade[ ] type." ;
VectorSelection := GradeSelection[ #, 1 ] & ;
BivectorSelection::usage = "BivectorSelection[ m ] selects the grade 2 (bivector) elements from the multivector m.  The selected result is represented internally as a grade[ ] type." ;
BivectorSelection := GradeSelection[ #, 2 ] & ;

binaryOperator[ f_, b_?bladeQ, m_grade ] := Total[ f[ b, # ] & /@ (GradeSelection[ m, # ] & /@ (Range[ 2+1 ] - 1)) ]
binaryOperator[ f_, m_grade, b_?bladeQ ] := Total[ f[ #, b ] & /@ (GradeSelection[ m, # ] & /@ (Range[ 2+1 ] - 1)) ]
binaryOperator[ f_, m1_grade, m2_grade ] := Total[ f[ # // First, # // Last ] & /@ (
    {GradeSelection[ m1, # ] & /@ (Range[ 2+1 ] - 1),
     GradeSelection[ m2, # ] & /@ (Range[ 2+1 ] - 1)} // Transpose) ]

(* Plus *)
grade /: (v1_?notGradeQ) + grade[ k_, v2_ ] := Scalar[ v1 ] + grade[ k, v2 ] ;
grade /: grade[ 0, v1_ ] + grade[ 0, v2_ ] := grade[ 0, v1 + v2 ] ;
grade /: grade[ 1, v1_ ] + grade[ 1, v2_ ] := grade[ 1, v1 + v2 ] ;
grade /: grade[ 2, v1_ ] + grade[ 2, v2_ ] := grade[ 2, v1 + v2 ] ;
grade /: grade[ _, v1_ ] + grade[ _, v2_ ] := grade[ -1, v1 + v2 ] ;

(* Times[ -1, _ ] *)
grade /: -grade[ k_, v_ ] := grade[ k, -v ] ;

(* Times *)
grade /: (v_?notGradeQ) grade[ k_, m_ ] := grade[ k, v m ] ;
grade /: grade[ 0, s_ ] grade[ k_, m_ ] := grade[ k, s.m ] ;

(* NonCommutativeMultiply *)
grade /: grade[ 0, s_ ] ** grade[ k_, m_ ] := grade[ k, s.m ] ;
grade /: grade[ k_, m_ ] ** grade[ 0, s_ ] := grade[ k, s.m ] ;
grade /: grade[ _, m1_ ] ** grade[ _, m2_ ] := grade[ -1, m1.m2 ] ;

(* Dot *)
grade /: (s_?notGradeQ).grade[ k_, m_ ] := grade[ k, s m ] ;
grade /: grade[ k_, m_ ].(s_?notGradeQ) := grade[ k, s m ] ;
grade /: grade[ 0, s_ ].grade[ k_, m_ ] := grade[ k, s m ] ;
grade /: grade[ k_, m_ ].grade[ 0, s_ ] := grade[ k, s m ] ;

grade /: (v1_?vectorQ).grade[ 1, v2_ ] := symmetric[ 0, v1, grade[ 1, v2 ] ] ;
grade /: (v_?vectorQ).grade[ 2, b_ ] := antisymmetric[ 1, v, grade[ 2, b ] ] ;
grade /: (b_?bivectorQ).grade[ 1, v_ ] := antisymmetric[ 1, b, grade[ 1, v ] ] ;
grade /: (b1_?bivectorQ).grade[ 2, b2_ ] := symmetric[ 0, b1, grade[ 2, b2 ] ] ;

grade /: grade[ g1_, m1_ ] . grade[ g2_, m2_ ]:= binaryOperator[ Dot, grade[ g1, m1 ], grade[ g2, m2 ] ] ;

(* == comparison operator *)
grade /: grade[ _, m1_ ] == grade[ _, m2_ ] := (m1 == m2) ;

grade[ _, {{0, 0}, {0, 0}} ] := 0

(*Define a custom wedge operator*)

grade /: grade[ 0, s_ ] \[Wedge] grade[ k_, v_ ] := grade[ k, s.v ] ;
grade /: grade[ k_, v_ ] \[Wedge] grade[ 0_, s_ ] := grade[ k, s.v ] ;
grade /: grade[ 1, v1_ ] \[Wedge] (v2_?vectorQ) := antisymmetric[ 2, grade[ 1, v1 ], v2 ] ;

(* Only e12 ^ scalar is non zero, and that is handled above *)
grade /: grade[ 2, _ ] \[Wedge] b_?bladeQ := 0 ;
grade /: b_?bladeQ \[Wedge] grade[ 2, _ ] := 0 ;

grade /: grade[ g1_, m1_ ] \[Wedge] grade[ g2_, m2_ ]:= binaryOperator[ Wedge, grade[ g1, m1 ], grade[ g2, m2 ] ] ;

(*Begin[ "`Private`" ]*)
pmagnitude::usage =
  "pmagnitude[ ].  select the 1,1 element from a pauli matrix assuming it represents \
a Scalar (i.e. scaled diagonal matrix)." ;
pmagnitude[ m_ ] := m[ [1, 1 ] ] ;
(*End[ "`Private`" ]*)

(* AngleBracket,single operand forms, enter with[ Esc ]<[ Esc ] \
v[ Esc ]>[ Esc ] *)
grade /: AngleBracket[ grade[ 0, s_ ] ] := pmagnitude[ s ]
grade /: AngleBracket[ grade[ 1, _ ]  ] := 0
grade /: AngleBracket[ grade[ 2, _ ]  ] := 0
grade /: AngleBracket[ grade[ _, m_ ] ] := ((pauliGradeSelect[ m, 0 ]) // pmagnitude)

ScalarValue::usage = "ScalarValue[ m ].  Same as AngleBracket[ m ], aka [ Esc ]<[ Esc ] m1 [ Esc ]>[ Esc ]." ;
ScalarValue[ m_grade ] := AngleBracket[ m ] ;

(* AngleBracket,two operand forms. *)

grade /: AngleBracket[ grade[ 0, s1_ ], grade[ 0, s2_ ] ] := (pmagnitude[ s1 ] pmagnitude[ s2 ]) ;
grade /: AngleBracket[ grade[ 0, s1_ ], grade[ -1, m_ ] ] := (pmagnitude[ s1 ] ((pauliGradeSelect[ m, 0 ]) // pmagnitude)) ;
grade /: AngleBracket[ grade[ -1, m_ ], grade[ 0, s1_ ] ] := (pmagnitude[ s1 ] ((pauliGradeSelect[ m, 0 ]) // pmagnitude)) ;

grade /: AngleBracket[ grade[ 0, s1_ ], grade[ _, _ ] ] := 0 ;
grade /: AngleBracket[ grade[ _, _ ], grade[ 0, s1_ ] ] := 0 ;

grade /: AngleBracket[ grade[ k1_, m1_ ], grade[ k2_, m2_ ] ] := (pauliGradeSelect[ m1.m2, 0 ] // pmagnitude) ;

ScalarProduct::usage = "ScalarProduct[ ].  Same as AngleBracket[ m1, m2 ], aka [ Esc ]<[ Esc ] m1, m2 [ Esc ]>[ Esc ]." ;
ScalarProduct[ m1_grade, m2_grade ] := AngleBracket[ m1, m2 ] ;

(*Begin[ "`Private`" ]*)
bold = Style[ #, Bold ] & ;
esub = Subscript[ bold[ "e" ], # ] & ;
displayMapping = {
   {Scalar[ 1 ], 1, 1, 1},
   {Vector[ 1, 1 ], esub[ 1 ], e[ 1 ], "\\mathbf{e}_1"},
   {Vector[ 1, 2 ], esub[ 2 ], e[ 2 ], "\\mathbf{e}_2"},
   {Bivector[ -1 ], esub[ "12" ], e[ 1 ]e[ 2 ], "\\mathbf{e}_{12}"}
} ;

GAdisplay[ v_grade, how_ ] :=
  Total[ (Times[ AngleBracket[ # // First, v ], #[ [how ] ] ]) & /@
    displayMapping ] ;
(*End[ "`Private`" ]*)

TraditionalForm[ m_grade ] := ((GAdisplay[ m, 2 ]) // TraditionalForm) ;
DisplayForm[ m_grade ] := GAdisplay[ m, 2 ] ;
Format[ m_grade ] := GAdisplay[ m, 2 ] ;
StandardForm[ m_grade ] := GAdisplay[ m, 3 ] ;

oneTeXForm[m_, blade_, disp_] := Module[{p},
  p = AngleBracket[blade, m];
  If[ p // PossibleZeroQ, 0,
   If[ p === 1, disp,
    Times[p // TeXForm, disp]]
   ]
  ]

TeXForm[m_grade] := Total[ oneTeXForm[m, # // First, #[[4]]] & /@ displayMapping ];

Protect[
   TraditionalForm,
   DisplayForm,
   StandardForm,
   Format,
   TeXForm,
   Scalar,
   Vector,
   Bivector,
   GradeSelection,
   ScalarSelection,
   VectorSelection,
   BivectorSelection,
   e,
   ScalarValue,
   ScalarProduct
]

EndPackage[ ]
