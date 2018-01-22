(* ::Package:: *)

(* copy this module to a directory in $Path.  Then invoke with <<GA30` *)
BeginPackage[ "GA30`" ]

(* Must reference any global symbol (or some of them) before Unprotecting it, since it may not have
   been loaded:

   http://mathematica.stackexchange.com/a/137007/10
 *)
{D, TraditionalForm, DisplayForm, StandardForm, Grad, Div, Curl, Format};

Unprotect[
   Bivector,
   BivectorSelection,
   Curl,
   D,
   DisplayForm,
   Div,
   Grad,
   GradeSelection,
   Scalar,
   ScalarProduct,
   ScalarSelection,
   ScalarValue,
   StandardForm,
   TeXForm,
   TraditionalForm,
   Trivector,
   TrivectorSelection,
   Vcurl,
   Vector,
   VectorSelection,
   complex,
   complexI,
   complexQ,
   conjugate,
   e,
   fMatrix,
   imag,
   matrixconj,
   matriximag,
   matrixreal,
   notComplexQ,
   real
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
   Trivector,
   TrivectorSelection,
   Vector,
   VectorSelection,
   antisymmetric,
   binaryOperator,
   bivectorQ,
   bladeQ,
   bold,
   complex,
   complexI,
   complexQ,
   conjugate,
   conjugateTranspose,
   directProduct,
   displayMapping,
   esub,
   fMatrix,
   grade,
   imag,
   matrixconj,
   matriximag,
   matrixreal,
   notComplexQ,
   oneTeXForm,
   pauliGradeSelect,
   pauliGradeSelect0,
   pauliGradeSelect01,
   pauliGradeSelect1,
   pauliGradeSelect2,
   pauliGradeSelect23,
   pauliGradeSelect3,
   pauliMatrix,
   pmagnitude,
   real,
   scalarQ,
   signedSymmetric,
   symmetric,
   trivectorQ,
   vectorQ
]

GA30::usage = "GA30: An implementation of Euclidean (CL(3,0)) Geometric Algebra.

Pauli matrices are used to represent the algebraic elements.  This provides an efficient and compact representation
of the entire algebraic space.

Internally, a multivector is represented by a pair (grade, pauli-representation).  The grade portion will be
obliterated when adding objects that have different grade, or multiplying vectors or bivectors.  When
it is available, certain operations can be optimized.  Comparison ignores the cached grade if it exists.

Elements of the algebra can be constructed with one of

   Scalar[ v ]
   Vector[ v, n ]
   Bivector[ v, n, m ]
   Trivector[ v ]

Example:

   m = Scalar[ Sin[ x ] ] + Vector[ Log[ z ], 3 ] + Trivector[ 7 ] ;
   m // StandardForm

> 7 e[ 123 ] + e[ 3 ] Log[ z ] + Sin[ x ]

A few operators are provided:
   ==         Compare two multivectors, ignoring the cached grade if any.
   m1 + m2
   m1 - m2
   - m
   st * vb    Scalars and trivectors can multiply vectors and bivectors in any order
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
   - TrivectorSelection
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
   - trivectorQ
   - bladeQ
   - gradeAnyQ
   - notGradeQ

TODO:

1) How to get better formatted output by default without using one of TraditionalForm, DisplayForm, StandardForm ?

2) Can a package have options (i.e. to define the name of the e[ ] operator used in StandardForm that represents a basis vector).

3) proper packaging stuff:  private for internals.
" ;
complex::usage =
  "complex.  A limited use complex number implementation to use internally in \
a Pauli or Dirac matrix basis representation, independent of any Complex" ;
norm::usage = "norm[ z ].  A Norm like function for complex[ ]" ;
complexQ::usage = "complexQ[ z ].  predicate pattern match for complex[ ]" ;
notComplexQ::usage = "notComplexQ[ z ].  predicate pattern match for !complex[ ]" ;
real::usage = "real[ z ].  Re[ z ] like function for complex[ ]" ;
imag::usage = "imag[ z ].  Im[ z ] like function for complex[ ]" ;
conjugate::usage = "conjugate[ z ].  Conjugate[ z ] like function for complex[ ]" ;
complexI::usage = "complexI.  I like unit imaginary for complex[ ]" ;
fMatrix::usage = "thread a function f over all the elements p in a list." ;
matrixreal::usage =
  "matrixreal.  method to apply real to all elements in matrix.  This is a hack.  \
Can probably set an attribute on the real function to do this." ;
matriximag::usage =
  "matriximag.  method to apply imag to all elements in matrix.  This is a hack.  \
Can probably set an attribute on the imag function to do this." ;
matrixconj::usage =
  "matrixconj.  method to apply conjugate to all elements in matrix.  This is a \
hack.  Can probably set an attribute on the conj function to do this." ;
pauliMatrix::usage =
  "pauliMatrix[ n ], n = 1,2,3.  PauliMatrix[ ] implemented with complex[ ], instead of Complex[ ]." ;
conjugateTranspose::usage =
  "conjugateTranspose[ ].  ConjugateTranspose[ ] like operation for pauliMatrix." ;
grade::usage = "grade.  (internal) An upvalue type that represents a CL(3,0) algebraic element as a pair {grade, v}, where v is a sum of products of Pauli matrices.  These matrices may be scaled by arbitrary numeric or symbolic factors." ;
Scalar::usage = "Scalar[ v ] constructs a scalar grade quantity with value v." ;
Vector::usage = "Vector[ v, n ], where n = {1,2,3} constructs a vector grade quantity with value v in direction n." ;
Bivector::usage = "Bivector[ v, n1, n2 ], where n1,n2 = {1,2,3} constructs a bivector grade quantity with value v in the plane n1,n2." ;
Trivector::usage = "Trivector[ v ] constructs a trivector (pseudoscalar) grade quantity scaled by v." ;
gradeQ::usage = "gradeQ[ m, n ] tests if the multivector m is of grade n.  n = -1 is used internally to represent values of more than one grade." ;
scalarQ::usage = "scalarQ[ m ] tests if the multivector m is of grade 0 (scalar)" ;
vectorQ::usage = "vectorQ[ m ] tests if the multivector m is of grade 1 (vector)" ;
bivectorQ::usage = "bivectorQ[ m ] tests if the multivector m is of grade 2 (bivector)" ;
trivectorQ::usage = "trivectorQ[ m ] tests if the multivector m is of grade 3 (trivector)" ;
bladeQ::usage = "bladeQ[ m ] tests if the multivector is of a single grade." ;
gradeAnyQ::usage = "gradeAnyQ[ ].  predicate pattern match for grade[ _ ]" ;
notGradeQ::usage = "notGradeQ[ ].  predicate pattern match for !grade[ ]" ;
GradeSelection::usage = "GradeSelection[ m, k ] selects the grade k elements from the multivector m.  The selected result is represented internally as a grade[ ] type (so scalar selection is not just a number)." ;
ScalarSelection::usage = "ScalarSelection[ m ] selects the grade 0 (scalar) elements from the multivector m.  The selected result is represented internally as a grade[ ] type (not just a number or an expression)." ;
VectorSelection::usage = "VectorSelection[ m ] selects the grade 1 (vector) elements from the multivector m.  The selected result is represented internally as a grade[ ] type." ;
BivectorSelection::usage = "BivectorSelection[ m ] selects the grade 2 (bivector) elements from the multivector m.  The selected result is represented internally as a grade[ ] type." ;
TrivectorSelection::usage = "TrivectorSelection[ m ] selects the grade 3 (trivector) element from the multivector m if it exists.  The selected result is represented internally as a grade[ ] type (not just an number or expression)." ;
pmagnitude::usage = "pmagnitude[ ].  select the 1,1 element from a pauli matrix assuming it represents \
either a Scalar, or a Trivector (i.e. scaled diagonal matrix)." ;
ScalarValue::usage = "ScalarValue[ m ].  Same as AngleBracket[ m ], aka [ Esc ]<[ Esc ] m1 [ Esc ]>[ Esc ]." ;
ScalarProduct::usage = "ScalarProduct[ ].  Same as AngleBracket[ m1, m2 ], aka [ Esc ]<[ Esc ] m1, m2 [ Esc ]>[ Esc ]." ;
(*Grad::usage = "grad[m,{x,y,z}] computes the vector product of the gradient with multivector m with respect to cartesian coordinates x,y,z..";*)
(*Div::usage = "div[m,{x,y,z}] of a grade k+1 blade m, computes < \[Del] m >_k, where the gradient is evaluated with respect to cartesian coordinates x,y,z." ;*)
(*Curl::usage = "Given a grade (k-1) blade m, curl[ m, {x, y, z} ] = < \[Del] m >_k, where the gradient is evaluated with respect to cartesian coordinates x,y,z." ;*)
Vcurl::usage = "Given a vector m, vcurl[m,{x,y,z}] computes the traditional vector valued curl of that vector with respect to cartesian coordinates x,y,z." ;

(* Begin Private Context *)
Begin["`Private`"]

complex /: complex[ r1_, i1_ ] + complex[ r2_, i2_ ] := complex[ r1 + r2, i1 + i2 ] ;
complex /: r1_ + complex[ r2_, i2_ ] := complex[ r1 + r2, i2 ] ;

complex /: -complex[ re_, im_ ] := complex[ -re, -im ] ;

complex /: complex[ re_ ] := re ;
complex /: complex[ re_, 0 ] := re ;

complex /: complex[ r1_, i1_ ] complex[ r2_, i2_ ] := complex[ r1 r2 - i1 i2, r1 i2 + r2 i1 ] ;

norm[ z_complex ] := ((z // First)^2 + (z // Last)^2) ;

(*special case this one to deal with the sort of products that are \
generated multiplying pauli matrices*)

complex /: Power[ z_complex, 2 ] := complex[ z ] complex[ z ] ;
complex /: Power[ z_complex, n_ ] :=
  Module[ {r = norm[ z ]^(n/2), theta = n ArcTan[ z // First, z // Last ]},
    r complex[ Cos[ theta ], Sin[ theta ] ] ] ;

complexQ[ z_complex ] := True ;
complexQ[ _ ] := False ;
notComplexQ[ v_ ] := Not[ complexQ[ v ] ] ;

complex /: (v_?notComplexQ) complex[ re_, im_ ] := complex[ v re, v im ] ;

real[ z_complex ] := (z // First) ;
imag[ z_complex ] := (z // Last) ;
real[ ex_ ] := ex ;
imag[ ex_ ] := 0 ;
conjugate[ z_complex ] := complex[ z // First, -z // Last ] ;
conjugate[ ex_ ] := ex ;

complexI := complex[ 0, 1 ] ;

fMatrix[ p_, f_ ] := (Function[ a, f@a, Listable ]@p)

matrixreal[ m_ ] := fMatrix[ m, real ] ;

matriximag[ m_ ] := fMatrix[ m, imag ] ;

matrixconj[ m_ ] := fMatrix[ m, conjugate ] ;

pauliMatrix[ 1 ] := PauliMatrix[ 1 ] ;
pauliMatrix[
   2 ] := (PauliMatrix[ 2 ] /. {Complex[ 0, 1 ] -> complexI,
     Complex[ 0, -1 ] -> -complexI}) ;
pauliMatrix[ 3 ] := PauliMatrix[ 3 ] ;
conjugateTranspose[ m_List ] := Transpose[ matrixconj[ m ] ] ;

TraditionalForm[ z_complex ] := (((z // real) + I (z // imag)) // TraditionalForm)
DisplayForm[ z_complex ] := (((z // real) + I (z // imag)) // DisplayForm)
StandardForm[ z_complex ] := (((z // real) + I (z // imag)) // StandardForm)

(* End of complex, and pauliMatrix section.  Define the basic CL(3,0) operations. *)

Scalar[ v_ ] := grade[ 0, v IdentityMatrix[ 2 ] ] ;
Vector[ v_, k_Integer /; k >= 1 && k <= 3 ] :=
  grade[ 1, v pauliMatrix[ k ] ] ;
Bivector[ v_, k_Integer /; k >= 1 && k <= 3, j_Integer /; j >= 1 && j <= 3 ] := grade[ 2, v pauliMatrix[ k ].pauliMatrix[ j ] ] ;
Trivector[ v_ ] := grade[ 3, complexI v IdentityMatrix[ 2 ] ] ;

gradeQ[ m_grade, n_Integer ] := ((m // First) == n)
scalarQ[ m_grade ] := gradeQ[ m, 0 ]
vectorQ[ m_grade ] := gradeQ[ m, 1 ]
bivectorQ[ m_grade ] := gradeQ[ m, 2 ]
trivectorQ[ m_grade ] := gradeQ[ m, 3 ]
bladeQ[ m_grade ] := ((m // First) >= 0)
gradeAnyQ[ m_grade ] := True
gradeAnyQ[ _ ] := False
notGradeQ[ v_ ] := Not[ gradeAnyQ[ v ] ]

directProduct[ t_, v1_, v2_ ] := grade[ t, (v1 // Last).(v2 // Last) ] ;
signedSymmetric[ t_, v1_, v2_, s_ ] :=
  Module[ {a = (v1 // Last), b = (v2 // Last)},
   grade[ t, (a.b + s b.a)/2 ] ] ;
symmetric[ t_, v1_, v2_ ] := signedSymmetric[ t, v1, v2, 1 ] ;
antisymmetric[ t_, v1_, v2_ ] := signedSymmetric[ t, v1, v2, -1 ] ;

(*These operator on just the Pauli matrix portions x of \
pauliGradeSelect[ ,x ]*)
pauliGradeSelect01 := ((# + (# // conjugateTranspose))/2) & ;
pauliGradeSelect23 := ((# - (# // conjugateTranspose))/2) & ;
pauliGradeSelect[ m_, 0 ] := IdentityMatrix[ 2 ] (m/2 // Tr // real // Simplify) ;
pauliGradeSelect[ m_, 1 ] := ((pauliGradeSelect01[ m ] - pauliGradeSelect[ m, 0 ]) // Simplify) ;
pauliGradeSelect[ m_, 2 ] := ((pauliGradeSelect23[ m ] - pauliGradeSelect[ m, 3 ]) // Simplify) ;
pauliGradeSelect[ m_, 3 ] := complexI IdentityMatrix[ 2 ] (m/2 // Tr // imag // Simplify) ;

pauliGradeSelect0 := pauliGradeSelect[ #, 0 ] & ;
pauliGradeSelect1 := pauliGradeSelect[ #, 1 ] & ;
pauliGradeSelect2 := pauliGradeSelect[ #, 2 ] & ;
pauliGradeSelect3 := pauliGradeSelect[ #, 3 ] & ;

GradeSelection[ m_?scalarQ, 0 ] := m ;
GradeSelection[ m_?vectorQ, 1 ] := m ;
GradeSelection[ m_?bivectorQ, 2 ] := m ;
GradeSelection[ m_?trivectorQ, 3 ] := m ;
GradeSelection[ m_, k_Integer /; k >= 0 && k <= 3 ] := grade[ k, pauliGradeSelect[ m // Last, k ] ] ;
ScalarSelection := GradeSelection[ #, 0 ] & ;
VectorSelection := GradeSelection[ #, 1 ] & ;
BivectorSelection := GradeSelection[ #, 2 ] & ;
TrivectorSelection := GradeSelection[ #, 3 ] & ;

binaryOperator[ f_, b_?bladeQ, m_grade ] := Total[ f[ b, # ] & /@ (GradeSelection[ m, # ] & /@ (Range[ 3+1 ] - 1)) ]
binaryOperator[ f_, m_grade, b_?bladeQ ] := Total[ f[ #, b ] & /@ (GradeSelection[ m, # ] & /@ (Range[ 3+1 ] - 1)) ]
binaryOperator[ f_, m1_grade, m2_grade ] := Total[ f[ # // First, # // Last ] & /@ (
    {GradeSelection[ m1, # ] & /@ (Range[ 3+1 ] - 1),
     GradeSelection[ m2, # ] & /@ (Range[ 3+1 ] - 1)} // Transpose) ]

(* Plus *)
grade /: (v1_?notGradeQ) + grade[ k_, v2_ ] := Scalar[ v1 ] + grade[ k, v2 ] ;
grade /: grade[ 0, v1_ ] + grade[ 0, v2_ ] := grade[ 0, v1 + v2 ] ;
grade /: grade[ 1, v1_ ] + grade[ 1, v2_ ] := grade[ 1, v1 + v2 ] ;
grade /: grade[ 2, v1_ ] + grade[ 2, v2_ ] := grade[ 2, v1 + v2 ] ;
grade /: grade[ 3, v1_ ] + grade[ 3, v2_ ] := grade[ 3, v1 + v2 ] ;
grade /: grade[ _, v1_ ] + grade[ _, v2_ ] := grade[ -1, v1 + v2 ] ;

(* Times[ -1, _ ] *)
grade /: -grade[ k_, v_ ] := grade[ k, -v ] ;

(* Times *)
grade /: (v_?notGradeQ) grade[ k_, m_ ] := grade[ k, v m ] ;
grade /: grade[ 0, s_ ] grade[ k_, m_ ] := grade[ k, s.m ] ;
grade /: grade[ 3, p_ ] grade[ 1, m_ ] := grade[ 2, p.m ] ;
grade /: grade[ 3, p_ ] grade[ 2, m_ ] := grade[ 1, p.m ] ;
grade /: grade[ 3, p_ ] grade[ 3, m_ ] := grade[ 0, p.m ] ;
grade /: grade[ 3, p_ ] grade[ _, m_ ] := grade[ -1, p.m ] ;

(* NonCommutativeMultiply *)
grade /: grade[ 0, s_ ] ** grade[ k_, m_ ] := grade[ k, s.m ] ;
grade /: grade[ k_, m_ ] ** grade[ 0, s_ ] := grade[ k, s.m ] ;
grade /: grade[ 3, s_ ] ** grade[ k_, m_ ] := grade[ 3, s ] grade[ k, m ] ;
grade /: grade[ k_, m_ ] ** grade[ 3, s_ ] := grade[ 3, s ] grade[ k, m ] ;
grade /: grade[ _, m1_ ] ** grade[ _, m2_ ] := grade[ -1, m1.m2 ] ;

(* Dot *)
grade /: (s_?notGradeQ).grade[ k_, m_ ] := grade[ k, s m ] ;
grade /: grade[ k_, m_ ].(s_?notGradeQ) := grade[ k, s m ] ;
grade /: grade[ 0, s_ ].grade[ k_, m_ ] := grade[ k, s m ] ;
grade /: grade[ k_, m_ ].grade[ 0, s_ ] := grade[ k, s m ] ;

grade /: (t_?trivectorQ).m_grade := t m ;
grade /: m_grade.(t_?trivectorQ) := t m ;

grade /: (v1_?vectorQ).grade[ 1, v2_ ] := symmetric[ 0, v1, grade[ 1, v2 ] ] ;
grade /: (v_?vectorQ).grade[ 2, b_ ] := antisymmetric[ 1, v, grade[ 2, b ] ] ;
grade /: (b_?bivectorQ).grade[ 1, v_ ] := antisymmetric[ 1, b, grade[ 1, v ] ] ;
grade /: (b1_?bivectorQ).grade[ 2, b2_ ] := symmetric[ 0, b1, grade[ 2, b2 ] ] ;

(* == comparison operator *)
grade /: grade[ _, m1_ ] == grade[ _, m2_ ] := (m1 == m2) ;

(* Dot ; handle dot products where one or more factors is a multivector.  *)
grade /: grade[ g1_, m1_ ] . grade[ g2_, m2_ ]:= binaryOperator[ Dot, grade[ g1, m1 ], grade[ g2, m2 ] ] ;

grade[ _, {{0, 0}, {0, 0}} ] := 0

(*Define a custom wedge operator*)

grade /: grade[ 0, s_ ] \[Wedge] grade[ k_, v_ ] := grade[ k, s.v ] ;
grade /: grade[ k_, v_ ] \[Wedge] grade[ 0_, s_ ] := grade[ k, s.v ] ;
grade /: grade[ 1, v1_ ] \[Wedge] (v2_?vectorQ) := antisymmetric[ 2, grade[ 1, v1 ], v2 ] ;

grade /: grade[ 1, v1_ ] \[Wedge] (v2_?bivectorQ) := symmetric[ 3, grade[ 1, v1 ], v2 ] ;
grade /: grade[ 2, v1_ ] \[Wedge] (v2_?vectorQ) := symmetric[ 3, grade[ 2, v1 ], v2 ] ;
grade /: grade[ 2, _ ] \[Wedge] (v2_?bivectorQ) := 0 ;

(* Only e123 ^ scalar is non zero, and that is handled above *)
grade /: grade[ 3, _ ] \[Wedge] b_?bladeQ := 0 ;
grade /: b_?bladeQ \[Wedge] grade[ 3, _ ] := 0 ;

grade /: grade[ g1_, m1_ ] \[Wedge] grade[ g2_, m2_ ]:= binaryOperator[ Wedge, grade[ g1, m1 ], grade[ g2, m2 ] ] ;

pmagnitude[ m_ ] := m[ [1, 1 ] ] ;

(* AngleBracket,single operand forms, enter with[ Esc ]<[ Esc ] \
v[ Esc ]>[ Esc ] *)
grade /: AngleBracket[ grade[ 0, s_ ] ] := pmagnitude[ s ]
grade /: AngleBracket[ grade[ 1, _ ] ] := 0
grade /: AngleBracket[ grade[ 2, _ ] ] := 0
grade /: AngleBracket[ grade[ 3, _ ] ] := 0
grade /: AngleBracket[ grade[ _, m_ ] ] := ((pauliGradeSelect[ m, 0 ]) // pmagnitude)

ScalarValue[ m_grade ] := AngleBracket[ m ] ;

(* AngleBracket,two operand forms. *)

grade /: AngleBracket[ grade[ 0, s1_ ], grade[ 0, s2_ ] ] := (pmagnitude[ s1 ] pmagnitude[ s2 ]) ;
grade /: AngleBracket[ grade[ 0, s1_ ], grade[ -1, m_ ] ] := (pmagnitude[ s1 ] ((pauliGradeSelect[ m, 0 ]) // pmagnitude)) ;
grade /: AngleBracket[ grade[ 0, s1_ ], grade[ _, _ ] ] := 0 ;
grade /: AngleBracket[ grade[ -1, m_ ], grade[ 0, s1_ ] ] := (pmagnitude[ s1 ] ((pauliGradeSelect[ m, 0 ]) // pmagnitude)) ;
grade /: AngleBracket[ grade[ _, _ ], grade[ 0, s1_ ] ] := 0 ;

grade /: AngleBracket[ grade[ 3, t1_ ], grade[ 3, t2_ ] ] := (pmagnitude[ t1 ] pmagnitude[ t2 ])
grade /: AngleBracket[ grade[ 3, t1_ ], grade[ -1, m_ ] ] := (pmagnitude[ t1 ] ((pauliGradeSelect[ m, 3 ]) // pmagnitude)) ;
grade /: AngleBracket[ grade[ 3, t1_ ], grade[ _, _ ] ] := 0 ;
grade /: AngleBracket[ grade[ -1, m_ ], grade[ 3, t1_ ] ] := (pmagnitude[ t1 ] ((pauliGradeSelect[ m, 3 ]) // pmagnitude)) ;
grade /: AngleBracket[ grade[ _, _ ], grade[ 3, t1_ ] ] := 0 ;

grade /: AngleBracket[ grade[ k1_, m1_ ], grade[ k2_, m2_ ] ] := (pauliGradeSelect[ m1.m2, 0 ] // pmagnitude) ;

ScalarProduct[ m1_grade, m2_grade ] := AngleBracket[ m1, m2 ] ;

bold = Style[ #, Bold ] & ;
esub = Subscript[ bold[ "e" ], # ] & ;
displayMapping = {
   {Scalar[ 1 ], 1, 1, 1},
   {Vector[ 1, 1 ], esub[ 1 ], e[ 1 ], "\\mathbf{e}_1"},
   {Vector[ 1, 2 ], esub[ 2 ], e[ 2 ], "\\mathbf{e}_2"},
   {Vector[ 1, 3 ], esub[ 3 ], e[ 3 ], "\\mathbf{e}_3"},
   {Bivector[ 1, 2, 1 ], esub[ "12" ], e[ 1 ]e[ 2 ], "\\mathbf{e}_{12}"},
   {Bivector[ 1, 3, 2 ], esub[ "23" ], e[ 2 ]e[ 3 ], "\\mathbf{e}_{23}"},
   {Bivector[ 1, 1, 3 ], esub[ "31" ], e[ 3 ]e[ 1 ], "\\mathbf{e}_{31}"},
   {Trivector[ -1 ], esub[ "123" ], e[ 1 ]e[ 2 ]e[ 3 ], "\\mathbf{e}_{123}"}
} ;

GAdisplay[ v_grade, how_ ] :=
  Total[ (Times[ (AngleBracket[ # // First, v ] (*// Simplify*)), #[ [how ] ] ]) & /@
    displayMapping ] ;

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

D[ m_grade, u_ ] := grade[ m // First,
   Map[
     complex[
         D[# // real // Simplify, u],
         D[# // imag // Simplify, u]
     ] &,
     m // Last,
     {2}
   ]
] ;

grade /: Grad[ grade[ k_, m_ ], u_List ] := ( ( Vector[1, #] ** D[ grade[k, m], u[[#]] ] ) & /@ Range[3] ) // Total ;

grade /: Div[ grade[ 1, m_], u_List ] := Grad[ grade[1, m], u ] // ScalarSelection ;
grade /: Div[ grade[ 2, m_], u_List ] := Grad[ grade[2, m], u ] // VectorSelection ;
grade /: Div[ grade[ 3, m_], u_List ] := Grad[ grade[3, m], u ] // BivectorSelection ;

grade /: Curl[ grade[ 1, m_], u_List ] := Grad[ grade[1, m], u ] // BivectorSelection ;
grade /: Curl[ grade[ 2, m_], u_List ] := Grad[ grade[2, m], u ] // TrivectorSelection ;
grade /: Curl[ grade[ 3, m_], u_List ] := 0

Vcurl[ m_?vectorQ, u_List ] := -Trivector[1] Curl[ m, u ] ;

(* End Private Context *)
End[]

Protect[
   Bivector,
   BivectorSelection,
   Curl,
   D,
   DisplayForm,
   Div,
   Grad,
   GradeSelection,
   Scalar,
   ScalarProduct,
   ScalarSelection,
   ScalarValue,
   StandardForm,
   TeXForm,
   TraditionalForm,
   Trivector,
   TrivectorSelection,
   Vcurl,
   Vector,
   VectorSelection,
   complex,
   complexI,
   complexQ,
   conjugate,
   e,
   fMatrix,
   imag,
   matrixconj,
   matriximag,
   matrixreal,
   notComplexQ,
   real
]

EndPackage[ ]
