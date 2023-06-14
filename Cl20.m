(* ::Package:: *)

(* copy this module to a directory in $Path.  Then invoke with <<Cl20` *)
BeginPackage[ "Cl20`" ]

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
   displayMapping,
   esub,
   grade,
   oneTeXForm,
   gradeSelect,
   pmagnitude,
   scalarQ,
   signedSymmetric,
   symmetric,
   vectorQ
]
(*
   gradeSelect0,
   gradeSelect1,
   gradeSelect2,
 *)

Cl20::usage = "Cl20: An implementation of Euclidean (CL(2,0)) Geometric Algebra.

A pair of complex numbers are used to represent the algebraic elements.  This provides an efficient and compact representation
of the entire algebraic space.

Internally, a multivector is represented by a triplet (grade, complex-odd, complex-even).  The grade portion will be
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
complex::usage =
  "complex.  A limited use complex number implementation, independent of the Mathematica Complex" ;
normsq::usage = "normsq[ z ].  A Norm like function for complex[ ]" ;
complexQ::usage = "complexQ[ z ].  predicate pattern match for complex[ ]" ;
notComplexQ::usage = "notComplexQ[ z ].  predicate pattern match for !complex[ ]" ;
real::usage = "real[ z ].  Re[ z ] like function for complex[ ]" ;
imag::usage = "imag[ z ].  Im[ z ] like function for complex[ ]" ;
conjugate::usage = "conjugate[ z ].  Conjugate[ z ] like function for complex[ ]" ;
complexI::usage = "complexI.  I like unit imaginary for complex[ ]" ;
grade::usage = "grade.  (internal) An upvalue type that represents a CL(2,0) algebraic element as a pair {grade, v}, where v is a sum of products of Pauli matrices.  These matrices may be scaled by arbitrary numeric or symbolic factors." ;
Scalar::usage = "Scalar[ v ] constructs a scalar grade quantity with value v." ;
Vector::usage = "Vector[ v, n ], where n = {1,2} constructs a vector grade quantity with value v in direction n." ;
Bivector::usage = "Bivector[ v ], constructs a bivector grade quantity with value v in the plane e1,e2." ;
gradeQ::usage = "gradeQ[ m, n ] tests if the multivector m is of grade n.  n = -1 is used internally to represent values of more than one grade." ;
scalarQ::usage = "scalarQ[ m ] tests if the multivector m is of grade 0 (scalar)" ;
vectorQ::usage = "vectorQ[ m ] tests if the multivector m is of grade 1 (vector)" ;
bivectorQ::usage = "bivectorQ[ m ] tests if the multivector m is of grade 2 (bivector)" ;
bladeQ::usage = "bladeQ[ m ] tests if the multivector is of a single grade." ;
gradeAnyQ::usage = "gradeAnyQ[ ].  predicate pattern match for grade[ _ ]" ;
notGradeQ::usage = "notGradeQ[ ].  predicate pattern match for !grade[ ]" ;
GradeSelection::usage = "GradeSelection[ m, k ] selects the grade k elements from the multivector m.  The selected result is represented internally as a grade[ ] type (so scalar selection is not just a number)." ;
ScalarSelection::usage = "ScalarSelection[ m ] selects the grade 0 (scalar) elements from the multivector m.  The selected result is represented internally as a grade[ ] type (not just a number or an expression)." ;
VectorSelection::usage = "VectorSelection[ m ] selects the grade 1 (vector) elements from the multivector m.  The selected result is represented internally as a grade[ ] type." ;
BivectorSelection::usage = "BivectorSelection[ m ] selects the grade 2 (bivector) elements from the multivector m.  The selected result is represented internally as a grade[ ] type." ;
pmagnitude::usage = "pmagnitude[ ].  select the 1,1 element from a pauli matrix assuming it represents \
a Scalar (i.e. scaled diagonal matrix)." ;
ScalarValue::usage = "ScalarValue[ m ].  Same as AngleBracket[ m ], aka [ Esc ]<[ Esc ] m1 [ Esc ]>[ Esc ]." ;
ScalarProduct::usage = "ScalarProduct[ ].  Same as AngleBracket[ m1, m2 ], aka [ Esc ]<[ Esc ] m1, m2 [ Esc ]>[ Esc ]." ;

(* Begin Private Context *)
Begin["`Private`"]

complex /: complex[ r1_, i1_ ] + complex[ r2_, i2_ ] := complex[ r1 + r2, i1 + i2 ] ;
complex /: r1_ + complex[ r2_, i2_ ] := complex[ r1 + r2, i2 ] ;

complex /: -complex[ re_, im_ ] := complex[ -re, -im ] ;

complex /: complex[ re_ ] := re ;
complex /: complex[ re_, 0 ] := re ;

complex /: complex[ r1_, i1_ ] complex[ r2_, i2_ ] := complex[ r1 r2 - i1 i2, r1 i2 + r2 i1 ] ;

normsq[ z_complex ] := ((z // First)^2 + (z // Last)^2) ;

(*special case this one to deal with the sort of products that are \
generated multiplying pauli matrices*)

complex /: Power[ z_complex, 2 ] := complex[ z ] complex[ z ] ;
complex /: Power[ z_complex, n_ ] :=
  Module[ {r = normsq[ z ]^(n/2), theta = n ArcTan[ z // First, z // Last ]},
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

complex0 := complex[ 0, 0 ] ;
complex1 := complex[ 1, 0 ] ;
complexI := complex[ 0, 1 ] ;

(* grade[g, even-grades, odd-grade] *)
Scalar[ v_ ] := grade[ 0, complex[v, 0], complex0 ] ;
Vector[ v_, k_Integer /; k == 1 ] := grade[ 1, complex0, complex[v, 0] ] ;
Vector[ v_, k_Integer /; k == 2 ] := grade[ 1, complex0, complex[0, v] ] ;
Bivector[ v_ ] := grade[ 2, complex[0, v], complex0 ] ;

gradeQ[ m_grade, n_Integer ] := ((m // First) == n)
scalarQ[ m_grade ] := gradeQ[ m, 0 ]
vectorQ[ m_grade ] := gradeQ[ m, 1 ]
bivectorQ[ m_grade ] := gradeQ[ m, 2 ]
bladeQ[ m_grade ] := ((m // First) >= 0)
gradeAnyQ[ m_grade ] := True
gradeAnyQ[ _ ] := False
notGradeQ[ v_ ] := Not[ gradeAnyQ[ v ] ]

gradeSelect[m_grade, 0] := grade[0, (m[[2]] // real), complex0]
gradeSelect[m_grade, 1] := grade[1, complex0, m[[3]]]
gradeSelect[m_grade, 2] := grade[2, complexI *(m[[2]] // imag), complex0]

GradeSelection[ m_?scalarQ, 0 ] := m ;
GradeSelection[ m_?vectorQ, 1 ] := m ;
GradeSelection[ m_?bivectorQ, 2 ] := m ;
GradeSelection[ m_, k_Integer /; k >= 0 && k <= 2 ] := gradeSelect[ m, k ] ;
ScalarSelection := GradeSelection[ #, 0 ] & ;
VectorSelection := GradeSelection[ #, 1 ] & ;
BivectorSelection := GradeSelection[ #, 2 ] & ;

binaryOperator[ f_, b_?bladeQ, m_grade ] := Total[ f[ b, # ] & /@ (GradeSelection[ m, # ] & /@ (Range[ 2+1 ] - 1)) ]
binaryOperator[ f_, m_grade, b_?bladeQ ] := Total[ f[ #, b ] & /@ (GradeSelection[ m, # ] & /@ (Range[ 2+1 ] - 1)) ]
binaryOperator[ f_, m1_grade, m2_grade ] := Total[ f[ # // First, # // Last ] & /@ (
    {GradeSelection[ m1, # ] & /@ (Range[ 2+1 ] - 1),
     GradeSelection[ m2, # ] & /@ (Range[ 2+1 ] - 1)} // Transpose) ]

(* Plus *)
grade /: (v1_?notGradeQ) + grade[k_, z1_, z2_] :=
  Scalar[v1] + grade[k, z1, z2];
grade /: grade[0, z1_, _] + grade[0, q1_, _] :=
  grade[0, z1 + q1, complex0];
grade /: grade[1, _, z2_] + grade[1, _, q2_] :=
  grade[1, complex0, z2 + q2];
grade /: grade[2, z1_, _] + grade[2, q1_, _] :=
  grade[2, z1 + q1, complex0];
grade /: grade[_, z1_, z2_] + grade[_, q1_, q2_] :=
  grade[-1, z1 + q1, z2 + q2];

(*
(* Times[ Scalar[], k ] *)
grade /: Power[grade[ 0, s_ ], k_] := Scalar[ Power[pmagnitude[ s ], k] ]
(* Vector inversion: Times[ Vector[], -1 ] *)
grade /: Power[grade[ 1, s_ ], -1] := grade[ 1, s ] * Power[ pmagnitude[s.s], -1]
*)

(* Times[ -1, _ ] *)
grade /: -grade[ k_, z1_, z2_ ] := grade[ k, -z1, -z2 ] ;

(* Scalar Times Blade *)
grade /: (v_?notGradeQ) grade[ k_, z1_, z2_ ] := grade[ k, v z1, v z2 ] ;
grade /: grade[ 0, s_, _ ] grade[ k_, q1_, q2_ ] := grade[ k, s q1, s q2 ] ;

(* NonCommutativeMultiply *)

(* vector products: M N \sim \Real(m_2 n_2^\conj) - i \Imag(m_2 n_2^\conj) *)
(* vector squared *)
grade /: grade[ 1, _, s_ ] ** grade[ 1, _, s_ ] := grade[ 0, s conjugate[s], complex0 ] ;
(* vector ** vector *)
grade /: grade[ 1, _, s_ ] ** grade[ 1, _, t_ ] := Module[{p}, p = s conjugate[t]; grade[ -1, complex[ p // real, -p // imag ], complex0 ] ] ;

(* scalar or bivector or 0,2 multivector products: M N \sim m_1 n_1 *)
(* scalar ** scalar *)
grade /: grade[ 0, s_, _ ] ** grade[ 0, t_, _ ] := grade[ 0, s t, complex0 ] ;
(* bivector ** bivector *)
grade /: grade[ 2, s_, _ ] ** grade[ 2, t_, _ ] := grade[ 0, s t, complex0 ] ;

(* scalar times multivector *)
grade /: grade[ 0, s_, _ ] ** grade[ k_, q1_, q2_ ] := grade[ k, s q1, s q2 ] ;
grade /: grade[ k_, z1_, z2_ ] ** grade[ 0, s_, _ ] := grade[ k, s z1, s z2 ] ;

(* special cases for R2 *)
(*  M_1 N_2 \sim \lr{ 0, n_{12} i m_2 }. *)
grade /: grade[ 1, _, m2_ ] ** grade[ 2, n1_, _ ] := grade[ 1, complex0, (n1//imag) complexI m2 ]
(* M_2 N_1 \sim \lr{ 0, - m_{12} i n_2 }. *)
grade /: grade[ 2, m1_, _ ] ** grade[ 1, _, n2_ ] := grade[ 1, complex0, -(m1//imag) complexI n2 ]

(*
   M = (m_1, m_2)
   N = (n_1, n_2)
   M N \sim \lr{ m_1 n_1 + \Real(m_2 n_2^\conj) - i \Imag(m_2 n_2^\conj), n_{11} m_2 + m_{11} n_2 + n_{12} i m_2 - m_{12} i n_2 }.
*)
grade /: grade[ _, m1_, m2_ ] ** grade[ _, n1_, n2_ ] :=
   grade[ -1,
          m1 n1 + real[m2 conjugate[n2]] - complexI imag[ m2 conjugate[n2]],
          real[n1] m2 + real[m1] n2 + imag[n1] complexI m2 - imag[m1] complexI n2 ]

signedSymmetric[v1_, v2_, s_] := (1/2) v1 ** v2 + (1/2) s v2 ** v1
symmetric[v1_, v2_] := signedSymmetric[v1, v2, 1];
antisymmetric[v1_, v2_] := signedSymmetric[v1, v2, -1];

(*Dot*)
grade /: (s_?notGradeQ) . grade[k_, z1_, z2_] :=
grade[k, s z1, s z2];
grade /: grade[k_, z1_, z2_] . (s_?notGradeQ) :=
grade[k, s z2, s z2];
grade /: grade[0, s_, _] . grade[k_, z1_, z2_] := grade[k, s z1, s z2];
grade /: grade[k_, z1_, z2_] . grade[0, s_, _] := grade[k, s z1, s s2];

grade /: (v_?vectorQ) . grade[1, _, z2_] := symmetric[v, grade[1, complex0, z2]];
grade /: (v_?vectorQ) . grade[2, z1_, _] := antisymmetric[v, grade[2, z1, complex0]];
grade /: (b_?bivectorQ) . grade[1, _, z2_] := antisymmetric[b, grade[1, complex0, z2]];
grade /: (b_?bivectorQ) . grade[2, z1_, _] := symmetric[b, grade[2, z1, complex0]];
grade /: grade[ g1_, z1_, z2_ ] . grade[ g2_, q1_, q2_ ]:= binaryOperator[ Dot, grade[ g1, z1, z2 ], grade[ g2, q1, q2 ] ] ;

(* == comparison operator *)
grade /: grade[ _, z1_, z2_ ] == grade[ _, q1_, q2_ ] := ((z1 == q1) && (z2 == q2));

grade[ _, complex0, complex0 ] := 0

(*Define a custom wedge operator*)

grade /: grade[ 0, s_, _ ] \[Wedge] grade[ k_, z1_, z2 ] := grade[ k, s z1, s z2 ] ;
grade /: grade[ k_, z1_, z2_ ] \[Wedge] grade[ 0_, s_ ] := grade[ k, s z1, s z2 ] ;
grade /: grade[ 1, _, v1_ ] \[Wedge] (v2_?vectorQ) := antisymmetric[ grade[ 1, complex0, v1 ], v2 ] ;

(* Only e12 ^ scalar is non zero, and that is handled above *)
grade /: grade[ 2, _, _ ] \[Wedge] b_?bladeQ := 0 ;
grade /: b_?bladeQ \[Wedge] grade[ 2, _, _ ] := 0 ;

grade /: grade[ g1_, z1_, z2_ ] \[Wedge] grade[ g2_, q1_, q2_ ]:= binaryOperator[ Wedge, grade[ g1, z1, z2 ], grade[ g2, q1, q2 ] ] ;

(*AngleBracket,single operand forms,enter with[Esc]<[Esc] \
v[Esc]>[Esc]*)
grade /: AngleBracket[grade[1, _, _]] := 0
grade /: AngleBracket[grade[2, _, _]] := 0
grade /: AngleBracket[grade[_, z1_, _]] := (z1 // real)

ScalarValue[m_grade] := AngleBracket[m];

(* AngleBracket,two operand forms. *)

(* \gpgrade{M N}{0} = \Real\lr{m_1 n_1 + m_2 n_2^\conj} *)
grade /: AngleBracket[ grade[  0, z1_, _ ], grade[  0, q1_, _ ] ] := real[ z1 q1 ]
grade /: AngleBracket[ grade[  0, z1_, _ ], grade[ -1, q1_, q2_ ] ] := real[ z1 q1 ]
grade /: AngleBracket[ grade[ -1, z1_, _ ], grade[  0, q1_, _ ] ] := real[ z1 q1 ]

grade /: AngleBracket[ grade[ 0, z1_, _ ], grade[ _, _, _ ] ] := 0 ;
grade /: AngleBracket[ grade[ _, _, _ ], grade[ 0, q1_, _ ] ] := 0 ;

grade /: AngleBracket[ grade[ k1_, z1_, z2_ ], grade[ k2_, q1_, q2_ ] ] := real[ z1 q1 + z2 conjugate[q2] ]

ScalarProduct[ m1_grade, m2_grade ] := AngleBracket[ m1, m2 ] ;

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

TraditionalForm[ m_grade ] := ((GAdisplay[ m, 2 ]) // TraditionalForm) ;
DisplayForm[ m_grade ] := GAdisplay[ m, 2 ] ;
Format[ m_grade ] := GAdisplay[ m, 2 ] ;
StandardForm[ m_grade ] := GAdisplay[ m, 3 ] ;

oneTeXForm[m_, blade_, disp_] := Module[{p},
  p = AngleBracket[blade, m];
  If[ p // PossibleZeroQ, 0,
   If[ p === 1, disp,
    Row[ {"(", p // TeXForm, ")(", disp, ")"} ]
    (* Times[p // TeXForm, disp] *)
     ]
   ]
  ]

TeXForm[m_grade] := Total[ oneTeXForm[m, # // First, #[[4]]] & /@ displayMapping ];

(* End Private Context *)
End[]

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
