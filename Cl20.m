(* ::Package:: *)

(* copy this module to a directory in $Path.  Then invoke with <<Cl20` *)
BeginPackage[ "Cl20`" ]

(* Must reference any global symbol (or some of them) before Unprotecting it, since it may not have
   been loaded:

   http://mathematica.stackexchange.com/a/137007/10
 *)
{TraditionalForm, DisplayForm, StandardForm, Format};

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
   multivector,
   oneTeXForm,
   gradeSelect,
   pmagnitude,
   scalarQ,
   signedSymmetric,
   symmetric,
   vectorQ
]

Cl20::usage = "Cl20: An implementation of Euclidean (CL(2,0)) Geometric Algebra.

A pair of Complex numbers are used to represent the algebraic elements.  This provides an efficient and compact representation
of the entire algebraic space.  For details see:

https://peeterjoot.com/archives/math2023/bicomplexCl20.pdf

with the motivational exploratory work leading to that here:

https://peeterjoot.com/archives/math2023/bicomplexGA20.v3.pdf

Internally, a multivector is represented by a triplet (multivector, Complex, Complex), where
 the first complex value encodes the even grades, and the second complex value encodes the odd (vector) multivector.
The numerical multivector portion will be
obliterated when adding objects that have different multivector, or multiplying vectors or bivectors.  When
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
   vb1 ** vb1 Vectors and bivectors when multiplied have to use the NonCommutativeMultiply operator, but any multivector object may also.
   m1 . m2    Dot product.  The functional form Dot[ m1, m2 ] may also be used.
   m1 ^ m2   Wedgeproduct.  Enter with m1 [ Esc ]^[ Esc ] m2.  The functional form Wedge[ m1, m2 ]
   <m>        Scalar selection.  Enter with [ Esc ]<[ Esc ] m [ Esc ]>[ Esc ].  The functional form ScalarValue[ m ] may also be used.  This returns the numeric (or expression) value of the scalar multivector of the multivector, and not a multivector[ ] object.
   <m1,m2>    Scalar product.  Enter with [ Esc ]<[ Esc ] m1,m2 [ Esc ]>[ Esc ].  The functional form ScalarProduct[ m1, m2 ] may also be used.  This returns the numeric (or expression) value of the scalar product of the multivectors, and not a multivector[ ] object.

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
   - multivectorQ
   - notMultivectorQ

TODO:

1) How to get better formatted output by default without using one of TraditionalForm, DisplayForm, StandardForm ?

2) Can a package have options (i.e. to define the name of the e[ ] operator used in StandardForm that represents a basis vector).

3) proper packaging stuff:  private for internals.
" ;
multivector::usage = "multivector.  (internal) An upvalue type that represents a CL(2,0) algebraic element as a triplet {grade, v, w}, where v encodes the even grades (using a Complex value), and w encodes the vector (odd) grade, also using a Complex value.  These Complex values may be scaled by arbitrary real numeric or symbolic factors." ;
Scalar::usage = "Scalar[ v ] constructs a scalar multivector quantity with value v." ;
Vector::usage = "Vector[ v, n ], where n = {1,2} constructs a vector multivector quantity with value v in direction n." ;
Bivector::usage = "Bivector[ v ], constructs a bivector multivector quantity with value v in the plane e1,e2." ;
gradeQ::usage = "gradeQ[ m, n ] tests if the multivector m is of grade n.  n = -1 is used internally to represent values of more than one grade." ;
scalarQ::usage = "scalarQ[ m ] tests if the multivector m is of grade 0 (scalar)" ;
vectorQ::usage = "vectorQ[ m ] tests if the multivector m is of grade 1 (vector)" ;
bivectorQ::usage = "bivectorQ[ m ] tests if the multivector m is of grade 2 (bivector)" ;
bladeQ::usage = "bladeQ[ m ] tests if the multivector is of a single grade." ;
multivectorQ::usage = "multivectorQ[ ].  predicate pattern match for multivector[ _ ]" ;
notMultivectorQ::usage = "notMultivectorQ[ ].  predicate pattern match for !multivector[ ]" ;
GradeSelection::usage = "GradeSelection[ m, k ] selects the grade k elements from the multivector m.  The selected result is represented internally as a multivector[ ] type (so scalar selection is not just a number)." ;
ScalarSelection::usage = "ScalarSelection[ m, asMv: True ] selects the multivector 0 (scalar) elements from the multivector m.  The selected result is represented internally as a multivector[ ] type (not just a number or an expression).  If asMv is true, the result will be returned as a multivector (multivector) object, not scalar." ;
VectorSelection::usage = "VectorSelection[ m, asMv_Boolean : True ] selects the multivector 1 (vector) elements from the multivector m.  The selected result is represented internally as a multivector[ ] type.  If asMv is False then the result will be converted to a List of coordinates." ;
BivectorSelection::usage = "BivectorSelection[ m, asMv ] selects the multivector 2 (bivector) elements from the multivector m.  If asMv is True, the selected result is represented internally as a multivector[ ] type (multivector), otherwise as a scalar.." ;
pmagnitude::usage = "pmagnitude[ ].  select the 1,1 element from a pauli matrix assuming it represents \
a Scalar (i.e. scaled diagonal matrix)." ;
ScalarValue::usage = "ScalarValue[ m ].  Same as AngleBracket[ m ], aka [ Esc ]<[ Esc ] m1 [ Esc ]>[ Esc ]." ;
ScalarProduct::usage = "ScalarProduct[ ].  Same as AngleBracket[ m1, m2 ], aka [ Esc ]<[ Esc ] m1, m2 [ Esc ]>[ Esc ]." ;

(* Begin Private Context *)
Begin["`Private`"]

(* multivector[g, even-grades, odd-multivector] *)
Scalar[ v_ ] := multivector[ 0, v, 0 ] ;
Vector[ v_, k_Integer /; k == 1 ] := multivector[ 1, 0, v ] ;
Vector[ v_, k_Integer /; k == 2 ] := multivector[ 1, 0, I v ] ;
Bivector[ v_ ] := multivector[ 2, I v, 0 ] ;

gradeQ[ m_multivector, n_Integer ] := ((m // First) == n)
scalarQ[ m_multivector ] := gradeQ[ m, 0 ]
vectorQ[ m_multivector ] := gradeQ[ m, 1 ]
bivectorQ[ m_multivector ] := gradeQ[ m, 2 ]
bladeQ[ m_multivector ] := ((m // First) >= 0)
multivectorQ[ m_multivector ] := True
multivectorQ[ _ ] := False
notMultivectorQ[ v_ ] := Not[ multivectorQ[ v ] ]

gradeSelect[m_multivector, 0] := multivector[0, (m[[2]] // Re), 0]
gradeSelect[m_multivector, 1] := multivector[1, 0, m[[3]]]
gradeSelect[m_multivector, 2] := multivector[2, I *(m[[2]] // Im), 0]

GradeSelection[ m_?scalarQ, 0 ] := m ;
GradeSelection[ m_?vectorQ, 1 ] := m ;
GradeSelection[ m_?bivectorQ, 2 ] := m ;
GradeSelection[ m_, k_Integer /; k >= 0 && k <= 2 ] := gradeSelect[ m, k ] ;
ScalarSelection[ v_multivector ] := GradeSelection[ v, 0 ] ;
ScalarSelection[ v_multivector, True ] := GradeSelection[ v, 0 ] ;
ScalarSelection[ v_multivector, False ] := (GradeSelection[ v, 0 ][[2]]) // Real ;
VectorSelection[ v_multivector ] := GradeSelection[ v, 1 ] ;
VectorSelection[ v_multivector, True ] := GradeSelection[ v, 1 ] ;
VectorSelection[ v_multivector, False ] := (GradeSelection[ v, 1 ] // Last) // ReIm;
BivectorSelection[ v_multivector ] := GradeSelection[ v, 2 ];
BivectorSelection[ v_multivector, True ] := GradeSelection[ v, 2 ];
BivectorSelection[ v_multivector, False ] := (GradeSelection[ v, 2 ][[2]]) // Imag;

binaryOperator[ f_, b_?bladeQ, m_multivector ] := Total[ f[ b, # ] & /@ (GradeSelection[ m, # ] & /@ (Range[ 2+1 ] - 1)) ]
binaryOperator[ f_, m_multivector, b_?bladeQ ] := Total[ f[ #, b ] & /@ (GradeSelection[ m, # ] & /@ (Range[ 2+1 ] - 1)) ]
binaryOperator[ f_, m1_multivector, m2_multivector ] := Total[ f[ # // First, # // Last ] & /@ (
    {GradeSelection[ m1, # ] & /@ (Range[ 2+1 ] - 1),
     GradeSelection[ m2, # ] & /@ (Range[ 2+1 ] - 1)} // Transpose) ]

(* Plus *)
multivector /: (v1_?notMultivectorQ) + multivector[k_, z1_, z2_] :=
  Scalar[v1] + multivector[k, z1, z2];
multivector /: multivector[0, z1_, _] + multivector[0, q1_, _] :=
  multivector[0, z1 + q1, 0];
multivector /: multivector[1, _, z2_] + multivector[1, _, q2_] :=
  multivector[1, 0, z2 + q2];
multivector /: multivector[2, z1_, _] + multivector[2, q1_, _] :=
  multivector[2, z1 + q1, 0];
multivector /: multivector[_, z1_, z2_] + multivector[_, q1_, q2_] :=
  multivector[-1, z1 + q1, z2 + q2];

(*
(* Times[ Scalar[], k ] *)
multivector /: Power[multivector[ 0, s_ ], k_] := Scalar[ Power[pmagnitude[ s ], k] ]
(* Vector inversion: Times[ Vector[], -1 ] *)
multivector /: Power[multivector[ 1, s_ ], -1] := multivector[ 1, s ] * Power[ pmagnitude[s.s], -1]
*)

(* Times[ -1, _ ] *)
multivector /: -multivector[ k_, z1_, z2_ ] := multivector[ k, -z1, -z2 ] ;

(* Scalar Times Blade *)
multivector /: (v_?notMultivectorQ) multivector[ k_, z1_, z2_ ] := multivector[ k, v z1, v z2 ] ;
multivector /: multivector[ 0, s_, _ ] multivector[ k_, q1_, q2_ ] := multivector[ k, s q1, s q2 ] ;

(* NonCommutativeMultiply *)

(* vector products: M N \sim \Real(m_2 n_2^\conj) - i \Imag(m_2 n_2^\conj) *)
(* vector squared *)
multivector /: multivector[ 1, _, s_ ] ** multivector[ 1, _, s_ ] := multivector[ 0, s Conjugate[s], 0 ] ;
(* vector ** vector *)
multivector /: multivector[ 1, _, s_ ] ** multivector[ 1, _, t_ ] := multivector[ -1, t Conjugate[s], 0 ] ;

(* scalar or bivector or 0,2 multivector products: M N \sim m_1 n_1 *)
(* scalar ** scalar *)
multivector /: multivector[ 0, s_, _ ] ** multivector[ 0, t_, _ ] := multivector[ 0, s t, 0 ] ;
(* bivector ** bivector *)
multivector /: multivector[ 2, s_, _ ] ** multivector[ 2, t_, _ ] := multivector[ 0, s t, 0 ] ;

(* scalar times multivector *)
multivector /: multivector[ 0, s_, _ ] ** multivector[ k_, q1_, q2_ ] := multivector[ k, s q1, s q2 ] ;
multivector /: multivector[ k_, z1_, z2_ ] ** multivector[ 0, s_, _ ] := multivector[ k, s z1, s z2 ] ;

(* special cases for R2 *)
(*  M_1 N_2 \sim \lr{ 0, n_{12} i m_2 }. *)
multivector /: multivector[ 1, _, m2_ ] ** multivector[ 2, n1_, _ ] := multivector[ 1, 0, (n1//Im) I m2 ]
(* M_2 N_1 \sim \lr{ 0, - m_{12} i n_2 }. *)
multivector /: multivector[ 2, m1_, _ ] ** multivector[ 1, _, n2_ ] := multivector[ 1, 0, -(m1//Im) I n2 ]

(*
   M = (m_1, m_2)
   N = (n_1, n_2)

   M N \sim \lr{ m_1 n_1 + \Real(m_2 n_2^\conj) - i \Imag(m_2 n_2^\conj), n_{11} m_2 + m_{11} n_2 + n_{12} i m_2 - m_{12} i n_2 }.
       =
            \lr{ m_1 n_1 + m_2^\conj n_2, n_{11} m_2 + m_{11} n_2 + n_{12} i m_2 - m_{12} i n_2 }.
       =
            \lr{ m_1 n_1 + m_2^\conj n_2, n_1 m_2 + m_1^\conj n_2 }.
*)
multivector /: multivector[ _, m1_, m2_ ] ** multivector[ _, n1_, n2_ ] :=
   multivector[ -1, m1 n1 + Conjugate[m2] n2, n1 m2 + Conjugate[m1] n2 ]

signedSymmetric[v1_, v2_, s_] := (1/2) v1 ** v2 + (1/2) s v2 ** v1
symmetric[v1_, v2_] := signedSymmetric[v1, v2, 1];
antisymmetric[v1_, v2_] := signedSymmetric[v1, v2, -1];

(*Dot*)
multivector /: (s_?notMultivectorQ) . multivector[k_, z1_, z2_] :=
multivector[k, s z1, s z2];
multivector /: multivector[k_, z1_, z2_] . (s_?notMultivectorQ) :=
multivector[k, s z2, s z2];
multivector /: multivector[0, s_, _] . multivector[k_, z1_, z2_] := multivector[k, s z1, s z2];
multivector /: multivector[k_, z1_, z2_] . multivector[0, s_, _] := multivector[k, s z1, s s2];

multivector /: (v_?vectorQ) . multivector[1, _, z2_] := symmetric[v, multivector[1, 0, z2]];
multivector /: (v_?vectorQ) . multivector[2, z1_, _] := antisymmetric[v, multivector[2, z1, 0]];
multivector /: (b_?bivectorQ) . multivector[1, _, z2_] := antisymmetric[b, multivector[1, 0, z2]];
multivector /: (b_?bivectorQ) . multivector[2, z1_, _] := symmetric[b, multivector[2, z1, 0]];
multivector /: multivector[ g1_, z1_, z2_ ] . multivector[ g2_, q1_, q2_ ]:= binaryOperator[ Dot, multivector[ g1, z1, z2 ], multivector[ g2, q1, q2 ] ] ;

(* == comparison operator *)
multivector /: multivector[ _, z1_, z2_ ] == multivector[ _, q1_, q2_ ] := ((z1 == q1) && (z2 == q2));

multivector[ _, 0, 0 ] := 0

(*Define a custom wedge operator*)

multivector /: multivector[ 0, s_, _ ] \[Wedge] multivector[ k_, z1_, z2 ] := multivector[ k, s z1, s z2 ] ;
multivector /: multivector[ k_, z1_, z2_ ] \[Wedge] multivector[ 0_, s_ ] := multivector[ k, s z1, s z2 ] ;
multivector /: multivector[ 1, _, v1_ ] \[Wedge] (v2_?vectorQ) := antisymmetric[ multivector[ 1, 0, v1 ], v2 ] ;

(* Only e12 ^ scalar is non zero, and that is handled above *)
multivector /: multivector[ 2, _, _ ] \[Wedge] b_?bladeQ := 0 ;
multivector /: b_?bladeQ \[Wedge] multivector[ 2, _, _ ] := 0 ;

multivector /: multivector[ g1_, z1_, z2_ ] \[Wedge] multivector[ g2_, q1_, q2_ ]:= binaryOperator[ Wedge, multivector[ g1, z1, z2 ], multivector[ g2, q1, q2 ] ] ;

(*AngleBracket,single operand forms,enter with[Esc]<[Esc] \
v[Esc]>[Esc]*)
multivector /: AngleBracket[multivector[1, _, _]] := 0
multivector /: AngleBracket[multivector[2, _, _]] := 0
multivector /: AngleBracket[multivector[_, z1_, _]] := (z1 // Re)

ScalarValue[m_multivector] := AngleBracket[m];

(* AngleBracket,two operand forms. *)

(* \gpgrade{M N}{0} = \Real\lr{m_1 n_1 + m_2 n_2^\conj} *)
multivector /: AngleBracket[ multivector[  0, z1_, _ ], multivector[  0, q1_, _ ] ] := Re[ z1 q1 ]
multivector /: AngleBracket[ multivector[  0, z1_, _ ], multivector[ -1, q1_, q2_ ] ] := Re[ z1 q1 ]
multivector /: AngleBracket[ multivector[ -1, z1_, _ ], multivector[  0, q1_, _ ] ] := Re[ z1 q1 ]

multivector /: AngleBracket[ multivector[ 0, z1_, _ ], multivector[ _, _, _ ] ] := 0 ;
multivector /: AngleBracket[ multivector[ _, _, _ ], multivector[ 0, q1_, _ ] ] := 0 ;

multivector /: AngleBracket[ multivector[ k1_, z1_, z2_ ], multivector[ k2_, q1_, q2_ ] ] := Re[ z1 q1 + z2 Conjugate[q2] ]

ScalarProduct[ m1_multivector, m2_multivector ] := AngleBracket[ m1, m2 ] ;

bold = Style[ #, Bold ] & ;
esub = Subscript[ bold[ "e" ], # ] & ;
displayMapping = {
   {Scalar[ 1 ], 1, 1, 1},
   {Vector[ 1, 1 ], esub[ 1 ], e[ 1 ], "\\mathbf{e}_1"},
   {Vector[ 1, 2 ], esub[ 2 ], e[ 2 ], "\\mathbf{e}_2"},
   {Bivector[ -1 ], esub[ "12" ], e[ 1 ]e[ 2 ], "\\mathbf{e}_{12}"}
} ;

GAdisplay[ v_multivector, how_ ] :=
  Total[ (Times[ AngleBracket[ # // First, v ], #[ [how ] ] ]) & /@
    displayMapping ] ;

TraditionalForm[ m_multivector ] := ((GAdisplay[ m, 2 ]) // TraditionalForm) ;
DisplayForm[ m_multivector ] := GAdisplay[ m, 2 ] ;
Format[ m_multivector ] := GAdisplay[ m, 2 ] ;
StandardForm[ m_multivector ] := GAdisplay[ m, 3 ] ;

oneTeXForm[m_, blade_, disp_] := Module[{p},
  p = AngleBracket[blade, m];
  If[ p // PossibleZeroQ, 0,
   If[ p === 1, disp,
    Row[ {"(", p // TeXForm, ")(", disp, ")"} ]
    (* Times[p // TeXForm, disp] *)
     ]
   ]
  ]

TeXForm[m_multivector] := Total[ oneTeXForm[m, # // First, #[[4]]] & /@ displayMapping ];

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
