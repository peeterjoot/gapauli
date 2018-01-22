(* ::Package:: *)

(* copy this module to a directory in $Path.  Then invoke with <<GA13` *)
BeginPackage[ "GA13`" ]

(* Must reference any global symbol (or some of them) before Unprotecting it,since it may not have been loaded:
   http://mathematica.stackexchange.com/a/137007/10
 *)
{D, TraditionalForm, DisplayForm, StandardForm, Format, Grad, Div, Curl};

Unprotect[
   Bivector,
   BivectorSelection,
   Curl,
   D,
   DisplayForm,
   Div,
   Format,
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
   bs,
   complex,
   complexI,
   complexQ,
   conjugate,
   conjugateTranspose,
   diracGammaMatrix,
   diracGradeSelect,
   diracGradeSelect0,
   diracGradeSelect1,
   diracGradeSelect2,
   diracGradeSelect3,
   diracGradeSelect4,
   diracPseudoscalar,
   directProduct,
   displayMapping,
   fMatrix,
   grade,
   gsub,
   imag,
   matrixconj,
   matriximag,
   matrixreal,
   notComplexQ,
   oneTeXForm,
   pairs,
   pauliMatrix,
   pmagnitude,
   real,
   scalarQ,
   signedSymmetric,
   smagnitude,
   symmetric,
   triplets,
   trivectorQ,
   trivectors,
   ts,
   vectorQ,
   vs
]

GA13::usage = "GA13: An implementation of a Minkowski (CL(1,3)) Geometric Algebra.

Dirac matrices are used to represent the algebraic elements.  This provides an fairly efficient and compact representation
of the entire algebraic space.  This representation unfortunately has a built in redundancy, since the complex
4x4 matrix has 32 degrees of freedom, while there are only 16 elements in the algebraic space.

Internally, a multivector is represented by a pair (grade, dirac-representation).  The grade portion will be
obliterated when adding objects that have different grade, or multiplying vectors or bivectors.  When
it is available, certain operations can be optimized.  Comparison ignores the cached grade if it exists.

Elements of the algebra can be constructed with one of

   Scalar[ v ]
   Vector[ v, n ]
   Bivector[ v, n, m ]
   Trivector[ v, n, m, o ]
   Quadvector[ v ]

Example:

   m = Scalar[ Sin[ x ] ] + Vector[ Log[ z ], 3 ] + Trivector[ 7, 0, 1, 3 ] ;
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
   - QuadvectorSelection
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
   - quadvectorQ
   - bladeQ
   - gradeAnyQ
   - notGradeQ

TODO:

1) How to get better formatted output by default without using one of TraditionalForm, DisplayForm, StandardForm ?

2) Can a package have options (i.e. to define the name of the e[ ] operator used in StandardForm that represents a basis vector).

3) proper packaging stuff:  private for internals.
" ;

complex::usage =
  "complex.  A limited use complex number implementation to use internally in a Dirac or Dirac matrix basis representation, independent of any Complex";
norm::usage = "norm[ z ].  A Norm like function for complex[ ]";
complexQ::usage =
  "complexQ[ z ].  predicate pattern match for complex[ ]";
notComplexQ::usage =
  "notComplexQ[ z ].  predicate pattern match for !complex[ ]";
real::usage = "real[ z ].  Re[ z ] like function for complex[ ]";
imag::usage = "imag[ z ].  Im[ z ] like function for complex[ ]";
conjugate::usage =
  "conjugate[ z ].  Conjugate[ z ] like function for complex[ ]";
complexI::usage = "complexI.  I like unit imaginary for complex[ ]";
fMatrix::usage =
  "thread a function f over all the elements p in a list.";
matrixreal::usage =
  "matrixreal.  method to apply real to all elements in matrix.  This is a hack.  Can probably set an attribute on the real function to do this.";
matriximag::usage =
  "matriximag.  method to apply imag to all elements in matrix.  This is a hack.  Can probably set an attribute on the imag function to do this.";
matrixconj::usage =
  "matrixconj.  method to apply conjugate to all elements in matrix.  This is a hack.  Can probably set an attribute on the conj function to do this.";
pauliMatrix::usage =
  "pauliMatrix[ n ], n = 1,2,3.  PauliMatrix[ ] implemented with complex[ ], instead of Complex[ ].";
diracGammaMatrix::usage =
  "diracGammaMatrix[ n ], n = 0,1,2,3.  This is like the DiracGammaMatrix[ ] mentioned in mathpages implemented with complex[ ], instead of Complex[ ].";
conjugateTranspose::usage =
  "conjugateTranspose[ ].  ConjugateTranspose[ ] like operation for diracGammaMatrix.";
grade::usage =
  "grade.  (internal) An upvalue type that represents a CL(1,3) algebraic element as a pair {grade, v}, where v is a sum of products of Dirac matrices.  These matrices may be scaled by arbitrary numeric or symbolic factors.";
Scalar::usage =
  "Scalar[ v ] constructs a scalar grade quantity with value v.";
Vector::usage =
  "Vector[ v, n ], where n = {0,1,2,3} constructs a vector grade quantity with value v in direction n.";
Bivector::usage =
  "Bivector[ v, n1, n2 ], where n1,n2 = {1,2,3} constructs a bivector grade quantity with value v in the plane n1,n2.";
Trivector::usage =
  "Trivector[ v, k, l, m ] constructs a trivector (pseudoscalar) grade quantity scaled by v.";
Quadvector::usage =
  "Quadvector[ v ] constructs a quadvector (pseudoscalar) grade quantity scaled by v.";
gradeQ::usage =
  "gradeQ[ m, n ] tests if the multivector m is of grade n.  n = -1 is used internally to represent values of more than one grade.";
scalarQ::usage =
  "scalarQ[ m ] tests if the multivector m is of grade 0 (scalar)";
vectorQ::usage =
  "vectorQ[ m ] tests if the multivector m is of grade 1 (vector)";
bivectorQ::usage =
  "bivectorQ[ m ] tests if the multivector m is of grade 2 (bivector)";
trivectorQ::usage =
  "trivectorQ[ m ] tests if the multivector m is of grade 3 (trivector)";
quadvectorQ::usage =
  "quadvectorQ[ m ] tests if the multivector m is of grade 4 (quadvector)";
bladeQ::usage =
  "bladeQ[ m ] tests if the multivector is of a single grade.";
gradeAnyQ::usage =
  "gradeAnyQ[ ].  predicate pattern match for grade[ _ ]";
notGradeQ::usage =
  "notGradeQ[ ].  predicate pattern match for !grade[ ]";
diracGradeSelect::usage = "diracGradeSelect[m, g] selects the grade g components of the multivector m, as represented by a Dirac matrix (internal).";
vs::usage = "vs[m, n], n in [0,3].  (internal) Select the n'th component of the vector as part of grade 1 selection.";
bs::usage = "bs[m, n, o], n,o in [0,3].  (internal) Select bivector components as part of grade 2 selection.";
ts::usage = "ts[m, n, o, p], n,o,p in [0,3].  (internal) Select trivector components as part of grade 3 selection.";
triplets::usage = "Unique indexes for trivector components (internal)";
trivectors::usage = "Set of all the trivector basis elements for a grade 3 blade, corresponding to indexes of the triplets variable. (internal)";
GradeSelection::usage =
  "GradeSelection[ m, k ] selects the grade k elements from the multivector m.  The selected result is represented internally as a grade[ ] type (so scalar selection is not just a number).";
ScalarSelection::usage =
  "ScalarSelection[ m ] selects the grade 0 (scalar) elements from the multivector m.  The selected result is represented internally as a grade[ ] type (not just a number or an expression).";
VectorSelection::usage =
  "VectorSelection[ m ] selects the grade 1 (vector) elements from the multivector m.  The selected result is represented internally as a grade[ ] type.";
BivectorSelection::usage =
  "BivectorSelection[ m ] selects the grade 2 (bivector) elements from the multivector m.  The selected result is represented internally as a grade[ ] type.";
TrivectorSelection::usage =
  "TrivectorSelection[ m ] selects the grade 3 (trivector) element from the multivector m if it exists.  The selected result is represented internally as a grade[ ] type (not just an number or expression).";
QuadvectorSelection::usage =
  "QuadvectorSelection[ m ] selects the grade 4 (trivector) element from the multivector m if it exists.  The selected result is represented internally as a grade[ ] type (not just an number or expression).";
smagnitude::usage =
  "smagnitude[ ].  select the 1,1 element from a dirac matrix assuming it represents a Scalar";
pmagnitude::usage =
  "pmagnitude[ ].  select the 1,3 element from a dirac matrix assuming it represents a scaled version of g_0123";
ScalarValue::usage =
  "ScalarValue[ m ].  Same as AngleBracket[ m ], aka [ Esc ]<[ Esc ] m1 [ Esc ]>[ Esc ].";
ScalarProduct::usage =
  "ScalarProduct[ ].  Same as AngleBracket[ m1, m2 ], aka [ Esc ]<[ Esc ] m1, m2 [ Esc ]>[ Esc ].";

(* Begin Private Context *) 
Begin["`Private`"] 

complex /: complex[r1_, i1_] + complex[r2_, i2_] := complex[r1 + r2, i1 + i2];
complex /: r1_ + complex[r2_, i2_] := complex[r1 + r2, i2];

complex /: -complex[re_, im_] := complex[-re, -im];

complex /: complex[re_] := re;
complex /: complex[re_, 0] := re;

complex /: complex[r1_, i1_] complex[r2_, i2_] := complex[r1 r2 - i1 i2, r1 i2 + r2 i1];

norm[z_complex] := ((z // First)^2 + (z // Last)^2);

(*special case this one to deal with the sort of products that are generated multiplying dirac matrices*)

complex /: Power[z_complex, 2] := complex[z] complex[z];
complex /: Power[z_complex, n_] :=
  Module[{r = norm[z]^(n/2), theta = n ArcTan[z // First, z // Last]},
    r complex[Cos[theta], Sin[theta]]];

complexQ[z_complex] := True;
complexQ[_] := False;
notComplexQ[v_] := Not[complexQ[v]];

complex /: (v_?notComplexQ) complex[re_, im_] := complex[v re, v im];

real[z_complex] := (z // First);
imag[z_complex] := (z // Last);
real[ex_] := ex;
imag[ex_] := 0;
conjugate[z_complex] := complex[z // First, -z // Last];
conjugate[ex_] := ex;

complexI := complex[0, 1];

fMatrix[p_, f_] := (Function[a, f@a, Listable]@p)

matrixreal[m_] := fMatrix[m, real];

matriximag[m_] := fMatrix[m, imag];

matrixconj[m_] := fMatrix[m, conjugate];


pauliMatrix[1] := PauliMatrix[1];
pauliMatrix[2] := (PauliMatrix[2] /. {Complex[0, 1] -> complexI, Complex[0, -1] -> -complexI});
pauliMatrix[3] := PauliMatrix[3];


diracGammaMatrix[0] = DiagonalMatrix[{1, 1, -1, -1}];
diracGammaMatrix[1] = ArrayFlatten[{{DiagonalMatrix[{0, 0}], pauliMatrix[1]}, {-pauliMatrix[1], DiagonalMatrix[{0, 0}]}}];
diracGammaMatrix[2] = ArrayFlatten[{{DiagonalMatrix[{0, 0}], pauliMatrix[2]}, {-pauliMatrix[2], DiagonalMatrix[{0, 0}]}}];
diracGammaMatrix[3] = ArrayFlatten[{{DiagonalMatrix[{0, 0}], pauliMatrix[3]}, {-pauliMatrix[3], DiagonalMatrix[{0, 0}]}}];

conjugateTranspose[m_List] := Transpose[matrixconj[m]];

TraditionalForm[ z_complex] := (((z // real) + I (z // imag)) // TraditionalForm)
DisplayForm[ z_complex] := (((z // real) + I (z // imag)) // DisplayForm)
StandardForm[ z_complex] := (((z // real) + I (z // imag)) // StandardForm)

(* End of complex,and diracGammaMatrix section.Define the basic CL(1,3) operations. *)



Scalar[v_] := grade[0, v IdentityMatrix[4]];
Vector[v_, k_Integer /; k >= 0 && k <= 3] := grade[1, v diracGammaMatrix[k]];
Bivector[v_, k_Integer /; k >= 0 && k <= 3, j_Integer /; j >= 0 && j <= 3] := grade[2, v diracGammaMatrix[k].diracGammaMatrix[j]];
Trivector[v_, k_, l_, m_] := grade[3, v diracGammaMatrix[k].diracGammaMatrix[l].diracGammaMatrix[ m]];
Quadvector[v_] := grade[4, v diracPseudoscalar];

gradeQ[m_grade, n_Integer] := ((m // First) == n)
scalarQ[m_grade] := gradeQ[m, 0]
vectorQ[m_grade] := gradeQ[m, 1]
bivectorQ[m_grade] := gradeQ[m, 2]
trivectorQ[m_grade] := gradeQ[m, 3]
quadvectorQ[m_grade] := gradeQ[m, 4]
bladeQ[m_grade] := ((m // First) >= 0)
gradeAnyQ[m_grade] := True
gradeAnyQ[_] := False
notGradeQ[v_] := Not[gradeAnyQ[v]]


directProduct[t_, v1_, v2_] := grade[t, (v1 // Last).(v2 // Last)];
signedSymmetric[t_, v1_, v2_, s_] := Module[{a = (v1 // Last), b = (v2 // Last)}, grade[t, (a.b + s b.a)/2]];
symmetric[t_, v1_, v2_] := signedSymmetric[t, v1, v2, 1];
antisymmetric[t_, v1_, v2_] := signedSymmetric[t, v1, v2, -1];

(* These operator on just the Dirac matrix portions x of diracGradeSelect[ ,x] *)
diracGradeSelect[m_, 0] := IdentityMatrix[4] (((m // Tr)/4) // Simplify);

vs[m_, 0] :=           ( m[[1, 1]] + m[[2, 2]] - m[[3, 3]] - m[[4, 4]])/4;
vs[m_, 1] :=           ( m[[1, 4]] + m[[2, 3]] - m[[3, 2]] - m[[4, 1]])/4;
vs[m_, 2] := -complexI (-m[[1, 4]] + m[[2, 3]] + m[[3, 2]] - m[[4, 1]])/4;
vs[m_, 3] :=           ( m[[1, 3]] - m[[2, 4]] - m[[3, 1]] + m[[4, 2]])/4;
diracGradeSelect[m_, 1] := (((vs[m, #] diracGammaMatrix[#]) & /@ (Range[4] - 1) ) // Total) ;

bs[m_, {1, 2}] := -complexI ( m[[4, 4]] - m[[3, 3]] + m[[2, 2]] - m[[1, 1]])/4;
bs[m_, {0, 3}] :=           ( m[[1, 3]] - m[[2, 4]] + m[[3, 1]] - m[[4, 2]])/4;
bs[m_, {0, 2}] := -complexI (-m[[1, 4]] + m[[2, 3]] - m[[3, 2]] + m[[4, 1]])/4;
bs[m_, {0, 1}] :=           ( m[[1, 4]] + m[[2, 3]] + m[[3, 2]] + m[[4, 1]])/4;
bs[m_, {1, 3}] :=           ( m[[1, 2]] - m[[2, 1]] + m[[3, 4]] - m[[4, 3]])/4;
bs[m_, {2, 3}] :=  complexI ( m[[1, 2]] + m[[2, 1]] + m[[3, 4]] + m[[4, 3]])/4;

pairs := {{1, 2}, {0, 3}, {0, 2}, {0, 1}, {1, 3}, {2, 3}};
diracGradeSelect[m_, 2] := (((bs[m, #] diracGammaMatrix[# // First] . diracGammaMatrix[# // Last]) & /@ pairs) // Total) // Simplify ;

ts[m_, {0, 1, 2}] := -complexI (-m[[1, 1]] + m[[2, 2]] + m[[3, 3]] - m[[4, 4]])/4;
ts[m_, {0, 1, 3}] :=           ( m[[1, 2]] - m[[2, 1]] - m[[3, 4]] + m[[4, 3]])/4;
ts[m_, {0, 2, 3}] := -complexI (-m[[1, 2]] - m[[2, 1]] + m[[3, 4]] + m[[4, 3]])/4;
ts[m_, {1, 2, 3}] := -complexI (-m[[1, 3]] - m[[2, 4]] + m[[3, 1]] + m[[4, 2]])/4;

triplets := {{0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3}}
trivectors := diracGammaMatrix[# // First].diracGammaMatrix[#[[2]]].diracGammaMatrix[# // Last] & /@ triplets;
diracGradeSelect[m_, 3] := ((ts[m, triplets[[#]]] trivectors[[#]] & /@ Range[4]) // Total) // Simplify;

diracPseudoscalar = diracGammaMatrix[0].diracGammaMatrix[1].diracGammaMatrix[2].diracGammaMatrix[3];
diracGradeSelect[m_, 4] := diracPseudoscalar (imag[-(m[[1, 3]] + m[[2, 4]] + m[[3, 1]] + m[[4, 2]] ) (1/4)] // Simplify);

diracGradeSelect0 := diracGradeSelect[#, 0] &;
diracGradeSelect1 := diracGradeSelect[#, 1] &;
diracGradeSelect2 := diracGradeSelect[#, 2] &;
diracGradeSelect3 := diracGradeSelect[#, 3] &;
diracGradeSelect4 := diracGradeSelect[#, 4] &;


GradeSelection[m_?scalarQ, 0] := m;
GradeSelection[m_?vectorQ, 1] := m;
GradeSelection[m_?bivectorQ, 2] := m;
GradeSelection[m_?trivectorQ, 3] := m;
GradeSelection[m_?quadvectorQ, 4] := m;
GradeSelection[m_, k_Integer /; k >= 0 && k <= 4] :=
  grade[k, diracGradeSelect[m // Last, k]];
ScalarSelection := GradeSelection[#, 0] &;
VectorSelection := GradeSelection[#, 1] &;
BivectorSelection := GradeSelection[#, 2] &;
TrivectorSelection := GradeSelection[#, 3] &;
QuadvectorSelection := GradeSelection[#, 4] &;


binaryOperator[f_, b_?bladeQ, m_grade] := Total[f[b, #] & /@ (GradeSelection[m, #] & /@ (Range[4 + 1] - 1))]
binaryOperator[f_, m_grade, b_?bladeQ] := Total[f[#, b] & /@ (GradeSelection[m, #] & /@ (Range[4 + 1] - 1))]
binaryOperator[f_, m1_grade, m2_grade] :=
 Total[f[# // First, # //
      Last] & /@ ({GradeSelection[m1, #] & /@ (Range[4 + 1] - 1),
      GradeSelection[m2, #] & /@ (Range[4 + 1] - 1)} // Transpose)]

(*Plus*)
grade /: (v1_?notGradeQ) + grade[k_, v2_] := Scalar[v1] + grade[k, v2];
grade /: grade[0, v1_] + grade[0, v2_] := grade[0, v1 + v2];
grade /: grade[1, v1_] + grade[1, v2_] := grade[1, v1 + v2];
grade /: grade[2, v1_] + grade[2, v2_] := grade[2, v1 + v2];
grade /: grade[3, v1_] + grade[3, v2_] := grade[3, v1 + v2];
grade /: grade[4, v1_] + grade[4, v2_] := grade[4, v1 + v2];
grade /: grade[_, v1_] + grade[_, v2_] := grade[-1, v1 + v2];

(*Times[-1,_]*)
grade /: -grade[k_, v_] := grade[k, -v];

(*Times*)
grade /: (v_?notGradeQ) grade[k_, m_] := grade[k, v m];
grade /: grade[0, s_] grade[k_, m_] := grade[k, s.m];


(*NonCommutativeMultiply*)

grade /: grade[0, s_] ** grade[k_, m_] := grade[k, s.m];
grade /: grade[k_, m_] ** grade[0, s_] := grade[k, s.m];

grade /: grade[4, p_] ** grade[1, m_] := grade[3, p.m];
grade /: grade[1, m_] ** grade[4, p_] := grade[3, m.p];

grade /: grade[4, p_] ** grade[2, m_] := grade[2, p.m];
grade /: grade[2, m_] ** grade[4, p_] := grade[2, m.p];

grade /: grade[4, p_] ** grade[3, m_] := grade[1, p.m];
grade /: grade[3, m_] ** grade[4, p_] := grade[1, m.p];

grade /: grade[4, p_] ** grade[4, m_] := grade[0, p.m];
grade /: grade[4, p_] ** grade[_, m_] := grade[-1, p.m];

grade /: grade[_, m1_] ** grade[_, m2_] := grade[-1, m1.m2];

(* \[Equal]comparison operator *)

grade /: grade[_, m1_] == grade[_, m2_] := (m1 == m2);


(*Dot

v1.v2=(1/2)(v1 v2+v2 v1);v1.v2=v2.v1
v.b=(1/2)(v b-b v);v.b=-b.v
v.t=(1/2)(v t+t v);v.t=t.v
v.q=(1/2)(v q-q v);v.q=-q.v
*)
grade /: (s_?notGradeQ).grade[k_, m_] := grade[k, s m];
grade /: grade[k_, m_].(s_?notGradeQ) := grade[k, s m];
grade /: grade[0, s_].grade[k_, m_] := grade[k, s m];
grade /: grade[k_, m_].grade[0, s_] := grade[k, s m];

grade /: (v1_?vectorQ).grade[1, v2_] := symmetric[0, v1, grade[1, v2]];
grade /: (v_?vectorQ).grade[2, b_] := antisymmetric[1, v, grade[2, b]];
grade /: (v_?vectorQ).grade[3, t_] := symmetric[2, v, grade[3, t]];
grade /: (v_?vectorQ).grade[4, q_] := antisymmetric[3, v, grade[3, q]];

grade /: (b_?bivectorQ).grade[1, v_] := antisymmetric[1, b, grade[1, v]];
grade /: (b_?bivectorQ).grade[2, b2_] := grade[0, (b // Last).b2 // diracGradeSelect0];
grade /: (b_?bivectorQ).grade[3, t_] := grade[1, (b // Last).t // diracGradeSelect1];
grade /: (b_?bivectorQ).grade[4, q_] := grade[1, (b // Last).q // diracGradeSelect2];
grade /: grade[3, t_] . (b_?bivectorQ) := grade[1, t. (b // Last) // diracGradeSelect1];
grade /: grade[4, q_] . (b_?bivectorQ) := grade[1, q. (b // Last) // diracGradeSelect2];

grade /: (t_?trivectorQ).grade[1, v_] := symmetric[2, t, grade[1, v]];
grade /: (t_?trivectorQ).grade[3, t2_] := grade[0, (t // Last).t2 // diracGradeSelect0];
grade /: (t_?trivectorQ).grade[4, q_] := grade[1, (t // Last).q // diracGradeSelect1];
grade /: grade[4, q_].(t_?trivectorQ) := grade[1, q. (t // Last) // diracGradeSelect1];

grade /: (q_?quadvectorQ).grade[1, v_] := antisymmetric[3, q, grade[1, v]];
grade /: (q_?quadvectorQ).m_grade := q ** m;
grade /: m_grade.(q_?quadvectorQ) := m ** q;

(*Dot;handle dot products where one or more factors is a multivector.*)

grade /: grade[g1_, m1_].grade[g2_, m2_] := binaryOperator[Dot, grade[g1, m1], grade[g2, m2]];

grade[_, ConstantArray[0, {4, 4}]] := 0

(* Define a custom wedge operator

v1^v2=(1/2)(v1 v2-v2 v1);v1^v2=-v1.v2
v^b=(1/2)(v b+b v);v^b=v.b
v^t=(1/2)(v t-t v);v^t=-v.t
v^q=(1/2)(v q+q v)=0;v^q=v.q
*)
grade /: grade[0, s_]\[Wedge]grade[k_, v_] := grade[k, s.v];
grade /: grade[k_, v_]\[Wedge]grade[0 _, s_] := grade[k, s.v];

grade /: grade[1, v1_]\[Wedge](v2_?vectorQ) := antisymmetric[2, grade[1, v1], v2];
grade /: grade[1, v_]\[Wedge](b_?bivectorQ) := symmetric[3, grade[1, v], b];
grade /: grade[1, v_]\[Wedge](t_?trivectorQ) := antisymmetric[4, grade[1, v], t];

grade /: grade[2, v1_]\[Wedge](v2_?vectorQ) := symmetric[3, grade[2, v1], v2];
grade /: grade[2, _]\[Wedge](v2_?bivectorQ) := 0;

(*Only e0123^scalar is non zero,and that is handled above*)

grade /: grade[4, _]\[Wedge]b_?bladeQ := 0;
grade /: b_?bladeQ\[Wedge]grade[4, _] := 0;

grade /: grade[g1_, m1_]\[Wedge]grade[g2_, m2_] := binaryOperator[Wedge, grade[g1, m1], grade[g2, m2]];


smagnitude[m_] := m[[1, 1]];
pmagnitude[m_] := m[[1, 3]] complexI;

(* AngleBracket,single operand forms,enter with[Esc]<[Esc] v[Esc]>[Esc] *)
grade /: AngleBracket[grade[0, s_]] := smagnitude[s]
grade /: AngleBracket[grade[1, _]] := 0
grade /: AngleBracket[grade[2, _]] := 0
grade /: AngleBracket[grade[3, _]] := 0
grade /: AngleBracket[grade[4, _]] := 0
grade /: AngleBracket[grade[_, m_]] := ((diracGradeSelect[m, 0]) // smagnitude)

ScalarValue[m_grade] := AngleBracket[m];

(* AngleBracket,two operand forms. *)

grade /: AngleBracket[grade[0, s1_], grade[0, s2_]] := (smagnitude[s1] smagnitude[s2]);
grade /: AngleBracket[grade[0, s1_], grade[-1, m_]] := (smagnitude[ s1] ((diracGradeSelect[m, 0]) // smagnitude));
grade /: AngleBracket[grade[-1, m_], grade[0, s1_]] := (smagnitude[s1] ((diracGradeSelect[m, 0]) // smagnitude));

grade /: AngleBracket[grade[0, s1_], grade[_, _]] := 0;
grade /: AngleBracket[grade[-1, m_], grade[0, s1_]] := (smagnitude[s1] ((diracGradeSelect[m, 0]) // smagnitude));
grade /: AngleBracket[grade[_, _], grade[0, s1_]] := 0;

grade /: AngleBracket[grade[4, q1_], grade[4, q2_]] := (pmagnitude[q1] pmagnitude[q2])
grade /: AngleBracket[grade[4, q1_], grade[-1, m_]] := (pmagnitude[ q1] ((diracGradeSelect[m, 4]) // pmagnitude));
grade /: AngleBracket[grade[4, q1_], grade[_, _]] := 0;

grade /: AngleBracket[grade[-1, m_], grade[4, q1_]] := (smagnitude[q1] ((diracGradeSelect[m, 4]) // pmagnitude));
grade /: AngleBracket[grade[_, _], grade[4, t1_]] := 0;

grade /: AngleBracket[grade[k1_, m1_], grade[k2_, m2_]] := (diracGradeSelect[m1.m2, 0] // smagnitude);



ScalarProduct[m1_grade, m2_grade] := AngleBracket[m1, m2];

gsub = Subscript["\[Gamma]", #] &;
displayMapping = {
   {Scalar[1], 1, 1, 1},
   {Vector[-1, 1], gsub[1], \[Gamma][1], "\\gamma_1"},
   {Vector[-1, 2], gsub[2], \[Gamma][2], "\\gamma_2"},
   {Vector[-1, 3], gsub[3], \[Gamma][3], "\\gamma_3"},
   {Vector[1, 0], gsub[3], \[Gamma][0], "\\gamma_0"},
   {Bivector[-1, 1, 2], gsub["12"], \[Gamma][1] \[Gamma][2], "\\gamma_{12}"},
   {Bivector[-1, 1, 3], gsub["13"], \[Gamma][1] \[Gamma][3], "\\gamma_{13}"},
   {Bivector[-1, 2, 3], gsub["23"], \[Gamma][2] \[Gamma][3], "\\gamma_{23}"},
   {Bivector[1, 1, 0], gsub["10"], \[Gamma][1] \[Gamma][0], "\\gamma_{10}"},
   {Bivector[1, 2, 0], gsub["20"], \[Gamma][2] \[Gamma][0], "\\gamma_{20}"},
   {Bivector[1, 3, 0], gsub["30"], \[Gamma][3] \[Gamma][0], "\\gamma_{30}"},
   {Trivector[1, 1, 2, 3], gsub["123"], \[Gamma][1] \[Gamma][2] \[Gamma][3], "\\gamma_{123}"},
   {Trivector[-1, 1, 2, 0], gsub["120"], \[Gamma][1] \[Gamma][2] \[Gamma][0], "\\gamma_{120}"},
   {Trivector[-1, 1, 3, 0], gsub["130"], \[Gamma][1] \[Gamma][3] \[Gamma][0], "\\gamma_{130}"},
   {Trivector[-1, 2, 3, 0], gsub["230"], \[Gamma][2] \[Gamma][3] \[Gamma][0], "\\gamma_{230}"},
   {Quadvector[1], gsub["0123"], \[Gamma][0] \[Gamma][1] \[Gamma][2] \[Gamma][3], "\\gamma_{0123}"}
};

GAdisplay[v_grade, how_] :=
  Total[(Times[(AngleBracket[# // First, v] (*//
        Simplify*)), #[[how]]]) & /@ displayMapping];

TraditionalForm[m_grade] := ((GAdisplay[m, 2]) // TraditionalForm);
DisplayForm[m_grade] := GAdisplay[m, 2];
Format[m_grade] := GAdisplay[m, 2];
StandardForm[m_grade] := GAdisplay[m, 3];


D[m_grade, u_] :=
  grade[m // First,
   Map[complex[D[# // real // Simplify, u],
      D[# // imag // Simplify, u]] &, m // Last, {2}]];


oneTeXForm[m_, blade_, disp_] := Module[{p},
  p = AngleBracket[blade, m];
  If[ p // PossibleZeroQ, 0,
   If[ p === 1, disp,
    Times[p // TeXForm, disp]]
   ]
  ]

TeXForm[m_grade] := Total[ oneTeXForm[m, # // First, #[[4]]] & /@ displayMapping ];

(* Grad::
usage="grad[m,{x,y,z}] computes the vector product of the gradient with multivector m with respect to coordinates ct,x,y,z..";
*)
grade /: Grad[grade[k_, m_], u_List] := ((Vector[1, #] ** D[grade[k, m], u[[#]]]) & /@ Range[4]) // Total;

(*Div::usage="div[m,{x,y,z}] of a grade k+1 blade m, computes < \[Del] m >_k, where the gradient is evaluated with respect to coordinates ct,x,y,z.";*)

grade /: Div[grade[1, m_], u_List] := Grad[grade[1, m], u] // ScalarSelection;
grade /: Div[grade[2, m_], u_List] := Grad[grade[2, m], u] // VectorSelection;
grade /: Div[grade[3, m_], u_List] := Grad[grade[3, m], u] // BivectorSelection;
grade /: Div[grade[4, m_], u_List] := Grad[grade[4, m], u] // TrivectorSelection;

(* Curl::
usage="Given a grade (k-1) blade m, curl[ m, {x, y, z} ] = < \[Del] m >_k, where the gradient is evaluated with respect to coordinates ct,x,y,z.";
*)

grade /: Curl[grade[1, m_], u_List] := Grad[grade[1, m], u] // BivectorSelection;
grade /: Curl[grade[2, m_], u_List] := Grad[grade[2, m], u] // TrivectorSelection;
grade /: Curl[grade[3, m_], u_List] := Grad[grade[3, m], u] // QuadvectorSelection;
grade /: Curl[grade[4, m_], u_List] := 0;

(* End Private Context *)
End[]

Protect[
   Bivector,
   BivectorSelection,
   Curl,
   D,
   DisplayForm,
   Div,
   Format,
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
   Vector,
   VectorSelection,
   \[Gamma],
   complex,
   complexI,
   complexQ,
   conjugate,
   fMatrix,
   imag,
   matrixconj,
   matriximag,
   matrixreal,
   notComplexQ,
   real
]

EndPackage[ ]
