(* ::Package:: *)

(* copy this module to a directory in $Path.  Then invoke with <<pauliGA` *)
BeginPackage["pauliGA`"]
pauliGA::usage = "An attempt to use the Pauli representation to compactly but efficiently implement 3D Euclidean Geometric Algebra operations.";

(*Begin["`Private`"]*)
complex::usage = 
  "A limited use complex number implementation to use internally in \
Pauli basis representation, independent of any Complex";
ClearAll[complex, complexQ, notComplexQ, real, imag, conjugate];

(*complex[0,0]:=0;*)

complex /: complex[r1_, i1_] + complex[r2_, i2_] := 
  complex[r1 + r2, i1 + i2];
complex /: r1_ + complex[r2_, i2_] := complex[r1 + r2, i2];

complex /: -complex[re_, im_] := complex[-re, -im];

complex /: complex[re_] := re;
complex /: complex[re_, 0] := re;

(*complex/:complex[re_]:=complex[re,0];*)

complex /: complex[r1_, i1_] complex[r2_, i2_] := 
  complex[r1 r2 - i1 i2, r1 i2 + r2 i1];

norm::usage = "Norm like function for complex[]";
norm[z_complex] := ((z // First)^2 + (z // Last)^2);

(*special case this one to deal with the sort of products that are \
generated multiplying pauli matrices*)

complex /: Power[z_complex, 2] := complex[z] complex[z];
complex /: Power[z_complex, n_] := 
  Module[{r = norm[z]^(n/2), theta = n ArcTan[z // First, z // Last]},
    r complex[Cos[theta], Sin[theta]]];

complexQ::usage = "predicate pattern match for complex[]";
complexQ[z_complex] := True;
complexQ[_] := False;
notComplexQ::usage = "predicate pattern match for !complex[]";
notComplexQ[v_] := Not[complexQ[v]];

complex /: (v_?notComplexQ) complex[re_, im_] := complex[v re, v im];

real::usage = "Re like function for complex[]";
real[z_complex] := (z // First);
imag::usage = "Im like function for complex[]";
imag[z_complex] := (z // Last);
real[ex_] := ex;
imag[ex_] := 0;
conjugate::usage = "Conjugate like function for complex[]";
conjugate[z_complex] := complex[z // First, -z // Last];
conjugate[ex_] := ex;

ClearAll[complexI, fMatrix, matrixreal, matriximag, matrixconj]
complexI::usage = "I like unit imaginary for complex[]";
complexI := complex[0, 1];

fMatrix[p_, f_] := (Function[a, f@a, Listable]@p)
matrixreal::usage = 
  "method to apply real to all elements in matrix.  This is a hack.  \
Can probably set an attribute on the real function to do this.";
matrixreal[m_] := fMatrix[m, real];
matriximag::usage = 
  "method to apply imag to all elements in matrix.  This is a hack.  \
Can probably set an attribute on the imag function to do this.";
matriximag[m_] := fMatrix[m, imag];
matrixconj::usage = 
  "method to apply conjugate to all elements in matrix.  This is a \
hack.  Can probably set an attribute on the conj function to do this.";
matrixconj[m_] := fMatrix[m, conjugate];

ClearAll[pauliMatrix, conjugateTranspose]
pauliMatrix::usage = 
  "PauliMatrix implemented with complex[], instead of Complex[].";
pauliMatrix[1] := PauliMatrix[1];
pauliMatrix[
   2] := (PauliMatrix[2] /. {Complex[0, 1] -> complexI, 
     Complex[0, -1] -> -complexI});
pauliMatrix[3] := PauliMatrix[3];
conjugateTranspose::usage = 
  "ConjugateTranspose like operation for pauliMatrix.";
conjugateTranspose[m_List] := Transpose[matrixconj[m]];
(*End["`Private`"]*)

Unprotect[TraditionalForm, DisplayForm, StandardForm];
TraditionalForm[
  z_complex] := (((z // real) + I (z // imag)) // TraditionalForm)
DisplayForm[
  z_complex] := (((z // real) + I (z // imag)) // DisplayForm)
StandardForm[
  z_complex] := (((z // real) + I (z // imag)) // StandardForm)
Protect[TraditionalForm, DisplayForm, StandardForm];

(* End of complex, and pauliMatrix section.  Define the basic CL(3,0) operations. *)

grade::usage = "(internal) An upvalue type that represents a CL(3,0) algebraic element as a pair {grade, v}, where v is a sum of products of Pauli matrices.  These matrices may be scaled by arbitrary numeric or symbolic factors." ;
ClearAll[Vector, Scalar, Bivector, Trivector, grade]
Scalar::usage = "Scalar[v] constructs a scalar grade quantity with value v." ;
Scalar[v_] := grade[0, v IdentityMatrix[2]];
Vector::usage = "Vector[v, n], where n = {1,2,3} constructs a vector grade quantity with value v in direction n." ;
Vector[v_, k_Integer /; k >= 1 && k <= 3] := 
  grade[1, v pauliMatrix[k]];
Bivector::usage = "Bivector[v, n1, n2], where n1,n2 = {1,2,3} constructs a bivector grade quantity with value v in the plane n1,n2." ;
Bivector[v_, k_Integer /; k >= 1 && k <= 3, 
   j_Integer /; j >= 1 && j <= 3] := 
  grade[2, v pauliMatrix[k].pauliMatrix[j]];
Bivector::usage = "Trivector[v] constructs a trivector (pseudoscalar) grade quantity scaled by v." ;
Trivector[v_] := grade[3, complexI v IdentityMatrix[2]];

(*Begin["`Private`"]*)
ClearAll[scalarQ, vectorQ, bivectorQ, trivectorQ, bladeQ]
gradeQ::usage = "gradeQ[m, n] tests if the multivector m is of grade n.  n = -1 is used internally to represent values of more than one grade.";
gradeQ[m_grade, n_Integer] := ((m // First) == n)
scalarQ::usage = "scalarQ[m] tests if the multivector m is of grade 0 (scalar)" ;
scalarQ[m_grade] := gradeQ[m, 0]
vectorQ::usage = "vectorQ[m] tests if the multivector m is of grade 1 (vector)" ;
vectorQ[m_grade] := gradeQ[m, 1]
bivectorQ::usage = "bivectorQ[m] tests if the multivector m is of grade 2 (bivector)" ;
bivectorQ[m_grade] := gradeQ[m, 2]
trivectorQ::usage = "trivectorQ[m] tests if the multivector m is of grade 3 (trivector)" ;
trivectorQ[m_grade] := gradeQ[m, 3]
bladeQ::usage = "bladeQ[m] tests if the multivector is of a single grade." ;
bladeQ[m_grade] := ((m // First) >= 0)

ClearAll[directProduct, signedSymmetric, symmetric, antisymmetric]

directProduct[t_, v1_, v2_] := grade[t, (v1 // Last).(v2 // Last)];
signedSymmetric[t_, v1_, v2_, s_] := 
  Module[{a = (v1 // Last), b = (v2 // Last)}, 
   grade[t, (a.b + s b.a)/2]];
symmetric[t_, v1_, v2_] := signedSymmetric[t, v1, v2, 1];
antisymmetric[t_, v1_, v2_] := signedSymmetric[t, v1, v2, -1];

(*These operator on just the Pauli matrix portions x of \
pauliGradeSelect[,x]*)
ClearAll[pauliGradeSelect01, \
pauliGradeSelect23, pauliGradeSelect]
pauliGradeSelect01 := ((# + (# // conjugateTranspose))/2) &;
pauliGradeSelect23 := ((# - (# // conjugateTranspose))/2) &;
pauliGradeSelect[m_, 0] := 
  IdentityMatrix[2] (m/2 // Tr // real // Simplify);
pauliGradeSelect[m_, 
   1] := ((pauliGradeSelect01[m] - pauliGradeSelect[m, 0]) // 
    Simplify);
pauliGradeSelect[m_, 
   2] := ((pauliGradeSelect23[m] - pauliGradeSelect[m, 3]) // 
    Simplify);
pauliGradeSelect[m_, 3] := 
  complexI IdentityMatrix[2] (m/2 // Tr // imag // Simplify);

ClearAll[pauliGradeSelect0, pauliGradeSelect1, pauliGradeSelect2, \
pauliGradeSelect3]
pauliGradeSelect0 := pauliGradeSelect[#, 0] &;
pauliGradeSelect1 := pauliGradeSelect[#, 1] &;
pauliGradeSelect2 := pauliGradeSelect[#, 2] &;
pauliGradeSelect3 := pauliGradeSelect[#, 3] &;
(*End["`Private`"]*)


(* Plus *)
grade /: grade[0, v1_] + grade[0, v2_] := grade[0, v1 + v2];
grade /: grade[1, v1_] + grade[1, v2_] := grade[1, v1 + v2];
grade /: grade[2, v1_] + grade[2, v2_] := grade[2, v1 + v2];
grade /: grade[3, v1_] + grade[3, v2_] := grade[3, v1 + v2];
grade /: grade[_, v1_] + grade[_, v2_] := grade[-1, v1 + v2];

(* Times[-1, _] *)
grade /: -grade[k_, v_] := grade[k, -v];

(* Times *)
grade /: grade[0, s_] grade[k_, m_] := grade[k, s*m];
grade /: grade[3, p_] grade[1, m_] := grade[2, p*m];
grade /: grade[3, p_] grade[2, m_] := grade[1, p*m];
grade /: grade[3, p_] grade[3, m_] := grade[0, p*m];
grade /: grade[3, p_] grade[_, m_] := grade[-1, p*m];

(* NonCommutativeMultiply *)

grade /: grade[_, m1_] ** grade[_, m2_] := grade[-1, m1.m2];

(* Dot *)
grade /: grade[0, s_].grade[k_, m_] := grade[k, s*m];
grade /: grade[k_, m_].grade[0, s_] := grade[k, s*m];

grade /: (t_?trivectorQ).m_grade := t m;
grade /: m_grade.(t_?trivectorQ) := t m;

grade /: (v1_?vectorQ).grade[1, v2_] := symmetric[0, v1, grade[1, v2]];
grade /: (v_?vectorQ).grade[2, b_] := antisymmetric[0, v, grade[2, b]];
grade /: (b_?bivectorQ).grade[1, v_] := 
  antisymmetric[0, b, grade[1, v]];
grade /: (b1_?bivectorQ).grade[2, b2_] := 
  symmetric[0, b1, grade[2, b2]];

ClearAll[GradeSelection, ScalarSelection, VectorSelection, \
BivectorSelection, TrivectorSelection]
GradeSelection[m_?scalarQ, 0] := m;
GradeSelection[m_?vectorQ, 1] := m;
GradeSelection[m_?bivectorQ, 2] := m;
GradeSelection[m_?trivectorQ, 3] := m;
GradeSelection[m_, k_Integer /; k >= 0 && k <= 3] := 
  grade[k, pauliGradeSelect[m // Last, k]];
ScalarSelection := GradeSelection[#, 0] &;
VectorSelection := GradeSelection[#, 1] &;
BivectorSelection := GradeSelection[#, 2] &;
TrivectorSelection := GradeSelection[#, 3] &;

(* Dot ; handle dot products where one or more factors is a \
multivector.
Fixme: there's probably a fancier "Mathematica" way to distribute the \
sums of all these dots without temporaries. *)
grade /: 
 grade[k_, m_] . v2_?bladeQ := 
 Module[{g0, g1, g2, g3}, {g0, g1, g2, g3} = 
   GradeSelection[grade[k, m], #] & /@ (Range[4] - 1);
  Total[{g0.v2, g1.v2, g2.v2, g3.v2}]]

grade /: (v1_?bladeQ ) . grade[k_, m_] := 
 Module[{g0, g1, g2, g3}, {g0, g1, g2, g3} = 
   GradeSelection[grade[k, m], #] & /@ (Range[4] - 1);
  Total[{v1.g0, v1.g1, v1.g2, v1.g3}] ]

grade /: m1_ . grade[k_, m_]  := 
 Module[{g0, g1, g2, g3}, {g0, g1, g2, g3} = 
   GradeSelection[grade[k, m], #] & /@ (Range[4] - 1);
  Total[{m1.g0, m1.g1, m1.g2, m1.g3}] ]

grade[_, {{0, 0}, {0, 0}}] := 0

(*Define a custom wedge operator*)

grade /: grade[0, s_]\[Wedge]grade[k_, v_] := grade[k, s.v];
grade /: grade[1, v1_]\[Wedge](v2_?vectorQ) := 
  antisymmetric[2, grade[1, v1], v2];

grade /: grade[1, v1_]\[Wedge](v2_?bivectorQ) := 
  symmetric[3, grade[1, v1], v2];
grade /: grade[2, v1_]\[Wedge](v2_?vectorQ) := 
  symmetric[3, grade[2, v1], v2];
grade /: grade[2, _]\[Wedge](v2_?bivectorQ) := 0;

grade /: grade[_, _]\[Wedge](v2_?trivectorQ) := 0;
grade /: grade[3, _]\[Wedge]m_ := 0;

ClearAll[pmagnitude]

(*Begin["`Private`"]*)
pmagnitude::usage = 
  "select the 1,1 element from a pauli matrix assuming it represents \
either a Scalar, or a Trivector (i.e. scaled diagonal matrix)." ;
pmagnitude[m_] := m[[1, 1]];
(*End["`Private`"]*)

(* AngleBracket,single operand forms, enter with[Esc]<[Esc] \
v[Esc]>[Esc] *)
grade /: AngleBracket[grade[0, s_]] := pmagnitude[s]
grade /: AngleBracket[grade[1, _]] := 0
grade /: AngleBracket[grade[2, _]] := 0
grade /: AngleBracket[grade[3, _]] := 0
grade /: AngleBracket[
  grade[_, m_]] := ((pauliGradeSelect[m, 0]) // pmagnitude)

(* AngleBracket,two operand forms. *)

grade /: AngleBracket[grade[0, s1_], 
   grade[0, s2_]] := (pmagnitude[s1] pmagnitude[s2]);
grade /: AngleBracket[grade[0, s1_], 
   grade[-1, 
    m_]] := (pmagnitude[
     s1] ((pauliGradeSelect[m, 0]) // pmagnitude));
grade /: AngleBracket[grade[0, s1_], grade[_, _]] := 0;
grade /: AngleBracket[grade[-1, m_], 
   grade[0, 
    s1_]] := (pmagnitude[s1] ((pauliGradeSelect[m, 0]) // pmagnitude));
grade /: AngleBracket[grade[_, _], grade[0, s1_]] := 0;

grade /: AngleBracket[grade[3, t1_], 
  grade[3, t2_]] := (pmagnitude[t1] pmagnitude[t2])
grade /: AngleBracket[grade[3, t1_], 
   grade[-1, 
    m_]] := (pmagnitude[
     t1] ((pauliGradeSelect[m, 3]) // pmagnitude));
grade /: AngleBracket[grade[3, t1_], grade[_, _]] := 0;
grade /: AngleBracket[grade[-1, m_], 
   grade[3, 
    t1_]] := (pmagnitude[t1] ((pauliGradeSelect[m, 3]) // pmagnitude));
grade /: AngleBracket[grade[_, _], grade[3, t1_]] := 0;

grade /: AngleBracket[grade[k1_, m1_], 
   grade[k2_, m2_]] := (pauliGradeSelect[m1.m2, 0] // pmagnitude);


(*Begin["`Private`"]*)
ClearAll[displayMapping, bold, esub, GAdisplay]
bold = Style[#, Bold] &;
esub = Subscript[bold["e"], #] &;
displayMapping = {{Scalar[1], 1}, {Vector[1, 1], 
    esub[1]}, {Vector[1, 2], esub[2]}, {Vector[1, 3], 
    esub[3]}, {Bivector[1, 2, 1], esub["12"]}, {Bivector[1, 3, 2], 
    esub["23"]}, {Bivector[1, 1, 3], esub["31"]}, {Trivector[-1], 
    esub["123"]}};


GAdisplay[v_grade] := 
  Total[(Times[AngleBracket[# // First, v], # // Last]) & /@ 
    displayMapping];
(*End["`Private`"]*)

Unprotect[TraditionalForm, DisplayForm];
TraditionalForm[m_grade] := ( m // GAdisplay) ;
DisplayForm[m_grade] := ( m // GAdisplay) ;
Protect[TraditionalForm, DisplayForm];


(*
Protect[ Scalar, Vector, Bivector, Trivector,
GradeSelection, ScalarSelection, VectorSelection, BivectorSelection, TrivectorSelection
]
*)

EndPackage[]
