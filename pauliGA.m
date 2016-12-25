(* ::Package:: *)

(* copy this module to a directory in $Path.  Then invoke with <<pauliGA` *)
BeginPackage["pauliGA`"]
pauliGA::usage = "An attempt to use the Pauli representation to compactly implement 3D Euclidean Geometric Algebra operations.";

(*Begin["`Private`"]*)
ClearAll[ scalarType, vectorType, bivectorType, trivectorType, multivectorType ]
{scalarType, vectorType, bivectorType, trivectorType, multivectorType} = Range[5];

complex::usage = "A limited use complex number implementation to use internally in Pauli basis representation, independent of any Complex" ;
ClearAll[ complex, complexQ, notComplexQ, real, imag, conjugate ];

(*complex[0,0] := 0;*)
complex /: complex[r1_, i1_] + complex[r2_, i2_] := complex[r1 + r2, i1 + i2 ] ;
complex /: r1_ + complex[r2_, i2_] := complex[r1 + r2, i2 ] ;

complex /: - complex[re_, im_] := complex[-re, -im];

complex /: complex[re_] := re ;
complex /: complex[re_, 0] := re ;

(*complex /: complex[re_] := complex[re, 0] ;*)

complex /: complex[r1_, i1_]  complex[r2_, i2_] := complex[r1 r2 - i1 i2, r1 i2 + r2 i1 ] ;

norm[z_complex] := ((z//First)^2 + (z//Last)^2) ;

(* special case this one to deal with the sort of products that are generated multiplying pauli matrices *)
complex /: Power[z_complex, 2] := complex[z] complex[z] ;
complex /: Power[z_complex, n_] :=
  Module[{r = norm[z]^(n/2), theta = n ArcTan[z // First, z // Last]},
    r complex[Cos[theta], Sin[theta]]] ;

complexQ[z_complex] := True ;
complexQ[_] := False ;
notComplexQ[v_] := Not[complexQ[v]] ;

complex /: (v_?notComplexQ ) complex[re_, im_] := complex[ v re, v im ];

real[z_complex] := (z // First) ;
imag[z_complex] := (z // Last) ;
real[ex_] := ex ;
imag[ex_] := 0 ;
conjugate[z_complex] := complex[z // First, -z // Last] ;
conjugate[ex_] := ex ;

ClearAll[ complexI, fMatrix, matrixreal, matriximag, matrixconj ]
complexI := complex[0, 1] ;

fMatrix[p_, f_] := (Function[a, f@a, Listable]@p)
matrixreal[m_] := fMatrix[m, real];
matriximag[m_] := fMatrix[m, imag];
matrixconj[m_] := fMatrix[m, conjugate];

ClearAll[pauliMatrix, conjugateTranspose]
pauliMatrix[1] := PauliMatrix[1] ;
pauliMatrix[2] := (PauliMatrix[2] /. {Complex[0, 1] -> complexI, Complex[0, -1] -> -complexI}) ;
pauliMatrix[3] := PauliMatrix[3] ;
conjugateTranspose[m_List] := Transpose[matrixconj[m]] ;

Unprotect[TraditionalForm, DisplayForm, StandardForm];
TraditionalForm[z_complex] := (((z // real) + I (z // imag)) // TraditionalForm)
DisplayForm[z_complex] := (((z // real) + I (z // imag)) // DisplayForm)
StandardForm[z_complex] := (((z // real) + I (z // imag)) // StandardForm)
Protect[TraditionalForm, DisplayForm, StandardForm];


ClearAll[vector, scalar, bivector, multivector, grade]
scalar[v_]                                                                 := grade[0, v IdentityMatrix[2] ] ;
vector[v_, k_Integer /; k >= 1 && k <= 3]                                  := grade[1, v pauliMatrix[k]] ;
bivector[v_, k_Integer /; k >= 1 && k <= 3, j_Integer /; j >= 1 && j <= 3] := grade[2, v pauliMatrix[k].pauliMatrix[j]];
trivector[v_]                                                              := grade[3, complexI v IdentityMatrix[2] ] ;

ClearAll[ scalarQ, vectorQ, bivectorQ, trivectorQ ]
gradeQ[m_grade, n_Integer] := ((m // First) == n)
scalarQ[m_grade]           := gradeQ[m, 0]
vectorQ[m_grade]           := gradeQ[m, 1]
bivectorQ[m_grade]         := gradeQ[m, 2]
trivectorQ[m_grade]        := gradeQ[m, 3]

ClearAll[ directProduct, signedSymmetric, symmetric, antisymmetric ]

directProduct[ t_, v1_, v2_ ] := grade[ t, (v1 // Last). (v2 // Last) ] ;
signedSymmetric[ t_, v1_, v2_, s_ ] := Module[ {a = (v1 // Last), b = (v2 // Last)}, grade[t, (a . b + s b . a)/2] ] ;
    symmetric[ t_, v1_, v2_ ] := signedSymmetric[ t, v1, v2, 1 ] ;
antisymmetric[ t_, v1_, v2_ ] := signedSymmetric[ t, v1, v2, -1 ] ;

(* These operator on just the Pauli matrix portions x of pauliGradeSelect[, x] *)
ClearAll[ pauliGradeSelect01, pauliGradeSelect23, pauliGradeSelect ]
pauliGradeSelect01 := ((# + (# // conjugateTranspose))/2) & ;
pauliGradeSelect23 := ((# - (# // conjugateTranspose))/2) & ;
pauliGradeSelect[m_, 0] := IdentityMatrix[2] ( m/2 // Tr // real // Simplify) ;
pauliGradeSelect[m_, 1] := ((pauliGradeSelect01[m] - pauliGradeSelect[m, 0]) // Simplify) ;
pauliGradeSelect[m_, 2] := ((pauliGradeSelect23[m] - pauliGradeSelect[m, 3]) // Simplify) ;
pauliGradeSelect[m_, 3] := complexI IdentityMatrix[2] ( m/2 // Tr // imag // Simplify) ;

ClearAll[ pauliGradeSelect0, pauliGradeSelect1, pauliGradeSelect2, pauliGradeSelect3 ]
pauliGradeSelect0 := pauliGradeSelect[#, 0] & ;
pauliGradeSelect1 := pauliGradeSelect[#, 1] & ;
pauliGradeSelect2 := pauliGradeSelect[#, 2] & ;
pauliGradeSelect3 := pauliGradeSelect[#, 3] & ;

(* Addition and Subtraction operations *)
grade /: grade[0, v1_] + grade[0, v2_] := grade[0, v1 + v2] ;
grade /: grade[1, v1_] + grade[1, v2_] := grade[1, v1 + v2] ;
grade /: grade[2, v1_] + grade[2, v2_] := grade[2, v1 + v2] ;
grade /: grade[3, v1_] + grade[3, v2_] := grade[3, v1 + v2] ;
grade /: grade[_, v1_] + grade[_, v2_] := grade[-1, v1 + v2] ;
grade /: - grade[k_, v_] := grade[k, -v ];

ClearAll[ GradeSelection, ScalarSelection, VectorSelection, BivectorSelection, TrivectorSelection ]
GradeSelection[m_?scalarQ, 0] := m ;
GradeSelection[m_?vectorQ, 1] := m ;
GradeSelection[m_?bivectorQ, 2] := m ;
GradeSelection[m_?trivectorQ, 3] := m ;
GradeSelection[m_, k_Integer /; k >= 1 && k <= 3] := grade[k, pauliGradeSelect[m // Last, k]] ;
ScalarSelection := GradeSelection[ #, 0 ] &;
VectorSelection := GradeSelection[ #, 1 ] &;
BivectorSelection := GradeSelection[ #, 2 ] &;
TrivectorSelection  := GradeSelection[ #, 3 ] &;

(* Unprotect functions that are used to define our rules *)
protected = Unprotect [Times, NonCommutativeMultiply] ;

Times[ s_?scalarQ, m_grade ]         := directProduct[ m // First, s, m  ] ;
Times[ t_?trivectorQ, m_?vectorQ]    := directProduct[ 2, t, m ] ;
Times[ t_?trivectorQ, m_?bivectorQ]  := directProduct[ 1, t, m ] ;
Times[ t_?trivectorQ, m_?trivectorQ] := directProduct[ 0, t, m ] ;
Times[ t_?trivectorQ, m_]            := directProduct[ -1, t, m ] ;

(* There is a built in multiply for non-commutative products:
   http://reference.wolfram.com/language/ref/NonCommutativeMultiply.html
 *)

ClearAll[ GeometricProduct ]
GeometricProduct[ m1_grade, m2_grade ]        := directProduct[ -1, m1, m2 ] ;

NonCommutativeMultiply[ v1_?vectorQ, v2_?vectorQ ] := GeometricProduct[ v1, v2 ] ;
NonCommutativeMultiply[ v1_?vectorQ, v2_?bivectorQ ] := GeometricProduct[ v1, v2 ] ;
NonCommutativeMultiply[ v1_?bivectorQ, v2_?vectorQ ] := GeometricProduct[ v1, v2 ] ;

Protect[Evaluate[protected]] ; (* Restore protection of the functions *)

ClearAll[ DotProduct, WedgeProduct ]

DotProduct[ s_?scalarQ, m_ ]                  := directProduct[ m // First, s, m  ] ;
DotProduct[ m_, s_?scalarQ ]                  := DotProduct[ s, m ] ;

DotProduct[ v_?vectorQ, b_?vectorQ ]          := symmetric[ 0, v, b  ] ;
DotProduct[ v_?vectorQ, b_?bivectorQ ]        := antisymmetric[ 0, v, b  ] ;
DotProduct[ b_?bivectorQ, v_?vectorQ ]        := -DotProduct[ v, b ] ;

DotProduct[ t_?trivectorQ, m_?vectorQ ]       := directProduct[ 2, t, m  ] ;
DotProduct[ t_?trivectorQ, m_?bivectorQ ]     := directProduct[ 1, t, m  ] ;
DotProduct[ t_?trivectorQ, m_?trivectorQ ]    := directProduct[ 0, t, m  ] ;
DotProduct[ m_, t_?trivectorQ ]               := DotProduct[ t, m ] ;

DotProduct[ b1_?bivectorQ, b2_?bivectorQ ]    := symmetric[ 0, b1, b2  ] ;

WedgeProduct[ s_?scalarQ, m_ ]                := directProduct[ m // First, s, m  ] ;
WedgeProduct[ m_, s_?scalarQ ]                := WedgeProduct[ s, m ] ;
WedgeProduct[ v_?vectorQ, b_?vectorQ ]        := antisymmetric[ 0, v, b  ] ;
WedgeProduct[ v_?vectorQ, b_?bivectorQ ]      := symmetric[ 0, v, b  ] ;
WedgeProduct[ b_?bivectorQ, v_?vectorQ ]      := WedgeProduct[ v, b ] ;
WedgeProduct[ b1_?bivectorQ, b2_?bivectorQ ]  := 0 ;
WedgeProduct[ t_?trivectorQ, m_ ] := 0 ;
WedgeProduct[ m_, t_?trivectorQ ] := 0 ;

ClearAll[ magnitude, scalarProduct ]

(* private *)
magnitude[ v_?scalarQ ] := (v // Last)[[1, 1]] ;

ScalarProduct[ s1_?scalarQ, s2_?scalarQ ]       := ( magnitude[s1] magnitude[s2] ) ;
ScalarProduct[ s_?scalarQ, v_?vectorQ ]         := 0
ScalarProduct[ s_?scalarQ, b_?bivectorQ ]       := 0
ScalarProduct[ s_?scalarQ, t_?trivectorQ ]      := 0

ScalarProduct[ v_?vectorQ, s_?scalarQ ]         := 0
ScalarProduct[ v1_?vectorQ, v2_?vectorQ ]       := ( GeometricProduct[ v1, v2 ] // magnitude ) ;
ScalarProduct[ v_?vectorQ, b_?bivectorQ ]       := 0
ScalarProduct[ v_?vectorQ, t_?trivectorQ ]      := 0

ScalarProduct[ b_?bivectorQ, s_?scalarQ ]       := 0
ScalarProduct[ b_?bivectorQ, v_?vectorQ ]       := 0
ScalarProduct[ b1_?bivectorQ, b2_?bivectorQ ]   := ( GeometricProduct[ b1, b2 ] // ScalarSelection // magnitude ) ;
ScalarProduct[ b_?bivectorQ, t_?trivectorQ ]    := 0

ScalarProduct[ t_?trivectorQ, s_?scalarQ ]      := 0
ScalarProduct[ t_?trivectorQ, v_?vectorQ ]      := 0
ScalarProduct[ t_?trivectorQ, b_?bivectorQ ]    := 0
ScalarProduct[ t1_?trivectorQ, t2_?trivectorQ ] := ( GeometricProduct[ t1, t2 ] // magnitude ) ;

ScalarProduct[ s_?scalarQ, m_grade]             := ( magnitude[s] (m // ScalarSelection // magnitude) ) ;
ScalarProduct[ v_?vectorQ, m_grade]             := ( GeometricProduct[v, (m // VectorSelection)] // magnitude ) ;
ScalarProduct[ b_?bivectorQ, m_grade]           := ( GeometricProduct[b, (m // BivectorSelection)] // magnitude ) ;
ScalarProduct[ t_?trivectorQ, m_grade]          := ( GeometricProduct[t, (m // TrivectorSelection)] // magnitude ) ;

ScalarProduct[ m_grade, s_?scalarQ ]             := ScalarProduct[ s, m ] ;
ScalarProduct[ m_grade, v_?vectorQ ]             := ScalarProduct[ v, m ] ;
ScalarProduct[ m_grade, b_?bivectorQ ]           := ScalarProduct[ b, m ] ;
ScalarProduct[ m_grade, t_?trivectorQ ]          := ScalarProduct[ t, m ] ;

(*
(*End[]*) (* private *)

ScalarMagnitude[ v_ ] := (GradeSelection[ v, 0 ] // magnitude) ;

(*
DotProduct[ v1_?pMultivectorQ, v2_?scalarQ ] := directproduct[ multivectorType, v1, v2 ] ;
DotProduct[ v1_?pMultivectorQ, v2_?pBladeQ ] := Module[{g0, g1, g2, g3},
   {g0, g1, g2, g3} = GradeSelection[v1, #] &/@ (Range[4]-1);
   gaPlus[
      gaPlus[
         gaPlus[
            DotProduct[ g0, v2 ],
            DotProduct[ g1, v2 ]
         ],
         DotProduct[ g2, v2 ]
         ],
      DotProduct[ g3, v2 ]
   ]
]
DotProduct[ v1_, v2_?pMultivectorQ ] := Module[{g0, g1, g2, g3},
   {g0, g1, g2, g3} = GradeSelection[v2, #] &/@ (Range[4]-1);
   gaPlus[
      gaPlus[
         gaPlus[
            DotProduct[ v1, g0 ],
            DotProduct[ v1, g1 ]
         ],
         DotProduct[ v1, g2 ]
         ],
      DotProduct[ v1, g3 ]
   ]
]
*)


(*Begin["`Private`"]*)
ClearAll[ displayMapping, bold, esub ]
bold = Style[ #, Bold] &;
esub = Subscript[ bold["e"], # ] & ;
displayMapping = {
   {Scalar[ 1 ], 1},
   {Vector[ 1, 1 ], esub[1] },
   {Vector[ 1, 2 ], esub[2] },
   {Vector[ 1, 3 ], esub[3] },
(* -1 : reversal for dot product selection *)
   {Bivector[ -1, 1, 2 ], esub["12"] },
   {Bivector[ -1, 2, 3 ], esub["23"] },
   {Bivector[ -1, 3, 1 ], esub["31"] },
   {Trivector[ -1 ], esub["123"] }
} ;
(*End[]  private *)

GAdisplay[v_] := Total[ (Times[DotProduct[# // First, v] // ScalarMagnitude, # // Last]) & /@ displayMapping] ;
*)

(*
Protect[ Scalar, Vector, Bivector, Trivector,
         GeometricProduct, DotProduct, WedgeProduct, GradeSelection, ScalarMagnitude,
         ScalarSelection, VectorSelection, BivectorSelection, PseudoScalarSelection ]
*)

EndPackage[]
