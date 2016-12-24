(* ::Package:: *)

(* copy this module to a directory in $Path.  Then invoke with <<pauliGA` *)
BeginPackage["pauliGA`"]
pauliGA::usage = "An attempt to use the Pauli representation to compactly implement 3D Euclidean Geometric Algebra operations.";

(*Begin["`Private`"]*)
ClearAll[ scalarType, vectorType, bivectorType, trivectorType, multivectorType ]
{scalarType, vectorType, bivectorType, trivectorType, multivectorType} = Range[5];

complex::usage = "A limited use complex number implementation to use internally in Pauli basis representation, independent of any Complex" ;
ClearAll[ complex, complexQ, notComplexQ, real, imag, conjugate ];

complex[0,0] := 0;
complex /: complex[r1_, i1_] + complex[r2_, i2_] := complex[r1 + r2, i1 + i2 ] ;
complex /: r1_ + complex[r2_, i2_] := complex[r1 + r2, i2 ] ;

complex /: - complex[re_, im_] := complex[-re, -im];

complex /: complex[re_] := complex[re, 0] ;

complex /: complex[re_] := complex[re, 0] ;

complex /: complex[r1_, i1_]  complex[r2_, i2_] := complex[r1 r2 - i1 i2, r1 i2 + r2 i1 ] ;

complex /: Power[z_complex, 2] := complex[z] complex[z] ;

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

(*
(* Unprotect functions that are used to define our rules *)
protected = Unprotect [NonCommutativeMultiply]

ClearAll[ directProduct, signedSymmetric, symmetric, antisymmetric ]

directProduct[ t_, v1_, v2_ ] := { t, (v1 // Last). (v2 // Last) } ;
signedSymmetric[ t_, v1_, v2_, s_ ] := Module[ {a = (v1 // Last), b = (v2 // Last)}, {scalarType, (a . b + s b . a)/2} ] ;
    symmetric[ t_, v1_, v2_ ] := signedSymmetric[ t, v1, v2, 1 ] ;
antisymmetric[ t_, v1_, v2_ ] := signedSymmetric[ t, v1, v2, -1 ] ;

ClearAll[ pScalarQ, pVectorQ, pBivectorQ, pTrivectorQ ]
pScalarQ[m : {scalarType, _}] := True ;
pScalarQ[m : {_Integer, _}] := False ;
pVectorQ[m : {vectorType, _}] := True ;
pVectorQ[m : {_Integer, _}] := False ;
pBivectorQ[m : {bivectorType, _}] := True ;
pBivectorQ[m : {_Integer, _}] := False ;
pTrivectorQ[m : {trivectorType, _}] := True ;
pTrivectorQ[m : {_Integer, _}] := False ;
pMultivectorQ[m : {multivectorType, _}] := True ;
pMultivectorQ[m : {_Integer, _}] := False ;
pBladeQ[m : {multivectorType, _}] := False ;
pBladeQ[m : {_Integer, _}] := True ;

ClearAll[ grade0, grade1, grade2, grade3, grade01, grade23 ]
grade0 := IdentityMatrix[2] ( #/2 // Tr // matrixreal // Simplify) &;
grade3 := complexI IdentityMatrix[2] ( #/2 // Tr // matriximag // Simplify) &;
grade01 := (((# + (# // conjugateTranspose))/2) // Simplify)  &;
grade23 := (((# - (# // conjugateTranspose))/2) // Simplify)  &;
grade1 := ((grade01[#] - grade0[#]) // Simplify) &;
grade2 := ((grade23[#] - grade3[#]) // Simplify) &;

ClearAll[ gaPlus, gaminus, magnitude ]

gaPlus[ v1_, v2_, f1_ : 1, f2_ : 1 ] := Module[
   {g1 = (v1 // First), g2 = (v2 // First), a = (v1 // Last), b = (v2 // Last)},
   { If[ g1 == g2, g1, multivectorType ],
   f1 a + f2 b} ]
gaMinus[ v1_, v2_ ] := gaPlus[ v1, v2, 1, -1 ] ;
gaNegate[ v_ ] := { v // First, - (v // Last) } ;

magnitude[ v_?pScalarQ ] := (v // Last)[[1, 1]] ;
(*End[]*) (* private *)

ClearAll[ Scalar, Vector, Bivector, Trivector ]

Scalar[v_] := {scalarType, v IdentityMatrix[2]};
Vector[v_, k_Integer /; k > 0 && k < 4] := {vectorType, v pauliMatrix[k] };
Bivector[v_, k_Integer /; k > 0 && k < 4, j_Integer  /; j > 0 && j < 4] := {bivectorType, v pauliMatrix[k] . pauliMatrix[j] };
Trivector[v_] := {trivectorType, v complexI IdentityMatrix[2]};

(* There is a built in multiply for non-commutative products:
   http://reference.wolfram.com/language/ref/NonCommutativeMultiply.html
 *)
ClearAll[ DotProduct, WedgeProduct, GeometricProduct, ScalarProduct ]
NonCommutativeMultiply[ v1_?pScalarQ, v2_?pScalarQ]    := directProduct[ scalarType, v1 , v2  ] ;
NonCommutativeMultiply[ v1_?pScalarQ, v2_?pVectorQ]    := directProduct[ vectorType, v1 , v2  ] ;
NonCommutativeMultiply[ v1_?pScalarQ, v2_?pBivectorQ]  := directProduct[ bivectorType, v1 , v2  ] ;
NonCommutativeMultiply[ v1_?pScalarQ, v2_?pTrivectorQ] := directProduct[ trivectorType, v1 , v2  ] ;
NonCommutativeMultiply[ v1_?pScalarQ, v2_?pMultivectorQ] := directProduct[ trivectorType, v1 , v2  ] ;
(* multiplication of (scalar, expressions) pairs: *)
NonCommutativeMultiply[ v1_?pScalarQ, v2_] := directProduct[ trivectorType, v1 v2 ] ;
NonCommutativeMultiply[ v1_, v2_?pScalarQ] := NonCommutativeMultiply[ v2, v1 ] ;

NonCommutativeMultiply[ v1_?pTrivectorQ, v2_?pScalarQ]    := directProduct[ trivectorType, v1 , v2  ] ;
NonCommutativeMultiply[ v1_?pTrivectorQ, v2_?pVectorQ]    := directProduct[ bivectorType, v1 , v2  ] ;
NonCommutativeMultiply[ v1_?pTrivectorQ, v2_?pBivectorQ]  := directProduct[ vectorType, v1 , v2  ] ;
NonCommutativeMultiply[ v1_?pTrivectorQ, v2_?pTrivectorQ] := directProduct[ scalarType, v1 , v2  ] ;
NonCommutativeMultiply[ v1_?pTrivectorQ, v2_?pMultivectorQ] := directProduct[ multivectorType, v1 , v2  ] ;
(* multiplication of (trivector, expression) pairs: *)
NonCommutativeMultiply[ v1_?pTrivectorQ, v2_] := directProduct[ trivectorType, v1 v2 ] ;
NonCommutativeMultiply[ v1_, v2_?pTrivectorQ] := NonCommutativeMultiply[ v2, v1 ] ;



DotProduct[ v1_?pScalarQ, v2_] := NonCommutativeMultiply[ v1, v2 ] ;
DotProduct[ v1_, v2_?pScalarQ] := NonCommutativeMultiply[ v1, v2 ] ;
WedgeProduct[ v1_?pScalarQ, v2_] := NonCommutativeMultiply[ v1, v2 ] ;
WedgeProduct[ v1_, v2_?pScalarQ] := NonCommutativeMultiply[ v1, v2 ] ;

DotProduct[ v1_?pTrivectorQ, v2_] := NonCommutativeMultiply[ v1, v2 ] ;
DotProduct[ v1_, v2_?pTrivectorQ] := NonCommutativeMultiply[ v1, v2 ] ;


DotProduct[ v1_?pVectorQ, v2_?pVectorQ]    :=     symmetric[ scalarType,   v1, v2 ] ;
DotProduct[ v1_?pVectorQ, v2_?pBivectorQ]  := antisymmetric[ vectorType,   v1, v2 ] ;
DotProduct[ v1_?pBivectorQ,  v2_?pVectorQ] := -DotProduct[ v2, v1 ] ;

DotProduct[ v1_?pBivectorQ, v2_?pBivectorQ]  := symmetric[ scalarType,     v1, v2 ] ;

(*
DotProduct[ v1_?pTrivectorQ,  v2_?pTrivectorQ] := directProduct[ scalarType, v1 , v2 ] ;
DotProduct[ v1_?pVectorQ, v2_?pTrivectorQ] :=     symmetric[ bivectorType, v1, v2 ] ;
DotProduct[ v1_?pTrivectorQ, v2_?pVectorQ] :=  DotProduct[ v2, v1 ] ;
DotProduct[ v1_?pBivectorQ, v2_?pTrivectorQ] := directProduct[ vectorType, v1, v2 ] ;
DotProduct[ v1_?pTrivectorQ,  v2_?pBivectorQ] := DotProduct[ v2, v1 ] ;

WedgeProduct[ v1_?pScalarQ, v2_?pScalarQ ]    := DotProduct[ v1, v2 ] ;
WedgeProduct[ v1_?pScalarQ, v2_?pVectorQ ]    := DotProduct[ v1, v2 ] ;
WedgeProduct[ v1_?pScalarQ, v2_?pBivectorQ ]  := DotProduct[ v1, v2 ] ;
WedgeProduct[ v1_?pScalarQ, v2_?pTrivectorQ ] := DotProduct[ v1, v2 ] ;
WedgeProduct[ v1_?pVectorQ,    v2_?pScalarQ ] := DotProduct[ v2, v1 ] ;
WedgeProduct[ v1_?pBivectorQ,  v2_?pScalarQ ] := DotProduct[ v2, v1 ] ;
WedgeProduct[ v1_?pTrivectorQ, v2_?pScalarQ ] := DotProduct[ v2, v1 ] ;

WedgeProduct[ v1_?pVectorQ, v2_?pVectorQ ]    := antisymmetric[ bivectorType,  v1, v2 ] ;
WedgeProduct[ v1_?pVectorQ, v2_?pBivectorQ ]  :=     symmetric[ trivectorType, v1, v2 ] ;
WedgeProduct[ v1_?pVectorQ, v2_?pTrivectorQ ] := Scalar[0] ;
WedgeProduct[ v1_?pBivectorQ,  v2_?pVectorQ ] := WedgeProduct[ v2, v1 ] ;
WedgeProduct[ v1_?pTrivectorQ, v2_?pVectorQ ] := Scalar[0] ;

WedgeProduct[ v1_?pBivectorQ, v2_?pBivectorQ ]  := Scalar[0] ;
WedgeProduct[ v1_?pBivectorQ, v2_?pTrivectorQ ] := Scalar[0] ;
WedgeProduct[ v1_?pTrivectorQ, v2_?pBivectorQ ] := Scalar[0] ;

WedgeProduct[ v1_?pTrivectorQ, v2_?pTrivectorQ ] := Scalar[0] ;
*)

GeometricProduct[ v1_, v2_ ] := directProduct[ multivectorType, v1, v2 ] ;


ClearAll[ GradeSelection, ScalarSelection, PseudoScalarSelection ]
GradeSelection[ v_, n_Integer /; n == 0 ] := Scalar[ (v // Last) // grade0 ] ;
GradeSelection[ v_, n_Integer /; n == 1 ] := Vector[ (v // Last) // grade1 ] ;
GradeSelection[ v_, n_Integer /; n == 2 ] := Bivector[ (v // Last) // grade2 ] ;
GradeSelection[ v_, n_Integer /; n == 3 ] := Trivector[ (v // Last) // grade3 ] ;
ScalarSelection[ v_ ] := GradeSelection[ v, 0 ] ;
VectorSelection[ v_ ] := GradeSelection[ v, 1 ] ;
BivectorSelection[ v_ ] := GradeSelection[ v, 2 ] ;
PseudoScalarSelection[ v_ ] := GradeSelection[ v, 3 ] ;

ScalarMagnitude[ v_ ] := (GradeSelection[ v, 0 ] // magnitude) ;

(*
DotProduct[ v1_?pMultivectorQ, v2_?pScalarQ ] := directproduct[ multivectorType, v1, v2 ] ;
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

Protect[Evaluate[protected]]   (* Restore protection of the functions *)

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
