(* ::Package:: *)

Begin["`Private`"]
ClearAll[ scalarType, vectorType, bivectorType, trivectorType, multivectorType ]
{scalarType, vectorType, bivectorType, trivectorType, multivectorType} = Range[5];

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
grade0 := IdentityMatrix[2] ( #/2 // Tr // Re // Simplify) &;
grade3 := I IdentityMatrix[2] ( #/2 // Tr // Im // Simplify) &;
grade01 := (((# + (# // ConjugateTranspose))/2) // Simplify)  &;
grade23 := (((# - (# // ConjugateTranspose))/2) // Simplify)  &;
grade1 := ((grade01[#] - grade0[#]) // Simplify) &;
grade2 := ((grade23[#] - grade3[#]) // Simplify) &;

ClearAll[ gaPlus, gaminus, magnitude ]

gaPlus[ v1_, v2_, f1_ : 1, f2_ : 1 ] := Module[ 
   {g1 = (v1 // First), g2 = (v2 // First), a = (v1 // Last), b = (v2 // Last)},
   { If[ g1 == g2, g1, multivectorType ],
   f1 a + f2 b} ]
gaMinus[ v1_, v2_ ] := gaPlus[ v1, v2, 1, -1 ] ;
gaNegate[ v_ ] := { v // First, - (v // Last) } ;

magnitude[ v_?pScalarQ ] := v[[1,1]] ;
End[]

(* copy this module to a directory in $Path.  Then invoke with <<peeters` *)
Begin["pauliGA`"]

ClearAll[ Scalar, Vector, Bivector, Trivector ]

Scalar[v_] := {scalarType, v IdentityMatrix[2]};
Vector[v_, k_Integer /; k > 0 && k < 4] := {vectorType, v PauliMatrix[k] };
Bivector[v_, k_Integer /; k > 0 && k < 4, j_Integer  /; j > 0 && j < 4] := {bivectorType, v PauliMatrix[k] . PauliMatrix[j] };
Trivector[v_] := {trivectorType, v I IdentityMatrix[2]};

ClearAll[ DotProduct, WedgeProduct, GeometricProduct ]
DotProduct[ v1_?pScalarQ, v2_?pScalarQ]    := directProduct[ scalarType, v1 , v2  ] ;
DotProduct[ v1_?pScalarQ, v2_?pVectorQ]    := directProduct[ vectorType, v1 , v2  ] ;
DotProduct[ v1_?pScalarQ, v2_?pBivectorQ]  := directProduct[ bivectorType, v1 , v2  ] ;
DotProduct[ v1_?pScalarQ, v2_?pTrivectorQ] := directProduct[ trivectorType, v1 , v2  ] ;
DotProduct[ v1_?pVectorQ,    v2_?pScalarQ] := DotProduct[ v2, v1 ] ;
DotProduct[ v1_?pBivectorQ,  v2_?pScalarQ] := DotProduct[ v2, v1 ] ;
DotProduct[ v1_?pTrivectorQ, v2_?pScalarQ] := DotProduct[ v2, v1 ] ;

DotProduct[ v1_?pVectorQ, v2_?pVectorQ]    :=     symmetric[ scalarType,   v1, v2 ] ;
DotProduct[ v1_?pVectorQ, v2_?pBivectorQ]  := antisymmetric[ vectorType,   v1, v2 ] ;
DotProduct[ v1_?pVectorQ, v2_?pTrivectorQ] :=     symmetric[ bivectorType, v1, v2 ] ;
DotProduct[ v1_?pBivectorQ,  v2_?pVectorQ] := -DotProduct[ v2, v1 ] ;
DotProduct[ v1_?pTrivectorQ, v2_?pVectorQ] :=  DotProduct[ v2, v1 ] ;

DotProduct[ v1_?pBivectorQ, v2_?pBivectorQ]  := symmetric[ scalarType,     v1, v2 ] ;
DotProduct[ v1_?pBivectorQ, v2_?pTrivectorQ] := directProduct[ vectorType, v1, v2 ] ;
DotProduct[ v1_?pTrivectorQ,  v2_?pBivectorQ] := DotProduct[ v2, v1 ] ;

DotProduct[ v1_?pTrivectorQ,  v2_?pTrivectorQ] := directProduct[ scalarType, v1 , v2 ] ;

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






(*
t1 = Scalar[1] ;
t2 = Vector[1, 2];
t3 = Bivector[1, 2, 1];
t4 = Trivector[1];
pScalarQ[t1]
pScalarQ[t2]
pScalarQ[t3]
pScalarQ[t4]

pVectorQ[t1]
pVectorQ[t2]
pVectorQ[t3]
pVectorQ[t4]

pBivectorQ[t1]
pBivectorQ[t2]
pBivectorQ[t3]
pBivectorQ[t4]

pTrivectorQ[t1]
pTrivectorQ[t2]
pTrivectorQ[t3]
pTrivectorQ[t4]
*)

(*Bivector[1, 5]
Vector[1, 0]
Vector[1, 1]
Vector[1, 2]
Vector[1, 3]
Vector[1, 4]
Bivector[1, 5, 4]
Bivector[1, 3, 1]*)

End[]
