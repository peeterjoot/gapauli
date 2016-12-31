(* ::Package:: *)

(* copy this module to a directory in $Path.  Then invoke with <<GA30` *)
BeginPackage[ "complex`" ]
complex::usage =
  "complex.  A limited use complex number implementation to use internally in \
a Pauli or Dirac matrix basis representation, independent of any Complex" ;
ClearAll[ complex, complexQ, notComplexQ, real, imag, conjugate ] ;

complex /: complex[ r1_, i1_ ] + complex[ r2_, i2_ ] := complex[ r1 + r2, i1 + i2 ] ;
complex /: r1_ + complex[ r2_, i2_ ] := complex[ r1 + r2, i2 ] ;

complex /: -complex[ re_, im_ ] := complex[ -re, -im ] ;

complex /: complex[ re_ ] := re ;
complex /: complex[ re_, 0 ] := re ;

complex /: complex[ r1_, i1_ ] complex[ r2_, i2_ ] := complex[ r1 r2 - i1 i2, r1 i2 + r2 i1 ] ;

norm::usage = "norm[ z ].  A Norm like function for complex[ ]" ;
norm[ z_complex ] := ((z // First)^2 + (z // Last)^2) ;

(*special case this one to deal with the sort of products that are \
generated multiplying pauli matrices*)

complex /: Power[ z_complex, 2 ] := complex[ z ] complex[ z ] ;
complex /: Power[ z_complex, n_ ] :=
  Module[ {r = norm[ z ]^(n/2), theta = n ArcTan[ z // First, z // Last ]},
    r complex[ Cos[ theta ], Sin[ theta ] ] ] ;

complexQ::usage = "complexQ[ z ].  predicate pattern match for complex[ ]" ;
complexQ[ z_complex ] := True ;
complexQ[ _ ] := False ;
notComplexQ::usage = "notComplexQ[ z ].  predicate pattern match for !complex[ ]" ;
notComplexQ[ v_ ] := Not[ complexQ[ v ] ] ;

complex /: (v_?notComplexQ) complex[ re_, im_ ] := complex[ v re, v im ] ;

real::usage = "real[ z ].  Re[ z ] like function for complex[ ]" ;
real[ z_complex ] := (z // First) ;
imag::usage = "imag[ z ].  Im[ z ] like function for complex[ ]" ;
imag[ z_complex ] := (z // Last) ;
real[ ex_ ] := ex ;
imag[ ex_ ] := 0 ;
conjugate::usage = "conjugate[ z ].  Conjugate[ z ] like function for complex[ ]" ;
conjugate[ z_complex ] := complex[ z // First, -z // Last ] ;
conjugate[ ex_ ] := ex ;

ClearAll[ complexI, fMatrix, matrixreal, matriximag, matrixconj ]
complexI::usage = "complexI.  I like unit imaginary for complex[ ]" ;
complexI := complex[ 0, 1 ] ;

fMatrix[ p_, f_ ] := (Function[ a, f@a, Listable ]@p)
matrixreal::usage =
  "matrixreal.  method to apply real to all elements in matrix.  This is a hack.  \
Can probably set an attribute on the real function to do this." ;
matrixreal[ m_ ] := fMatrix[ m, real ] ;
matriximag::usage =
  "matriximag.  method to apply imag to all elements in matrix.  This is a hack.  \
Can probably set an attribute on the imag function to do this." ;
matriximag[ m_ ] := fMatrix[ m, imag ] ;
matrixconj::usage =
  "matrixconj.  method to apply conjugate to all elements in matrix.  This is a \
hack.  Can probably set an attribute on the conj function to do this." ;
matrixconj[ m_ ] := fMatrix[ m, conjugate ] ;

Protect[ complex, complexQ, notComplexQ, real, imag, conjugate, complexI, fMatrix, matrixreal, matriximag, matrixconj ] ;

EndPackage[ ]
