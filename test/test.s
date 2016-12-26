{ Assert[GradeSelection[#, 0] == s0], Assert[ScalarSelection[#] == s0] } & /@ { s0, m01, m02, m03, m013, m012, m023 } ;
{ Assert[GradeSelection[#, 0] == 0], Assert[ScalarSelection[#] == 0] } & /@ { e1, e2, e3, b23, b31, b12, t123, m12, m13, m23, m123 } ;

{ Assert[GradeSelection[#, 1] == e1], Assert[VectorSelection[#] == e1] } & /@ { e1, m01, m12, m13, m012, m013, m123 } ;
{ Assert[GradeSelection[#, 1] == 0], Assert[VectorSelection[#] == 0] } & /@ { s0, b23, b31, b12, t123, m02, m03, m23, m023 };

{ Assert[GradeSelection[#, 2] == b23], Assert[BivectorSelection[#] == b23] } & /@ { b23, m02, m12, m23, m023, m123, m012 };
{ Assert[GradeSelection[#, 2] == 0], Assert[BivectorSelection[#] == 0] } & /@ { s0, e1, e2, e3, t123, m01, m03, m13, m013 };

{ Assert[GradeSelection[#, 3] == t123], Assert[TrivectorSelection[#] == t123] } & /@ { t123, m03, m13, m23, m013, m023, m123 } ;
{ Assert[GradeSelection[#, 3] == 0], Assert[TrivectorSelection[#] == 0] } & /@ { s0, e1, e2, e3, b23, b31, b12, m01, m02, m12, m012 } ;


