{ Assert[bladeQ[#]] } &/@ { s0, e1, e2, e3, b23, b31, b12, t123 } ;
{ Assert[!bladeQ[#]] } &/@ { m01, m02, m03, m12, m13, m23, m012, m013, m023, m123 } ;
{ Assert[gradeAnyQ[#]] } &/@ { s0, e1, e2, e3, b23, b31, b12, t123, m01, m02, m03, m12, m13, m23, m012, m013, m023, m123 } ;
{ Assert[!gradeAnyQ[#]] } &/@ { 1, Sin[x], Exp[ I theta] } ;
{ Assert[!notGradeQ[#]] } &/@ { s0, e1, e2, e3, b23, b31, b12, t123, m01, m02, m03, m12, m13, m23, m012, m013, m023, m123 } ;
{ Assert[notGradeQ[#]] } &/@ { 1, Sin[x], Exp[ I theta] } ;
{ Assert[gradeQ[#, 0]], Assert[scalarQ[#]] } &/@ { s0 } ;
{ Assert[!gradeQ[#, 0]], Assert[!scalarQ[#]] } &/@ { e1, e2, e3, b23, b31, b12, t123, m01, m02, m03, m12, m13, m23, m012, m013, m023, m123 } ;
{ Assert[gradeQ[#, 1]], Assert[vectorQ[#]] } &/@ { e1, e2, e3 } ;
{ Assert[!gradeQ[#, 1]], Assert[!vectorQ[#]] } &/@ { s0, b23, b31, b12, t123, m01, m02, m03, m12, m13, m23, m012, m013, m023, m123 } ;
{ Assert[gradeQ[#, 2]], Assert[bivectorQ[#]] } &/@ { b23, b31, b12 } ;
{ Assert[!gradeQ[#, 2]], Assert[!bivectorQ[#]] } &/@ { s0, e1, e2, e3, t123, m01, m02, m03, m12, m13, m23, m012, m013, m023, m123 } ;
{ Assert[gradeQ[#, 3]], Assert[trivectorQ[#]] } &/@ { t123 } ;
{ Assert[!gradeQ[#, 3]], Assert[!trivectorQ[#]] } &/@ { s0, e1, e2, e3, b23, b31, b12, m01, m02, m03, m12, m13, m23, m012, m013, m023, m123 } ;
{ Assert[gradeQ[#, -1]] } &/@ { m01, m02, m03, m12, m13, m23, m012, m013, m023, m123 } ;
{ Assert[!gradeQ[#, -1]] } &/@ { s0, e1, e2, e3, b23, b31, b12, t123 } ;

