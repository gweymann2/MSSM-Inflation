(* ::Package:: *)

(* ::Package:: *)
(**)


" LLe sans les conventions de suspect3 "

g[i_][mu_]:= gmu0[i]/Sqrt[1 - b[i] gmu0[i]^2 Log[mu/mu0]];

M[i_][mu_]:= Mmu0[i] (g[i][mu]/gmu0[i])^2;

mphi2[mu_]:= mphi2mu0 + (M[2][mu0]^2 - M[2][mu]^2)  + 1/11 (M[1][mu0]^2 - M[1][mu]^2);
\:200b
A6[mu_]:= A6mu0 + 6(M[2][mu0] - M[2][mu]) + 6/11 (M[1][mu0] - M[1][mu]);
\:200b
lambda6[mu_]:= lambda6mu0 (g[2][mu0]/g[2][mu])^(6) (g[1][mu0]/g[1][mu])^(6/11); 
\:200b
lzero={mu D[g[1][mu], mu] - 11/16/Pi^2 g[1][mu]^3,
       mu D[g[2][mu], mu] - 1/16/Pi^2 g[2][mu]^3,
       D[M[1][mu]/g[1][mu]^2, mu],
       D[M[2][mu]/g[2][mu]^2, mu],
       mu D[mphi2[mu], mu] + 1/6/Pi^2 (3/2 M[2][mu]^2 g[2][mu]^2 +
                                       3/2 M[1][mu]^2 g[1][mu]^2),
       mu D[A6[mu], mu] + 1/2/Pi^2 (3/2 M[2][mu] g[2][mu]^2 +
                                    3/2 M[1][mu] g[1][mu]^2),
       mu D[lambda6[mu], mu] + 1/4/Pi^2 lambda6[mu] (3/2 g[2][mu]^2 +
                                       3/2   g[1][mu]^2)} //. {b[1] -> 11/8/Pi^2, b[2] -> 1/8/Pi^2};

Simplify[lzero]


(* ::Package:: *)
(**)


"lle avec les bonnes conventions:"

g[i_][mu_]:= gmu0[i]/Sqrt[1 - b[i] gmu0[i]^2 Log[mu/mu0]];
\:200b
M[i_][mu_]:= Mmu0[i] (g[i][mu]/gmu0[i])^2;
\:200b
mphi2[mu_]:= mphi2mu0 + (M[2][mu0]^2 - M[2][mu]^2)  + 1/11 (M[1][mu0]^2 - M[1][mu]^2);
\:200b
A6[mu_]:= A6mu0 - 6(M[2][mu0] - M[2][mu]) - 6/11 (M[1][mu0] - M[1][mu]);
\:200b
lambda6[mu_]:= lambda6mu0 (g[2][mu0]/g[2][mu])^(6) (g[1][mu0]/g[1][mu])^(6/11); 

lzero={mu D[g[1][mu], mu] - 33/5/16/Pi^2 g[1][mu]^3,
       mu D[g[2][mu], mu] - 1/16/Pi^2 g[2][mu]^3,
       D[M[1][mu]/g[1][mu]^2, mu],
       D[M[2][mu]/g[2][mu]^2, mu],
       mu D[mphi2[mu], mu] + 1/6/Pi^2 (3/2 M[2][mu]^2 g[2][mu]^2 +
                                       9/10 M[1][mu]^2 g[1][mu]^2),
       mu D[A6[mu], mu] - 1/2/Pi^2 (3/2 M[2][mu] g[2][mu]^2 +
                                    9/10 M[1][mu] g[1][mu]^2),
       mu D[lambda6[mu], mu] + 1/4/Pi^2 lambda6[mu] (3/2 g[2][mu]^2 +
                                       9/10   g[1][mu]^2)} //. {b[1] -> 33/40/Pi^2, b[2] -> 1/8/Pi^2};
\:200b
Simplify[lzero]



\:200b "udd avec les bonnes conventions"
\:200b\:200b
g[i_][mu_]:= gmu0[i]/Sqrt[1 - b[i] gmu0[i]^2 Log[mu/mu0]];
\:200b
M[i_][mu_]:= Mmu0[i] (g[i][mu]/gmu0[i])^2;
\:200b
mphi2[mu_]:= mphi2mu0 + (-8)/9*(M[3][mu0]^2 - M[3][mu]^2)  + 4/99 (M[1][mu0]^2 - M[1][mu]^2);
\:200b
A6[mu_]:= A6mu0 - 6*(-8)/9*(M[3][mu0] - M[3][mu]) - 6*4/99 (M[1][mu0] - M[1][mu]);
\:200b
lambda6[mu_]:= lambda6mu0 (g[3][mu0]/g[3][mu])^(6*(-8)/9) (g[1][mu0]/g[1][mu])^(6*4/99); 

lzero={mu D[g[1][mu], mu] - 33/5/16/Pi^2 g[1][mu]^3,
       mu D[g[3][mu], mu] +3/16/Pi^2 g[3][mu]^3,
       D[M[1][mu]/g[1][mu]^2, mu],
       D[M[3][mu]/g[3][mu]^2, mu],
       mu D[mphi2[mu], mu] + 1/6/Pi^2 (4 M[3][mu]^2 g[3][mu]^2 +
                                       2/5 M[1][mu]^2 g[1][mu]^2),
       mu D[A6[mu], mu] - 1/2/Pi^2 (4 M[3][mu] g[3][mu]^2 +
                                    2/5 M[1][mu] g[1][mu]^2),
       mu D[lambda6[mu], mu] + 1/4/Pi^2 lambda6[mu] (4 g[3][mu]^2 +
                                       2/5   g[1][mu]^2)} //. {b[1] -> 33/40/Pi^2, b[3] -> -3/8/Pi^2};
\:200b
Simplify[lzero]




