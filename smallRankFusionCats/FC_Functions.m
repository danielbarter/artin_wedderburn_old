(* ::Package:: *)

polarForm=Expand[#/.z_?NumericQ:>Abs[z] Exp[I Arg[z]]]&;


fusionproduct[a_,b_]/;MemberQ[obs,a]&&MemberQ[obs,b]:=fusionproduct[a,b]=DeleteCases[Table[If[V[a,b][c]=!=0,c],{c,obs}],Null]

dim[a_]/;MemberQ[obs,a]:=dim[a]=Module[{d0,dimeqs,dimsol,da},
dimeqs=Join[
Outer[d0[#1]d0[#2]==Sum[V[#1,#2][x]d0[x],{x,fusionproduct[#1,#2]}]&,obs,obs]//Flatten
,d0[#]>=1&/@obs];
da:=d0[a]//.Solve[dimeqs][[1]]//RootReduce//FullSimplify//ToRadicals;
Return[da]];

dtot[]:=dtot[]=Sqrt[Total[dim[#]^2&/@obs]]//RootReduce//FullSimplify//ToRadicals
dtot[X_]:=dtot[X]=Sqrt[Total[dim[#]^2&/@X]]//RootReduce//FullSimplify//ToRadicals

dual[a_]/;MemberQ[obs,a]:=dual[a]=Module[{v=Abs[V[a,#][unit]]&/@obs},obs[[Position[v,1][[1,1]]]]]
OverBar[a_]/;MemberQ[obs,a]:=dual[a]

dim[{A__}]:=Sum[dim[a],{a,{A}}]


F[a_,b_,c_][d_][e_,f_]:=0

F[a_,b_,u_][c_][c_,b_]/;(u===unit&&MemberQ[fusionproduct[a,b],c]):=1;
F[a_,u_,b_][c_][a_,b_]/;(u===unit&&MemberQ[fusionproduct[a,b],c]):=1;
F[u_,a_,b_][c_][a_,c_]/;(u===unit&&MemberQ[fusionproduct[a,b],c]):=1;

F[a_,b_,c_][d_][e_,f_]/;MemberQ[fusionproduct[a,b],e]&&MemberQ[fusionproduct[b,c],f]&&MemberQ[Intersection[fusionproduct[e,c],fusionproduct[a,f]],d]:=F[a,b,c][d][e,f]=
X0[a,b,c][d][e,f]//.rep0


pentagons[a_,b_,c_,d_]:=pentagons[a,b,c,d]=Table[
F[f,c,d][e][g,x]F[a,b,x][e][f,y]==
Sum[F[a,b,c][g][f,z]F[a,z,d][e][g,y]F[b,c,d][y][z,x],{z,fusionproduct[b,c]}]
,{f,fusionproduct[a,b]},{g,fusionproduct[f,c]},{x,fusionproduct[c,d]},{y,fusionproduct[b,x]},{e,Intersection[fusionproduct[g,d],fusionproduct[a,y]]}
]//Flatten//DeleteCases[True]

allpentagons:=allpentagons=Table[pentagons[a,b,c,d],{a,obs},{b,obs},{c,obs},{d,obs}]//Flatten//DeleteCases[True]//DeleteDuplicates

nicepentagons:=nicepentagons=Module[
{a,b,c,d,pents={}},
Do[If[Length[fusionproduct[b,c]]==1,
AppendTo[pents,Table[pentagons[a,b,c,d],{a,obs},{d,obs}]//Flatten]
],{b,obs},{c,obs}];
Return[pents//Flatten//DeleteDuplicates];
]

unitaryFL:=unitaryFL=Table[
Sum[F[a,b,c][d][e,f]Conjugate[F[a,b,c][d][g,f]],{f,fusionproduct[b,c]}]==If[e===g,1,0]
,{a,obs},{b,obs},{c,obs},{e,fusionproduct[a,b]},{g,fusionproduct[a,b]},{d,fusionproduct[e,c]}]//Flatten//DeleteCases[True]//DeleteDuplicates

unitaryFR:=unitaryFR=Table[
Sum[F[a,b,c][d][f,e]Conjugate[F[a,b,c][d][f,g]],{f,fusionproduct[a,b]}]==If[e===g,1,0]
,{a,obs},{b,obs},{c,obs},{e,fusionproduct[b,c]},{g,fusionproduct[b,c]},{d,fusionproduct[a,e]}]//Flatten//DeleteCases[True]//DeleteDuplicates

unitaryF:=unitaryF=Join[unitaryFL,unitaryFR]//DeleteDuplicates;
Fs:=Fs=Table[F[a,b,c][d][e,f],{a,obs},{b,obs},{c,obs},{e,fusionproduct[a,b]},{f,fusionproduct[b,c]},{d,Intersection[fusionproduct[a,f],fusionproduct[e,c]]}]//Flatten

newF[a_,b_,c_][d_][e_,f_]/;MemberQ[fusionproduct[a,b],e]&&MemberQ[fusionproduct[b,c],f]&&MemberQ[Intersection[fusionproduct[e,c],fusionproduct[a,f]],d]:=newF[a,b,c][d][e,f]=
F[a,b,c][d][e,f]U[b,c][f]U[a,f][d]Conjugate[U[a,b][e]U[e,c][d]]

\[Kappa][x_]:=F[x,dual[x],x][x][unit,unit]dim[x]//FullSimplify

new\[Kappa][x_]:=newF[x,dual[x],x][x][unit,unit]dim[x]//FullSimplify

cleanupFs[ns_]:=Module[{relabF0,relabR0,P,T,M=$Assumptions,Q},
relabF0=Table[X0[a,b,c][d][e,f]->FullSimplify[F[a,b,c][d][e,f]//.ns],{a,obs},{b,obs},{c,obs},{e,fusionproduct[a,b]},{f,fusionproduct[b,c]},{d,Intersection[fusionproduct[e,c],fusionproduct[a,f]]}]//Flatten;
P=DeleteCases[relabF0,_?(#[[1]]===#[[2]]&)];
Q=Join[P];
Return[Q];
]

reguageF[ns_]:=Module[{relabF0,relabR0,P,T,M=$Assumptions,Q},
relabF0=Table[X0[a,b,c][d][e,f]->Simplify[newF[a,b,c][d][e,f]//.ns//RootReduce]//Refine,{a,obs},{b,obs},{c,obs},{e,fusionproduct[a,b]},{f,fusionproduct[b,c]},{d,Intersection[fusionproduct[e,c],fusionproduct[a,f]]}]//Flatten;
P=DeleteCases[relabF0,_?(#[[1]]===#[[2]]&)];

Q=Join[P];
Return[Q];
]


A[a_,b_][c_]:=Conjugate[F[dual[a],a,b][b][unit,c]]Sqrt[(dim[a]dim[b])/dim[c]]
A[c_][a_,b_]:=F[dual[a],a,b][b][unit,c]Sqrt[(dim[a]dim[b])/dim[c]]

B[a_,b_][c_]:=F[a,b,dual[b]][a][c,unit]Sqrt[(dim[a]dim[b])/dim[c]]
B[c_][a_,b_]:=Conjugate[F[a,b,dual[b]][a][c,unit]]Sqrt[(dim[a]dim[b])/dim[c]]

symmetry1:=symmetry1=
Table[
F[a,b,c][d][e,f]==1/Abs[\[Kappa][a]]^2 Sqrt[(dim[e]dim[f])/(dim[b]dim[d])]A[d][a,f]A[a,b][e]Conjugate[F[dual[a],e,c][f][b,d]],
{a,obs},{b,obs},{c,obs},{e,fusionproduct[a,b]},{f,fusionproduct[b,c]},{d,Intersection[fusionproduct[a,f],fusionproduct[e,c]]}]//Flatten(*//RootReduce*)//DeleteCases[True]//DeleteDuplicates;

symmetry2:=symmetry2=
Table[
F[a,b,c][d][e,f]==1/Conjugate[\[Kappa][b]] Sqrt[(dim[e]dim[f])/(dim[a]dim[c])]A[f][b,c]B[a,b][e]Conjugate[F[e,dual[b],f][d][a,c]],
{a,obs},{b,obs},{c,obs},{e,fusionproduct[a,b]},{f,fusionproduct[b,c]},{d,Intersection[fusionproduct[a,f],fusionproduct[e,c]]}]//Flatten(*//RootReduce*)//DeleteCases[True]//DeleteDuplicates;

symmetry3:=symmetry3=
Table[
F[a,b,c][d][e,f]==1/Abs[\[Kappa][c]]^2 Sqrt[(dim[e]dim[f])/(dim[b]dim[d])]B[f][b,c]B[e,c][d]Conjugate[F[a,f,dual[c]][e][d,b]],
{a,obs},{b,obs},{c,obs},{e,fusionproduct[a,b]},{f,fusionproduct[b,c]},{d,Intersection[fusionproduct[a,f],fusionproduct[e,c]]}]//Flatten(*//RootReduce*)//DeleteCases[True]//DeleteDuplicates;

pivotalstar:=Table[
F[a,b,c][d][e,f]==1/(Conjugate[\[Kappa][a]]\[Kappa][c]) A[d][a,f]A[dual[a],e][b]B[f][b,c]B[d,dual[c]][e]F[dual[a],d,dual[c]][b][f,e],{a,obs},{b,obs},{c,obs},{e,fusionproduct[a,b]},{f,fusionproduct[b,c]},{d,Intersection[fusionproduct[a,f],fusionproduct[e,c]]}]//Flatten//DeleteDuplicates//DeleteCases[True];


F[a_,b_][c_,d_][e_,f_]:=1/\[Kappa][d] B[e][c,d]F[a,b,dual[d]][c][e,f]B[b,dual[d]][f]


U[a_,u_][b_]/;(u===unit):=1;
U[u_,a_][b_]/;(u===unit):=1;


newFs:=newFs=Table[newF[a,b,c][d][e,f],{a,obs},{b,obs},{c,obs},{e,fusionproduct[a,b]},{f,fusionproduct[b,c]},{d,Intersection[fusionproduct[a,f],fusionproduct[e,c]]}]//Flatten

unitaryU:=unitaryU=Table[
Abs[U[a,b][c]]==1
,{a,obs},{b,obs},{c,fusionproduct[a,b]}
]//Flatten//DeleteCases[True]//DeleteDuplicates


Unprotect[Abs,Conjugate,Power,Times,Sqrt,Im,Re];
Abs[\[Phi][a_,b_,c_][d_][e_,f_]]:=1;
Conjugate[\[Phi][a_,b_,c_][d_][e_,f_]]:=1/\[Phi][a,b,c][d][e,f];

Abs[\[Zeta][a_,b_,c_][d_][e_,f_]]:=1;
Conjugate[\[Zeta][a_,b_,c_][d_][e_,f_]]:=\[Zeta][a,b,c][d][e,f];
\[Zeta][a_,b_,c_][d_][e_,f_]^m_Integer/;(m<0||m>1):=\[Zeta][a,b,c][d][e,f]^Mod[m,2]

Abs[\[CurlyPhi][x_Integer][a_,b_,c_][d_][e_,f_]]:=1
Conjugate[\[CurlyPhi][x_Integer][a_,b_,c_][d_][e_,f_]]:=\[CurlyPhi][x][a,b,c][d][e,f]^-1
\[CurlyPhi][x_Integer][a_,b_,c_][d_][e_,f_]^m_Integer/;(m<0||m>=x):=\[CurlyPhi][x][a,b,c][d][e,f]^Mod[m,x]

Conjugate[Re0[a_,b_,c_][d_][e_,f_]]:=Re0[a,b,c][d][e,f];
Re[Re0[a_,b_,c_][d_][e_,f_]]:=Re0[a,b,c][d][e,f];
Im[Re0[a_,b_,c_][d_][e_,f_]]:=0;
Abs[Re0[a_,b_,c_][d_][e_,f_]]:=Sign[Re0[a,b,c][d][e,f]];

Conjugate[P0[a_,b_,c_][d_][e_,f_]]:=P0[a,b,c][d][e,f];
Re[P0[a_,b_,c_][d_][e_,f_]]:=P0[a,b,c][d][e,f];
Im[P0[a_,b_,c_][d_][e_,f_]]:=0;
Abs[P0[a_,b_,c_][d_][e_,f_]]:=P0[a,b,c][d][e,f];

Conjugate[U[a_,b_][c_]]/;V[a,b][c]==1:=1/U[a,b][c]
Protect[Abs,Conjugate,Power,Times,Sqrt];
