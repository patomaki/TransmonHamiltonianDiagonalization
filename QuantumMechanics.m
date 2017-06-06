(* ::Package:: *)

BeginPackage["QuantumMechanics`"];
IdM::usage="Identity matrix of size N; alias for built-in IdentityMatrix.\nInput: N (positive integer)";
eList::usage="List of (orthonormal) basis vectors of size Nquanta. The definition is based on that of Id.\nInput: Nquana (positive integer)";
e::usage="A basis vector of size Nquanta, with ith nonzero element. Returns a zero vector of size Nquanta for i > Nquanta.\nInput: Nquanta (positive integer); i (positive integer)";
RaiseM::usage="Raising operator, acting on states of size Nquanta, in matrix form.\nInput: Nquanta (positive integer)";
LowerM::usage="Lowering operator, acting on states of size Nquanta, in matrix form.\nInput: Nquanta (positive integer)";
KetBraM::usage="The matrix |i><j|, which couples the ket i and bra j, which are basis vectors of size Nquanta. Returns a zero matrix for i or j > Nquanta.\nInput: Nquanta; i; j";
NumberM::usage="Number operator Sum_m m |m><m| for number states m=0,1,...,Nquanta in matrix form. \nInput: Nquanta";
Begin["Private`"];
IdM[N_]:=IdentityMatrix[N];
eList[Nquanta_]:=IdM[Nquanta];
e[Nquanta_,i_]:=Module[{},
	If[i>Nquanta||i<1,
	Print["Invalid index for e.\n"];
	Table[0,{j,1,Nquanta}],
	eList[Nquanta][[i]]
	]
];
RaiseM[Nquanta_]:=Table[If[i==(j+1),Sqrt[j],0],{i,1,Nquanta},{j,1,Nquanta}];
LowerM[Nquanta_]:=Table[If[j==(i+1),Sqrt[i],0],{i,1,Nquanta},{j,1,Nquanta}];
KetBraM[Nquanta_,i_,j_]:=Module[{},
	If[i>Nquanta||j>Nquanta||i<1||j<1,
	Print["Invalid index for KetBraM.\n"];
	ConstantArray[0,{Nquanta,Nquanta}],
	ArrayFlatten[TensorProduct[e[Nquanta,i],e[Nquanta,j]]]
	]
];
NumberM[Nquanta_]:=Sum[m*KetBraM[m,m],{m,1,Nquanta}];
End[];
EndPackage[];
