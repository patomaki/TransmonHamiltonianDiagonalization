(* ::Package:: *)

BeginPackage["QuantumMechanics`"];
IdM::usage="Identity matrix of size N; alias for built-in IdentityMatrix.\nInput: N (positive integer)";
MProduct::usage="Returns ArrayFlatten[TensorProduct[A,B]].";
eList::usage="List of (orthonormal) basis vectors of size Nquanta. The definition is based on that of Id.\nInput: Nquana (positive integer)";
e::usage="A basis vector of size Nquanta, with ith nonzero element. Returns a zero vector of size Nquanta for i > Nquanta.\nInput: Nquanta (positive integer); i (positive integer)";
RaiseM::usage="Raising operator, acting on states of size Nquanta, in matrix form.\nInput: Nquanta (positive integer)";
LowerM::usage="Lowering operator, acting on states of size Nquanta, in matrix form.\nInput: Nquanta (positive integer)";
KetBraM::usage="The matrix |i><j|, which couples the ket i and bra j, which are basis vectors of size Nquanta. Returns a zero matrix for i or j > Nquanta.\nInput: Nquanta; i; j";
NumberM::usage="Number operator Sum_m m |m><m| for number states m=0,1,...,Nquanta in matrix form. \nInput: Nquanta";
QubitStateList::usage="Returns a list of 1-qubit states for a list of indices.";
QubitStateVector::usage="Returns a vector representation of a basis vector labeled by 1-qubit indices (indexing starts from 1) indexList";
FindQubitStateIndices::usage="Returns a list of 1-qubit indices (indexing starts from 1) that correspond to a basis vector (in the chosen representation).";
FindVectorIndex::usage="";
MapBasisElement::usage="Re-represent as basis vector with permuted indices.";
MapMatrix::usage="Re-represent a matrix in a basis with permuted indices.";
Begin["Private`"];
IdM[N_]:=IdentityMatrix[N];
MProduct[A_List,B_List]:=ArrayFlatten[TensorProduct[A,B]];
eList[Nquanta_]:=IdM[Nquanta];
e[Nquanta_,i_]:=Module[{},
	If[i>Nquanta||i<1,
	Print["Invalid index for e.\n"];
	Table[0,{i,1,Nquanta}],
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
QubitStateList[indexList_]:=Module[{eList,eVector},
  eList={};
  Do[
    If[(indexList[[k]]!=1)&&(indexList[[k]]!=2),
      Print["Warning: incorrect index in argument number ",k,"indexList[[k]]:",indexList[[k]]]
    ];
    eList=Append[eList,If[indexList[[k]]==1,{1,0},If[indexList[[k]]==2,{0,1},{0,0}]]];
  ,{k,1,Length[indexList]}];
  eList
];
QubitStateVector[indexList_]:=Module[{eList,eVector},
  eList=QubitStateList[indexList];
  eVector=1;
  Do[
    eVector=TensorProduct[eVector,eList[[k]]]
  ,{k,1,Length[indexList]}];
  Flatten[eVector]
];
FindQubitStateIndices[e_List]:=Module[{tuples,indexList},
  tuples=Tuples[{1,2},Log[2,Length[e]]];
  indexList={};
  Do[
	If[e==QubitStateVector[tuples[[i]]],
      indexList=tuples[[i]];
      Break[]
    ];
  ,{i,1,Length[tuples]}];
  indexList
];
FindVectorIndex[indexList_]:=Module[{index,vector},
  index=0;
  vector=QubitStateVector[indexList];
  Do[
    If[vector[[k]]==1,
      index=k;
      Break[]
    ]
  ,{k,1,Length[vector]}];
  index
];
MapBasisElement[e_List,perm_]:=Module[{indices,vector},
(*perms is the vector with the places of the new indices, i.e. a
  permutation of {1,2,...,Nqubits}*)
  indices=FindQubitStateIndices[e];
  On[Assert];
  Assert[Length[perm]==Length[indices]];
  Off[Assert];
  vector=QubitStateVector[Permute[indices,perm]];
  vector
];
MapMatrix[M_List,perm_]:=Module[{newM},
  newM=Table[0,{i,1,Dimensions[M][[1]]},{j,1,Dimensions[M][[2]]}];
  On[Assert];
  Assert[Dimensions[M][[1]]==Dimensions[M][[2]]];
  Off[Assert];
  Do[
    Do[
      newM=newM+M[[i,j]]*ArrayFlatten[TensorProduct[MapBasisElement[e[Dimensions[M][[1]],i],perm],MapBasisElement[e[Dimensions[M][[1]],j],perm]]];
    ,{i,1,Dimensions[M][[1]]}];
  ,{j,1,Dimensions[M][[2]]}];
  newM
]
End[];
EndPackage[];
