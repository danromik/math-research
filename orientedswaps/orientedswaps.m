(* OrientedSwaps: A companion Mathematica package to the paper 										*)
(* "Sorting networks, staircase shape Young tableaux and last passage percolation"					*)
(* by Elia Bisi, Fabio Deelan Cunden, Shane Gibbons and Dan Romik									*)
(* Version 1.0 (October 17, 2019)																	*)


(**************************************************************************************************
Instructions:
The goal of this package is to verify the conjectural identity between the two generating functions 
F_n(x_1,...,x_{n-1}) and G_n(x_1,...,x_{n-1}) defined in the paper, for small values of n.

To verify the identity, run the following command: 

genfuns = VerifyStaircaseSYTSortingNetworksIdentity[n]; 

where n is 3, 4, 5 or 6. The calculations for n=3,4,5 are almost instantaneous.
The calculation for n=6 will take around 10-20 minutes on a modern machine.

The output variable "genfuns" will contain a list of the pairs of components 
of the two generating functions, indexed by permutations of order n-1. 
Feel free to inspect it and try to understand why the identity is true.
**************************************************************************************************)



(* Typing OrientedSwapsHelp[] on the Mathematica command line will print a help message *)
OrientedSwapsHelp[] := Print[
ToString[Style["\nOrientedSwaps (Version 1.0, October 17, 2019) \n\n",Bold],StandardForm]<>
"OrientedSwaps is a companion Mathematica package to the 2019 paper \n"<>"\"Sorting networks, staircase shape Young tableaux and last passage percolation\",\n"<>
"by Elia Bisi, Fabio Deelan Cunden, Shane Gibbons and Dan Romik \n\n"<>
"The goal of this package is to verify the conjectural identity between the two generating functions \n" <>
ToString[Subscript[Style["F",Italic],Style["n",Italic]],StandardForm]<>"("<>
ToString[Subscript[Style["x",Italic],"1"],StandardForm]<>", . . . , "<>
ToString[Subscript[Style["x",Italic],"n-1"],StandardForm]<>")"
<>" and "
ToString[Subscript[Style["G",Italic],Style["n",Italic]],StandardForm]<>"("<>
ToString[Subscript[Style["x",Italic],"1"],StandardForm]<>", . . . , "<>
ToString[Subscript[Style["x",Italic],"n-1"],StandardForm]<>")"
<> " defined in the paper, for small values of "<>
ToString[Style["n",Italic],StandardForm]<>".\n\n"<>
"To verify the identity, run the following command: \n\n"<>
ToString[Style["genfuns = VerifyStaircaseSYTSortingNetworksIdentity[n]; ",Blue],StandardForm]<>"\n\n"
"where "<>ToString[Style["n",Italic],StandardForm]<>" is 3, 4, 5 or 6. "<>
"The calculations for "<>ToString[Style["n",Italic],StandardForm]<>"=3,4,5 are almost instantaneous.\n"<>
"The calculation for "<>ToString[Style["n",Italic],StandardForm]<>"=6 will take around 10-20 minutes on a modern machine.\n\n"<>
"The output variable "<>ToString[Style["genfuns", Blue],StandardForm]<>" will contain a list of the pairs of components \n"<>
"of the two generating functions, indexed by permutations of order "<>ToString[Style["n",Italic],StandardForm]<>"-1. \n"<>
"Feel free to inspect it and try to understand why the identity is true."
];




(* The section below contains code to generate and manipulate sorting networks *)


(* ForwardSwaps[perm_] accepts a permutation and returns the positions of adjacent swaps *)
(* that increase the number of inversions *)
ForwardSwaps[perm_]:=Select[Table[If[perm[[k+1]]>perm[[k]],k,0],{k,1,Length[perm]-1}],#!=0&];

(* BackwardSwaps[perm_] accepts a permutation and returns the positions of adjacent swaps *)
(* that decrease the number of inversions *)
BackwardSwaps[perm_]:=Select[Table[If[perm[[k+1]]<perm[[k]],k,0],{k,1,Length[perm]-1}],#!=0&];

(* ReducedPaths[perm_] accepts a permutation and returns a list of all minimal-length paths *)
(* connecting the identity permutation to perm with adjacent transpositions as steps *)
ReducedPaths[perm_]:=Module[{paths,backswaps,prevperm},
	If[perm==Range[Length[perm]],{{perm}},
		backswaps=BackwardSwaps[perm];
		paths={};
		Do[ prevperm=perm;
			prevperm[[backswaps[[k]]]]=perm[[backswaps[[k]]+1]];
			prevperm[[backswaps[[k]]+1]]=perm[[backswaps[[k]]]];
			paths=Join[paths,Map[Append[#,perm]&,ReducedPaths[prevperm]]];
					,{k,1,Length[backswaps]}];
		paths
	]
];

(* SortingNetworks[n_] returns a list of all sorting networks of order n, that is, minimal-length *)
(* paths connecting the identity permutation to the reverse permutation *)
SortingNetworks[n_]:=ReducedPaths[Range[n, 1, -1]];

(* ReducedPathOutdegrees[path_] takes a path of permutations and returns the list of outdegrees *)
(* (that is, numbers of forward swaps) at each step of the path except the last one *)
ReducedPathOutdegrees[path_]:=Map[Length[ForwardSwaps[#]]&,Most[path]];


(* LastSwapTime[path_,k_] returns for a sorting network path and an integer k the k-th last swap time *)
(* associated with the sorting network *)
LastSwapTime[path_,k_]:=-1+Min[Flatten[Position[Map[Map[Function[x,If[x<=k,1,0]],#]&,path],Table[If[j<=Length[path[[1]]]-k,0,1],{j,1,Length[path[[1]]]}]]]];

(* VectorOfLastSwapTimes[path_] returns a list of all the last swap times associated with the *)
(* sorting network *)
VectorOfLastSwapTimes[path_]:=Table[LastSwapTime[path,k],{k,1,Length[path[[1]]]-1}];

(* LastSwapTimesPermutation[path_] returns the permutation encoding the ordering of the last swap times *)
LastSwapTimesPermutation[path_]:=InversePermutation[Ordering[VectorOfLastSwapTimes[path]]];






(* The next section contains code to generate and manipulate standard Young tableaux *)

Needs["Combinatorica`"];  (* WARNING: DEPRECATED PACKAGE - improve code so it doesn't depend on it *)


(* StaircaseShapeYoungTableaux[n_] returns a list of all staircase shape standard Young *)
(* tableaux of order n [defined as the Young diagram (n-1,n-2,...,1)] *)
(* TODO: write an implementation that doesn't depend on a deprecated package *)
StaircaseShapeYoungTableaux[n_]:=Tableaux[Range[n-1,1,-1]];

NumberOfStaircaseShapeSYT[n_] := (n(n-1)/2)!/Product[(2k-1)^(n-k),{k,1,n-1}];


(* SYTToPath[tableau_] takes a Young tableau and converts it to a path on the Young lattice *)
(* connecting the empty diagram with the diagram of shape shape(tableau) *)
SYTToPath[tableau_]:=Module[{diagrams,cellsequence,n},
	n=Apply[Plus,Map[Length,tableau]];
	diagrams={{}};
	cellsequence=Table[Position[tableau,k][[1]],{k,1,n}];
	Do[
		AppendTo[diagrams,diagrams[[k]]];
		If[cellsequence[[k]][[1]]>Length[diagrams[[k+1]]],
			AppendTo[diagrams[[k+1]],1],
			diagrams[[k+1]][[cellsequence[[k]][[1]]]]++
		]
		,{k,1,n}
	];
	diagrams
];

(* StaircaseShapeSYTCorners[tableau] returns the list of corner entries *)
(* of the tableau  *)
StaircaseShapeSYTCorners[tableau_]:=Reverse[Map[Last,tableau]];

(* StaircaseShapeSYTCornersPermutation[tableau] returns the permutation encoding the ordering of the corner entries *)
StaircaseShapeSYTCornersPermutation[tableau_]:=InversePermutation[Ordering[StaircaseShapeSYTCorners[tableau]]];


(* StaircaseSubdiagramOutdegree[diagram_,n_] returns the number of boxes that may be added *)
(* to diagram to create a larger diagram that is still a subdiagram of the staircase shape of *)
(* order n *)
StaircaseSubdiagramOutdegree[diagram_,n_]:=If[diagram=={},1,
If[diagram[[1]]<n-1,1,0]+
Sum[If[diagram[[j]]<diagram[[j-1]]&&diagram[[j]]<n-j,1,0],{j,2,Length[diagram]}]+
If[Length[diagram]<n-1,1,0]];

(* StaircaseShapeSYTOutDegrees[tableau_] returns the sequence of outdegrees associated with the subdiagrams *)
(* in the Young diagram growth path encoded by a staircase shape SYT tableau *)
StaircaseShapeSYTOutDegrees[tableau_]:= Module[{path,n},
    n = Length[tableau[[1]]]+1;
    path = SYTToPath[tableau];
    Map[StaircaseSubdiagramOutdegree[#,n]&,Most[path]]
];






(* This final section has code to compute the generating functions F_n and G_n and verify the identity between them *)

(* A helper function used by the functions GeneratingFactorForStaircaseShapeSYT[...] and *)
(* GeneratingFactorForSortingNetwork[...]  defined below *)
GeneratingFactorFromParameters[cornerslastswapvector_, degreeseq_, varname_]:=
Module[{sortedvector},
	sortedvector = Sort[cornerslastswapvector];
	1/Product[
		Product[Subscript[varname, k] + degreeseq[[j]], {j, If[k==1,1,sortedvector[[k-1]]+1], sortedvector[[k]]}]
			,{k,1,Length[sortedvector]}]
];

(* GeneratingFactorForStaircaseShapeSYT[tableau, varname] returns the generating factor associated with a *)
(* staircase shape Young tableau, in the variable varname *)
GeneratingFactorForStaircaseShapeSYT[tableau_, varname_]:=
	GeneratingFactorFromParameters[StaircaseShapeSYTCorners[tableau],StaircaseShapeSYTOutDegrees[tableau],varname];

(* GeneratingFactorForSortingNetwork[sortingnet, varname] returns the generating factor associated with a *)
(* sorting network, in the variable varname *)
GeneratingFactorForSortingNetwork[sortingnet_, varname_]:=
	GeneratingFactorFromParameters[VectorOfLastSwapTimes[sortingnet],ReducedPathOutdegrees[sortingnet],varname];
	
(* GeneratingFunctionForAllStaircaseShapeSYT[n, varname] returns the vector-valued generating function *)
(* (denoted F_n(...) in the paper), in the variable varname,  defined in terms of a summation *)
(* over staircase shape Young tableaux of order n *)
(* The list is of length (n-1)!, with element corresponding to the component of the vector associated *)
(* with a particular permutation of order n-1 *)
GeneratingFunctionForAllStaircaseShapeSYT[n_, varname_]:=Module[
		{tableaux, cornerpermutations, generatingfactors, allpermutations},
	allpermutations = Permutations[Range[n-1]];
	tableaux = StaircaseShapeYoungTableaux[n];
	cornerpermutations = Map[StaircaseShapeSYTCornersPermutation, tableaux];
	generatingfactors = Map[GeneratingFactorForStaircaseShapeSYT[#, varname]&, tableaux];
	Table[Total[Map[Part[generatingfactors,#[[1]]]&,Position[cornerpermutations,allpermutations[[k]]]]] ,{k,1,Length[allpermutations]}]
];

(* GeneratingFunctionForAllStaircaseShapeSYT[n, varname] returns the vector-valued generating function *)
(* (denoted G_n(...) in the paper), in the variable varname,  defined in terms of a summation *)
(* over staircase shape Young tableaux of order n *)
(* The list is of length (n-1)!, with element corresponding to the component of the vector associated *)
(* with a particular permutation of order n-1 *)
GeneratingFunctionForAllSortingNetworks[n_, varname_]:=Module[
		{sortingnets, lastswappermutations, generatingfactors, allpermutations},
	allpermutations = Permutations[Range[n-1]];
	sortingnets = SortingNetworks[n];
	lastswappermutations = Map[LastSwapTimesPermutation, sortingnets];
	generatingfactors = Map[GeneratingFactorForSortingNetwork[#, varname]&, sortingnets];
	Table[Total[Map[Part[generatingfactors,#[[1]]]&,Position[lastswappermutations,allpermutations[[k]]]]] ,{k,1,Length[allpermutations]}]
];

(* GeneratingFunctionsForAllStaircaseShapeSYTAndSortingNetworks[n, varname] returns a list of (n-1)! triples *)
(* { p, f, g } where p is a permutation of order n-1, f is the component of the generating function F_n(...) associated *)
(* with the permutation p, and f is the component of the generating function F_n(...) associated with the same permutation *)
GeneratingFunctionsForAllStaircaseShapeSYTAndSortingNetworks[n_, varname_]:=Module[{allperms, genfun1, genfun2},
	allperms = Permutations[Range[n-1]];
	genfun1 = GeneratingFunctionForAllStaircaseShapeSYT[n, varname];
	genfun2 = GeneratingFunctionForAllSortingNetworks[n, varname];
	Table[{allperms[[k]], genfun1[[k]], genfun2[[k]]}, {k, 1, Length[allperms]}]
];

(* VerifyStaircaseSYTSortingNetworksIdentity[n] verifies the combinatorial identity between the generating *)
(* functions F_n(...) and G_n(...) (Conjecture 2 in the paper) *)
(* Valid values for n are 3, 4, 5, and 6. (Values higher than 6 are technically allowed but the computation *)
(* will take an impractically long amount of time.) *)
(* The function prints its output regarding the correctness of the identity, and returns the output of *)
(* GeneratingFunctionsForAllStaircaseShapeSYTAndSortingNetworks[n,varname], used in the calculation, as its return value *)
VerifyStaircaseSYTSortingNetworksIdentity[n_]:=Module[{twogeneratingfunctions, isidentitytrue, formattedoutput},
	twogeneratingfunctions = GeneratingFunctionsForAllStaircaseShapeSYTAndSortingNetworks[n, x];
	isidentitytrue = Apply[And,Map[Simplify[#[[2]]-#[[3]]==0]&, twogeneratingfunctions]];
	formattedoutput = "The conjecture is " <> If[isidentitytrue,"TRUE","FALSE"] <> " for n="<>ToString[n]<>" (summmed over "<>ToString[NumberOfStaircaseShapeSYT[n]]<>" staircase shape tableaux and sorting networks)";
	Print[formattedoutput];
	Return[twogeneratingfunctions];
];



Print["OrientedSwaps.m package (Version 1.0) loaded. Type: 'OrientedSwapsHelp[]' for help."];
