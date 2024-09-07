(* ::Package:: *)

BeginPackage["Rational`Interpolation`"]
w::Usage="w[i_,xList_] -> Calculates orthogonal polynomial weights"
CalcInitialWeights::Usage="CalcInitialWeights[xi_, fi_] -> Calculates the initial weights for the TTR algorithm"
TTR::Usage="TTR[xi_, fi_] -> Performs the TTR alogorithm to calcuate matrix needed for coefficients matrix"
CoefficientsMatrix::Usage="CoefficientsMatrix[Q_, xi_,alpha_, beta_] -> Calculates the denominator coefficients matrix"
UpdateCoefficientValues::Usage="UpdateCoefficientValues[A_,Q_,P_,fi_ ] -> Calculates the numerator coefficients matrix using 1/fi"
GetPolynomial::Usage="GetPolynomial[B_, A_, degdenom_, degnumer_] -> Constructs rational polynomail from B,A to specific degrees"
yPointsFromFunction::Usage="yPointsFromFunction[points_, f_] -> Calculates the y-points if user provides function instead of y-points"
ConvertEquation::Usage="ConvertEquation[eq_, var_] -> Creates the function needed for yPointsFromFunction to work"
PerformRationalInterpolation::Usage="PerformRationalInterpolation[xPoints_, yPoints_, DegreeNumer_,DegreeDenom_] -> Peforms the rational interpolation"
Plotpq::Usage="Plotpq[p,q,xi] -> Plot the results of p/q"
PlotCompare::Usage="PlotCompare[p_,q_ , xi_, func_:None] -> Plots results. If function was provided, it compares approximation to the results"

(*Orthogonal Polynomial Weights*)
w[i_,xList_]:=(
tempx=Drop[xList,{i}];
Times@@(xList[[i]]-tempx)
)
(*Initial Three-term reccurance weights*)
CalcInitialWeights[xi_, fi_]:= (
nx=Length[xi];
\[Beta]={0};R = {ConstantArray[1,nx]};
wi=w[#,xi]&/@Table[i,{i,1,nx}];
ci=fi/wi; \[Gamma]={ci . xi};\[Theta]={Total[ci]};
If[\[Theta][[1]]==0, (Print["This isn't going to work, symmetric y points!"]; Abort[]),\[Alpha]={\[Gamma][[-1]]/\[Theta][[-1]]};];
res00 = xi-\[Alpha][[1]];(* Calculating Row 2*)
AppendTo[R,res00];
AppendTo[\[Gamma],ci . (xi*res00^2)];
AppendTo[\[Theta],ci . (res00^2)];
AppendTo[\[Alpha],\[Gamma][[-1]]/\[Theta][[-1]]];
AppendTo[\[Beta],\[Theta][[2]]/\[Theta][[1]]];
{wi,ci,\[Gamma],\[Theta],\[Alpha],\[Beta]}
)
(*TTR Algorithm*)
TTR[xi_, fi_]:=(
n=Length[xi];
{wi,ci,\[Gamma],\[Theta],\[Alpha],\[Beta]}=CalcInitialWeights[xi,fi];
For[j=2, j<n, j++,
res00=(xi-\[Alpha][[j]])*R[[j,All]]-\[Beta][[j]]*R[[j-1,All]];
AppendTo[R,res00];
AppendTo[\[Gamma],ci . (xi*res00^2)];
AppendTo[\[Theta],ci . (res00^2)];
If[\[Theta][[-1]]==0, (Break[]; Print["Procedure Cut Short, we only calculated ",j," rows instead of ",n])];AppendTo[\[Alpha], \[Gamma][[-1]]/\[Theta][[-1]] ];AppendTo[\[Beta], \[Theta][[-1]]/\[Theta][[j]]];];
If[Cases[\[Theta], _Symbol]!={},(Print["Procedure Cut Short, we only calculated some rows instead of ",n]; Abort[])];
{R, \[Alpha],\[Beta]}
)
(*Coefficeint Matricies*)
CoefficientsMatrix[Q_, xi_,alpha_, beta_]:=(
(*Calculates the coefficient matrix for a polynomial, each row being a different degree*)
n=Length[xi];
R =DiagonalMatrix[Table[1, n]];R[[2,1]]=-1 alpha[[1]];
For[j=2, j<n, j++,
For[k=1, k<j+1, k++,
If[k==1,R[[j+1, k]] = - alpha[[j]] R[[j,k]] - beta[[j]] R[[j-1, k]],
R[[j+1, k]] = R[[j, k-1]] - alpha[[j]] R[[j,k]] - beta[[j]] R[[j-1, k]]
];]];
Return[R];
)
UpdateCoefficientValues[A_,Q_,P_,fi_ ]:=(
(*Function to rescale the coefficents for polynomial q when doing a rational interpolation of r_{m,n} = p_{m}(x)/q_{n}(x)*)
n=Length[fi]+1;tempA = A;
For[j=1, j<n,j++,
For[k=1, k<j+1, k++,
tempA[[j,k]]= (fi[[1]])*(Q[[n-j, 1]] / P[[j,1]])*A[[j, k]];
]];
tempA
)
(*Get Polynomial*)
GetPolynomial[B_, A_, degdenom_, degnumer_]:=(
(*p-Numerator, q-denominator*)
n00=Length[B];
p00 = x^Range[0,n00-1] . A[[degnumer+1, All]];
q00 = x^Range[0,n00-1] . B[[degdenom+1, All]];
{p00,q00}
)
(*Equation manipulation functions*)
yPointsFromFunction[points_, f_]:=(
f&[x]/.x->points
)
ConvertEquation[eq_, var_]:=(
eq00 =eq /.var -> #;
eq00
)
(*Actually perform the interpolation*)
PerformRationalInterpolation[xPoints_, yPoints_, DegreeNumer_,DegreeDenom_]:=(
x00=xPoints;y00=yPoints;var00=DeleteCases[DeleteDuplicates@Cases[yPoints,_Symbol,Infinity],E];
If[Length[x00]==0 Or Length[y00]==0, (Print["Please provide points to interpolate"]; Abort[])];
If[Length[var00]>1, (Print["Only provide an equation with one variable please"]; Abort[])];
If[
Head @ yPoints===Head@{},AppendTo[var00,x] ,y00=Chop[yPointsFromFunction[x00, ConvertEquation[y00, var00[[1]]]]]
];
If[Cases[y00, 0]!={}, (Print["y points contain 0! This is not supported... Abort!"]; Abort[])];
{Qx00,a00,b00}=TTR[x00, y00];
Bx00=CoefficientsMatrix[Qx00,x00,a00,b00];
{Px00,a200,b200}=TTR[x00,1/y00];
Ax00=UpdateCoefficientValues[CoefficientsMatrix[Px00,x00,a200,b200],Qx00,Px00, y00];
If[Dimensions[Qx00][[1]]<=TakeLargest[{DegreeNumer,DegreeDenom},1][[1]], (Print["TTR has not extrapolated to high enough degree.. Abort!"]; Abort[])];
{pm,qn} = GetPolynomial[Bx00, Ax00, DegreeNumer, DegreeDenom]/.x->var00[[1]];
{pm, qn}
)
Plotpq[p_, q_, xi_]:=(
Plot[p/q,{x,xi[[1]],xi[[-1]]},PlotLegends->{"p/q"}]
)
PlotCompare[func_,p_, q_, xi_]:=(
Plot[{func,p/q},{x,xi[[1]],xi[[-1]]},PlotLegends->{"f","p/q"}]
)

EndPackage[]












