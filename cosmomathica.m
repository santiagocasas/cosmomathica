(* ::Package:: *)

(************************************************************************)
(* This file was generated automatically by the Mathematica front end.  *)
(* It contains Initialization cells from a Notebook file, which         *)
(* typically will have the same name as this file except ending in      *)
(* ".nb" instead of ".m".                                               *)
(*                                                                      *)
(* This file is intended to be loaded into the Mathematica kernel using *)
(* the package loading commands Get or Needs.  Doing so is equivalent   *)
(* to using the Evaluate Initialization Cells menu command in the front *)
(* end.                                                                 *)
(*                                                                      *)
(* DO NOT EDIT THIS FILE.  This entire file is regenerated              *)
(* automatically each time the parent Notebook file is saved in the     *)
(* Mathematica front end.  Any changes you make to this file will be    *)
(* overwritten.                                                         *)
(************************************************************************)



BeginPackage["Cosmology`"]


$location=DirectoryName[$InputFileName];


Cosmology::usage="This list of rules defines the background cosmology."


SetCosmology::usage="Modify the fiducial cosmology by passing the new options to this function."


Fiducial::usage="Fiducial[X] returns the fiducial value of parameter X."


Redshift::usage="Redshift[a] returns the redshift as a function of the scale factor, i.e. z=1/a+1."


Hubble::usage="Hubble[z] returns the Hubble parameter at a given redshift z."


ScaleFactor::usage="ScaleFactor[t] returns a at cosmic time t."


DimensionlessHubble::usage="DimensionlessHubble[z] returns the Hubble parameter at a given redshift z divided by the Hubble constant today."


w::usage="w[z] returns the equation of state at redshift z using the CPL parameterization."


GrowthFunction::usage="GrowthFunction[z] returns the growth function at redshift z, parameterized by the growth index."


GrowthFactor::usage="GrowthFactor[z] returns the growth factor at redshift z, parameterized by the growth index."


OmegaM::usage="OmegaM[z] returns the relative matter density at redshift z."


Distance::usage="Distance[z] accepts the Option DistanceType with the possible values Comoving, TransverseComoving, AngularDiameter, Luminosity or LookBack. Distance[z1, z2] gives the distance between two galaxies at redshifts z1 and z2 along the line of sight."


TransferFunction::usage="TransferFunction[k,path] returns a list of five values using the fitting formulae by Eisenstein&Hu(1999). The values are: The full transfer function, the baryon transfer function, the CDM transfer function, the no-wiggles transfer function, and the no-baryon transfer function. The path to the MathLink executable needs to be passed by path."


LinearPS::usage="LinearPS[k,z,path] returns a list of five values using the fitting formulae by Eisenstein&Hu(1999). The values are: The full linear power spectrum, the baryon power spectrum, the CDM power spectrum, the no-wiggles power spectrum, and the no-baryon power spectrum."


PowerSpectrum::usage="PowerSpectrum[k,z] returns the power spectrum at wavenumber k and redshift z. The algorithm can be specified with PSType."


NonlinearPowerSpectrumTRG::usage="NonlinearPowerSpectrumTRG[k,LinearPS,Omega21,Omega22] returns the power spectrum at wavenumber k and redshift z using the background functions Omega21 and Omega22 as specified by the TRG package. LinearPS[k] is a one-parameter function representing the linear power spectrum."


NonlinearPowerSpectrumHalofit::usage="NonlinearPowerSpectrumHalofit[k,LinearPS] returns the power spectrum at wavenumber k and redshift z using the Halofit algorithm. LinearPS[k] is a one-parameter function representing the linear power spectrum."


(*Options*)
OmegaB::usage="Baryon density"
OmegaC::usage="Cold dark matter density"
OmegaL::usage="Dark energy density"
OmegaTot::usage="Total density"
OmegaK::usage="Curvature 'density'"
h::usage="Reduced Hubble parameter"
H0::usage="Hubble parameter"
omegaM::usage="Reduced matter density"
omegaL::usage="Reduced dark energy density"
omegaB::usage="Reduced baryon density"
omegaC::usage="Reduced cold dark matter density"
omegaK::usage="Reduced curvature 'density'"
w0::usage="Equation of state today"
w1::usage="Equation of state running parameter"
gamma::usage="Growth index"
sigma8::usage="Power spectrum amplitude"
ns::usage="Spectral index"
Tcmb::usage="Temperature of the CMB"
DistanceType::usage="Distance type"
Comoving::usage="Comoving distance"
TransverseComoving::usage="Transeverse comoving distance"
AngularDiameter::usage="Angular diameter distance"
Luminosity::usage="Luminosity distance"
TransferType::usage="Type of transfer function (FitFull etc)"
FitFull::usage="Full transfer function from the Eisenstein&Hu1999 fitting function"
FitBaryon::usage="Baryon transfer function from the Eisenstein&Hu1999 fitting function"
FitCDM::usage="CDM transfer function from the Eisenstein&Hu1999 fitting function"
FitNoWiggles::usage="'No wiggle' transfer function (wiggles removed, but baryon suppression still there) from the Eisenstein&Hu1999 fitting function"
FitZeroBaryon::usage="Transfer function in Universe with no baryons from the Eisenstein&Hu1999 fitting function"
PSType::usage="Specifies the algorithm used to compute the power spectrum"


Begin["`Private`"]


zmax=100;


SetCosmology[o:OptionsPattern[]]:=Module[{clist},
clist=List[o];
(DownValues[#]={Last@DownValues[#]})&/@{odesol, antideriv1, antideriv2,transferfitint, normalization,nonlineartrgint,nonlinearhaloint};
Unprotect[Cosmology];
Do[Cosmology[[First@First@Position[Cosmology,opt[[1]],{2}]]]=opt,{opt,clist}];
Protect[Cosmology];
];


Unprotect[Cosmology];
Cosmology={OmegaM->.25,
OmegaB->.05,
OmegaC->OmegaM-OmegaB,
OmegaL->.75,
OmegaTot->OmegaM+OmegaL,
OmegaK->1-OmegaTot,
h->.7,
H0->h *1./3000,(*Units: c/Mpc*)
omegaM->OmegaM h^2,
omegaL->OmegaL h^2,
omegaB->OmegaB h^2,
omegaC->OmegaC h^2,
omegaK->OmegaK h^2,
w0->-1,
w1->0,
gamma->.545,(*growth index*)
sigma8->.8,
ns->.9,
Tcmb->2.728};
Protect[Cosmology];
cosmoopts:=Cosmology~Join~{DistanceType->"Comoving",PSType->"EH",TransferType->FitFull};


OV[x_,opts:OptionsPattern[]]:=(OptionValue@x/.If[List@opts!={},List@opts,1->1]//.cosmoopts);
Options[OV]:=cosmoopts;


Fiducial[x_,opts:OptionsPattern[]]:=OV[x,opts];
Options[Fiducial]:=cosmoopts;


Redshift=Compile[{{a,_Real}},1/a+1];


w[z_,opts:OptionsPattern[]]:=OV[w0,opts]+OV[w1,opts]*z/(1+z);
Options[w]:=cosmoopts;
fw[z_,opts:OptionsPattern[]]=Integrate[(1+w[z,opts])/(1+zx),{zx,0,z},GenerateConditions->False];


(* H in units of H_0 *)
DimensionlessHubble[z_?NumericQ,opts:OptionsPattern[]]:=Sqrt[OV[OmegaM,opts]*(1+z)^3+OV[OmegaK,opts]*(1+z)^2+OV[OmegaL,opts]*Exp[3*fw[z,opts]]];
Options[DimensionlessHubble]:=cosmoopts;


Hubble[z_?NumericQ,opts:OptionsPattern[]]:=OV[H0,opts]DimensionlessHubble[z,opts];
Options[Hubble]:=cosmoopts;



(*TODO: Cover closed Universe; make sure a=0 at some point*)
odesol[opts:OptionsPattern[]]:=odesol[opts]=Quiet@NDSolve[{a'[t]/a[t]== Hubble[1/a[t]-1,opts],a[1]==1},a,{t,-2/Hubble[0,opts],1/Hubble[0,opts]}][[1,1,2]];
ScaleFactor[t_?NumericQ,opts:OptionsPattern[]]:=odesol[opts][t];
Options[ScaleFactor]:=cosmoopts;
Options[odesol]:=cosmoopts;


OmegaM[z_?NumericQ,opts:OptionsPattern[]]:=OV[OmegaM,opts]*(1+z)^3/DimensionlessHubble[z,opts]^2;
Options[OmegaM]:=cosmoopts;


sinn=Compile[{{x,_Real},{ok,_Real}},Evaluate@If[ok==0,x,Sinh[Sqrt[ok]x]/Sqrt[ok]]];



antideriv1[opts:OptionsPattern[]]:=antideriv1[opts]=NDSolve[{y'[zx]==1/Hubble[zx,opts],y[0]==0},y,{zx,0,zmax},MaxStepSize->.02][[1,1,2]];
antideriv2[opts:OptionsPattern[]]:=antideriv2[opts]=NDSolve[{y'[zx]==1/Hubble[zx,opts]/(1+z),y[0]==0},y,{zx,0,zmax},MaxStepSize->.02][[1,1,2]];
Distance[z_?NumericQ,opts:OptionsPattern[]]:=Switch[OptionValue@DistanceType,
"Comoving",antideriv1[opts][z],
"TransverseComoving",sinn[Distance[z,DistanceType->Comoving,opts],OV[OmegaK,opts]],
"AngularDiameter",Distance[z,DistanceType->TransverseComoving,opts]/(1+z),
"Luminosity",Distance[z,DistanceType->TransverseComoving,opts]*(1+z),
"LookBack",antideriv2[opts][z]
];
Distance[z1_?NumericQ,z2_?NumericQ,opts:OptionsPattern[]]:=Switch[OptionValue@DistanceType,
"Comoving",Distance[z2,opts]-Distance[z1,opts],
"TransverseComoving",Null,
"AngularDiameter",Null,
"Luminosity",Null,(*TODO*)
"LookBack",Distance[z2,opts]-Distance[z1,opts]];
Options[Distance]:=cosmoopts;


(*growth factor*)

GrowthFactor[z_?NumericQ,opts:OptionsPattern[]]:=OmegaM[z,opts]^OV[gamma,opts];
Options[GrowthFactor]:=cosmoopts;
antideriv3[opts:OptionsPattern[]]:=antideriv3[opts]=NDSolve[{y'[zx]==-GrowthFactor[zx,opts]/(1+zx),y[0]==0},y,{zx,0,zmax}][[1,1,2]];
Options[antideriv3]:=cosmoopts;
GrowthFunction[z_?NumericQ,opts:OptionsPattern[]]:=Exp[antideriv3[opts][z]];
Options[GrowthFunction]:=cosmoopts;


(*Transfer function*)

transferfitint[opts:OptionsPattern[]]:= transferfitint[opts]=Module[{result,link,fitonek},
link=Install[$location<>"/tf_fit_link"];
Global`TFSetParameters[OV[omegaM,opts],OV[OmegaB,opts]/OV[OmegaC,opts],OV[Tcmb,opts]];
fitonek[k_]:={Global`TFFitOneK[k],LinkRead[link],LinkRead[link],
Global`TFNoWiggles[OV[OmegaC,opts],OV[OmegaB,opts]/OV[OmegaC,opts],OV[h,opts],OV[Tcmb,opts],k/OV[h,opts]],
Global`TFZeroBaryon[OV[OmegaC,opts],OV[h,opts],OV[Tcmb,opts],k/OV[h,opts]]}; 
result=Interpolation@Table[{10^lk,fitonek[10^lk]},{lk,-6,4,.01}];
Uninstall[link];
result
];
TransferFunction[k_?NumericQ,opts:OptionsPattern[]]:=transferfitint[opts][k][[Switch[OV[TransferType,opts],FitFull,1,FitBaryon,2,FitCDM,3,FitNoWiggles,4,FitZeroBaryon,5]]];
Options[transferfitint]:=cosmoopts;
Options[TransferFunction]:=cosmoopts;


(*linear ps*)
(*See Amendola&Tsujikawa 4.212*)
(*LinearPS[k_,z_,opts:OptionsPattern[]]:=TransferFunction[k,opts]^2 * (k/OV[H0,opts])^OV[ns,opts]* GrowthFunction[z,opts]^2/OV[H0,opts]^3*2 Pi^2*(3.2*^-10);*)

tophatFT[x_]:=3/x^3(Sin[x]-x Cos[x]);
normalization[opts:OptionsPattern[]]:=normalization[opts]=OV[sigma8,opts]/Sqrt@NIntegrate[Exp[3 lk]/(2Pi^2) (TransferFunction[Exp[lk]*OV[h,opts],opts]^2 * Exp[lk*OV[ns,opts]])tophatFT[8Exp[lk]]^2,{lk,Log[10^-6],Log[10^4]},PrecisionGoal->5];
LinearPS[k_,z_,opts:OptionsPattern[]]:=TransferFunction[k*OV[h,opts],opts]^2 * k^OV[ns,opts]* GrowthFunction[z,opts]^2*normalization[opts]^2;Options[LinearPS]:=cosmoopts;
Options[normalization]:=cosmoopts;


(*Umberella function for the power spectrum*)
PowerSpectrum[k_?NumericQ,z_?NumericQ,opts:OptionsPattern[]]:=Module[{},
LinearPS[k,z,opts]
];
Options[PowerSpectrum]:=cosmoopts;


Get[$location<>"halofit.m"];
nonlinearhaloint[LinearPS_,oM_,oL_]:=nonlinearhaloint[LinearPS,oM,oL]=Module[{result},
HaloFitInit[LinearPS];
result=Table[{k,HaloFitPS[k,oM,oL,LinearPS]},{k,10^Range[-4,3,.01]}];
Interpolation@result
];
NonlinearPowerSpectrumHalofit[k_,LinearPS_,opts:OptionsPattern[]]:=nonlinearhaloint[LinearPS,OV[OmegaM,opts],OV[OmegaL,opts]][k];
Options[NonlinearPowerSpectrumHalofit]:=cosmoopts;


Get[$location<>"psnonlinear.m"];
Needs["NumericalCalculus`"];
nonlineartrgint[LinearPS_,Omega21_,Omega22_,dlogGdz_,locH0_,locns_]:=nonlineartrgint[LinearPS,Omega21,Omega22,dlogGdz,locH0,locns]=Module[{pstable,klist,result,zz},
pstable=Transpose@Table[{k,LinearPS[k]},{k,10^Range[-4,3,.01]}];
result=Transpose@TRGPowerSpectrum[pstable,35,0,dlogGdz,locH0,locns,Omega21,Omega22];
Interpolation[result[[All,{1,2}]]]
];
NonlinearPowerSpectrumTRG[k_,LinearPS_,Omega21_,Omega22_,opts:OptionsPattern[]]:=nonlineartrgint[LinearPS,Omega21,Omega22,Abs@ND[Log[GrowthFunction[zz,opts]],zz,35],OV[H0,opts],OV[ns,opts]][k];
Options[NonlinearPowerSpectrumTRG]:=cosmoopts;


End[ ]
EndPackage[ ]