(* ::Package:: *)

(*Cosmomathica version 0.9, September 2017. By Santiago Casas, Adrian Vollmer, Institute for Theoretical Physics, University Heidelberg.*)


BeginPackage["cosmomathica`interface`"]


Transfer::usage="Transfer[OmegaM, fBaryon, Tcmb, h] provides an interface to Eisenstein & Hu's fitting formula for the transfer function (high baryon case). It takes the total matter density \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(M\)]\), the fraction of baryons \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(b\)]\)/\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(M\)]\), the CMB temperature and the dimensionless Hubble constant h as input, and returns the sound horizon, the wavenumber \!\(\*SubscriptBox[\(k\), \(peak\)]\) where the spectrum has a maximum, the transfer function for CDM, baryons, both, or with no baryons or no wiggles.";

TFPower::usage="TFPower[OmegaM, OmegaB, OmegaMN, DegenNu, OmegaL, h, z] provides an interface to Eisenstein & Hu's fitting formula for the transfer function (mixed matter case). It takes the total matter density \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(M\)]\), the baryon density \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(b\)]\), the massive neutrino density \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(\[Nu]\)]\), the integer number of degenerate massive neutrino species, the dark energy density \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(L\)]\), the dimensionless Hubble constant h, and the redshift z as input, and return the transfer function at the given redshift.";

Halofit::usage="Halofit[OmegaM, OmegaL, gammaShape, sigma8, ns, betaP, z0] provides an interface to the halofit algorithm by Robert E. Smith et al. (reimplemented in C by Martin Kilbinger). It takes the total matter density \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(M\)]\), the vacuum energy density \!\(\*SubscriptBox[\(\[CapitalOmega]\), \(L\)]\), a shape factor, \!\(\*SubscriptBox[\(\[Sigma]\), \(8\)]\), \!\(\*SubscriptBox[\(n\), \(s\)]\), \!\(\*SubscriptBox[\(\[Beta]\), \(p\)]\), and a fixed redshift \!\(\*SubscriptBox[\(z\), \(0\)]\) as input, and returns the nonlinear matter power spectrum (computed in three ways: linearly with the BBKS alogorithm, nonlinear with the algorithm by Peacock and Dobbs, or nonlinearly with the Halofit algorithm) at 20 different values of the scale factor and the convergence power spectrum in tabulated form.";

HaloFitCorrection::usage="HaloFitCorrection[LinearPS, OmegaM, OmegaL] applies nonlinear corrections to a linear power spectrum given in tabulated form {{\!\(\*SubscriptBox[\(k\), \(1\)]\),P(\!\(\*SubscriptBox[\(k\), \(1\)]\))},...}."

CosmicEmu::usage="CosmicEmu[omegaM, omegaB, sigma8, ns, w] provides an interface to the CosmicEmulator by Earl Lawrence. It takes \!\(\*SubscriptBox[\(\[Omega]\), \(M\)]\), \!\(\*SubscriptBox[\(\[Omega]\), \(b\)]\), \!\(\*SubscriptBox[\(\[Sigma]\), \(8\)]\), \!\(\*SubscriptBox[\(n\), \(s\)]\), and the equation of state w, and returns the nonlinear matter power spectrum at five different redshifts as well as z, H, d (all at last scattering), and the sound horizon.";

FrankenEmu::usage="FrankenEmu[omegaM, omegaB, h, sigma8, ns, w] provides an interface to FrankenEmu by Earl Lawrence. It takes \!\(\*SubscriptBox[\(\[Omega]\), \(M\)]\), \!\(\*SubscriptBox[\(\[Omega]\), \(b\)]\), h, \!\(\*SubscriptBox[\(\[Sigma]\), \(8\)]\), \!\(\*SubscriptBox[\(n\), \(s\)]\), and the equation of state w, and returns the nonlinear matter power spectrum at five different redshifts as well as z, H, d (all at last scattering), and the sound horizon. The Hubble parameter h can also be omitted, in which case it will be determined from the CMB just as in CosmicEmu. Additional cosmological parameters are only returned if h is missing.";

CAMB::usage="CAMB[OmegaC, OmegaB, h] :
provides an interface to CAMB by Antony Lewis and Anthony Challinor. It takes a few parameters as well as a number of options as input, and returns various cosmological quantities. The distinction between parameters and options is in principle arbitrary. However, since some physical parameters are often assumed to take on a default value, they are being interpreted as an option here. To see the default options, type `Options[CAMB]`.";

Class::usage="Class[\"parameter1\" -> \"value1\",...] runs Class. You need to specifiy all Class parameters in the form of options. Cosmomathica then writes all options into an ini-file, where each line contains \"parameter = value\". This file is fed to Class and the output is returned."; 




MGflag::usage="An option for CAMB/MGCAMB: 
          MG_flag = 0 :  default GR
          MG_flag = 1 :  pure MG models
          MG_flag = 2 :  alternative MG models
          MG_flag = 3 :  QSA models
"
pureMGflag::usage="An option for CAMB/MGCAMB: 
pure_MG_flag = 1 : mu, gamma parametrization
pure_MG_flag = 2 : mu, sigma parametrization
pure_MG_flag = 3 : Q, R  parametrization
"
altMGflag::usage="An option for CAMB/MGCAMB: alt_MG_flag = 1 : Linder Gamma parametrization"

QSAflag::usage="An option for CAMB/MGCAMB: 
QSA_flag = 1 : f(R)
QSA_flag = 2 : Symmetron      
QSA_flag = 3 : Dilaton
QSA_flag = 4 : Hu-Sawicki f(R)"

mugammapar::usage="An option for CAMB/MGCAMB
mugamma_par = 1 : BZ parametrization
mugamma_par = 2 : Planck parametrization
"
musigmapar::usage="An option for CAMB/MGCAMB:
musigma_par = 1 : DES parametrization
"
QRpar::usage="An option for CAMB/MGCAMB:
QR_par = 1 : (Q,R)(arXiv:1002.4197 )
QR_par = 2 : (Q0,R0,s)(arXiv:1002.4197 )"

DEmodel::usage="An option for CAMB/MGCAMB:
DE_model = 0 : LCDM
DE_model = 1 : wCDM        
DE_model = 2 : (w0,wa)CDM  
DE_model = 3 : user defined
"

GRtrans::usage="Parameter options for MGCAMB"
B1::usage="Parameter options for MGCAMB"
lambda12::usage="Parameter options for MGCAMB"
B2::usage="Parameter options for MGCAMB"
lambda22::usage="Parameter options for MGCAMB"
ss::usage="Parameter options for MGCAMB"
E11::usage="Parameter options for MGCAMB"
E22::usage="Parameter options for MGCAMB"
mu0::usage="Parameter options for MGCAMB"
sigma0::usage="Parameter options for MGCAMB"
MGQfix::usage="Parameter options for MGCAMB"
MGRfix::usage="Parameter options for MGCAMB"       
Qnot::usage="Parameter options for MGCAMB"
Rnot::usage="Parameter options for MGCAMB"
sss::usage="Parameter options for MGCAMB"
Lindergamma::usage="Parameter options for MGCAMB"
betastar::usage="Parameter options for MGCAMB"
astar::usage="Parameter options for MGCAMB"
xistar::usage="Parameter options for MGCAMB"
beta0::usage="Parameter options for MGCAMB"
xi0::usage="Parameter options for MGCAMB"
DilS::usage="Parameter options for MGCAMB"
DilR::usage="Parameter options for MGCAMB"
FR0::usage="Parameter options for MGCAMB" 
FRn::usage="Parameter options for MGCAMB"

w0ppf::usage="An option for CAMB";
wappf::usage="An option for CAMB";
Tcmb::usage="An option for CAMB";
OmegaNu::usage="An option for CAMB";
OmegaK0::usage="An option for CAMB";
YHe::usage="An option for CAMB";
MasslessNeutrinos::usage="An option for CAMB";
MassiveNeutrinos::usage="An option for CAMB";
NuMassDegeneracies::usage="An option for CAMB";
NuMassFractions::usage="An option for CAMB";
ScalarInitialCondition::usage="An option for CAMB";
NonLinear::usage="An option for CAMB";
HalofitVersion::usage="An option for CAMB. Integer number that chooses the Halofit version used in the non-linear matter power spectrum.
1: Original, 2: Bird, 3: Peacock, 4: Takahashi, 5: Mead, 6: Halo model, 7: Casarini,
8: Mead2015, 9: Winther-f(R)";
WantCMB::usage="An option for CAMB";
WantTransfer::usage="An option for CAMB";
WantCls::usage="An option for CAMB";
ScalarSpectralIndex::usage="An option for CAMB";
ScalarRunning::usage="An option for CAMB";
TensorSpectralIndex::usage="An option for CAMB";
RatioScalarTensorAmplitudes::usage="An option for CAMB";
ScalarPowerAmplitude::usage="An option for CAMB";
PivotScalar::usage="An option for CAMB";
PivotTensor::usage="An option for CAMB";
DoReionization::usage="An option for CAMB";
UseOpticalDepth::usage="An option for CAMB";
OpticalDepth::usage="An option for CAMB";
ReionizationRedshift::usage="An option for CAMB";
ReionizationFraction::usage="An option for CAMB";
ReionizationDeltaRedshift::usage="An option for CAMB";
AccuracyBoost::usage="An option for CAMB";
TransferHighPrecision::usage="An option for CAMB";
TransferInterpMatterPower::usage="An option for CAMB";
WantScalars::usage="An option for CAMB";
WantVectors::usage="An option for CAMB";
WantTensors::usage="An option for CAMB";
WantZstar::usage="An option for CAMB";
WantZdrag::usage="An option for CAMB";
OutputNormalization::usage="An option for CAMB";
MaxEll::usage="An option for CAMB";
MaxEtaK::usage="An option for CAMB";
MaxEtaKTensor::usage="An option for CAMB";
MaxEllTensor::usage="An option for CAMB";
TransferKmax::usage="An option for CAMB";
TransferKperLogInt::usage="An option for CAMB";
TransferRedshifts::usage="An option for CAMB";
AccuratePolarization::usage="An option for CAMB";
AccurateReionization::usage="An option for CAMB";
AccurateBB::usage="An option for CAMB";
DoLensing::usage="An option for CAMB";
OnlyTransfers::usage="An option for CAMB";
DerivedParameters::usage="An option for CAMB";
MassiveNuMethod::usage="An option for CAMB";
debugIntFloats::usage="An option for CAMB";



CAMB::WrongOptionsFloat="List of floats passed to CAMB contains non-numeric values. Check arguments and passed options."
CAMB::WrongOptionsInteger="List of integers passed to CAMB contains non-integer values. Check arguments and passed options."
CAMB::Eigenstates="NuMassEigenstates and NuMassFractions must have the same length  (can be zero).";
CAMB::InvalidOption="Option `1` is '`2`', but must be one of the following: `3`";
CAMB::Lists="The following options need to be non-empty lists of the same length: `1`";
CAMB::Error="CAMB exited with the following error code: `1`";
Interface::LinkBroken="`1` crashed. See if there is anything useful on stdout.";
Interface::OutsideBounds="Parameter out of bounds. `5` requires `3` <= `1` <= `4`, but you have `1`=`2`.";
Interface::NotInstalled="`1` appears to be unavailable on your system.";
Interface::MathLinkFail="The MathLink failed. Check if there is anything useful on stdout.";
Class::Error="`1` => `2`";


Begin["`Private`"]


Needs["NumericalCalculus`"]


$location=DirectoryName[$InputFileName];

$lightspeedc=299792.458;

bool2int[b_]:=If[b,1,0];
validatestring[val_,name_,poss_]:=If[!MemberQ[poss,val],Message[CAMB::InvalidOption,name,val,StringJoin@@Riffle[poss,", "]];Abort[]];
validatelimits[val_,name_,lower_,upper_,module_]:=If[!(lower<=val<=upper),Message[Interface::OutsideBounds,name,val,lower,upper,module];Abort[]];
validatelists[list_]:=If[1!=Length@Union[Length/@list],Message[CAMB::Lists,list];Abort[]];
validateresult[x_,name_]:=Switch[x,$Failed,Message[Interface::LinkBroken,name];Abort[];False,Null,Message[Interface::NotInstalled,name];Abort[];False,_,True];


reshape[list_,dimensions_]:=First[Fold[Partition[#1,#2]&,Flatten[list],Reverse[dimensions]]]


verbosestring="background_verbose = 1
thermodynamics_verbose = 1
perturbations_verbose = 1
bessels_verbose = 1
transfer_verbose = 1
primordial_verbose = 1
spectra_verbose = 1
nonlinear_verbose = 1
lensing_verbose = 1
output_verbose = 1";

Class[options:OptionsPattern[]]:=Module[{link,inifile,tempdir,h,result,limits,klen,taulen,unflatten},
tempdir=CreateDirectory[];
inifile=tempdir<>"/mathematica.ini";
Export[inifile,Table[o[[1]]<>" = "<>o[[2]],{o,{options}}]~Join~{"root = "<>tempdir<>"/out_",verbosestring},"Text"];

link=Install[$location<>"ext/math_link"];
result=Global`MlClass[inifile];
validateresult[result,"Class"];
Uninstall[link];
DeleteDirectory[tempdir,DeleteContents->True];

If[And@@StringQ/@result,Message[Class::Error,result[[1]],result[[2]]];Return[$Failed];Abort[]];

limits=Select[First@result,#!=0&]/.{-1->"kvalues",-2->"tauvalues",-3->"sigma8",-4->"transfer",-5->"PSlinear",-6->"PSnonlinear",-7->"background"};
limits=Partition[limits,2];

(*reformat the numbers*)
result=Table[
With[{count=If[i==1,0,Total@limits[[;;i-1,2]]]},
Class[limits[[i,1]]]->result[[2,count+1;;count+limits[[i,2]]]]],
{i,Length@limits}];

klen=Length[Class["kvalues"]/.result];
taulen=Length[Class["tauvalues"]/.result];
h=First[Class["background"]/.result];


result=If[MemberQ[{"kvalues","tauvalues","PSlinear","PSnonlinear"},#[[1,1]]],#[[1]]->Exp@#[[2]],#]&/@result;
result=If[#[[1,1]]=="transfer",#[[1]]->Partition[#[[2]],Length[#[[2]]]/klen],#]&/@result;
result=If[#[[1,1]]=="kvalues",#[[1]]->#[[2]]/h,#]&/@result;
result=If[MemberQ[{"PSlinear","PSnonlinear"},#[[1,1]]],#[[1]]->Partition[#[[2]],klen],#]&/@result;
result
];


CAMB[OmegaC_?NumericQ,OmegaB_?NumericQ,h_?NumericQ,opts:OptionsPattern[]]:=Module[{j,
  link,result,resultfloat,resultint,floats,ints,initialcond,
  nonlinear,massivenu,limits,check,getDimensions,dimensions,array,redshifts},

getDimensions[list_]:=Module[{i,r},
i=1;r={};
While[i<Length@list,AppendTo[r,list[[i+1;;i+list[[i]]]]];i=i+1+list[[i]]];
r];


(*some parameters must be within certain limits*)
limits={{ReionizationFraction,0,1.5},
  {OpticalDepth,0,.9},
  {h,"h",.2,1.},
  {Tcmb,2.7,2.8},
  {YHe,.2,.8},
  {HalofitVersion, 1,9},
  {MasslessNeutrinos,0,3.1},
(*MassiveNeutrinos?*)
{OmegaB,"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(B\)]\)",.001/h^2,1./h^2},
{OmegaC ,"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(C\)]\)",0./h^2,3./h^2}
};


validatelists[OptionValue@{ScalarSpectralIndex,ScalarRunning,TensorSpectralIndex,RatioScalarTensorAmplitudes,ScalarPowerAmplitude}];
validatelists[OptionValue@{NuMassDegeneracies,NuMassFractions}];
validatelimits[Sequence@@If[NumericQ[#[[1]]],#,{OptionValue[#[[1]]],ToString[#[[1]]]}~Join~Rest@#],"CAMB"]&/@limits;

initialcond={"vector","adiabatic","iso_CDM","iso_baryon","iso_neutrino","iso_neutrino_vel"};
nonlinear={"none","pk","lens","both"};
massivenu={"int","trunc","approx","best"};

validatestring[OptionValue[ScalarInitialCondition],"ScalarInitialCondition",initialcond];
validatestring[OptionValue[NonLinear],"NonLinear",nonlinear];
validatestring[OptionValue[MassiveNuMethod],"MassiveNuMethod",massivenu];

tkmaxh=OptionValue[TransferKmax]*h;  (*rescale passed k_max with h, so that output of CAMB is correctly in h/Mpc*)

floatsMG=OptionValue[#]&/@{GRtrans,B1,lambda12,B2,lambda22,ss,E11,E22,mu0,sigma0,MGQfix,MGRfix,Qnot,Rnot,sss,Lindergamma,betastar,astar,xistar,beta0,xi0,DilS,DilR,FR0,FRn};

floats=Flatten@{OmegaC,OmegaB,h*100,OptionValue[#]&/@{OmegaNu,OmegaK0,w0ppf,wappf,Tcmb,YHe},
                floatsMG,
                OptionValue[#]&/@{MasslessNeutrinos,NuMassDegeneracies,
                NuMassFractions,ScalarSpectralIndex,ScalarRunning,TensorSpectralIndex,RatioScalarTensorAmplitudes,
                ScalarPowerAmplitude,PivotScalar,PivotTensor,
                OpticalDepth,ReionizationRedshift,ReionizationFraction,
                ReionizationDeltaRedshift,AccuracyBoost,MaxEtaK,MaxEtaKTensor}, tkmaxh,
                Reverse@Sort@OptionValue@TransferRedshifts};

intsMG=OptionValue[#]&/@{MGflag,pureMGflag,altMGflag,QSAflag,mugammapar,musigmapar,QRpar,DEmodel};
(*NuMassDegeneracies ignored at the moment in camb_wrapper*)
ints=Flatten@{intsMG, OptionValue[MassiveNeutrinos],Length@OptionValue[NuMassFractions],
  Position[initialcond,OptionValue[ScalarInitialCondition]][[1,1]]-1,
  Position[nonlinear,OptionValue[NonLinear]][[1,1]]-1,
  OptionValue[HalofitVersion],
  Length@OptionValue@ScalarSpectralIndex,
  bool2int/@OptionValue@{DoReionization,UseOpticalDepth,TransferHighPrecision,
    TransferInterpMatterPower,WantCMB,WantTransfer,
    WantCls,WantScalars,WantVectors,WantTensors,WantZstar, WantZdrag},
  OptionValue[#]&/@{OutputNormalization,MaxEll,MaxEllTensor,TransferKperLogInt},Length@OptionValue@TransferRedshifts,
  bool2int/@OptionValue@{AccuratePolarization,AccurateReionization,AccurateBB,DoLensing,OnlyTransfers,DerivedParameters},
  Position[massivenu,OptionValue[MassiveNuMethod]][[1,1]]-1};

If[OptionValue[debugIntFloats]==True,
  Print["list of ints and floats passed to camb_wrapper: "];
  Print[N@(ints),N@(floats)];
];

If[AllTrue[floats, NumericQ]==False,
  Message[CAMB::WrongOptionsFloat];
  Abort[]
];
If[AllTrue[ints, IntegerQ]==False,
 Message[CAMB::WrongOptionsInteger];
  Abort[]
];

SetDirectory[$location<>"ext/mgcamb"];
link=Install[$location<>"ext/math_link"];
result=Global`CAMBrun[N/@floats,ints];
ResetDirectory[];
If[!validateresult[result,"CAMB"],Return[$Failed];Abort[]];
Uninstall[link];

{resultfloat,resultint}=result;
If[resultint[[1]]!=0,Message[CAMB::Error,resultint[[1]]];Return[$Failed];Abort[]];
resultint=Drop[resultint,1];
dimensions=getDimensions@resultint;

Do[array=Take[resultfloat,Times@@d];
resultfloat=Drop[resultfloat,Times@@d];
AppendTo[resultfloat,reshape[array,d]],{d,dimensions}];

j=2;
DeleteCases[{
CAMB["age"]->resultfloat[[1,1]],
CAMB["zstar"]->resultfloat[[1,2]],
CAMB["rstar"]->resultfloat[[1,3]],
CAMB["100thetastar"]->resultfloat[[1,4]],
CAMB["zdrag"]->resultfloat[[1,5]],
CAMB["rdrag"]->resultfloat[[1,6]],
CAMB["kD"]->resultfloat[[1,7]],
CAMB["100thetaD"]->resultfloat[[1,8]],
CAMB["zEQ"]->resultfloat[[1,9]],
CAMB["100thetaEQ"]->resultfloat[[1,10]],
CAMB["H(z)"]->Flatten[resultfloat[[-3]]],
CAMB["sigma8(z)"]->Flatten[resultfloat[[-2]]],
CAMB["f(z)"]->Flatten[resultfloat[[-1]]/((resultfloat[[-2]])^2)], (* sigma2_vdelta_8 / sigma_8 ^2 = f(z)*)
CAMB["D+(z)"]->Flatten[resultfloat[[-2]]]/(Flatten[resultfloat[[-2]]][[-1]]),
If[OptionValue[WantScalars],CAMB["CLscalar"]->Transpose@First@Transpose[resultfloat[[j++]]]],
If[OptionValue[WantVectors],CAMB["CLvector"]->Transpose@First@Transpose[resultfloat[[j++]]]],
If[OptionValue[WantTensors],CAMB["CLtensor"]->Transpose@First@Transpose[resultfloat[[j++]]]],
If[OptionValue[WantTransfer],CAMB["redshifts"]->(redshifts=resultfloat[[-4]])],
If[OptionValue[WantTransfer],Sequence@@Flatten@{
CAMB["PSlinear"]->Table[Exp@Transpose@{resultfloat[[j]],resultfloat[[j+1,All,k]]},{k,Length@redshifts}],
CAMB["PSnonlinear"]->Table[Exp@Transpose@{resultfloat[[j]],resultfloat[[j+2,All,k]]},{k,Length@redshifts}
]}~Join~{CAMB["Transfer"]->Table[Transpose[resultfloat[[j+3,All,All,k]]],{k,Length@redshifts}]}
](*,CAMB["ints"]->resultint,CAMB["floats"]->resultfloat*)
},Null]
];
Options[CAMB]={
MGflag->0, pureMGflag->1, altMGflag->1, QSAflag->1, mugammapar->2, musigmapar->1, QRpar->1, DEmodel->1,
GRtrans->0.001,
B1->0.0,
lambda12->0.0,
B2->0.0,
lambda22->0.0,
ss->4.0,
E11->0.0,
E22->0.0,
mu0->-1.0,
sigma0->0.0,
MGQfix->0.0,
MGRfix->0.0,
Qnot->1.0,
Rnot->1.0,
sss->0.0,
Lindergamma->0.545,
betastar->0.0,
astar->0.5,
xistar->0.001,
beta0->1.0,
xi0->0.0001,
DilS->0.24,
DilR->1.0,
FR0->0.0001,
FRn->1.0, 
w0ppf->-1.0, wappf->0., Tcmb->2.7255,
  OmegaNu->0, OmegaK0->0., YHe->.24,MasslessNeutrinos->3.046,
  MassiveNeutrinos->0,NuMassDegeneracies->{0},NuMassFractions->{1},
  ScalarInitialCondition->"adiabatic",NonLinear->"pk", HalofitVersion->4, WantCMB->True,WantTransfer->True,WantCls->False,
  WantScalars->False,WantVectors->False,WantTensors->False,
  ScalarSpectralIndex->{.96},ScalarRunning->{0},
  TensorSpectralIndex->{0},RatioScalarTensorAmplitudes->{1},ScalarPowerAmplitude->{2.1*^-9},PivotScalar->.05,PivotTensor->.05,
  DoReionization->True,UseOpticalDepth->False,OpticalDepth->0.,ReionizationRedshift->10.,ReionizationFraction->1.,ReionizationDeltaRedshift->.5,
  AccuracyBoost->2,
  TransferHighPrecision->True, TransferInterpMatterPower->False,
  WantZstar->True, WantZdrag->True,OutputNormalization->1,MaxEll->1500,MaxEtaK->4000.,MaxEtaKTensor->800.,
  MaxEllTensor->400,TransferKmax->5,TransferKperLogInt->50,TransferRedshifts->{0.},
  AccuratePolarization->True,AccurateReionization->False,AccurateBB->False,DoLensing->False,
  OnlyTransfers->False,DerivedParameters->True,MassiveNuMethod->"best",
  debugIntFloats->False};

Options[validatestring]=Options[CAMB];

(*Transfer function*)
Transfer[OmegaM_?NumericQ,fBaryon_?NumericQ,Tcmb_?NumericQ,h_?NumericQ]:=Module[{result,link,krange,fitonek,horizon,peak,OmegaC},
validatelimits[fBaryon,"fBaryon",.0001,1,"Transfer"];
OmegaC=OmegaM-fBaryon*OmegaM;
link=Install[$location<>"ext/math_link"];

(*These random factors of h are ridiculous. Read the comments in the source (tf.c) very carefully*)
Global`TFSetParameters[N@OmegaM*h^2,N@fBaryon,N@Tcmb];
horizon=Global`TFSoundHorizon[N@OmegaM,N@fBaryon,N@h];
peak=Global`TFkPeak[N@OmegaM,N@fBaryon,N@h];
If[!validateresult[horizon,"transfer"],Return[$Failed];Abort[]];

fitonek[k_]:=Join[Global`TFFitOneK[k*N@h],{
Global`TFNoWiggles[N@OmegaM,N@fBaryon,N@h,N@Tcmb,k],
Global`TFZeroBaryon[N@OmegaM,N@h,N@Tcmb,k]}]; 
krange=10^Range[-6.,4.,.01];
result=Transpose[fitonek/@krange];
Uninstall[link];

{Transfer["soundhorizon"]->horizon,
Transfer["peak"]->peak,
Transfer["kvalues"]->krange,
Transfer["full"]->result[[1]],
Transfer["baryon"]->result[[2]],
Transfer["cdm"]->result[[3]],
Transfer["nowiggles"]->result[[4]],
Transfer["zerobaryons"]->result[[5]]}
];


TFPower[OmegaM_,OmegaB_,OmegaH_,Degen_,OmegaL_,h_,z_]:=Module[{result,link,krange},
krange=10^Range[-6.,4.,.01];
link=Install[$location<>"ext/math_link"];
Global`TFSetCosmology[N@OmegaM,N@OmegaB,N@OmegaH,Degen,N@OmegaL,N@h,N@z];
result=Global`TFOneK/@krange;
Uninstall[link];
Transpose[{krange,result}]
];


(**Halofit Correction formulas**)
deltaLsq[k_?NumericQ]:=linearPS[k]k^3/(2 Pi^2);
sigmaSq[R_?NumericQ] :=NIntegrate[deltaLsq[Exp@lk]Exp[-Exp[2lk] R^2],{lk,Log[kmin],Log[kmax]},PrecisionGoal->5,MaxRecursion->12];

(*if omega_m < 1 and lambda = 0*)
f1a=omegaM^-.0732;
f2a=omegaM^-.1423;
f3a=omegaM^.0725;

(*if omega_m + lambda = 1*)
f1b=omegaM^-.0307;
f2b=omegaM^-.0585;
f3b=omegaM^.0743;

(*interpolation of the two*)
f1=(1-omegaM-omegaL)/(1-omegaM)f1a+omegaL/(1-omegaM) f1b;
f2=(1-omegaM-omegaL)/(1-omegaM)f2a+omegaL/(1-omegaM)f2b;
f3=(1-omegaM-omegaL)/(1-omegaM)f3a+omegaL/(1-omegaM)f3b;

an:=10^(1.5222+2.8553 neff+2.3706 neff^2 +.9903neff^3+.2250neff^4-0.6038curve);
bn:=10^(-.5642+.5864neff+.5716neff^2-1.5474curve);
cn:=10^(.3698+2.0404neff+.8161neff^2+.5869curve);
mun:=0;
nun:=10^(5.2105+3.6902neff);
gamman:=.1971-.0843neff+.8460curve;
alphan:=Abs[6.0835+1.3373neff-.1959neff^2-5.5274curve];
betan:=2.0379-0.7354neff+.3157neff^2+1.249neff^3+.3980neff^4-.1682curve;

deltaHsqPrime[k_]:=an (k/ksig)^(3f1)/(1+bn (k/ksig)^f2+(cn f3 k/ksig)^(3-gamman));

deltaHsq[k_]:=deltaHsqPrime[k]/(1+ ksig/k (mun+ nun ksig/k));

f[y_]:=y/4+y^2/8;

deltaQsq[k_]:=deltaLsq[k]*(1+deltaLsq[k])^betan/(1+alphan deltaLsq[k])Exp[-f[k/ksig]];

deltaNLsq[k_]:=deltaQsq[k]+deltaHsq[k];
nonlinearPS[k_]:=deltaNLsq[k]*2Pi^2/k^3;

HaloFitCorrection[LinearPS_,OmegaM_,OmegaL_]:=Block[{linearPS=Interpolation[LinearPS],omegaM=OmegaM,omegaL=OmegaL,kmin=LinearPS[[1,1]],kmax=LinearPS[[-1,1]]},
ksig=Exp@x/.FindRoot[sigmaSq[1/Exp@x]==1,{x,-3,3},Method->"Secant"];
neff=-3+2NIntegrate[deltaLsq[Exp@lk]Exp[2lk]/ksig^2Exp[-Exp[2lk]/ksig^2],{lk,Log[kmin],Log[kmax]},PrecisionGoal->5,MaxRecursion->12];
curve=(3+neff)^2+4NIntegrate[deltaLsq[Exp@lk](Exp[2lk]/ksig^2-Exp[4lk]/ksig^4)Exp[-Exp[2lk]/ksig^2],{lk,Log[kmin],Log[kmax]},PrecisionGoal->5,MaxRecursion->12];
(*neff=-3-ND[Log@sigmaSq[Exp@lR],lR,Log@ksig];
curve=-ND[Log@sigmaSq[Exp@lR],{lR,2},Log@ksig];*)
Table[{k,nonlinearPS[k]},{k,LinearPS[[All,1]]}]
];





Halofit[OmegaM_?NumericQ,OmegaL_?NumericQ,gammaShape_?NumericQ,sigma8_?NumericQ,ns_?NumericQ,betaP_?NumericQ,z0_?NumericQ]:=Module[{link,Tf={},Kappa={},unitsfact,arange,krange,ellrange,labels,limits,parameters,check},
link=Install[$location<>"ext/math_link"];

arange=Most[10^Range[-2,0,.1]]~Join~{.99999};
krange=10^Range[-4,4,.1];
unitsfact=$lightspeedc/100;
krange=krange*unitsfact;(*halofit uses units H0/c Mpc^-1*)
ellrange=10^Range[-2,6,.1];

labels={"\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(M\)]\)","\!\(\*SubscriptBox[\(\[CapitalOmega]\), \(b\)]\)","\!\(\*SubscriptBox[\(\[Sigma]\), \(8\)]\)","\!\(\*SubscriptBox[\(n\), \(s\)]\)","Gamma","\!\(\*SubscriptBox[\(\[Beta]\), \(p\)]\)","\!\(\*SubscriptBox[\(z\), \(0\)]\)"};
limits={{.1,1.5},{.1,1.5},{.1,1.5},{.7,1.3},{0.05,.5},{1.0,2.0},{.2,1.5}};
(*these are soft limits as given by the authors of halofit+. Results may not be reliable if parameters are outside these bounds*)
parameters={OmegaM,OmegaL,sigma8,ns,gammaShape,betaP,z0};

check=(#[[2,1]]<=#[[1]]<=#[[2,2]])&/@Transpose[{parameters,limits}];
Do[If[!check[[i]],Message[Interface::OutsideBounds,labels[[i]],parameters[[i]],limits[[i,1]],limits[[i,2]],"Halofit"]],{i,Length@check}];
If[!And@@check,Abort[]];


Do[
Global`HFSetParameters[N@OmegaM,N@OmegaL,N@gammaShape,N@sigma8,N@ns,N@betaP,N@z0,i];
AppendTo[Tf,Table[Global`HFGetPkNL[a,k],{a,arange},{k,krange}]];
AppendTo[Kappa,Table[Global`HFGetKappa[ell],{ell,ellrange}]],
{i,0,2}];
validateresult[Global`HFGetKappa[10.],"halofit"];
Uninstall[link];

arange[[-1]]=1.;
(*Just return the raw numbers*)
{Halofit["avalues"]->arange,
Halofit["kvalues"]->krange/(unitsfact),
Halofit["ellvalues"]->ellrange,
Halofit["kappaBBKS"]->Kappa[[1]],Halofit["BBKS"]->Tf[[1]]*((unitsfact)^3),
Halofit["kappaPD96"]->Kappa[[2]],Halofit["PD96"]->Tf[[2]]*((unitsfact)^3),
Halofit["kappaHalofit"]->Kappa[[3]],Halofit["Halofit"]->Tf[[3]]*((unitsfact)^3)}
];


CosmicEmu[omegaM_?NumericQ,omegaB_?NumericQ,sigma8_?NumericQ,ns_?NumericQ,w_?NumericQ]:=Module[{link,result,labels,limits,parameters,check,arange},

labels={"\!\(\*SubscriptBox[\(\[Omega]\), \(M\)]\)","\!\(\*SubscriptBox[\(\[Omega]\), \(b\)]\)","\!\(\*SubscriptBox[\(\[Sigma]\), \(8\)]\)","\!\(\*SubscriptBox[\(n\), \(s\)]\)","w"};
limits={{.12,.155},{.0214,.0235},{.61,.9},{.85,1.05},{-1.3,-.7}};
(*these are hard limits as given by the authors of the cosmic emulator - the program will crash if any parameter is outside its bounds*)
parameters={omegaM,omegaB,sigma8,ns,w};
arange=Range[.5,1.,.1];

check=(#[[2,1]]<=#[[1]]<=#[[2,2]])&/@Transpose[{parameters,limits}];
Do[If[!check[[i]],Message[Interface::OutsideBounds,labels[[i]],parameters[[i]],limits[[i,1]],limits[[i,2]],"CosmicEmu"]],{i,Length@check}];
If[!And@@check,Abort[]];

link=Install[$location<>"ext/math_link"];
result=Table[{Transpose@Partition[#[[1]],Length@#[[1]]/2],#[[2]]}&@Global`CEGetPkNL[N@omegaM,N@omegaB,N@ns,N@sigma8,N@w,1/a-1],{a,arange}];
 (*CosmicEmu only does these five redshifts, everything else is interpolated*)
validateresult[(result[[All,2]])[[1,4]],"CosmicEmu"];
Uninstall[link];

(*Just return the raw numbers*)
{CosmicEmu["zvalues"]->Table[1/a-1,{a,arange}],
CosmicEmu["pk"]->result[[All,1]],
CosmicEmu["soundhorizon"]->(result[[All,2]])[[1,1]],
CosmicEmu["zlss"]->(result[[All,2]])[[1,2]],
CosmicEmu["dlss"]->(result[[All,2]])[[1,3]],
CosmicEmu["hubblecmb"]->(result[[All,2]])[[1,4]]}
];


FrankenEmu[omegaM_?NumericQ,omegaB_?NumericQ,h_?NumericQ,sigma8_?NumericQ,ns_?NumericQ,w_?NumericQ]:=Module[
{link,result,labels,limits,parameters,check, arange},

arange=Range[.2,1.,.1];
labels={"\!\(\*SubscriptBox[\(\[Omega]\), \(M\)]\)","\!\(\*SubscriptBox[\(\[Omega]\), \(b\)]\)","\!\(\*SubscriptBox[\(\[Sigma]\), \(8\)]\)","\!\(\*SubscriptBox[\(n\), \(s\)]\)","w","h"};
limits={{.12,.155},{.0215,.0235},{.6,.9},{.85,1.05},{-1.3,-.7},{.55,.85}};
(*these are hard limits as given by the authors of the cosmic emulator - the program will crash if any parameter is outside its bounds*)
parameters={omegaM,omegaB,sigma8,ns,w,If[h<0,.7,h]};

check=(#[[2,1]]<=#[[1]]<=#[[2,2]])&/@Transpose[{parameters,limits}];
Do[If[!check[[i]],Message[Interface::OutsideBounds,labels[[i]],parameters[[i]],limits[[i,1]],limits[[i,2]],"FrankenEmu"]],{i,Length@check}];
If[!And@@check,Abort[]];

link=Install[$location<>"ext/math_link2"];
result=Table[{Transpose@Partition[#[[1]],Length@#[[1]]/2],#[[2]]}&@Global`FrankenCEGetPkNL[N@omegaM,N@omegaB,N@h,N@ns,N@sigma8,N@w,1/a-1],{a,arange}];
 (*CosmicEmu only does these five redshifts, everything else is interpolated*)
validateresult[(result[[1,1]]),"FrankenEmu"];
Uninstall[link];

(*Just return the raw numbers*)
{FrankenEmu["zvalues"]->Table[1/a-1,{a,arange}],
FrankenEmu["pk"]->result[[All,1]]}~Join~If[h<0,
{FrankenEmu["soundhorizon"]->(result[[All,2]])[[1,1]],
FrankenEmu["zlss"]->(result[[All,2]])[[1,2]],
FrankenEmu["dlss"]->(result[[All,2]])[[1,3]],
FrankenEmu["hubblecmb"]->(result[[All,2]])[[1,4]]},{}]
];
FrankenEmu[omegaM_?NumericQ,omegaB_?NumericQ,sigma8_?NumericQ,ns_?NumericQ,w_?NumericQ]:=FrankenEmu[omegaM,omegaB,-1.,sigma8,ns,w];



End[ ]
EndPackage[ ]



