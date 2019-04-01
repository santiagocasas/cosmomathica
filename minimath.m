location=If[aa=DirectoryName@$InputFileName; aa=="", NotebookDirectory[], aa] 
SetDirectory[location<>"/ext/mgcamb"]; 
link=Install[location<>"/ext/math_link"];

intslcdm={0,1,1,1,4,1,1,0,1,0,1,1,1,1,0,0,1,1,1,1500,400,50,2,1,0,0,0,0,1,3};

floatslcdm={0.25,0.05,67.,0,0.,-1.,0.,2.7255,0.24,3.046,0,1,0.96,0,0,1,2.1*10^-9,0.05,0.05,0.,10.,1.,0.5,2,4000.,800.,3.35,1.5,0.};

intsMG={0,1,1,1,2,1,1,1,0,1,1,1,4,1,1,0,1,0,1,1,1,1,0,0,1,1,1,1500,400,50,2,1,0,0,0,0,1,3};

floatsMG={0.25,0.05,67.,0,0.,-1.,0.,2.7255,0.24,0.001,0.,0.,0.,0.,4.,0.,0.,-1.,0.,0.,0.,1.,1.,0.,0.545,0.,0.5,0.001,1.,0.0001,0.24,1.,0.0001,1.,3.046,0,1,0.96,0,0,1,2.1*10^-9,0.05,0.05,0.,10.,1.,0.5,2,4000.,800.,3.35,1.5,0.};

cambcode="mgcamb";

Which[cambcode=="camb",
floats = N/@floatslcdm; (*the N is just there to convert exact irrationals (like Pi) to floats*)
ints = intslcdm;
,
cambcode=="mgcamb",
floats = N/@floatsMG; (*the N is just there to convert exact irrationals (like Pi) to floats*)
ints = intsMG;
];


result=Global`CAMBrun[floats,ints];

Print["Dimensions of result[[1]] :"];
Print[Dimensions@result[[1]] ];
Print["Numeric result[[1]] ?"];
Print[VectorQ[(result[[1]]), NumericQ]];

 
