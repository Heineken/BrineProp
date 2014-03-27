within BrineProp.Examples;
model InitializationTest
  // Parameter n_g_norm_start in Baseprops zu Variable machen
package Medium = BrineProp.Brine_5salts_TwoPhase_3gas;
  Medium.BaseProperties props(phase=0);
  SI.Density d= props.d;  /**/

  Partial_Units.Molality[:] b_l=massFractionsToMolalities(props.X_l, Medium.MM_vec);

initial equation
// props.T = 273.16+60;
equation
//  (h_francke,val2) =  SpecificEnthalpies.specificEnthalpy_pTX_liq_Francke_cp(props.p,props.T,props.X);
//  h_salt_100[2:end]={0,0,0,0};

 props.n_g_norm_start={1,1,1,1}*time;

  props.p = 200e5 "1.01325e5";
//  props.p = 1.01325e5;

//  props.T = 273.16+80+time;

 props.T = 125+273.16;
/*  props.p = 9.13e5;*/

//props.Xi = {     0.08214,   0.0053126,     0.11052,   0*0.0011094,   0*0.0019676,  0.00018083,  0.00074452,  6.661e-005}     "Elvira 9/11";
//    props.Xi = {0.083945671051201,0.00253479771131107,0.122842299461699,0*0.000612116692496665,0*0.00214041137028575,  0.00016883,  0.00073459, 6.5652e-005}     "Elvira 2-2013 1.1775g/ml";
    props.Xi = {0.0839077010751,0.00253365118988,0.122786737978,0,0,7.2426359111e-05,0.000689505657647,6.14906384726e-05}
    "Elvira 2-2013 1.1775g/ml V2";

algorithm
//  print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
//  print(Modelica.Math.Matrices.toString(transpose([props.Xi])));

  annotation (experiment(
      StartTime=0.01,
      StopTime=0.99,
      __Dymola_NumberOfIntervals=50),
      __Dymola_experimentSetupOutput);
end InitializationTest;
