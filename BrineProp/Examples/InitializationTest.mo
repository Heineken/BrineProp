within BrineProp.Examples;
model InitializationTest "Test influnce of start vector of VLE algorithm"
  /* In order to run this you have to 
   make n_g_norm_start in Baseprops a variable
   in BrineProp.PartialBrine_ngas_Newton.BaseProperties
  */
package Medium = BrineProp.Brine_5salts_TwoPhase_3gas;
  Medium.BaseProperties props(phase=0);
  SI.Density d= props.d;
  SI.MassFraction x= props.x;

initial equation
// props.T = 273.16+60;
equation

 props.n_g_norm_start={1,1,1,1}*time;

 props.T = 125+273.16;
  props.p = 9.13e5;

  props.Xi = {0.0839077010751,0.00253365118988,0.122786737978,0,0,7.2426359111e-05,0.000689505657647,6.14906384726e-05}
    "Feldbusch 2-2013 1.1775g/ml V2";
  annotation (experiment(
      StartTime=0.01,
      StopTime=0.99,
      __Dymola_NumberOfIntervals=50),
      __Dymola_experimentSetupOutput);
end InitializationTest;
