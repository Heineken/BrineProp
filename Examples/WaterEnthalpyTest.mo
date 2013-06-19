within BrineProp.Examples;
model WaterEnthalpyTest
//package Medium = Brine_5salts_noGas;
package Medium = Modelica.Media.Water.WaterIF97_pT;
  Medium.BaseProperties props1;
  Medium.BaseProperties props2;

  Modelica.SIunits.Density d= props1.d;
  Modelica.SIunits.Density d2= (props2.p-props1.p)/(props2.h-props1.h)
    "v=dh/dp";
equation
  props1.p = 1.0e5;
  props2.p = 1.1e5;
  props1.T = 273.16+60;
  props2.T = props1.T;

//  props.Xi = {     0.08214,   0.0053126,     0.11052,   0*0.0011094,   0*0.0019676} "Elvira 9/11";
//   props.Xi = {    0,   0,     0,   0,  0};

algorithm
  Modelica.Utilities.Streams.print(String(d2)+"="+String(d)+"?");
end WaterEnthalpyTest;
