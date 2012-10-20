within BrineProp.Examples;
model BrineProps1Phase
package Medium = Brine_5salts_noGas;
  Medium.BaseProperties props;

  Modelica.SIunits.Density d= props.d;
  Modelica.SIunits.SpecificEnthalpy h= props.h;
  Modelica.SIunits.Temp_C T_C= props.T-273.15;
  Modelica.SIunits.DynamicViscosity eta=Medium.dynamicViscosity(props.state);

equation
  props.p = 15e5;
//  props.h = 379778;
//  props.p = (10+time)*1.01325e5 "STP";
 props.T = 273.16+60+time*80;
/*
  props.p = 9.13e5;
  props.T = 273.16+99.61;
 */

//3-GAS
   props.Xi = {    0.081109,   0.0047275,     0.12519,   0*0.0013225,  0*0.0023842}
    "Entsprechend STR04/16 bei GG mit d_l=1199.48 kg/m³ - X_g stimmt";
//  props.Xi = {6*Salt_Data.M_NaCl/(1+6*Salt_Data.M_NaCl),0,0,0,0, 0,0,0} "NaCl, KCl, CaCl2, MgCl2, SrCl2, CO2, N2, CH4";

algorithm
//  Modelica.Utilities.Streams.print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  Modelica.Utilities.Streams.print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
end BrineProps1Phase;
