within BrineProp.Examples;
model BrineProps1Phase
package Medium = Brine_5salts_noGas;
//package Medium = Modelica.Media.Water.WaterIF97_pT;
  Medium.BaseProperties props;

  Modelica.SIunits.Density d= props.d;
  Modelica.SIunits.SpecificEnthalpy h= props.h;
  Modelica.SIunits.Temp_C T_C= props.T-273.15;
  Modelica.SIunits.DynamicViscosity eta=Medium.dynamicViscosity(props.state);
  Modelica.SIunits.ThermalConductivity lambda=Medium.thermalConductivity(props.state);
//  Modelica.SIunits.DynamicViscosity etaC=BrineProp.Viscosities.dynamicViscosity_Zhang_pTXd(props.p,props.T,props.X,props.d,Medium.MM_vec);
  Modelica.SIunits.SpecificHeatCapacity c_p_brine= Medium.specificHeatCapacityCp(props.state);
  Modelica.SIunits.SpecificHeatCapacity c_p_brine2=(Medium.specificEnthalpy_pTX(props.p,props.T+0.1,
                                                                                                   props.X)-Medium.specificEnthalpy_pTX(props.p,props.T-0.1,
                                                                                                    props.X))/0.2;
//  Modelica.SIunits.SurfaceTension sigma= Medium.surfaceTension_T(props.T);
/*  Real n_K=0*2.2846;
  Real n_Ca=0*.05/2;
  Real n_Na=5.998 "n_Ca*1.0001";*/
equation
  props.p = 10e5 "1.01325e5";
//  props.h = 379778;
//  props.p = (10+time)*1.01325e5 "STP";
// props.T = 273.16+60;
 props.T = 273.15+time;
/*
  props.p = 9.13e5;
  props.T = 273.16+99.61;
 */

//3-GAS
//   props.Xi = {    0.081109,   0.0047275,     0.12519,   0*0.0013225,  0*0.0023842}     "Entsprechend STR04/16 bei GG mit d_l=1199.48 kg/m³ - X_g stimmt";
//  props.Xi = {0.0828264747936192,0.00250100280471612,0.12120451826216,0.000603955715249515,0.00211187457541221};
  props.Xi = {     0.08214,   0.0053126,     0.11052,   0*0.0011094,   0*0.0019676}
    "Elvira 9/11";
//   props.Xi = {    0,   0,     0,   0,  0};
//  props.Xi = {6*SaltData.M_NaCl/(1+6*SaltData.M_NaCl),0,0,0,0} "NaCl, KCl, CaCl2, MgCl2, SrCl2";
//  props.Xi = {n_Na*SaltData.M_NaCl,n_K*SaltData.M_KCl,n_Ca*SaltData.M_CaCl2,0,0}/(1+n_Na*SaltData.M_NaCl+n_K*SaltData.M_KCl+n_Ca*SaltData.M_CaCl2)     "NaCl, KCl, CaCl2, MgCl2, SrCl2";

algorithm
//  Modelica.Utilities.Streams.print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  Modelica.Utilities.Streams.print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
end BrineProps1Phase;
