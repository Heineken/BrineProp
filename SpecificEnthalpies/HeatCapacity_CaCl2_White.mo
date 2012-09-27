within BrineProp.SpecificEnthalpies;
function HeatCapacity_CaCl2_White
//2D-fit Reproduction of measurements of heat capacity of KCl solution
//  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input BrineProp.Partial_Units.Molality b "n_KCl/m_H2O";
//  output Modelica.SIunits.SpecificHeatCapacity cp=1 "=cp_by_cpWater*cp_Water";
//  output Real cp_by_cpWater;
  output Real int_cp_by_cpWater;
  //Parameters of MATLAB 2D-Fit
protected
  Real P00= 0.7399;
  Real P10= 0.001297;
  Real P01=-0.1279;
  Real P20=-1.564E-06;
  Real P11=-0.00006349;
  Real P02= 0.0202;

//  Modelica.SIunits.SpecificHeatCapacity cp_Water =  Modelica.Media.Water.IF97_Utilities.cp_pT(p, T);
algorithm
//cp_by_cpWater :=p00 + p10*T + p01*b + p20*T^2 + p11*T*b + p02*b^2;
int_cp_by_cpWater:=(P00 + P01*b + P10/2*T + P20/3*T^2 + P11/2*T*b + P02*b^2)*T
    "cp ratio integrated over T";

//cp = cp_by_cpWater*cp_Water;
//Modelica.Utilities.Streams.print("Brine.specificEnthalpy_pTX_Francke: "+String(p*1e-5)+"bar."+String(T)+"°C->"+String(h)+" J/kg");
end HeatCapacity_CaCl2_White;
