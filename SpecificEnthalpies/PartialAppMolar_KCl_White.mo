within BrineProp.SpecificEnthalpies;
partial function PartialAppMolar_KCl_White
//2D-fit Reproduction of measurements of heat capacity of KCl solution
//  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input BrineProp.Partial_Units.Molality mola "n_KCl/m_H2O";
//  output Modelica.SIunits.SpecificHeatCapacity cp=1 "=cp_by_cpWater*cp_Water";
  //Parameters of MATLAB 2D-Fit
protected
  Real b =   0.09818;
  Real c =    -1.244;
  Real k =    -327.9;
  Real l = -1.31e+05;
  Real m =     628.8;

/*  Real a= -10.63;
  Real b= 57.96;
  Real c= 1.789;
  Real d= -56.53;
  Real e= 168.9;
  Real f= -274.6;
  Real g= -57.32;
  Real h= 102;
  Real i= -179.4;

  BrineProp.Partial_Units.Molality b_mean=1.188;
  BrineProp.Partial_Units.Molality b_std=1.103;
  Modelica.SIunits.Temp_K T_mean=475.1;
  Modelica.SIunits.Temp_K T_std=103.5;

  Real bn= (mola-b_mean)/b_std "normalized & centered";
  Real Tn= (T-T_mean)/T_std "normalized & centered";*/
//  Modelica.SIunits.SpecificHeatCapacity cp_Water =  Modelica.Media.Water.IF97_Utilities.cp_pT(p, T);
end PartialAppMolar_KCl_White;
