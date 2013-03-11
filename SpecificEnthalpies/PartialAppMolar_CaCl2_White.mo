within BrineProp.SpecificEnthalpies;
partial function PartialAppMolar_CaCl2_White
//2D-fit Reproduction of measurements of heat capacity of KCl solution
//  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input BrineProp.Partial_Units.Molality mola "n_KCl/m_H2O";
//  output Modelica.SIunits.SpecificHeatCapacity cp=1 "=cp_by_cpWater*cp_Water";
  //Parameters of MATLAB 2D-Fit
protected
  Real b =   -0.001977;
  Real c =     -0.9958;
  Real k =        1373;
  Real l =   6.736e+06;
  Real m =         628;

/*
  Real a= -50.08;
  Real b= 145.9;
  Real c= 146.3;
  Real d= -89.3;
  Real e= 227.8;
  Real f= -210.7;
  Real g= -101.9;
  Real h= 116.6;
  Real i= -213.8;

  BrineProp.Partial_Units.Molality b_mean=0.9803;
  BrineProp.Partial_Units.Molality b_std=1.047;
  Modelica.SIunits.Temp_K T_mean=437.7;
  Modelica.SIunits.Temp_K T_std=104.5;

  Real bn= (mola-b_mean)/b_std "normalized & centered";
  Real Tn= (T-T_mean)/T_std "normalized & centered";*/
//  Modelica.SIunits.SpecificHeatCapacity cp_Water =  Modelica.Media.Water.IF97_Utilities.cp_pT(p, T);
end PartialAppMolar_CaCl2_White;
