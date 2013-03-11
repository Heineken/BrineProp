within BrineProp.SpecificEnthalpies;
function specificHeatCapacity_pTX_Driesner
  "cp calculation according to Driesner 2007 et al: 0-1000°C; 0.1-500MPa (doi:10.1016/j.gca.2007.05.026)"
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input Modelica.SIunits.MassFraction X_NaCl "mass fraction m_NaCl/m_Sol";
  output Modelica.SIunits.SpecificHeatCapacity cp;
protected
  Modelica.SIunits.Temp_K T_Scale_h;
  Real q2;
algorithm
  (T_Scale_h,q2):=T_Scale_h_Driesner(p,T,X_NaCl);
  cp := Modelica.Media.Water.IF97_Utilities.cp_pT(p,T_Scale_h)*q2;
end specificHeatCapacity_pTX_Driesner;
