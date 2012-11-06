within BrineProp.SpecificEnthalpies;
function specificEnthalpy_pTX_Driesner
  "enthalpy calculation according to Driesner 2007 et al: 0-1000°C; 0.1-500MPa (doi:10.1016/j.gca.2007.05.026)"
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input Modelica.SIunits.MassFraction X_NaCl "mass fraction m_NaCl/m_Sol";
  output Modelica.SIunits.SpecificEnthalpy h;/**/
algorithm

//  h := Modelica.Media.Water.WaterIF97_base.specificEnthalpy_pT(p, Modelica.SIunits.Conversions.from_degC(T_Scale_h));
  h := Modelica.Media.Water.WaterIF97_base.specificEnthalpy_pT(p, T_Scale_h_Driesner(p,T,X_NaCl));

end specificEnthalpy_pTX_Driesner;
