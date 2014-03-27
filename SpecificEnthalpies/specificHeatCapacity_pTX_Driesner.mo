within BrineProp.SpecificEnthalpies;
function specificHeatCapacity_pTX_Driesner
  "cp calculation according to Driesner 2007 et al: 0-1000°C; 0.1-500MPa (doi:10.1016/j.gca.2007.05.026)"
  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X_NaCl "mass fraction m_NaCl/m_Sol";
  output SI.SpecificHeatCapacity cp;
protected
  SI.Temp_K T_Scale_h;
  Real q_2;
algorithm
  (T_Scale_h,q_2):=T_Scale_h_Driesner(p,T,X_NaCl);
  cp := Modelica.Media.Water.IF97_Utilities.cp_pT(p,T_Scale_h)*q_2;
//  print("T_Scale_h "+String(T_Scale_h_Driesner(p,T,X_NaCl))+" K. q_2="+String(q_2)+" ->"+String(cp)+" J/kKg");

end specificHeatCapacity_pTX_Driesner;
