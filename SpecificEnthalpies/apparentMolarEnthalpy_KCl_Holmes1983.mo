within BrineProp.SpecificEnthalpies;
function apparentMolarEnthalpy_KCl_Holmes1983
  "enthalpy calculation according to H. F. Holmes, R. E. Mesmer 1983: 0-250°C; 0.1-MPa (doi:10.1021/j100230a030)"
//  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
//  input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
  output Modelica.SIunits.SpecificEnthalpy h_app "apparent molar enthalpy";
protected
  Real[ 5] q={1.5774e5, -1037.86, 2.7739, -0.00284332, -686};
algorithm
  h_app := q[1] + q[2]*T + q[3]*T^2 + q[4]*T^3 + q[5] * ln(T-270);

//  Modelica.Utilities.Streams.print("Brine_Driesner.specificEnthalpy_pTX: "+String(p*1e-5)+"bar."+String(T_Scale_h)+"°C->"+String(h)+" J/kg");
end apparentMolarEnthalpy_KCl_Holmes1983;
