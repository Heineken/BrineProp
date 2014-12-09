within BrineProp.SpecificEnthalpies;
function apparentMolarEnthalpy_KCl_Holmes1983
  "enthalpy calculation according to H. F. Holmes, R. E. Mesmer 1983: 0-250degC; <6mol/kg (doi:10.1021/j100230a030)"
//  input SI.Pressure p;
  input SI.Temp_K T;
//  input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
  output SI.SpecificEnthalpy h_app "apparent molar enthalpy";
protected
  Real[ 5] q={1.5774e5, -1037.86, 2.7739, -0.00284332, -686};
  SI.Temp_C T_C = SI.Conversions.to_degC(T);
  constant SI.Temp_C T_min=273.15;
  constant SI.Temp_C T_max=250+273.15;
algorithm
  if AssertLevel>0 then
    assert(ignoreLimitSalt_T[3] or (T_C>=T_min and T_C<=T_max),"\nTemperature is out of validity range: T=" + String(T) + " mol/kg.\nTo ignore set ignoreLimitSalt_T[3]=true",aLevel);
  end if;

  h_app := q[1] + q[2]*T + q[3]*T^2 + q[4]*T^3 + q[5] * log(T-270);

//  print("Brine_Driesner.specificEnthalpy_pTX: "+String(p*1e-5)+"bar."+String(T_Scale_h)+"degC->"+String(h)+" J/kg");
end apparentMolarEnthalpy_KCl_Holmes1983;
