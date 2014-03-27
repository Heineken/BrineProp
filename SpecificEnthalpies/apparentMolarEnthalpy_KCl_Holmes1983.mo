within BrineProp.SpecificEnthalpies;
function apparentMolarEnthalpy_KCl_Holmes1983
  "enthalpy calculation according to H. F. Holmes, R. E. Mesmer 1983: 0-250°C; <6mol/kg (doi:10.1021/j100230a030)"
//  input SI.Pressure p;
  input SI.Temp_K T;
//  input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
  output SI.SpecificEnthalpy h_app "apparent molar enthalpy";
protected
  Real[ 5] q={1.5774e5, -1037.86, 2.7739, -0.00284332, -686};
  String msg = "";
algorithm
  if outOfRangeMode>0 then
    if not (T>=273.15 and T<=250+273.15) then
      msg :="Temperature is " + String(T_C) + "°C, but must be between " +
        String(T_min) + "°C and " + String(T_max) + "°C (BrineProp.SpecificEnthalpies.T_Scale_h_Driesner)";
    end if;
    if msg<>"" then
      if outOfRangeMode==1 then
      print(msg);
      elseif outOfRangeMode==2 then
       assert(true, msg);
      end if;
    end if;
  end if;

  h_app := q[1] + q[2]*T + q[3]*T^2 + q[4]*T^3 + q[5] * ln(T-270);

//  print("Brine_Driesner.specificEnthalpy_pTX: "+String(p*1e-5)+"bar."+String(T_Scale_h)+"°C->"+String(h)+" J/kg");
end apparentMolarEnthalpy_KCl_Holmes1983;
