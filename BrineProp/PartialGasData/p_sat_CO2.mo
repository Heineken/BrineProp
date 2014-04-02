within BrineProp.PartialGasData;
function p_sat_CO2
  "calculates saturation pressure, polynom derived from EES calculations"
  input SI.Temp_K T;
  output SI.Pressure p_sat;
algorithm
  assert(T<305,"Temperature above critical Temperature (304.1 K)");
  p_sat :=1178.4*T^2 - 555378*T + 7E+07;
end p_sat_CO2;
