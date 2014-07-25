within BrineProp.GasData;
partial function partial_degassingPressure_pTX
  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X[:] "mass fractions m_x/m_Sol";
  input SI.MolarMass MM_vec[:] "molar masses of components";
  output SI.Pressure p_gas;
end partial_degassingPressure_pTX;
