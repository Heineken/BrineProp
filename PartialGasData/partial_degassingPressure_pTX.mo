within BrineProp.PartialGasData;
partial function partial_degassingPressure_pTX
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_Sol";
  input Modelica.SIunits.MolarMass MM_vec[:] "molar masses of components";
  output Modelica.SIunits.Pressure p_gas;
end partial_degassingPressure_pTX;
