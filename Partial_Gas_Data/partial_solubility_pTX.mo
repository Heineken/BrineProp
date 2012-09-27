within BrineProp.Partial_Gas_Data;
partial function partial_solubility_pTX
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_Sol";
  input Modelica.SIunits.MolarMass MM_vec[:] "molar masses of components";
  input Modelica.SIunits.Pressure p_gas;
//  output Modelica.SIunits.MassFraction c_gas "gas concentration in kg_gas/kg_H2O";
  output Modelica.SIunits.MassFraction X_gas
    "gas concentration in kg_gas/kg_fluid";
protected
  Partial_Units.Molality solu "gas solubility";
//algorithm
//    Modelica.Utilities.Streams.print("mola("+String(X_gas)+","+String(T-273.16)+")=->k="+String(X_gas/max(1,p_gas))+" (partial_solubility_pTX)");
end partial_solubility_pTX;
