within BrineProp.GasData;
partial function partial_solubility_pTX
  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X[:] "mass fractions m_x/m_Sol";
  input SI.MolarMass MM_vec[:] "molar masses of components";
  input SI.Pressure p_gas;
  input Boolean ignoreTlimit "activated by temperature_phX"; //can't use default value here, because it breaks solveOneNonlinearEquation in degassingPressure_xxx
//  output SI.MassFraction c_gas "gas concentration in kg_gas/kg_H2O";
//  input PartialFlags flags;
  output SI.MassFraction X_gas "gas concentration in kg_gas/kg_fluid";
protected
  Types.Molality solu "gas solubility";
//algorithm
  //    print("mola("+String(X_gas)+","+String(T-273.16)+")=->k="+String(X_gas/max(1,p_gas))+" (partial_solubility_pTX)");
end partial_solubility_pTX;
