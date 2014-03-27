within BrineProp.PartialGasData;
partial function partial_solutionEnthalpy
  "template for calculation of solution enthalpy"
  input SI.Temp_K T;
  output SI.SpecificEnthalpy Delta_h_solution;
protected
  SI.Pressure p_H2O;

end partial_solutionEnthalpy;
