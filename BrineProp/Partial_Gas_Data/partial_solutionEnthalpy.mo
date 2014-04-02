within BrineProp.Partial_Gas_Data;
partial function partial_solutionEnthalpy
  "template for calculation of solution enthalpy"
  input Modelica.SIunits.Temp_K T;
  output Modelica.SIunits.SpecificEnthalpy Delta_h_solution;
protected
  Modelica.SIunits.Pressure p_H2O;

end partial_solutionEnthalpy;
