within BrineProp.SpecificEnthalpies;
function appMolarEnthalpy_KCl_White
//2D-fit Reproduction of measurements of heat capacity of KCl solution
extends PartialAppMolar_KCl_White;
  output Partial_Units.PartialMolarEnthalpy H_app_mol;
protected
  Modelica.SIunits.Temperature T0=293.16
    "Temperature at which HeatOfSolution is taken";
algorithm
//  H_app_mol := HeatOfSolution_KCl_Sanahuja1986(T) + (a + b*bn + c*Tn/2 + d*bn^2 + e*bn*Tn/2 + f*Tn^2/3 + g*bn^2*Tn/2 + h*bn*Tn^2/3 + i*Tn^3/4)*T_std*Tn;
  H_app_mol:= HeatOfSolution_KCl_Sanahuja1986(T0) + (mola.^b+c).*(k*(T-T0)+l*log((m-T)/(m-T0)));
end appMolarEnthalpy_KCl_White;
