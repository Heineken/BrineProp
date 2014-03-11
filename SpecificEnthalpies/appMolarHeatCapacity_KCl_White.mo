within BrineProp.SpecificEnthalpies;
function appMolarHeatCapacity_KCl_White
//2D-fit Reproduction of measurements of heat capacity of KCl solution
extends PartialAppMolar_KCl_White;
  output Partial_Units.PartialMolarHeatCapacity Cp_app_mol;
algorithm
//  Cp_app_mol:=a + b*bn + c*Tn + d*bn^2 + e*bn*Tn + f*Tn^2 + g*bn^2*Tn + h*bn*Tn^2 + i*Tn^3;
    Cp_app_mol:=(mola^b+c)*(k-l*(m-T)^(-1));
//   Modelica.Utilities.Streams.print("Cp_app_mol_KCl= "+String(Cp_app_mol)+"J/kg");
end appMolarHeatCapacity_KCl_White;
