within BrineProp.SpecificEnthalpies;
function appMolarHeatCapacity_CaCl2_White
//2D-fit Reproduction of measurements of heat capacity of KCl solution
  extends PartialAppMolar_CaCl2_White;
  output Partial_Units.PartialMolarHeatCapacity Cp_app_mol;
algorithm
//    Cp_app_mol:=a + b*bn + c*Tn + d*bn^2 + e*bn*Tn + f*Tn^2 + g*bn^2*Tn + h*bn* Tn^2 + i*Tn^3;
    Cp_app_mol:=(mola.^b+c).*(k-l*(m-T).^(-1));
end appMolarHeatCapacity_CaCl2_White;
