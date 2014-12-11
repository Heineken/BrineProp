within BrineProp.SpecificEnthalpies;
function appMolarHeatCapacity_CaCl2_White
//2D-fit Reproduction of measurements of heat capacity of KCl solution
  extends PartialAppMolar_CaCl2_White;
  output Types.PartialMolarHeatCapacity Cp_app_mol;
algorithm
  if AssertLevel>0 then
    assert(ignoreLimitSalt_T[iCaCl2] or (T_min<T and T<T_max),"Temperature T=" + String(T-273.15) + " C out of validity range ["+String(T_min-273.15)+"..."+String(T_max-273.15)+"]C.\nTo ignore set ignoreLimitSalt_T["+String(iCaCl2)+"]=true",aLevel);
  end if;

  Cp_app_mol:=(mola^b+c)*(k-l*(m-T)^(-1));
//   print("Cp_app_mol_CaCl2= "+String(Cp_app_mol)+"J/kg");
end appMolarHeatCapacity_CaCl2_White;
