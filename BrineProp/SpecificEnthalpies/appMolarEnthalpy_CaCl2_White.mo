within BrineProp.SpecificEnthalpies;
function appMolarEnthalpy_CaCl2_White
//2D-fit Reproduction of measurements of heat capacity of CaCl2 solution
  extends PartialAppMolar_CaCl2_White;
  output Types.PartialMolarEnthalpy H_app_mol;
protected
  constant SI.MolarInternalEnergy Delta_h_solution_CaCl2 = 81850
    "[J/mol_CaCl2] @ 298.15K Sinke1985 http://dx.doi.org/10.1016/0021-9614(85)90083-7";
  SI.Temperature T0=298.15 "Temperature at which HeatOfSolution is taken";
algorithm
  if AssertLevel>0 then
    assert(ignoreLimitSalt_T[iCaCl2] or (T_min<T and T<T_max),"\nTemperature T=" + String(T-273.15) + " C out of validity range ["+String(T_min-273.15)+"..."+String(T_max-273.15)+"]C.\nTo ignore set ignoreLimitSalt_T["+String(iCaCl2)+"]=true\n",aLevel);
  end if;
//  H_app_mol := Delta_h_solution_CaCl2 + (a + b*bn + c*Tn/2 + d*bn^2 + e*bn*Tn/2 + f*Tn^2/3 + g*bn^2*Tn/2 + h*bn*Tn^2/3 + i*Tn^3/4)*T_std*Tn;
//  H_app_mol:= Delta_h_solution_CaCl2 + (mola.^b+c).*(k*T+l*log(m-T));
  H_app_mol:= Delta_h_solution_CaCl2 + (mola^b+c)*(k*(T-T0)+l*log((m-T)/(m-T0)));

end appMolarEnthalpy_CaCl2_White;
