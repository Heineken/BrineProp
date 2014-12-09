within BrineProp.SpecificEnthalpies;
function appMolarEnthalpy_KCl_White
//2D-fit Reproduction of measurements of heat capacity of KCl solution
  extends PartialAppMolar_KCl_White;
  output Types.PartialMolarEnthalpy H_app_mol;
protected
  SI.Temperature T0=293.16 "Temperature at which HeatOfSolution is taken";
algorithm
  if AssertLevel>0 then
    assert(ignoreLimitSalt_T[iKCl] or (T_min<T and T<T_max),"\nTemperature T=" + String(T-273.15) + " C out of validity range ["+String(T_min-273.15)+"..."+String(T_max-273.15)+"]C.\nTo ignore set ignoreLimitSalt_T["+String(iKCl)+"]=true",aLevel);
  end if;

//  H_app_mol := HeatOfSolution_KCl_Sanahuja1986(T) + (a + b*bn + c*Tn/2 + d*bn^2 + e*bn*Tn/2 + f*Tn^2/3 + g*bn^2*Tn/2 + h*bn*Tn^2/3 + i*Tn^3/4)*T_std*Tn;
  H_app_mol:= HeatOfSolution_KCl_Sanahuja1986(T0) + (mola^b+c)*(k*(T-T0)+l*log((m-T)/(m-T0)));
end appMolarEnthalpy_KCl_White;
