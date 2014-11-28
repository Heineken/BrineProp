within BrineProp;
partial package PartialFlags
  "Flags to deactivate error messages or p/T-range limits"
  //extended in packages, flags can be set on Medium instantiation
//  record Flags
//    extends Modelica.Icons.Record;
    constant Boolean debugmode = false "print messages in functions";
    constant Integer outOfRangeMode=2
    "when out of validity range: 0-do nothing, 1-show warnings, 2-throw error";
    constant Boolean ignoreLimitN2_T=false;
    constant Boolean ignoreLimitN2_p=false;
    // constant Boolean[5] ignoreLimitSalt_visc={false,false,false,false,false};
    constant Boolean[5] ignoreLimitSalt_p={false,false,false,false,false}
    "ignore pressure limits            {NaCl, KCl, CaCl2, MgCl2, SrCl2}";
    constant Boolean[5] ignoreLimitSalt_T={false,false,false,false,false}
    "ignore temperature limits         {NaCl, KCl, CaCl2, MgCl2, SrCl2}";
    constant Boolean[5] ignoreLimitSalt_b={false,false,false,false,false}
    "ignore salinity limits            {NaCl, KCl, CaCl2, MgCl2, SrCl2}";
    annotation (Documentation(info="<html></html>"));
//  end Flags;
    constant Boolean ignoreLimitInh_KCl_Tmin=false
    "ignore Tmin in appMolarEnthalpy_KCl_White and appMolarHeatCapacity_KCl_White";
    constant Boolean ignoreLimitInh_CaCl2_Tmin=false
    "ignore Tmin in appMolarEnthalpy_CaCl2_White and appMolarHeatCapacity_CaCl2_White";

//  constant Flags flags;
end PartialFlags;
