within BrineProp;
partial package PartialFlags
  "Flags to deactivate error messages or p/T-range limits"
  //extended in packages, flags can be set on Medium instantiation
//  record Flags
//    extends Modelica.Icons.Record;
    constant Boolean debugmode = false "print messages in functions";
//    replaceable type Level=enumeration(none, warning, error)"Substances in Fluid";
    constant Integer AssertLevel= 2
    "when out of validity range: 0-do nothing, 1-show warnings, 2-throw error";
    constant AssertionLevel aLevel = if AssertLevel==1 then AssertionLevel.warning else AssertionLevel.error;

    //SALT stuff
    // constant Boolean[5] ignoreLimitSalt_visc={false,false,false,false,false};
    constant Boolean[5] ignoreLimitSalt_p={false,false,false,false,false}
    "ignore pressure limits            {NaCl, KCl, CaCl2, MgCl2, SrCl2}";
    constant Boolean[5] ignoreLimitSalt_T={false,false,false,false,false}
    "ignore temperature limits         {NaCl, KCl, CaCl2, MgCl2, SrCl2}";
    constant Boolean[5] ignoreLimitSalt_b={false,false,false,false,false}
    "ignore salinity limits            {NaCl, KCl, CaCl2, MgCl2, SrCl2}";

/*    constant Integer iNaCl "order in vectors";
    constant Integer iKCl=0 "order in vectors";
    constant Integer iCaCl2=0 "order in vectors";
    constant Integer iMgCl2=0 "order in vectors";
    constant Integer iSrCl2=0 "order in vectors";
*/
    //GAS stuff
    constant Boolean ignoreLimitN2_T=false;
    constant Boolean ignoreLimitN2_p=false;
    constant Boolean ignoreLimitN2_b=false;
    constant Boolean ignoreLimitCO2_T=false;
    constant Boolean ignoreLimitCO2_p=false;
    constant Boolean ignoreLimitCH4_T=false;
    constant Boolean ignoreLimitCH4_p=false;
    constant Boolean ignoreLimitCH4_b=false;
    constant Boolean ignoreNoCompositionInBrineGas=false;

//  end Flags;
//    constant Boolean ignoreLimitInh_KCl_Tmin=false "ignore Tmin in appMolarEnthalpy_KCl_White and appMolarHeatCapacity_KCl_White";
//    constant Boolean ignoreLimitInh_CaCl2_Tmin=false "ignore Tmin in appMolarEnthalpy_CaCl2_White and appMolarHeatCapacity_CaCl2_White";

//  constant Flags flags;
end PartialFlags;
