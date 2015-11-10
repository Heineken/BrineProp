within BrineProp;
package PartialFlags
    constant Boolean debugmode = false "print messages in functions";
//    replaceable type Level=enumeration(none, warning, error)"Substances in Fluid";
    constant Integer AssertLevel= 2
    "when out of validity range: 0-do nothing, 1-show warnings, 2-throw error";
    constant AssertionLevel aLevel = if AssertLevel==1 then AssertionLevel.warning else AssertionLevel.error;

    //SALT stuff
    constant Boolean[:] ignoreLimitSalt_p=fill(false,nX_salt_)
    "ignore pressure limits {NaCl, KCl, CaCl2, MgCl2, SrCl2}";
    constant Boolean[:] ignoreLimitSalt_T=fill(false,nX_salt_)
    "ignore temperature limits         {NaCl, KCl, CaCl2, MgCl2, SrCl2}";
    constant Boolean[:] ignoreLimitSalt_b=fill(false,nX_salt_)
    "ignore salinity limits            {NaCl, KCl, CaCl2, MgCl2, SrCl2}";

    constant Boolean[:] ignoreLimitSalt_soluCO2=fill(false,nX_salt_)
    "ignore salinity limits in CO2 solubility {NaCl, KCl, CaCl2, MgCl2, SrCl2}";
    constant Boolean[:] ignoreLimitSalt_soluN2=fill(false,nX_salt_)
    "ignore salinity limits in N2 solubility {NaCl, KCl, CaCl2, MgCl2, SrCl2}";
    constant Boolean[:] ignoreLimitSalt_soluCH4=fill(false,nX_salt_)
    "ignore salinity limits in CH4 solubility {NaCl, KCl, CaCl2, MgCl2, SrCl2}";

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

protected
 constant Integer nX_salt_;
end PartialFlags;
