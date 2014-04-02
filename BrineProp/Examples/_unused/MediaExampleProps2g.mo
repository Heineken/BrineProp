within BrineProp.Examples._unused;
model MediaExampleProps2g
//package Medium = Brine_Duan_Multi_TwoPhase_2gas;
//package Medium2 = Brine_Duan_Multi_TwoPhase_2gas(final explicitVars="ph");
package Medium = Media.Brine.Brine_Duan_Multi_TwoPhaseSep_2gas (final
        explicitVars = "pT");
  Medium.BaseProperties props;

 package Medium2 = Media.Brine.Brine_Duan_Multi_TwoPhaseSep_2gas (final
        explicitVars = "ph");
  Medium2.BaseProperties props2; /**/
equation
  props.p = 235e5;
  props.T = 130+273.15;
/*  props.T_gas = 293.15;*/
//  props.h = 513166;
//  props.h_phases = {3.21918e6, 8.41059e4};
 props.Xi = {0,0,0,0,0,1e-4,0}
    "NaCl, c_KCl, c_CaCl2, c_MgCl2, c_SrCl2, CO2, N2";

/*  props2.Xi = props.Xi "NaCl, c_KCl, c_CaCl2, c_MgCl2, c_SrCl2, N2, CO2";
  props2.p = props.p;
  props2.h = props.h;
*/
//  props.Xi = {0,0,0,0,0,3.076e-5} "NaCl, c_KCl, c_CaCl2, c_MgCl2, c_SrCl2, CO2/N2";

  props2.p = 215e5;
//  props2.h_phases*{props.q,1-props.q} = -(props.h_phases-Modelica.Constants.g_n*200*{1,1})*{props2.q,1-props2.q};
  props2.h = props.h;
  props2.Xi = props.Xi;
  props2.T = props2.T_gas; /**/

end MediaExampleProps2g;
