within BrineProp.Partial_Gas_Data;
function degassingPressure_CH4_Duan2006
  "calculates degassing pressure from concentration of dissolved gas"
  extends partial_degassingPressure_pTX;
/*protected 
    Modelica.SIunits.MassFraction solu_soll=X[end-3];*/

algorithm
//  Modelica.Utilities.Streams.print("p_sat_CH4("+String(X[end-1])+") (degassingPressure_CH4_Duan2006)");
/*  p_gas := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
      function solubility_CH4_Duan2006_res(p=p,T=T,X=X,MM_vec=MM_vec),
      0,
      1e8,
      1e-8);
*/

  p_gas := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
      function solubility_res(
        solufun=function solubility_CH4_pTX_Duan2006(),p=p,T=T,X=X,MM_vec=MM_vec,
        c_gas=X[end-1]),
      0,
      1e10,
      1e-8);

//  Modelica.Utilities.Streams.print("p_sat_CH4("+String(X[end-1])+")="+String(p_gas)+" (degassingPressure_CH4_Duan2006)");

end degassingPressure_CH4_Duan2006;
