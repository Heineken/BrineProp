within BrineProp.Partial_Gas_Data;
function degassingPressure_CO2_Duan2003
  "calculates degassing pressure from concentration of dissolved gas"
  extends partial_degassingPressure_pTX;
/*protected 
    Modelica.SIunits.MassFraction solu_soll=X[end-3];*/

algorithm
/*  p_gas := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
      function solubility_CO2_Duan2003_res(p=p,T=T,X=X,MM_vec=MM_vec),
      0,
      1e8,
      1e-8);
*/
  p_gas := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
      function solubility_res(
        solufun=function solubility_CO2_pTX_Duan2003(),p=p,T=T,X=X,MM_vec=MM_vec,
        c_gas=X[end-3]),
      0,
      1e10,
      1e-8);

//  Modelica.Utilities.Streams.print("p_sat_CO2("+String(X[end-3])+")="+String(p_gas)+" (degassingPressure_CO2_pTX)");

/*    while abs(solu-solu_soll)>1e-8 loop
    p_gas:=p_gas_neu;
    solu:=solubility_CO2_pTX_Duan2003(
        p=p,
        T=T,
        X=X,
        MM=MM_vec,
        p_gas=p_gas);
    z:=z + 1;
    assert(z<1000," Reached maximum number of iterations for solution equilibrium calculation. (quality_pTX)");
    Modelica.Utilities.Streams.print("p_gas="+String(p_gas_neu)+"->"+String(abs(p_gas-p_gas_neu)));
  end while;*/
end degassingPressure_CO2_Duan2003;
