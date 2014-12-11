within BrineProp.BrineGas3Gas;
function extends specificHeatCapacityCp_pTX
  "calculation of specific heat capacities of gas mixture"
  import Modelica.Media.IdealGases.Common.SingleGasNasa;
  import Modelica.Media.IdealGases.SingleGases;

  import Modelica.Media.Water;

protected
    SingleGases.H2O.ThermodynamicState state=SingleGases.H2O.ThermodynamicState(p=0,T=T);
    SI.SpecificHeatCapacity cp_CO2=SingleGases.CO2.specificHeatCapacityCp(state);
    SI.SpecificHeatCapacity cp_N2=SingleGases.N2.specificHeatCapacityCp(state);
    SI.SpecificHeatCapacity cp_CH4=SingleGases.CH4.specificHeatCapacityCp(state);
    SI.SpecificHeatCapacity cp_H2O=Water.IF97_Utilities.cp_pT(min(p,Water.IF97_Utilities.BaseIF97.Basic.psat(T)-1),T=T)
    "below psat -> gaseous";

    SI.SpecificHeatCapacity cp_vec[:]={cp_CO2,cp_N2,cp_CH4,cp_H2O};

    /*  SI.SpecificHeatCapacity cp_vec[:]={
    SingleGasNasa.cp_T(data=SingleGasesData.CO2,T=T),
    SingleGasNasa.cp_T(data=SingleGasesData.N2,T=T),
    SingleGasNasa.cp_T(data=SingleGasesData.CH4,T=T),
   Water.IF97_Utilities.cp_pT(min(p,Water.IF97_Utilities.BaseIF97.Basic.psat(T)-1),T=T)} 
*/
algorithm
  if debugmode then
    print("Running specificHeatCapacityCp_pTX("+String(p/1e5)+" bar,"+String(T-273.15)+" degC, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
  end if;

  if not ignoreNoCompositionInBrineGas and not min(X)>0 then
    print("No gas composition, assuming water vapour.(BrineProp.BrineGas_3Gas.specificHeatCapacityCp_pTX)");
  end if;

/*  if waterSaturated then
    cp := cp_vec * waterSaturatedComposition_pTX(p,T,X[end - nX+1:end]);
  else */
//    cp := cp_vec * X[end - nX+1:end];
  cp := cp_vec * cat(1,X[1:end-1],{if min(X)>0 then X[end] else 1});
    //  end if;

/*  print("cp_CO2: "+String(cp_vec[1])+" J/kg");
  print("cp_N2: "+String(cp_vec[2])+" J/kg");
  print("cp_CH4: "+String(cp_vec[3])+" J/kg");
  print("cp_H2O: "+String(cp_vec[4])+" J/kg"); */

end specificHeatCapacityCp_pTX;
