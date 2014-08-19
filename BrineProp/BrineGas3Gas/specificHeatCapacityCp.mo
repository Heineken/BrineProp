within BrineProp.BrineGas3Gas;
function extends specificHeatCapacityCp
  "water-saturated heat capacity of gas phase"
algorithm
     cp := specificHeatCapacityCp_pTX(
        p=state.p,
        T=state.T,
        X= if waterSaturated then
     waterSaturatedComposition_pTX(state.p,state.T,state.X)
  else state.X);
//  else state.X[end - nX + 1:end]);

end specificHeatCapacityCp;
