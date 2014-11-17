within BrineProp.Vacuum;
function extends specificHeatCapacityCp
  "water-saturated heat capacity of gas phase"
algorithm
     cp := specificHeatCapacityCp_pTX(
        p=state.p,
        T=state.T,
        X= state.X);

end specificHeatCapacityCp;
