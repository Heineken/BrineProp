within BrineProp.Vacuum;
function extends specificHeatCapacityCp_pTX
  "calculation of specific heat capacities of gas mixture"
algorithm
  if debugmode then
    print("Running specificHeatCapacityCp_pTX("+String(p/1e5)+" bar,"+String(T-273.15)+" degC, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
  end if;

  cp := 0;

end specificHeatCapacityCp_pTX;
