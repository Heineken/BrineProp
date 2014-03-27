within BrineProp.PartialGasData;
partial function partial_fugacity_pTX
  input SI.Pressure p "in Pa";
  input SI.Temp_K T "in K";
  output Real phi "fugacity coefficient";
protected
    Partial_Units.Pressure_bar p_bar=SI.Conversions.to_bar(p);
end partial_fugacity_pTX;
