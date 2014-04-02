within BrineProp.Partial_Gas_Data;
partial function partial_fugacity_pTX
  input Modelica.SIunits.Pressure p "in Pa";
  input Modelica.SIunits.Temp_K T "in K";
  output Real phi "fugacity coefficient";
protected
    Partial_Units.Pressure_bar p_bar=Modelica.SIunits.Conversions.to_bar(p);
end partial_fugacity_pTX;
