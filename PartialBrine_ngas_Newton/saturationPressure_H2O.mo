within BrineProp.PartialBrine_ngas_Newton;
function saturationPressure_H2O
  extends Modelica.Icons.Function;
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_Sol";
  input Modelica.SIunits.MolarMass MM[:] "molar masses of components";
  input Real nM[:] "molar masses of components";
  output Modelica.SIunits.Pressure p_sat;
  output Modelica.SIunits.Pressure p_H2O=0;
protected
   BrineProp.Partial_Units.Molality ionMoleFractions[nX];
algorithm
  if debugmode then
    Modelica.Utilities.Streams.print("Running saturationPressure_H2O("+String(p/1e5)+" bar,"+String(T-273.15)+" °C, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
//    Modelica.Utilities.Streams.print("p_H2O("+String(T)+")="+String(p_sat/1e5)+" bar (PartialBrine_Multi_TwoPhase_ngas.saturationPressure_H2O)");
  end if;
//  printVector(nM);
  assert(max(X)-1<=1e-8 and min(X)>=-1e-8, "X ="+String(max(X)-1)+" out of range [0...1] = "+Modelica.Math.Matrices.toString(transpose([X]))+" (saturationPressure_H2O())");
  if X[end]>0 then
    ionMoleFractions:=massFractionsToMoleFractions(X, MM).*nM;
    ionMoleFractions:=ionMoleFractions/sum(ionMoleFractions);
    p_H2O:=Modelica.Media.Water.WaterIF97_base.saturationPressure(T);
    p_sat:= p_H2O * ionMoleFractions[end];
  else
    p_sat:=10*p;
  end if;
// Modelica.Utilities.Streams.print("p_H2O="+String(p_H2O));
end saturationPressure_H2O;
