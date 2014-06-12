within BrineProp.PartialBrine_MultiSalt_2Phase_MultiGas;
function saturationPressure_H2O "brine water vapour pressure"
  extends Modelica.Icons.Function;
  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X[:] "mass fractions m_x/m_Sol";
  input SI.MolarMass MM[:] "molar masses of components";
  input Real nM[:] "number of atoms per molecule of components";
  output SI.Pressure p_sat;
  output SI.Pressure p_H2O=0 "pure water vapour pressure";
protected
   BrineProp.Partial_Units.Molality ionMoleFractions[nX];
algorithm
  if debugmode then
    print("Running saturationPressure_H2O("+String(p/1e5)+" bar,"+String(T-273.15)+" °C, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
//    print("p_H2O("+String(T)+")="+String(p_sat/1e5)+" bar (PartialBrine_Multi_TwoPhase_ngas.saturationPressure_H2O)");
  end if;
//  printVector(nM);
  assert(max(X)-1<=1e-8, "X ="+String(max(X))+" out of range [0...1] = "+Modelica.Math.Matrices.toString(transpose([X]))+" (saturationPressure_H2O())");
  assert(min(X)>=-1e-8, "X ="+String(min(X))+" out of range [0...1] = "+Modelica.Math.Matrices.toString(transpose([X]))+" (saturationPressure_H2O())");
  if X[end]>0 then
    ionMoleFractions:=massFractionsToMoleFractions(X, MM).*nM;
    ionMoleFractions:=ionMoleFractions/sum(ionMoleFractions) "normalize";
    p_H2O:=Modelica.Media.Water.WaterIF97_base.saturationPressure(T);
    p_sat:= p_H2O * ionMoleFractions[end];
  else
    p_sat:=10*p;
  end if;
// print("p_H2O="+String(p_H2O));
end saturationPressure_H2O;
