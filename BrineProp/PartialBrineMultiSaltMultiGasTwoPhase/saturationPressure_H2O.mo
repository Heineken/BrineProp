within BrineProp.PartialBrineMultiSaltMultiGasTwoPhase;
function saturationPressure_H2O "brine water vapour pressure"
  import BrineProp;
  extends Modelica.Icons.Function;
  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X[:] "mass fractions m_x/m_Sol";
  input SI.MolarMass MM[:] "molar masses of components";
  input Integer nM[:] "number of atoms per molecule of components";
  output SI.Pressure p_sat;
  output SI.Pressure p_H2O=0 "pure water vapour pressure";
protected
   SI.MoleFraction ionMoleFractions[nX];
algorithm
  if debugmode then
    print("Running saturationPressure_H2O("+String(p/1e5)+" bar,"+String(T-273.15)+" C, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
//    print("p_H2O("+String(T)+")="+String(p_sat/1e5)+" bar (PartialBrine_Multi_TwoPhase_ngas.saturationPressure_H2O)");
  end if;
//  printVector(nM);
  assert(max(X)-1<=1e-8, "X ="+String(max(X))+" out of range [0...1] = "+Modelica.Math.Matrices.toString(transpose([X])));
  assert(min(X)>=-1e-8, "X ="+String(min(X))+" out of range [0...1] = "+Modelica.Math.Matrices.toString(transpose([X])));
  if X[end]>0 then
//    ionMoleFractions:=BrineProp.massToMoleFractions(X, MM) .* nM;
    ionMoleFractions:=Utilities.massToMoleFractions(X, MM) .* nM;
    ionMoleFractions:=ionMoleFractions/sum(ionMoleFractions) "normalize";
    p_H2O:=Modelica.Media.Water.WaterIF97_pT.saturationPressure(T);
//    p_H2O:=Modelica.Media.Water.WaterIF97_pT.saturationPressure(if noEvent(T>300) then T else 300);
    p_sat:= p_H2O * ionMoleFractions[end];
  else
    p_sat:=10*p;
  end if;
// print("p_H2O="+String(p_H2O));
end saturationPressure_H2O;
