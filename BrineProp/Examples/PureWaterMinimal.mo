within BrineProp.Examples;
model PureWaterMinimal "Minimal example for pure water"
/*  package Medium = BrineProp.Brine_5salts_TwoPhase_3gas "specify medium";
  Real[Medium.nXi] Xi = fill(0,Medium.nXi)
    "GrSk brine (Feldbusch 2-2013 1.1775g/ml V2)";*/

  package Medium = BrineProp.Water_MixtureTwoPhase_pT;
    Real[Medium.nXi] Xi= fill(0,0);

  Medium.BaseProperties props;
equation
  //specify thermodynamic state
//  props.p = 1.01325e5;
  props.p = 8e5;
//  props.T = 220+273.15-time;
//  props.T = 20+273.15;
  props.h = 8.5e5;
  //specify brine composition
  props.Xi = Xi;
end PureWaterMinimal;
