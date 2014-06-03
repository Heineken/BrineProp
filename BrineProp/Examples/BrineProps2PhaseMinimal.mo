within BrineProp.Examples;
model BrineProps2PhaseMinimal
  "Minimal example for 2-phase brine property model"
  package Medium = BrineProp.Brine_5salts_TwoPhase_3gas "specify medium";

  //DEFINE BRINE COMPOSITION (NaCl, KCl, CaCl2, MgCl2, SrCl2, CO2, N2, CH4)
  Real[Medium.nXi] Xi = {0.0839077010751,0.00253365118988,0.122786737978,0,0,7.2426359111e-05,0.000689505657647,6.14906384726e-05}
    "GrSk brine (Feldbusch 2-2013 1.1775g/ml V2)";

/*  package Medium = BrineProp.Water_MixtureTwoPhase_pT;
    Real[Medium.nXi] Xi= fill(0,0);*/

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
end BrineProps2PhaseMinimal;
