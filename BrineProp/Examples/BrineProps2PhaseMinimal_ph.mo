within BrineProp.Examples;
model BrineProps2PhaseMinimal_ph
  "Minimal inversion (h(p,T) -> T(p,h)) example for 2-phase brine property model"
  //needs "Advanced.PedanticModelica:=false" to run

//SPECIFY MEDIUM and COMPOSITION

//   package Medium = BrineProp.Brine5salts3gas(ignoreLimitN2_T=true);
   package Medium = BrineProp.Brine5salts3gas;

//DEFINE BRINE COMPOSITION (NaCl, KCl, CaCl2, MgCl2, SrCl2, CO2, N2, CH4)
  Real[Medium.nXi] Xi = {0.0839077010751,0.00253365118988,0.122786737978,0,0,7.2426359111e-05,0.000689505657647,6.14906384726e-05}
    "GrSk brine (Feldbusch 2-2013 1.1775g/ml V2)";

/*  package Medium = BrineProp.Water_MixtureTwoPhase_pT;
    Real[Medium.nXi] Xi= fill(0,0);*/

  Medium.BaseProperties props;
  Real GVF = props.GVF "<- PLOT ME!";
equation
  //SPECIFY THERMODYNAMIC STATE
/*  //degassing by heating starting at STP
  props.p = 1.01325e5;
  props.T = 20+273.15+time;
*/
  //degassing by decompression starting at reservoir conditions
  props.p = 455*1e5;
  props.h = 541025;

  //specify brine composition
  props.Xi = Xi;
  annotation (experiment(StopTime=80, __Dymola_NumberOfIntervals=100),
      __Dymola_experimentSetupOutput);
end BrineProps2PhaseMinimal_ph;
