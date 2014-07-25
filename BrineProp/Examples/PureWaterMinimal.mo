within BrineProp.Examples;
model PureWaterMinimal "Minimal example for pure water"
//specify medium
  /*  package Medium = BrineProp.Brine_5salts_TwoPhase_3gas;
  Real[Medium.nXi] Xi = fill(0,Medium.nXi) "Brine with zero salinity";*/

  package Medium = BrineProp.WaterMixtureTwoPhase_pT "IAPWS water in wrapper";
    Real[Medium.nXi] Xi= fill(0,0);

  Medium.BaseProperties props;
equation
  //specify thermodynamic state
  props.p = 8e5;
  props.h = 8.5e5;

/*  props.p = 1.01325e5;
  props.T = 20+273.15;
  */

  //specify brine composition
  props.Xi = Xi;
end PureWaterMinimal;
