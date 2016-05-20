within BrineProp.Examples;
model BrineProps2PhaseMinimal_state
  "Using the state record instead of Baseproperties"
  //needs "Advanced.PedanticModelica:=false" to run

package Medium = BrineProp.Brine5salts3gas(ignoreLimitN2_T=true);
//  package Medium = Modelica.Media.Water.WaterIF97_ph;
  SI.Temperature T= 273.15+145;
  SI.Pressure p= (10+time)*1.01325e5;
  SI.MassFraction[Medium.nXi] Xi= {0.083945671051201,0.00253479771131107,0.122842299461699,0*0.000612116692496665,0*0.00214041137028575,  0.00016883,  0.00073459, 6.5652e-005}
    "Feldbusch 2-2013 1.1775g/ml";
  SI.MassFraction[:] X=cat(1,Xi,{1-sum(Xi)});
  Medium.ThermodynamicState state=Medium.setState_pTX(p,T,X);

  annotation (experiment(StopTime=100, __Dymola_NumberOfIntervals=100),
      __Dymola_experimentSetupOutput);
end BrineProps2PhaseMinimal_state;
