within BrineProp.Examples;
model PureWaterFlashing "gradual evaporation by pressure reduction"
//package Medium = Modelica.Media.Water.WaterIF97_ph;
//package Medium = REFPROPMedium(final substanceNames={"water"},  final explicitVars = "pT");
//package Medium = MediaTwoPhaseMixture.Water_MixtureTwoPhase_pT;
  package Medium = BrineProp.Brine5salts3gas "specify medium";

//  Modelica.SIunits.Pressure psat=Medium.saturationPressure(props.T);
//  Real x= Medium.vapourQuality(props.state);
  Real x;
  SI.Pressure p = 2e5-time*1e5;
//  Medium.BaseProperties props;
equation
    x=Medium.vapourQuality(Medium.setState_phX(p, 8.5e5,fill(0,Medium.nX)));
/*    props.p = 2e5-time*1e5;
    props.h = 8.5e5;
    props.Xi = fill(0,Medium.nXi);*/
  annotation (experiment(StopTime=1.95, __Dymola_NumberOfIntervals=20),
      __Dymola_experimentSetupOutput,
    Documentation(info="<html>
    <p>Pure water flashing only shows the transient gas fraction when setState_phX is called directly, due to the discontinuity in the p-T state space.</p>
    <p>The property model is temperature based on p-T functions, therefore cannot reproduce it.</p>
</html>"),
    __Dymola_Commands(file="Resources/Scripts/PureWaterFlashing.mos"
        "Run and plot"));
end PureWaterFlashing;
