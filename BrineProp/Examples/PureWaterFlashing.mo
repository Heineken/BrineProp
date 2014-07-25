within BrineProp.Examples;
model PureWaterFlashing "gradual evaporation by pressure reduction"
//package Medium = Modelica.Media.Water.WaterIF97_ph;
//package Medium = REFPROPMedium(final substanceNames={"water"},  final explicitVars = "pT");
//package Medium = MediaTwoPhaseMixture.Water_MixtureTwoPhase_pT;
  package Medium = BrineProp.Brine5salts3gas "specify medium";

  Medium.BaseProperties props;
//  Modelica.SIunits.Pressure psat=Medium.saturationPressure(props.T);
//  Real x= Medium.vapourQuality(props.state);
  Real x;
equation
  x=  Medium.vapourQuality(Medium.setState_phX( 2e5-time*1e5,8.5e5,fill(0,Medium.nX)));
    props.p = 2e5-time*1e5;
    props.h = 8.5e5;
    props.Xi = fill(0,Medium.nXi);
  annotation (experiment(StopTime=1.95, __Dymola_NumberOfIntervals=20),
      __Dymola_experimentSetupOutput,
    Documentation(info="<html>
<p>Pure water flashing only shows the correct gas fraction when setState_phX is called directly.</p>
<p>The property model is temperature based, therefore cannot reproduce it.</p>
</html>"));
end PureWaterFlashing;
