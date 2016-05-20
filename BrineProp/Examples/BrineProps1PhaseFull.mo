within BrineProp.Examples;
model BrineProps1PhaseFull
  "Example for 1-phase brine property model demonstrating all features"
  //needs "Advanced.PedanticModelica:=false" to run

  //SPECIFY MEDIUM
  package Medium = Brine5salts(AssertLevel=2,ignoreLimitSalt_T={false,true,true,false,false});
//package Medium = Modelica.Media.Water.WaterIF97_pT;

  Medium.BaseProperties props;

  SI.Density d= props.d "density";
  SI.SpecificEnthalpy h= props.h "specific enthalpy";
  SI.DynamicViscosity eta=Medium.dynamicViscosity(props.state)
    "dynamic viscosity";
  SI.ThermalConductivity lambda=Medium.thermalConductivity(props.state)
    "thermal conductivity";
  SI.SpecificHeatCapacity c_p_brine= Medium.specificHeatCapacityCp(props.state)
    "specific heat capacity";
  SI.TemperatureDifference dT= 0.1
    "temperature intervall for differential quotient";
  SI.SpecificHeatCapacity c_p_brine2=(Medium.specificEnthalpy_pTX(props.p,props.T+dT/2,props.X)
                                     -Medium.specificEnthalpy_pTX(props.p,props.T-dT/2,props.X))/dT
    "specific heat capacity from enthalpy derivative";

  SI.Pressure dp= 0.1e5 "temperature intervall for differential quotient";
  Real isothermalThrottlingCoefficient=(Medium.specificEnthalpy_pTX(props.p+dp/2,props.T,props.X)
                                       -Medium.specificEnthalpy_pTX(props.p-dp/2,props.T,props.X))/dp
    "isothermalThrottlingCoefficient";

  Real beta=Medium.isobaricExpansionCoefficient(props.state)
    "isobaric expansion coefficient";
    //  Real beta2=props.d*(1/props.d-1/Medium.density_pTX(props.p,props.T-1,props.X)) "isobaric expansion coefficient from density - as calculated in model";
  SI.LinearTemperatureCoefficient beta2=(1-props.d/Medium.density_pTX(props.p,props.T-1,props.X))
    "isobaric expansion coefficient from density - central differential coefficient";

  Real dtdp_s = beta/c_p_brine*props.T/props.d "?";
  //  SI.SurfaceTension sigma= Medium.surfaceTension_T(props.T);

  Real region = Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.region_pT(props.p, props.T)
    "IAPWS region";

  SI.MassConcentration TDS= sum(props.Xi) * props.d "total dissolved solids";
equation
//SPECIFY MEDIUM COMPOSITION {NaCl, KCl, CaCl2, MgCl2, SrCl2}
//   props.Xi = {    0,   0,   0,   0,  0} "pure water";
//  props.Xi = {6*SaltData.M_NaCl/(1+6*SaltData.M_NaCl),0,0,0,0} "6-molar NaCl solution";
//  props.Xi = {n_Na*SaltData.M_NaCl,n_K*SaltData.M_KCl,n_Ca*SaltData.M_CaCl2,0,0}/(1+n_Na*SaltData.M_NaCl+n_K*SaltData.M_KCl+n_Ca*SaltData.M_CaCl2) "specify molalities above";
    props.Xi = {0.0839077010751,0.00253365118988,0.122786737978,0,0}
    "GrSk brine composition (Feldbusch 2-2013 1.1775g/ml V2)";

//SPECIFY THERMODYNAMIC STATE
  props.p = 10*1.01325e5;
 props.T = 273.15+20;
 // props.T = 273.15+50+time "transient - DOES NOT WORK";
end BrineProps1PhaseFull;
