within BrineProp.Examples;
model BrineProps1PhaseFull
  "Example for 1-phase brine property model using all properties"

  //SPECIFY MEDIUM
  package Medium = Brine_5salts;
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
  SI.SpecificHeatCapacity dh_dT_p=(Medium.specificEnthalpy_pTX(props.p,props.T+0.1,props.X)
                                     -Medium.specificEnthalpy_pTX(props.p,props.T-0.1,props.X))
                                      /0.2
    "specific heat capacity from enthalpy derivative";
  Real isothermalThrottlingCoefficient=(Medium.specificEnthalpy_pTX(props.p+10,props.T,props.X)
                                       -Medium.specificEnthalpy_pTX(props.p-10,props.T,props.X))/20
    "isothermalThrottlingCoefficient";

  Real beta=Medium.isobaricExpansionCoefficient(props.state)
    "isobaric expansion coefficient";
  Real beta2=props.d*(1/props.d-1/Medium.density_pTX(props.p,props.T-1,props.X))
    "isobaric expansion coefficient from density";

  Real dtdp_s = beta/c_p_brine*props.T/props.d "?";
//  SI.SurfaceTension sigma= Medium.surfaceTension_T(props.T);
  Real region = Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.region_pT(props.p, props.T)
    "IAPWS region";

  SI.MassConcentration TDS= sum(props.Xi) * props.d "total dissolved solids";
    parameter Real n_K=1 "potassium molality";
    parameter Real n_Ca=2 "calcium molality";
    parameter Real n_Na=2 "natrium molality";
equation
//SPECIFY MEDIUM COMPOSITION {NaCl, KCl, CaCl2, MgCl2, SrCl2}
//   props.Xi = {    0,   0,   0,   0,  0} "pure water";
//  props.Xi = {6*SaltData.M_NaCl/(1+6*SaltData.M_NaCl),0,0,0,0} "6-molar NaCl solution";
//  props.Xi = {n_Na*SaltData.M_NaCl,n_K*SaltData.M_KCl,n_Ca*SaltData.M_CaCl2,0,0}/(1+n_Na*SaltData.M_NaCl+n_K*SaltData.M_KCl+n_Ca*SaltData.M_CaCl2) "specify molalities above";
    props.Xi = {0.0839077010751,0.00253365118988,0.122786737978,0,0}
    "GrSk brine composition (Feldbusch 2-2013 1.1775g/ml V2)";

//SPECIFY THERMODYNAMIC STATE
  props.p = 10e5;
 props.T = 273.15+50 "+time";
// props.T = 273.15+50+time "transient";

end BrineProps1PhaseFull;
