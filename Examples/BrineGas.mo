within BrineProp.Examples;
model BrineGas
package Medium = BrineProp.BrineGas_3Gas;
//package Medium = Modelica.Media.Air.SimpleAir;
//package Medium = PartialBrineGas;
  Medium.BaseProperties props;

/*  package Medium = BrineProp.Brine_5salts_TwoPhase_3gas;
  Medium.BaseProperties props;

  SI.Density d= props.d;
  SI.SpecificEnthalpy h= props.h;
  SI.Temp_C T_C= props.T-273.15;

  SI.ThermalConductivity lambda=Medium.thermalConductivity(props.state);
*/
  SI.DynamicViscosity eta=Medium.dynamicViscosity(props.state);
  //  SI.DynamicViscosity etaC=BrineProp.Viscosities.dynamicViscosity_Zhang_pTXd(props.p,props.T,props.X,props.d,Medium.MM_vec);
  SI.SpecificHeatCapacity c_p= Medium.specificHeatCapacityCp_pTX(props.p,props.T,props.X);
  SI.SpecificHeatCapacity c_p_sat= Medium.specificHeatCapacityCp(props.state)
    "c_p of gas phase after VLE";

//  SI.Density d_gas= BrineProp.BrineGas_3Gas.density_pTX_unsued(props.p,props.T, props.X[end - Medium.nX_gas:end]);

//  SI.Density d_gas1= BrineGas.density(st);
//  SI.Density d_gas2= BrineGas.density_pTX(props.p,props.T, BrineProp.BrineGas_3Gas.waterSaturatedComposition_pTX(props.p,props.T,props.X[end - Medium.nX_gas:end]));
//  SI.Density d_gas3= BrineGas.density_pTX(props.p,props.T, props.X[end - Medium.nX_gas:end]);

  /*
//  SI.SpecificHeatCapacity c_p_gas1= BrineProp.BrineGas_3Gas.specificHeatCapacityCp_pTX(props.p,props.T, props.X) 
    "cp of gas part";

  //SI.SpecificHeatCapacity c_p_gas1= BrineGas_3Gas.specificHeatCapacityCp(BrineGas_3Gas.ThermodynamicState(props.p,props.T, props.X));
  SI.SpecificHeatCapacity c_p_MoistGas= MoistAir.specificHeatCapacityCp(st);

  SI.SpecificHeatCapacity c_p_gasphase= Medium.specificHeatCapacityCp_gas_pTX_unused(props.p,props.T, props.state.X_g) 
    "cp of gas part";

  SI.SpecificHeatCapacity c_p_DryGas=MoistAir.specificHeatCapacityCp_pTX(props.p,props.T, props.state.X_g);

  SI.SpecificHeatCapacity c_p_watervapour=
  Modelica.Media.Water.WaterIF97_base.specificHeatCapacityCp(
    Modelica.Media.Water.WaterIF97_base.setState_pTX(p=Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.psat(props.T)-1,
    T=props.T));

 SI.MassFraction[:] X_sat=BrineProp.BrineGas_3Gas.waterSaturatedComposition_pTX(props.p,props.T,props.state.X_g);
MoistAir.ThermodynamicState st=MoistAir.ThermodynamicState(props.p,props.T,props.state.X_g);

//  SI.DynamicViscosity eta_gas = MoistAir.dynamicViscosity(MoistAir.ThermodynamicState(p=0,T=293.15,X={0}));
*/
equation
  props.p = 10*1e5;
  props.T = 293.15;
  props.Xi={0.1,0.1,0.1};
//  gasprops.h = 2e6;
/*  props.p = 2.00138*1e5;
  props.T = 352.089;
// props.Xi={0.8};
 props.Xi={0.0892134,0.0051995,0.137699,0,0,0.000639217,7.0356e-6,1.5464e-5};
 */
algorithm
//  print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
//    print("sum(X_sat)="+String(sum(X_sat)));

end BrineGas;
