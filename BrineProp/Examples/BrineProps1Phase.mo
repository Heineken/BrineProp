within BrineProp.Examples;
model BrineProps1Phase
   package Medium = Brine_5salts;
//package Medium = Modelica.Media.Water.WaterIF97_pT;
  Medium.BaseProperties props;

  SI.Density d= props.d;
  SI.SpecificEnthalpy h= props.h;
  SI.Temp_C T_C= props.T-273.15;
//  SI.DynamicViscosity eta=Medium.dynamicViscosity(props.state);
  SI.ThermalConductivity lambda=Medium.thermalConductivity(props.state);
//  SI.DynamicViscosity etaC=BrineProp.Viscosities.dynamicViscosity_Zhang_pTXd(props.p,props.T,props.X,props.d,Medium.MM_vec);
  SI.SpecificHeatCapacity c_p_brine= Medium.specificHeatCapacityCp(props.state);
  SI.SpecificHeatCapacity c_p_brine2=(Medium.specificEnthalpy_pTX(props.p,props.T+0.1,props.X)
                                     -Medium.specificEnthalpy_pTX(props.p,props.T-0.1,props.X))
                                      /0.2;
  Real isothermalThrottlingCoefficient=(Medium.specificEnthalpy_pTX(props.p+10,props.T,props.X)
                                       -Medium.specificEnthalpy_pTX(props.p-10,props.T,props.X))/20;

  Real beta2=props.d*(1/props.d-1/Medium.density_pTX(props.p,props.T-1,props.X));
  Real beta=Medium.isobaricExpansionCoefficient(props.state);

  Real dtdp_s = beta/c_p_brine*props.T/props.d;
//  SI.SurfaceTension sigma= Medium.surfaceTension_T(props.T);
/*  Real n_K=0*2.2846;
  Real n_Ca=0*.05/2;
  Real n_Na=5.998 "n_Ca*1.0001";*/
  Real region = Modelica.Media.Water.IF97_Utilities.BaseIF97.Regions.region_pT(props.p, props.T);
equation
//  props.p = 2.21e6 "1.01325e5";
  props.p = 100e5;

//  props.T = 273.16+80+time;

 props.T = 245+273.15;

// props.T = 273.15+50+time;
/*
  props.p = 9.13e5;
  props.T = 273.16+99.61;
 */

//3-GAS
//   props.Xi = {    0.081109,   0.0047275,     0.12519,   0*0.0013225,  0*0.0023842}     "Entsprechend STR04/16 bei GG mit d_l=1199.48 kg/m³ - X_g stimmt";
//  props.Xi = {0.0828264747936192,0.00250100280471612,0.12120451826216,0.000603955715249515,0.00211187457541221};
//   props.Xi[:]={0.0830158933516516,0.00536930693677626,0.111697141833139,0.00112120951982728,0.00198855862851682, 0.00018083,  0.00074452,  6.661e-005} "Elvira 9/11";

//  props.Xi = {     0.08214,   0.0053126,     0.11052,   0*0.0011094,   0*0.0019676}     "Elvira 9/11";
//   props.Xi = {    0,   0,     0,   0,  0};
//  props.Xi = {6*SaltData.M_NaCl/(1+6*SaltData.M_NaCl),0,0,0,0} "NaCl, KCl, CaCl2, MgCl2, SrCl2";
//  props.Xi = {n_Na*SaltData.M_NaCl,n_K*SaltData.M_KCl,n_Ca*SaltData.M_CaCl2,0,0}/(1+n_Na*SaltData.M_NaCl+n_K*SaltData.M_KCl+n_Ca*SaltData.M_CaCl2)     "NaCl, KCl, CaCl2, MgCl2, SrCl2";
    props.Xi = {0.0839077010751,0.00253365118988,0.122786737978,0,0}
    "Elvira 2-2013 1.1775g/ml V2";

algorithm
//  print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
//  print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");
end BrineProps1Phase;
