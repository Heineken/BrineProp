within BrineProp.Examples;
model Degassing
package Medium = BrineProp.Brine_5salts_TwoPhase_3gas;
  Medium.BaseProperties props;
  SI.Density d= props.d;  /**/
SI.DynamicViscosity eta_l = Medium.dynamicViscosity_liq(props.state);
SI.DynamicViscosity eta_g = Medium.dynamicViscosity_gas(props.state);

//  Real c_CO2=Partial_Gas_Data.solubility_CO2_pTX_Duan2003(props.p,props.T,props.X,Medium.MM_vec,props.p);
//  Real c_N2=Partial_Gas_Data.solubility_N2_pTX_Duan2006(props.p,props.T,props.X,Medium.MM_vec,props.p);
//  Real c=Partial_Gas_Data.solubility_CH4_pTX_Duan2006(props.p,props.T,props.X,Medium.MM_vec,props.p_gas[2]);

/*  Real k_CH4b=Partial_Gas_Data.HenryCoefficients_CH4(props.T); 
*/
/*    SI.Pressure p_degas_CO2=Medium.degassingPressure_CO2_Duan2003(props.p,props.T,props.X,Medium.MM_vec);
  SI.Pressure p_degas_N2=Medium.degassingPressure_N2_Duan2006(props.p,props.T,props.X,Medium.MM_vec);
  SI.Pressure p_degas_CH4=Medium.degassingPressure_CH4_Duan2006(props.p,props.T,props.X,Medium.MM_vec);
  */
/*
  Real[:] y=cat(1,props.p_gas,{props.p_H2O})/props.p 
    "volume fraction of gas phase";
//  Real[:] yy=y[2:3]./fill(1-y[4],2) "volume fraction of gas phase w/o H2O";

/*  Real[:] y_l=if not max(props.X_l[6:8])>0 then fill(0,Medium.nX_gas) else props.X_l[6:8]./Medium.MM_gas / sum(props.X_l[6:8]./Medium.MM_gas) 
    "mol fraction of dissolved gases";
  Real V_l = sum(props.X_l[6:8]./Medium.MM_gas)*22.4/props.X_l[end] 
    "Liter of dissolved gas per kg_brine would have after complete degassing at standard conditions";
  Real ratioGasLiquid = V_l*props.d_l/1000;
*/
//  Real V_l = props.X_l[8]/Medium.MM_gas[3]/22.4*1000/props.X_l[end];
/*
  Real[:] m=y[2:3]./fill(1-y[4],2) "mass fraction of gas in gas phase";

  Real m_t = sum(props.X[6:8]) "Total mass fraction of gases in fluid";
  Real m_g = if not max(props.X_l[6:8])>0 then 0 else sum(X_g[6:8])*props.x/sum(props.X[6:8]) 
    "Fraction of gas masses in gas phase";
*/

equation
  props.p = (455-time)*1e5;
  props.T = 112.164+273.15;

//  props.Xi = {6*SaltData.M_NaCl/(1+6*SaltData.M_NaCl),0,0,0,0, 0,0,0} "NaCl, KCl, CaCl2, MgCl2, SrCl2, CO2, N2, CH4";
//  props.Xi = {0,SaltData.M_KCl/(1+SaltData.M_KCl),0,0,0,0,0,0} "NaCl, KCl, CaCl2, MgCl2, SrCl2, CO2, N2, CH4";
//  props.Xi = {0,0,SaltData.M_CaCl2/(1+SaltData.M_CaCl2),0,0,0,0,0} "NaCl, KCl, CaCl2, MgCl2, SrCl2, CO2, N2, CH4";
    props.Xi = {0.083945671051201,0.00253479771131107,0.122842299461699,0*0.000612116692496665,0*0.00214041137028575,  0.00016883,  0.00073459, 6.5652e-005}
    "Elvira 2-2013 1.1775g/ml";

  annotation (experiment(StopTime=455), __Dymola_experimentSetupOutput);
end Degassing;
