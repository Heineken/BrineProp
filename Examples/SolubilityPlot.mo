within Brine.Examples;
model SolubilityPlot
  package Medium = Brine_Duan_Multi_TwoPhase_ngas_3;
    Modelica.SIunits.Pressure p=1.013e5;
    Modelica.SIunits.Temperature T=373.11;
    Modelica.SIunits.Pressure p_gas=1e5*time;
    Real solu=Partial_Gas_Data.solubility_N2_pTX_Duan2006(p,T,
     {0.089190167,0.005198142,0.137663206,0*0.001453819,0*0.002621571, 7.85e-4,  5.73e-5,6.98e-5,.767036}, Medium.MM_vec,
      p_gas);
//      Medium.BaseProperties props;
equation
/*  props.p=p;
  props.T=T;
  props.Xi = {1*Medium.Salt_data.M_NaCl,0,0,0,0, 1e-3,1e-3,1e-3};
  solu=solubility_N2_pTX_Duan2006(props.p,props.T,
     props.X, Medium.MM_vec,
      1e5*time);*/
end SolubilityPlot;
