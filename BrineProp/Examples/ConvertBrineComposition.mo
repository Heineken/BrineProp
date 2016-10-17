within BrineProp.Examples;
model ConvertBrineComposition
  "convert volume related composition into mass fraction g/l -> kg/kg"
//  Finds brine composition as mass fraction from infos in g/l by iterative density calculation

  package Medium = BrineProp.Brine5salts3gas(AssertLevel=0);

  constant SI.MolarMass M_cation[:] = {22.98977,39.0983,40.078,24.305,87.62}/1000 "kg/mol";
  constant SI.MolarMass M_Cl = 35.453/1000;
  constant Integer n_cation[:] = {1,1,2,2,2} "chlorine ions per cation";
/*
  // From cation molarities mmol/l Na+,K+,Ca2+,Mg2+,Sr2+
  constant SI.Concentration b_cation[:]={1691.386,40.037,1303.355,0*7.570,0*15.899} 
    "Feldbusch 02/2013";
  constant SI.MassConcentration  c_cation[:]=b_cation.*M_cation;
  constant SI.MassConcentration  c[:]= c_cation+M_Cl*b_cation.*n_cation 
    "g/l NaCl,KCl,CaCl2,MgCl2,SrCl2";
*/

// From cation concentrations g/l Na+,K+,Ca2+,Mg2+,Sr2+
//  constant SI.MassConcentration  c_cation[:]={38.88456464,1.565384615,52.23585483,0*0.184,0*1.393076923} "Feldbusch 02/2013";
  constant SI.MassConcentration  c_cation[:]={112,16.8,13.2,0*4.98,0*0.419}/1000
    "GTN 02/2013";
  constant SI.MassConcentration  c[:]= c_cation+M_Cl*c_cation./M_cation.*n_cation
    "g/l NaCl,KCl,CaCl2,MgCl2,SrCl2";

/*
//From Chloride concentrations g/l NaCl,KCl,CaCl2,MgCl2,SrCl2
//  SI.MassConcentration  c[:]={97.331,5.673,150.229,1.587,2.861} "STR04/16";
    constant SI.MassConcentration c[:]={98.84925634,2.984821797,144.6515323,0*0.720790948,0*2.520416712} "Feldbusch 02/2013";
*/

  //volume fractions of gas phase CO2,N2,CH4
  //  Real[:] gasVolumeFractions={0.03825,0.8285,0.12825} "STR04/16";
/*  Real[:] gasVolumeFractions={0.0935,78.5,13.5} "Feldbusch 02/2013";
Real gasLiquidRatio=0.8475 "GVF/(1-GVF) at STP";*/

  Real[:] gasVolumeFractions={0.418,0.202,0.326} "GTN";
  Real gasLiquidRatio=1 "GVF/(1-GVF) at STP GTN";

  Medium.BaseProperties props;
  SI.Density d_l(start=1200)= props.state.d_l;

  Real gasVolumePerLiquidVolume[:]= gasVolumeFractions*gasLiquidRatio;
  SI.Density[:] rho_gas = props.p*Medium.MM_gas/(Modelica.Constants.R*props.T);
  Real gasMassPerLiquidMass[:]= gasVolumePerLiquidVolume.*rho_gas/props.state.d_l;
  SI.MassConcentration[:] X_g_gas(each start=0)= gasMassPerLiquidMass/(1+sum(gasMassPerLiquidMass))
    "X_g.q gas mass in gas phase per fluid mass";
  // SI.MassConcentration[:] X_g_gas(each start=0)={max(0,g/(1+sum(g)))  for  g in  gasMassPerLiquidMass};
  SI.Pressure[:] PartialPressures(each start=1e4)=(props.p-props.p_H2O)*gasVolumeFractions;
  SI.MassConcentration[:] X_l_gas(each start=0) = {
    Medium.solubility_CO2_pTX_Duan2006(
        props.p,
        props.T,
        props.X,
        Medium.MM_vec,
        PartialPressures[1]),
    Medium.solubility_N2_pTX_Mao2006(
        props.p,
        props.T,
        props.X,
        Medium.MM_vec,
        PartialPressures[2]),
    Medium.solubility_CH4_pTX_Duan2006(
        props.p,
        props.T,
        props.X,
        Medium.MM_vec,
        PartialPressures[3])}/(1 - props.x)
    "gas mass in liquid phase per fluid mass";

  Real V_l = sum(props.X_l[i_gas]./Medium.MM_gas)*22.4/props.X_l[Medium.nX]
    "Volume dissolved gas would have at standard conditions / OM workaround: Medium.nX instead of end";
  SI.MassConcentration TDS= sum(props.X_l[1:Medium.nX_salt])*props.d_l;
  Real[:] y_l=if not max(props.X_l[i_gas])>0 then fill(0,Medium.nX_gas) else props.X_l[i_gas]./Medium.MM_gas / sum(props.X_l[i_gas]./Medium.MM_gas)
    "mol fraction of dissolved gases";
  Real[:] y=props.p_gas/props.p "volume fraction of gas phase";
  SI.MassFraction[:] X_g=props.X-props.X_l*(1-props.x)
    "mass in gas phase per kg fluid";

    /*    SI.Pressure p_degas_CO2=Medium.degassingPressure_CO2_Duan2003(props.p,props.T,props.X,Medium.MM_vec);
  SI.Pressure p_degas_N2=Medium.degassingPressure_N2_Duan2006(props.p,props.T,props.X,Medium.MM_vec);
  SI.Pressure p_degas_CH4=Medium.degassingPressure_CH4_Duan2006(props.p,props.T,props.X,Medium.MM_vec);
  */
protected
  Integer i_gas[:]=(Medium.nX-Medium.nX_gas):(Medium.nX-1);
equation
  props.p = 1.01325e5 "STP";
  props.T = 273.16+23 "STP";
  props.Xi = cat(1,c/d_l/0.999190406582784,X_g_gas+X_l_gas) "props.state.d_l";

algorithm
//  print("rho_l="+String(props.state.d_l)+" kg/m^3, TDS = " + String(TDS) + " g/l");
  print("X=" + Modelica.Math.Matrices.toString({props.X}) + " ");

//  print("sum(X_l)="+String(sum(props.state.X_l)-1)+"");

  annotation (experiment(__Dymola_NumberOfIntervals=1),
      __Dymola_experimentSetupOutput);
end ConvertBrineComposition;
