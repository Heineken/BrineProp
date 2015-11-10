within BrineProp.GasData;
function solubility_CO2_pTX_Duan2006 "CO2 solubility in aqueous saltsolutions"
/*  Zhenhao Duan et al. (2006) An improved model for the calculation of CO2 solubility in aqueous
solutions containing Na+,K+,Ca2+,Mg2+,Cl-, and SO4_2-. Marine Chemistry 98:131-139. 
  fugacity from doi:10.1016/j.marchem.2005.09.001*/
  extends partial_solubility_pTX;
/*  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X[:] "mass fractions m_x/m_Sol";
  input SI.MolarMass MM[:] "molar masses of components";
  input SI.Pressure p_gas;
 output SI.MassFraction c_gas "gas concentration in kg_gas/kg_H2O";*/
protected
  Types.Molality solu "CO2 solubility in mol_CO2/kg H2O";
  Real[:] mu_l0_CO2_RT_c = { 28.9447706,
                        -0.0354581768,
                     -4770.67077,
                         1.02782768e-5,
                        33.8126098,
                         9.04037140e-3,
                        -1.14934031e-3,
                        -0.307405726,
                        -0.0907301486,
                        9.32713393e-4,
                        0};

  Real[:] lambda_CO2_Na_c = {-0.411370585,
                            6.07632013e-4,
                           97.5347708,
                            0,
                            0,
                            0,
                            0,
                           -0.0237622469,
                            0.0170656236,
                            0,
                            1.41335834e-5};

  Real[:] zeta_CO2_NaCl_c = {3.36389723e-4,
                          -1.98298980e-5,
                          0,
                          0,
                          0,
                          0,
                          0,
                          2.12220830e-3,
                          -5.24873303e-3,
                          0,
                          0};

  SI.Pressure p_H2O = Modelica.Media.Water.WaterIF97_pT.saturationPressure(T)
    "SPEEDUP: get from outside";
//  SI.Pressure p_H2O = p_sat_H2O_Duan2003(T);
//  Partial_Units.Pressure_bar p_bar=SI.Conversions.to_bar(p);
//  Real y = p_gas/p "(p-p_H2O)/p mole fraction of CO2 in vapor phase";
  Real phi;
  Real mu_l0_CO2_RT;
  Real lambda_CO2_Na;
  Real zeta_CO2_NaCl;

  //constant
  Types.Molality molalities[size(X, 1)]=Utilities.massFractionsToMolalities(X,MM_vec)
    "TODO neglecting CO2?";
  Types.Molality m_Na=if iNaCl>0 then molalities[iNaCl] else 0;
  Types.Molality m_K=if iKCl>0 then molalities[iKCl] else 0;
  Types.Molality m_Ca=if iCaCl2>0 then molalities[iCaCl2] else 0;
  Types.Molality m_Mg=if iMgCl2>0 then molalities[iMgCl2] else 0;
  Types.Molality m_Sr=if iSrCl2>0 then molalities[iSrCl2] else 0;
  Types.Molality m_SO4=0;
  Types.Molality m_Cl=m_Na + m_K + 2*m_Mg + 2*m_Ca;
algorithm
  if debugmode then
    print("Running solubility_CO2_pTX_Duan2006("+String(p)+","+String(T)+","+String(X[end-3])+","+String(p_gas)+")");
  end if;
  if not p_gas>0 then
    X_gas:=0;
  else
  if AssertLevel>0 then
     assert(ignoreTlimit or ignoreLimitCO2_T or (273<T and T<573),"Temperature out of validity range: T=" + String(T) + "K.\nTo ignore set ignoreLimitCO2_T=true",aLevel);
     assert(ignoreLimitCO2_p or (0<p and p<2000e5),"Pressure out of validity range p=" + String(p/1e5) + " bar.\nTo ignore set ignoreLimitCO2_p=true",aLevel);
     assert(m_Sr==0 or ignoreLimitSalt_soluCO2[iSrCl2],  "The SrCl2 content is not considered here.",aLevel);
     assert(max(molalities[1:end-1])<=4.5,"Maximum salt molality (="+String(max(molalities))+") too high",aLevel);
  end if;
  //equ. 9
    phi :=fugacity_CO2_Duan2006(p_gas+p_H2O, T);
    mu_l0_CO2_RT :=Par_CO2_Duan2003(p_gas+p_H2O,T,mu_l0_CO2_RT_c);
    lambda_CO2_Na :=Par_CO2_Duan2003(p_gas+p_H2O,T,lambda_CO2_Na_c);
    zeta_CO2_NaCl :=Par_CO2_Duan2003(p_gas+p_H2O,T,zeta_CO2_NaCl_c);

    solu :=  phi*SI.Conversions.to_bar(p_gas)* exp(-mu_l0_CO2_RT
            -2*lambda_CO2_Na*(m_Na + m_K + 2*m_Ca + 2*m_Mg)
            -zeta_CO2_NaCl*m_Cl*(m_Na + m_K + m_Mg + m_Ca)
            +0.07*m_SO4);

//    solu := max(0, solu) "algorithm can return negative values";
    X_gas :=solu*M_CO2*X[end] "molality->mass fraction";
  end if;
//  c_gas:=solu*M_CO2 "kg_gas / kg_H2O";
//    print("mola_CO2(p_gas="+String(p_gas)+",T="+String(T)+")=X_gas="+String(X_gas)+"->k="+String(X_gas/max(1,p_gas))+" (solubility_CO2_pTX_Duan2003)");
//    print("solu(p="+String(p)+",T="+String(T)+",y="+String(y)+")="+String(solu)+" (solubility_CO2_pTX_Duan2003)");
end solubility_CO2_pTX_Duan2006;
