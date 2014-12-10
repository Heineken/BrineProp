within BrineProp.GasData;
function solubility_N2_pTX_Duan2006 "solubility calculation of N2 in seawater Mao&Duan(2006)
  Shide Mao and Zhenhao Duan (2006) A thermodynamic model for calculating nitrogen solubility, gas phase composition and density of the H2O-N2-NaCl system. Fluid Phase Equilibria, 248 (2): 103-114
  273-400 K, 1-600 bar and 0-6 mol/kg
  http://dx.doi.org/10.1016/j.fluid.2006.07.020
  http://www.geochem-model.org/wp-content/uploads/2009/09/46-FPE_248_103.pdf"
//redeclare function extends solubility_N2_pTX
  extends partial_solubility_pTX;
/*  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X[:] "mass fractions m_x/m_H2O";
  input SI.MolarMass MM[:] "molar masses of components";
  input SI.Pressure p_gas;
  output SI.MassFraction c_gas "gas concentration in kg_gas/kg_H2O";
*/
protected
  Real[:] mu_l0_N2_RT_c = { -0.23093813E+02,
                             0.56048525E-01,
                             0.98808898E+04,
                            -0.51091621E-04,
                            -0.13220298E+07,
                            -0.49542866E-03,
                             0.12698747E-05,
                             0.51411144E+00,
                            -0.64733978E-04};

  Real[:] lambda_N2_Na_c = {-0.24434074E+01,
                            0.36351795E-02,
                            0.44747364E+03,
                            0,
                            0,
                           -0.13711527E-04,
                            0,
                            0,
                            0.71037217E-05};

  Real[:] xi_N2_NaCl_c = {-0.58071053E-02,
                          0,
                          0,
                          0,
                          0,
                          0,
                          0,
                          0,
                          0};

  SI.MolarMass M_H2O = MM_vec[end];
  Types.Molality molalities[size(X, 1)]=
      Utilities.massFractionsToMolalities(X,MM_vec);
  Types.Molality m_Cl=molalities[iNaCl] + molalities[iKCl] + 2*molalities[iMgCl2]
       + 2*molalities[iCaCl2];
  Types.Molality m_Na=molalities[iNaCl];
  Types.Molality m_K=molalities[iKCl];
  Types.Molality m_Ca=molalities[iCaCl2];
  Types.Molality m_Mg=molalities[iMgCl2];
  Types.Molality m_SO4=0 "TODO";

  SI.Pressure p_H2O=Modelica.Media.Water.WaterIF97_pT.saturationPressure(T);
//  Pressure_bar p_H2O_bar=SI.Conversions.to_bar(p_sat_H2O_Duan2003(T));
//  Partial_Units.Pressure_bar p_bar=SI.Conversions.to_bar(p);
  SI.MassFraction X_NaCl = molalities[iNaCl]*M_H2O
    "mole fraction of NaCl in liquid phase";
  SI.MolarVolume v_l_H2O=M_H2O/Modelica.Media.Water.WaterIF97_pT.density_pT(p,T);
  Real phi_H2O = fugacity_H2O_Duan2006N2(p,T);
  final constant Real R(final unit="bar.cm3/(mol.K)") = 83.14472
    "Molar gas constant";
/*  Real y_H20 = (1-2*X_NaCl) * p_H2O/(phi_H2O*p) * exp(v_l_H2O*(p-p_H2O)/(Modelica.Constants.R*T)) 
    "equ. 4 (gamma_H2O=1";
  Real y_N2 = 1-y_H20 "mole fraction of CO2 in vapor phase";*/
//  Real y_N2 = p_gas/p "mole fraction of N2 in vapor phase";
  Real phi_N2;
  Real mu_l0_N2_RT;
  Real lambda_N2_Na;
  Real xi_N2_NaCl;
algorithm
// print("mola_N2("+String(p_gas)+","+String(T-273.16)+") (solubility_N2_pTX_Duan2006)");
  if not p_gas>0 then
    X_gas:=0;
  else

    if AssertLevel>0 then
     assert(ignoreLimitN2_T or (273<T and T<400), "\nTemperature out of validity range: T=" + String(T - 273.15) + ".\nTo ignore set ignoreLimitN2_T=true",aLevel);
     assert(ignoreLimitN2_p or (1e5<p and p<600e5),"\nPressure out of validity rangep=" + String(p/1e5) + " bar.\nTo ignore set ignoreLimitN2_p=true",aLevel);
     assert(ignoreLimitN2_b or molalities[iNaCl]<6,"\nMolality out of validity range: mola[NaCl]=" + String(molalities[iNaCl]) + " mol/kg.\nTo ignore set ignoreLimitN2_b=true",aLevel);
    end if;

    phi_N2 :=fugacity_N2_Duan2006(p_gas+p_H2O, T);
    mu_l0_N2_RT :=Par_N2_Duan2006(p_gas+p_H2O,T,mu_l0_N2_RT_c);
    lambda_N2_Na :=Par_N2_Duan2006(p_gas+p_H2O,T,lambda_N2_Na_c);
    xi_N2_NaCl :=Par_N2_Duan2006(p_gas+p_H2O,T,xi_N2_NaCl_c);

  //equ. 9
//    solu := y_N2*p_bar*phi_N2 * exp(-mu_l0_N2_RT - 2*lambda_N2_Na*(m_Na + m_K + 2*m_Ca + 2*m_Mg) - xi_N2_NaCl*(m_Cl+2*m_SO4)*(m_Na + m_K + 2*m_Ca + 2*m_Mg) - 4*0.0371*m_SO4);
    solu := SI.Conversions.to_bar(p_gas)*phi_N2 * exp(-mu_l0_N2_RT - 2*lambda_N2_Na*(m_Na + m_K + 2*m_Ca + 2*m_Mg) - xi_N2_NaCl*(m_Cl+2*m_SO4)*(m_Na + m_K + 2*m_Ca + 2*m_Mg) - 4*0.0371*m_SO4);

//    solu := max(0, solu) "algorithm can return negative values";
//  solu := p_H2O;
//  solu := 0;
    X_gas :=solu*M_N2*X[end] "molality->mass fraction";
  end if;
  //    print("mola_N2("+String(p_gas)+","+String(T-273.16)+","+String(X[1])+")="+String(c_gas)+" (solubility_N2_pTX_Duan2006)");

//  c_gas:=solu*M_N2 "kg_gas / kg_H2O";
//    print("mola_N2("+String(p_gas)+","+String(T-273.16)+")="+String(X_gas)+"->k="+String(X_gas/max(1,p_gas))+" (solubility_N2_pTX_Duan2006)");
end solubility_N2_pTX_Duan2006;
