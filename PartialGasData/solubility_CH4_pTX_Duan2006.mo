within BrineProp.PartialGasData;
function solubility_CH4_pTX_Duan2006 "Duan ZH, Mao SD. (2006) A thermodynamic model for calculating methane solubility, density and gas phase composition of methane-bearing aqueous fluids from 273 to 523 K and from 1 to 2000 bar. Geochimica et Cosmochimica Acta, 70 (13): 3369-3386. 
  http://geochem-model.org/Publications/43-GCA_2006_3369.pdf
  http://dx.doi.org/10.1016/j.gca.2006.03.018TODO Umrechnung andere Salz in NaCl"
  extends partial_solubility_pTX;
//  output SI.MassFraction c_gas "gas concentration in kg_gas/kg_H2O";

protected
  Real[:] mu_l0_CH4_RT_c = { 0.83143711E+01,
                            -0.72772168E-03,
                             0.21489858E+04,
                            -0.14019672E-04,
                            -0.66743449E+06,
                             0.76985890E-02,
                            -0.50253331E-05,
                            -0.30092013E+01,
                             0.48468502E+03,
                             0};

  Real[:] lambda_CH4_Na_c = {-0.81222036E+00,
                              0.10635172E-02,
                              0.18894036E+03,
                              0,
                              0,
                              0.44105635E-04,
                              0,
                              0,
                              0,
                             -0.46797718E-10};

  Real[:] xi_CH4_NaCl_c = {-0.29903571E-02,
                              0,
                              0,
                              0,
                              0,
                              0,
                              0,
                              0,
                              0,
                              0};

  SI.MolarMass M_H2O = MM_vec[end];
  Partial_Units.Molality molalities[size(X,1)];
//  Partial_Units.Molality molefractions[size(X,1)];
  Partial_Units.Molality m_Cl;
  Partial_Units.Molality m_Na;
  Partial_Units.Molality m_K;
  Partial_Units.Molality m_Ca;
  Partial_Units.Molality m_Mg;
  Partial_Units.Molality m_SO4;
//  Real X_NaCl = molalities[NaCl]*M_H2O "mole fraction of NaCl in liquid phase";

//  SI.Pressure p_H2O = Modelica.Media.Water.WaterIF97_base.saturationPressure(T);
  SI.Pressure p_H2O = p_sat_H2O_Duan2003(T) "TODO mit übergeben";
//  Partial_Units.Pressure_bar p_bar=SI.Conversions.to_bar(p);
  SI.MolarVolume v_l_H2O=M_H2O/Modelica.Media.Water.WaterIF97_base.density_pT(p,T);
  Real phi_H2O = fugacity_H2O_Duan2006CH4(p,T);
//  Real y_H20 "mole fraction of H2O in vapor phase";
//  Real y_CH4 =p_gas/p;
  Real phi_CH4;
  Real mu_l0_CH4_RT;
  Real lambda_CH4_Na;
  Real xi_CH4_NaCl;
algorithm
  if not p_gas>0 then
    X_gas:=0;
  else
   if outOfRangeMode==1 then
     if (273>T or T>273+250) then
        print("T="+String(T)+" K, but  CH4 solubility  calculation is only valid for temperatures between 0 and 250°C (Partial_Gas_Data.solubility_CH4_pTX_Duan1992)");
     end if;
     if (p<1e5 or p>2000e5) then
        print("p="+String(p/1e5)+" bar, but CH4 fugacity calculation only valid for pressures between 1 and 1600 bar (Partial_Gas_Data.solubility_CH4_pTX_Duan1992)");
     end if;
   elseif outOfRangeMode==2 then
     assert(273.15<=T and T<=273+250, "T="+String(T-273.15)+"°C, but CH4 solubility calculation is only valid for temperatures between 0 and 250°C (solubility_CH4_pTX_Duan1992)");
     assert(p<=1600e5, "p="+String(p/1e5)+"bar, but CH4 fugacity calculation only valid for pressures up to 1600 bar (solubility_CH4_pTX_Duan1992)");
   end if;

  //  (molefractions,molalities):=massFractionsToMoleFractions(X, MM);
    molalities:=massFractionsToMolalities(X, MM_vec);
  // print("molefractions[NaCl]="+String(molefractions[NaCl])+" (Partial_Gas_Data.solubility_CH4_pTX_Duan1992)");
    m_Cl :=molalities[NaCl] + molalities[KCl] + 2*molalities[MgCl2] + 2*
      molalities[CaCl2];
    m_Na :=molalities[NaCl];
    m_K :=molalities[KCl];
    m_Ca :=molalities[CaCl2];
    m_Mg :=molalities[MgCl2];
    m_SO4 :=0 "TODO";

  /*  y_H20 :=(1 - 2*molefractions[NaCl])*p_H2O/(phi_H2O*p)*exp(v_l_H2O*(p -
      p_H2O)/(Modelica.Constants.R*T)) "equ. 5 (x_H2O=1-2x_NaCl)";
    y_CH4 :=1 - y_H20 "mole fraction of CH4 in vapor phase"; Das ist wie im Paper*/
    //y_CH4 = (p-p_H2O)/p "mole fraction of CH4 in vapor phase TODO:neglecting other phases?" Das ist wie CO2_Duan2003;

    phi_CH4 :=fugacity_CH4_Duan1992(p_gas+p_H2O, T)
      "TODO: Hier p_H2O hinzugefügt so wie bei N2 und CO2. Is das richtig?";
    mu_l0_CH4_RT :=Par_CH4_Duan2006(p_gas+p_H2O,T,mu_l0_CH4_RT_c);
    lambda_CH4_Na :=Par_CH4_Duan2006(p_gas+p_H2O,T,lambda_CH4_Na_c);
    xi_CH4_NaCl :=Par_CH4_Duan2006(p_gas+p_H2O,T,xi_CH4_NaCl_c);

    //equ. 10
  //    solu := y_CH4*phi_CH4*p_bar*exp(-mu_l0_CH4_RT - 2*lambda_CH4_Na*(m_Na + m_K + 2*m_Ca + 2*m_Mg) - xi_CH4_NaCl*(m_Na + m_K + 2*m_Ca + 2*m_Mg)*(m_Cl + 2*m_SO4) - 4*0.0332*m_SO4);
      solu := SI.Conversions.to_bar(p_gas)*phi_CH4*exp(-mu_l0_CH4_RT - 2*lambda_CH4_Na*(m_Na + m_K + 2*m_Ca + 2*m_Mg) - xi_CH4_NaCl*(m_Na + m_K + 2*m_Ca + 2*m_Mg)*(m_Cl + 2*m_SO4) - 4*0.0332*m_SO4);

  //  solu := max(0, solu) "algorithm can return negative values";
  //  solu := p_H2O;
  //  solu := 0;
    X_gas :=solu*M_CH4*X[end] "molality->mass fraction";
  end if;
//  c_gas:=solu*M_CH4 "kg_gas / kg_H2O";
//      print("mola_CH4("+String(p_gas)+","+String(T-273.16)+")="+String(X_gas)+"->k="+String(X_gas/max(1,p_gas))+" (solubility_CH4_pTX_Duan2006)");
end solubility_CH4_pTX_Duan2006;
