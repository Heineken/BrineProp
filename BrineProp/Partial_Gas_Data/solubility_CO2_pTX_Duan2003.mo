within BrineProp.Partial_Gas_Data;
function solubility_CO2_pTX_Duan2003 "solubility calculation of CO2 in seawater Duan,Sun(2003)
  Zhenhao Duan and Rui Sun (2003) An improved model calculating CO2 solubility in pure water and aqueous NaCl solutions from 273 to 533 K and from 0 to 2000bar. Chemical Geology, 193:253-271. 
  http://dx.doi.org/10.1016/S0009-2541(02)00263-2
  http://www.geochem-model.org/wp-content/uploads/2009/09/28-CG2003-193.pdf
  fugacity from doi:10.1016/j.marchem.2005.09.001"
  extends partial_solubility_pTX;
/*  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_Sol";
  input Modelica.SIunits.MolarMass MM[:] "molar masses of components";
  output Modelica.SIunits.MassFraction c_gas "gas concentration in kg_gas/kg_H2O";*/
protected
  Partial_Units.Molality solu "CO2 solubility in mol_CO2/kg H2O";
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

  Modelica.SIunits.Pressure p_H2O = Modelica.Media.Water.WaterIF97_base.saturationPressure(T);
//  Modelica.SIunits.Pressure p_H2O = p_sat_H2O_Duan2003(T);
  Partial_Units.Pressure_bar p_bar=Modelica.SIunits.Conversions.to_bar(p);
  Real y_CO2 = p_gas/p "(p-p_H2O)/p mole fraction of CO2 in vapor phase";
  Real phi_CO2 = fugacity_CO2_Duan2006(p,T);
  Real mu_l0_CO2_RT = Par_CO2_Duan2003(p,T,mu_l0_CO2_RT_c);
  Real lambda_CO2_Na = Par_CO2_Duan2003(p,T,lambda_CO2_Na_c);
  Real zeta_CO2_NaCl = Par_CO2_Duan2003(p,T,zeta_CO2_NaCl_c);

 /*  constant Modelica.SIunits.MolarMass M_NaCl = 0.058443 "[kg/mol]";
   constant Modelica.SIunits.MolarMass M_KCl = 0.074551 "[kg/mol]";
   constant Modelica.SIunits.MolarMass M_CaCl2 = 0.1109840 "[kg/mol]";
   constant Modelica.SIunits.MolarMass M_MgCl2 = 0.095236 "[kg/mol]";
   constant Modelica.SIunits.MolarMass M_SrCl2 = 0.158536 "[kg/mol]";
   constant Modelica.SIunits.MolarMass M_H2O = 0.018015 "[kg/mol]";
  constant Real[:] MM = {
    M_NaCl,
    M_KCl,
    M_CaCl2,
    M_MgCl2,
    M_SrCl2,
    M_H2O};
*/
  //constant
  Partial_Units.Molality molalities[size(X,1)]=massFractionsToMolalities(X,MM_vec)
    "TODO neglecting CO2?";
  Partial_Units.Molality m_Cl = molalities[NaCl] + molalities[KCl] + 2*molalities[MgCl2] + 2*molalities[CaCl2];
  Partial_Units.Molality m_Na = molalities[NaCl];
  Partial_Units.Molality m_K = molalities[KCl];
  Partial_Units.Molality m_Ca = molalities[CaCl2];
  Partial_Units.Molality m_Mg = molalities[MgCl2];
  Partial_Units.Molality m_SO4 = 0;

algorithm
  if outOfRangeMode==1 then
    if T<273 or T>533 then
      Modelica.Utilities.Streams.print("T="+String(T)+"K, but CO2 solubility calculation is only valid for temperatures between 0 and 260°C (Partial_Gas_Data.solubility_CO2_pTX_Duan2003)");
    end if;
   if (p<0 or p>2000e5) then
      Modelica.Utilities.Streams.print("p="+String(p/1e5)+" bar, but CO2 fugacity calculation only valid for pressures between 0 and 2000 bar (Partial_Gas_Data.solubility_CO2_pTX_Duan2003)");
   end if;
  elseif outOfRangeMode==2 then
    assert(273<=T and T<=533, "T="+String(T)+"K, but CO2 solubility calculation is only valid for temperatures between 0 and 260°C");
    assert(p<=2000e5, "p="+String(p/1e5)+"bar, but CO2 fugacity calculation only valid for pressures between 0 and 2000 bar");
  end if;

  //equ. 9
    solu := y_CO2*p_bar*phi_CO2* exp(- mu_l0_CO2_RT - 2 * lambda_CO2_Na*(m_Na + m_K + 2*m_Ca + 2*m_Mg) - zeta_CO2_NaCl*m_Cl*(m_Na + m_K + m_Ca + m_Mg) + 0.07*m_SO4);
//    solu := max(0, solu) "algorithm can return negative values";
//  solu := p_H2O;
//  c_gas:=solu*M_CO2 "kg_gas / kg_H2O";
  X_gas :=solu*M_CO2*X[end];
//    Modelica.Utilities.Streams.print("mola_CO2("+String(X_gas)+","+String(T-273.16)+")=->k="+String(X_gas/max(1,p_gas))+" (solubility_CO2_pTX_Duan2003)");
end solubility_CO2_pTX_Duan2003;
