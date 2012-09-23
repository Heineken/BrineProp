within Brine;
partial package Partial_Gas_Data "Contains molar masses of gases"

  /*replaceable record GasConstants
    extends Modelica.Icons.Record;
    Modelica.SIunits.MolarMass M_salt "Molar Mass in kg/mol";
    annotation(Documentation(info="<html></html>"));
  end GasConstants;*/
//  constant Real[:] MM_gas;

  constant Modelica.SIunits.MolarMass M_CO2 = Modelica.Media.IdealGases.SingleGases.CO2.data.MM
    "0.0440095 [kg/mol]";
  constant Integer nM_CO2 = 1 "number of ions per molecule";
   constant Modelica.SIunits.MolarMass M_N2 = Modelica.Media.IdealGases.SingleGases.N2.data.MM
    "0.0280134 [kg/mol]";
  constant Integer nM_N2 = 1 "number of ions per molecule";
   constant Modelica.SIunits.MolarMass M_CH4 = Modelica.Media.IdealGases.SingleGases.CH4.data.MM
    "0.01604246 [kg/mol]";
  constant Integer nM_CH4 = 1 "number of ions per molecule";

/*    constant Modelica.SIunits.MolarMass M_CO2 = 0.0440095 "[kg/mol]";
   constant Modelica.SIunits.MolarMass M_N2 = 0.0280134 "[kg/mol]";
 */
//TODO braucht man die Variablen direkt oder nur indirekt via MM(CO2)

/* passiert in PartialBrine_Multi_TwoPhase_xgas 
 if gasNames{1} == "carbondioxide" then
      constant Real[:] MM_gas = {M_CO2};
  elseif gasNames{1} == "nitrogen" then
      constant Real[:] MM_gas = {M_N2};
  end if;
  
 constant Real[:] MM_gas = {
    M_CO2,M_N2};*/

  function polyval "Calculates polynomial value"
    input Real c[:];
    input Real x;
    output Real y;
  algorithm
      y := 0;
      for j in 0:size(c, 1)-1 loop
        y := y + c[end-j]*x^j;
      end for;
  end polyval;

  replaceable function massFractionsToMoleFractions
    "Return mole_i/sum(mole_i) from mass fractions X"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.MassFraction X[:] "Mass fractions of mixture";
    input Modelica.SIunits.MolarMass MM[:] "molar masses of components";
    output Partial_Units.Molality molefractions[size(X, 1)] "Molalities";
    output Partial_Units.Molality molalities[size(X, 1)]
      "Molalities moles/m_H2O";
  //  Real invMMX[size(X, 1)] "inverses of molar weights";
  //  Modelica.SIunits.MolarMass Mmix "molar mass of mixture";
  //  Real n[size(X, 1)] "number of Moles";
  protected
    Real n_total;
  algorithm
   assert(size(X, 1)==size(MM, 1), "Inconsistent vectors for mass fraction("+String(size(X, 1))+") and molar masses("+String(size(MM, 1))+")");
  // Modelica.Utilities.Streams.print(String(size(X,1))+" "+String(X[end]));
  //  printVector(MM);
    for i in 1:size(X, 1) loop
  // Modelica.Utilities.Streams.print("MMX["+String(i)+"]="+String(MMX[i]));
      molalities[i] := if X[end]>0 then X[i]/(MM[i]*X[end]) else -1;
  //    n[i] := X[i]/MMX[i];
    end for;
    n_total :=sum(molalities);
    for i in 1:size(X, 1) loop
      molefractions[i] := molalities[i]/n_total;
    end for;
    annotation(smoothOrder=5);
  end massFractionsToMoleFractions;

  partial function partial_fugacity_pTX
    input Modelica.SIunits.Pressure p "in Pa";
    input Modelica.SIunits.Temp_K T "in K";
    output Real phi "fugacity coefficient";
  protected
      Partial_Units.Pressure_bar p_bar=Modelica.SIunits.Conversions.to_bar(p);
  end partial_fugacity_pTX;

  function p_sat_CO2
    "calculates saturation pressure, polynom derived from EES calculations"
    input Modelica.SIunits.Temp_K T;
    output Modelica.SIunits.Pressure p_sat;
  algorithm
    assert(T<305,"Temperature above critical Temperature (304.1 K)");
    p_sat :=1178.4*T^2 - 555378*T + 7E+07;
  end p_sat_CO2;

  partial function partial_solubility_pTX
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_Sol";
    input Modelica.SIunits.MolarMass MM_vec[:] "molar masses of components";
    input Modelica.SIunits.Pressure p_gas;
  //  output Modelica.SIunits.MassFraction c_gas "gas concentration in kg_gas/kg_H2O";
    output Modelica.SIunits.MassFraction X_gas
      "gas concentration in kg_gas/kg_fluid";
  protected
    Partial_Units.Molality solu "gas solubility";
  //algorithm
  //    Modelica.Utilities.Streams.print("mola("+String(X_gas)+","+String(T-273.16)+")=->k="+String(X_gas/max(1,p_gas))+" (partial_solubility_pTX)");
  end partial_solubility_pTX;

  function solubility_res
  extends Modelica.Math.Nonlinear.Interfaces.partialScalarFunction;
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_Sol";
    input Modelica.SIunits.MolarMass MM_vec[:] "molar masses of components";
    /**/
    input Partial_Gas_Data.partial_solubility_pTX solufun;
  //  input Modelica.Icons.Function solufun;
    input Modelica.SIunits.MassFraction c_gas;
  //  input Modelica.SIunits.Pressure p_gas=u;
  protected
    Modelica.SIunits.MassFraction solu;
  algorithm
      y:=c_gas-solufun(p=p,T=T,X=X,MM_vec=MM_vec,p_gas=u) "*X[end]";
  end solubility_res;

  partial function partial_degassingPressure_pTX
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_Sol";
    input Modelica.SIunits.MolarMass MM_vec[:] "molar masses of components";
    output Modelica.SIunits.Pressure p_gas;
  end partial_degassingPressure_pTX;

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

  function p_sat_H2O_Duan2003
    "calculates saturation pressure of water, with the equation given in Duan2003"
    /*An improved model calculating CO2 solubility in pure water and aqueous NaCl solutions from 273 to 533 K and from 0 to 2000 bar"*/
    input Modelica.SIunits.Temp_K T;
    output Modelica.SIunits.Pressure p_sat;
  protected
    Real[:] c={-38.640844,
                 5.894842,
                59.876516,
                26.654627,
                10.637097};
    Modelica.SIunits.Pressure P_c = 220.85e5;
    Modelica.SIunits.Temperature T_c = 647.29;
    Real t=(T-T_c)/T_c;
  algorithm
  //  assert(T<305,"Temperature above critical Temperature (304.1 K)");
    p_sat :=(P_c*T/T_c)*(1+c[1]*(-t)^1.9+c[2]*t+c[3]*t^2+c[4]*t^3+c[5]*t^4);
  end p_sat_H2O_Duan2003;

  function Par_CO2_Duan2003 "Duan,Sun(2003)"
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input Real[:] c;
    output Real Par;
  protected
    Partial_Units.Pressure_bar p_bar=Modelica.SIunits.Conversions.to_bar(p);
  algorithm
    //eq. 7
    Par :=  c[1] + c[2]*T + c[3]/T + c[4]*T^2 + c[5]/(630-T) + c[6]*p_bar + c[7]*p_bar*Modelica.Math.log(T) + c[8]*p_bar/T + c[9]*p_bar/(630-T) + c[10]*p_bar^2/(630-T)^2 + c[11]*T*Modelica.Math.log(p_bar);
  end Par_CO2_Duan2003;

  function fugacity_CO2_Duan2006
    "Calculation of fugacity coefficient according to (Duan 2006)"
    extends partial_fugacity_pTX;
  /*  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  output Real phi;*/

  protected
    Partial_Units.Pressure_bar P_1;
    Real c[15];
  algorithm
    if outOfRangeMode==1 then
      if T<273 or T>573 then
        Modelica.Utilities.Streams.print("T="+String(T)+"K, but CO2 solubility calculation is only valid for temperatures between 0 and 260°C (Partial_Gas_Data.fugacity_CO2_Duan2006)");
      end if;
     if (p<0 or p>2000e5) then
        Modelica.Utilities.Streams.print("p="+String(p/1e5)+" bar, but CO2 fugacity calculation only valid for pressures between 0 and 2000 bar (Partial_Gas_Data.fugacity_CO2_Duan2006)");
     end if;
    elseif outOfRangeMode==2 then
      assert(T>273 and T<573, "T="+String(T-273.15)+"°C out of range(0...300°C) for CO2 fugacity calculation (fugacity_CO2_Duan2006)");
      assert(p<2000e5, "p="+String(p/1e5)+" bar out of range for CO2 fugacity calculation (fugacity_CO2_Duan2006)");
    end if;

      if T<305 then
        P_1 := Modelica.SIunits.Conversions.to_bar(p_sat_CO2(T));
      elseif T<405 then
        P_1 := 75 + (T-305)*1.25;
      else
        P_1 := 200;
      end if;

    if  p_bar<P_1 then
      //1 273<T<573 and p_bar<P_1
        c[:] := {1,
                 1.0,
                 4.7586835E-3,
                -3.3569963E-6,
                 0
                -1.3179396,
                -3.8389101E-6,
                 0,
                 2.2815104E-3,
                 0, 0, 0, 0, 0, 0, 0};
    elseif T>435 then
            //6 T>435 and p_bar>P_1
            c[:] :={-1.5693490E-1,
                        4.4621407E-4,
                       -9.1080591E-7,
                         0,
                         0,
                         1.0647399E-7,
                         2.4273357E-10,
                         0,
                         3.5874255E-1,
                         6.3319710E-5,
                      -249.89661,
                         0,
                         0,
                       888.76800,
                        -6.6348003E-7};

    elseif p_bar<1000 then
      //P_1<p_bar<1000 and 273<T<435
      if T<340 then
       //2 273<T<340 and P_1<p_bar<1000
             c[:] := {-7.1734882E-1,
                      1.5985379E-4,
                     -4.9286471E-7,
                      0,
                      0,
                     -2.7855285E-7,
                      1.1877015E-9,
                      0,
                      0,
                      0,
                      0,
                    -96.539512,
                      4.4774938E-1,
                    101.81078,
                      5.3783879E-6};
      else
        //T>340
            //4 340<T<435 and P_1<p_bar<1000
            c[:] := {5.0383896,
                    -4.4257744E-3,
                     0,
                     1.9572733,
                     0,
                     2.4223436E-6,
                     0,
                    -9.3796135E-4,
                    -1.5026030,
                     3.0272240E-3,
                    -31.377342,
                    -12.847063,
                     0,
                     0,
                    -1.5056648E-5};
      end if;

        else
      //p_bar>1000 bar and 273<T<435
          if T<340 then
            //3 273<T<340 and p_bar>1000
            c[:] :={-6.5129019E-2,
                    -2.1429977E-4,
                    -1.1444930E-6,
                     0.0,
                     0.0,
                    -1.1558081E-7,
                     1.1952370E-9,
                     0.0,
                     0.0,
                     0.0,
                     0.0,
                    -221.34306,
                     0.0,
                     71.820393,
                     6.6089246E-6};
          else
           //T>340"
            //5 340<T<435 and p_bar>1000
            c[:] := {-16.063152,
                      -2.7057990E-3,
                       0,
                       1.4119239E-1,
                       0,
                       8.1132965E-7,
                       0,
                      -1.1453082E-4,
                       2.3895671,
                       5.0527457E-4,
                     -17.763460,
                     985.92232,
                       0,
                       0,
                      -5.4965256E-7};
          end if;
    end if;

    phi := c[1] + (c[2] + c[3]*T + c[4]/T + c[5]/(T-150))*p_bar + (c[6] + c[7]*T + c[8]/T)*p_bar^2 + (c[9] + c[10]*T + c[11]/T)*log(p_bar) +(c[12]+c[13]*T)/p_bar + c[14]/T + c[15]*T^2;
  //  PowerPlant.Components.PipeStuff.print_msg(phi,"phi_CO2=");
  end fugacity_CO2_Duan2006;

  function degassingPressure_CO2_Duan2003
    "calculates degassing pressure from concentration of dissolved gas"
    extends partial_degassingPressure_pTX;
  /*protected 
    Modelica.SIunits.MassFraction solu_soll=X[end-3];*/

  algorithm
  /*  p_gas := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
      function solubility_CO2_Duan2003_res(p=p,T=T,X=X,MM_vec=MM_vec),
      0,
      1e8,
      1e-8);
*/
    p_gas := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
        function solubility_res(
          solufun=function solubility_CO2_pTX_Duan2003(),p=p,T=T,X=X,MM_vec=MM_vec,
          c_gas=X[end-3]),
        0,
        1e10,
        1e-8);

  //  Modelica.Utilities.Streams.print("p_sat_CO2("+String(X[end-3])+")="+String(p_gas)+" (degassingPressure_CO2_pTX)");

  /*    while abs(solu-solu_soll)>1e-8 loop
    p_gas:=p_gas_neu;
    solu:=solubility_CO2_pTX_Duan2003(
        p=p,
        T=T,
        X=X,
        MM=MM_vec,
        p_gas=p_gas);
    z:=z + 1;
    assert(z<1000," Reached maximum number of iterations for solution equilibrium calculation. (quality_pTX)");
    Modelica.Utilities.Streams.print("p_gas="+String(p_gas_neu)+"->"+String(abs(p_gas-p_gas_neu)));
  end while;*/
  end degassingPressure_CO2_Duan2003;

  function solubility_N2_pTX_Duan2006 "solubility calculation of N2 in seawater Mao&Duan(2006)
  Shide Mao and Zhenhao Duan (2006) A thermodynamic model for calculating nitrogen solubility, gas phase composition and density of the H2O-N2-NaCl system. Fluid Phase Equilibria, 248 (2): 103-114
  273-400 K, 1-600 bar and 0-6 mol/kg
  http://dx.doi.org/10.1016/j.fluid.2006.07.020
  http://www.geochem-model.org/wp-content/uploads/2009/09/46-FPE_248_103.pdf"
  //redeclare function extends solubility_N2_pTX
    extends partial_solubility_pTX;

  /*  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_H2O";
  input Modelica.SIunits.MolarMass MM[:] "molar masses of components";
  output Modelica.SIunits.MassFraction c_gas "gas concentration in kg_gas/kg_H2O";
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

    Modelica.SIunits.MolarMass M_H2O = MM_vec[end];
    Partial_Units.Molality molalities[size(X,1)]=massFractionsToMolalities(X,MM_vec);
    Partial_Units.Molality m_Cl = molalities[NaCl] + molalities[KCl] + 2*molalities[MgCl2] + 2*molalities[CaCl2];
    Partial_Units.Molality m_Na = molalities[NaCl];
    Partial_Units.Molality m_K = molalities[KCl];
    Partial_Units.Molality m_Ca = molalities[CaCl2];
    Partial_Units.Molality m_Mg = molalities[MgCl2];
    Partial_Units.Molality m_SO4 = 0 "TODO";

    Modelica.SIunits.Pressure p_H2O=Modelica.Media.Water.WaterIF97_base.saturationPressure(T);
  //  Pressure_bar p_H2O_bar=Modelica.SIunits.Conversions.to_bar(p_sat_H2O_Duan2003(T));
    Partial_Units.Pressure_bar p_bar=Modelica.SIunits.Conversions.to_bar(p);
    Modelica.SIunits.MassFraction X_NaCl = molalities[NaCl]*M_H2O
      "mole fraction of NaCl in liquid phase";
    Modelica.SIunits.MolarVolume v_l_H2O=M_H2O/Modelica.Media.Water.WaterIF97_base.density_pT(p,T);
    Real phi_H2O = fugacity_H2O_Duan2006(p,T);
    final constant Real R(final unit="bar.cm3/(mol.K)") = 83.14472
      "Molar gas constant";
  /*  Real y_H20 = (1-2*X_NaCl) * p_H2O/(phi_H2O*p) * exp(v_l_H2O*(p-p_H2O)/(Modelica.Constants.R*T)) 
    "equ. 4 (gamma_H2O=1";
  Real y_N2 = 1-y_H20 "mole fraction of CO2 in vapor phase";*/
    Real y_N2 = p_gas/p "mole fraction of N2 in vapor phase";
    Real phi_N2 = fugacity_N2_Duan2006(p,T);
    Real mu_l0_N2_RT = Par_N2_Duan2006(p,T,mu_l0_N2_RT_c);
    Real lambda_N2_Na = Par_N2_Duan2006(p,T,lambda_N2_Na_c);
    Real xi_N2_NaCl = Par_N2_Duan2006(p,T,xi_N2_NaCl_c);
  algorithm
  // Modelica.Utilities.Streams.print("mola_N2("+String(p_gas)+","+String(T-273.16)+") (solubility_N2_pTX_Duan2006)");

   if outOfRangeMode==1 then
     if not ignoreLimitN2_T and (273>T or T>400) then
        Modelica.Utilities.Streams.print("T="+String(T)+" K, but N2 solubility calculation is only valid for temperatures between 0 and 127°C (Partial_Gas_Data.solubility_N2_pTX_Duan2006)");
     end if;
     if (p<1e5 or p>600e5) then
        Modelica.Utilities.Streams.print("p="+String(p/1e5)+" bar, but N2 solubility calculation only valid for pressures between 1 and 600 bar (Partial_Gas_Data.solubility_N2_pTX_Duan2006)");
     end if;
     if molalities[NaCl]>6 then
       Modelica.Utilities.Streams.print("mola[NaCl]="+String(molalities[NaCl])+" mol/kg, but N2 solubility calculation only valid for salinities up to 6 mol/kg (Partial_Gas_Data.solubility_N2_pTX_Duan2006)");
     end if;
   elseif outOfRangeMode==2 then
       assert(ignoreLimitN2_T or (273<=T and T<=400), "T="+String(T)+" K, but N2 solubility calculation is only valid for temperatures between 0 and 127°C");
       assert(p>=1e5 and p<=600e5, "p="+String(p/1e5)+" bar, but N2 solubility calculation only valid for pressures between 1 and 600 bar");
       assert(molalities[NaCl]<6, "mola[NaCl]="+String(molalities[NaCl])+" mol/kg, but N2 solubility calculation only valid for salinities up to 6 mol/kg");
   end if;

    //equ. 9
      solu := y_N2*p_bar*phi_N2 * exp(-mu_l0_N2_RT - 2*lambda_N2_Na*(m_Na + m_K + 2*m_Ca + 2*m_Mg) - xi_N2_NaCl*(m_Cl+2*m_SO4)*(m_Na + m_K + 2*m_Ca + 2*m_Mg) - 4*0.0371*m_SO4);
  //    solu := max(0, solu) "algorithm can return negative values";
  //  solu := p_H2O;
  //  solu := 0;
  //  c_gas:=solu*M_N2 "kg_gas / kg_H2O";
    X_gas :=solu*M_N2*X[end];
    //    Modelica.Utilities.Streams.print("mola_N2("+String(p_gas)+","+String(T-273.16)+","+String(X[1])+")="+String(c_gas)+" (solubility_N2_pTX_Duan2006)");
  //    Modelica.Utilities.Streams.print("mola_N2("+String(p_gas)+","+String(T-273.16)+")="+String(c_gas)+"->k="+String(c_gas/max(1,p_gas))+" (solubility_N2_pTX_Duan2006)");
  end solubility_N2_pTX_Duan2006;

  function solubility_N2_pTX_Harting "..."
  //from Harting1981, considering only NaCl
    extends partial_solubility_pTX;

  /*  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_H2O";
  input Modelica.SIunits.MolarMass MM[:] "molar masses of components";
  output Modelica.SIunits.MassFraction c_gas "gas concentration in kg_gas/kg_H2O";
*/
  protected
    Molality molalities[size(X,1)]=massFractionsToMolalities(X,MM_vec);
    Modelica.SIunits.Temp_C T_C = Modelica.SIunits.Conversions.to_degC(T);
    Real L_0=.252 "N2-Löslichkeit in H2O bei 25 atm 75°C";
    Real L_rel_p "pressure influence";
    Real L_rel_c "salinity influence";
    Real L_rel_T "Temperature influence";
  //  Real c = X[1]/Salt_Data.M_NaCl/X[end];
    Real c = sum(molalities[1:2])+sum(molalities[3:5])*1.8;
    Real p_atm = p_gas/101325;
  algorithm
  // Modelica.Utilities.Streams.print("mola_N2("+String(p_gas)+","+String(T-273.16)+") (solubility_N2_pTX_Duan2006)");

  //  assert(273<=T and T<=400, "T="+String(T)+" K, but N2 solubility calculation is only valid for temperatures between 0 and 127°C");
   if 273>T or  T>400 then
      Modelica.Utilities.Streams.print("T="+String(T)+" K, but N2 solubility calculation is only valid for temperatures between 0 and 127°C (Partial_Gas_Data.solubility_N2_pTX_Duan2006())");
    end if; /**/
    assert(p>=1e5 and p<=600e5, "p="+String(p/1e5)+" bar, but N2 solubility calculation only valid for pressures between 1 and 600 bar");
  //  assert(molalities[NaCl]<6, "mola[NaCl]="+String(molalities[NaCl])+" mol/kg, but N2 solubility calculation only valid for salinities up to 6 mol/kg");
    if molalities[NaCl]>6 then
      Modelica.Utilities.Streams.print("mola[NaCl]="+String(molalities[NaCl])+" mol/kg, but N2 solubility calculation only valid for salinities up to 6 mol/kg");
    end if;

  //page 19
    L_rel_p :=(0.2569*p_atm - 2.471e-4*p_atm^2 + 1.617e-7*p_atm^3)/100*15.9474650632155 "N2";

    L_rel_c :=exp(-.315*c +.01452 *c^2) "S. 15";

    L_rel_T := -0.0000003493*(T_C-75)^3 + 0.0001054*(T_C-75)^2 - 0.000293*(T_C-75) + 1.01
      "fitted with Excel";
    solu :=L_0*L_rel_p*L_rel_c*L_rel_T "l/kg_H2O";

    X_gas :=solu/22.4*M_N2*X[end];
  //      Modelica.Utilities.Streams.print("mola_N2("+String(p_gas)+"Pa,"+String(T-273.16)+"°C,"+String(molalities[1])+")="+String(solu)+" (solubility_N2_pTX_Harting)");
  //    Modelica.Utilities.Streams.print("mola_N2("+String(p_gas)+","+String(T-273.16)+")="+String(c_gas)+"->k="+String(c_gas/max(1,p_gas))+" (solubility_N2_pTX_Duan2006)");
  end solubility_N2_pTX_Harting;

  function Par_N2_Duan2006 "Mao,Duan(2006)"
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input Real[:] c;
    output Real Par;
  protected
    Partial_Units.Pressure_bar p_bar=Modelica.SIunits.Conversions.to_bar(p);
  algorithm
    //eq. 7
    Par :=  c[1] + c[2]*T + c[3]/T + c[4]*T^2 + c[5]/T^2 + c[6]*p_bar + c[7]*p_bar*T + c[8]*p_bar/T + c[9]*p_bar^2/T;
  end Par_N2_Duan2006;

  function fugacity_N2_Duan2006 "Nullstellensuche mit EOS aus Duan2006"
    //doi:10.1016/j.?uid.2006.07.020
    //Shide Mao, Zhenhao Duan:A thermodynamic model for calculating nitrogen solubility, gas phase composition and density of the N2?H2O?NaCl system
    extends partial_fugacity_pTX;
  protected
    Modelica.SIunits.SpecificVolume V_neu=.024 "Startwert";
    Modelica.SIunits.SpecificVolume V=0;
    Real a[:]= {3.75504388E-02,
               -1.08730273E+04,
                1.10964861E+06,
                5.41589372E-04,
                1.12094559E+02,
               -5.92191393E+03,
                4.37200027E-06,
                4.95790731E-01,
               -1.64902948E+02,
               -7.07442825E-08,
                9.65727297E-03,
                4.87945175E-01,
                1.62257402E+04,
                8.99000000E-03};
  //  Real beta= a[14];
  //  Real gamma = a[15];
    Real sigma=3.63;
    Modelica.SIunits.Temperature epsilon=101;
    Real T_m = 154*T/epsilon;
    Real B = a[1]+a[2]/T_m^2+a[3]/T_m^3;
    Real C = a[4]+a[5]/T_m^2+a[6]/T_m^3;
    Real D = a[7]+a[8]/T_m^2+a[9]/T_m^3;
    Real E = a[10]+a[11]/T_m^2+a[12]/T_m^3;
    Real F = a[13]/T_m^3;

    Real P_m = 3.0626*sigma^3*p_bar/epsilon;
    Real G;
    Real ln_phi;
    Real Z;
    Real V_m;
    Integer z=0 "only a counter to avoid getting caught in the iteration loop";
    Real d=.7 "dampening factor 0=no dampening, 1=no progress";
  algorithm
  /*  V := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
      function fugacity_N2_Duan2006_res(p=p,T=T,a=a,B=B,C=C,D=D,E=E,F=F),
      1e-6,
      Modelica.Constants.R*T/p,
      1e-6);*/

  //iterative solution
  while abs(V-V_neu)>1e-8 loop
  //    V:=V_neu;
      V:=if z<5 then V_neu else (1-d)*V_neu+d*V;
      V_m :=V*1e6/(1e3*(sigma/3.691)^3);
      Z := 1 + B/V_m + C/V_m^2 + D/V_m^4 + E/V_m^5 + F/V_m^2*(1 + a[14]/V_m^2)*exp(-a[14]/V_m^2);
      V_neu :=Z/p*Modelica.Constants.R*T;
  //    Modelica.Utilities.Streams.print("V("+String(z)+")="+String(V_neu));
      z:=z + 1;
      assert(z<1000," Reached maximum number of iterations for fugacity calculation.(fugacity_N2_Duan2006)");
    end while;

    V_m := 1e3*V/(sigma/3.691)^3 "m³/mol -> dm³/mol";
  //  V_m := V/(Modelica.Constants.R*T_c/P_c);
  //  Z := 1 + B/V_m + C/V_m^2 + D/V_m^4 + E/V_m^5 + F/V_m^2*(1 + a[14]/V_m^2)*exp(-a[14]/V_m^2);
    Z := 1 + (a[1]+a[2]/T_m^2+a[3]/T_m^3)/V_m
           + (a[4]+a[5]/T_m^2+a[6]/T_m^3)/V_m^2
           + (a[7]+a[8]/T_m^2+a[9]/T_m^3)/V_m^4
           + (a[10]+a[11]/T_m^2+a[12]/T_m^3)/V_m^5
           + a[13]/T_m^3/V_m^2*(1 + a[14]/V_m^2)*exp(-a[14]/V_m^2);
    /*  ln_phi := Z-1-log(Z) + B/V_r + C/(2*V_r^2) + D/(4*V_r^4) + E/(5*V_r^5) + G;
  phi := exp(ln_phi) "fugacity coefficient";*/
    G :=a[13]/T_m^3/(2*a[14])*(2 - (2 + a[14]/V_m^2)*exp(-a[14]/V_m^2));
    phi := exp(Z-1 + B/V_m + C/(2*V_m^2) + D/(4*V_m^4) + E/(5*V_m^5) + G)/Z
      "fugacity coefficient";
  /*  phi:=exp(Z-1 + (a[1]+a[2]/T_m^2+a[3]/T_m^3)/V_m 
         + (a[4]+a[5]/T_m^2+a[6]/T_m^3)/(2*V_m^2)
         + (a[7]+a[8]/T_m^2+a[9]/T_m^3)/(4*V_m^4)
         + (a[10]+a[11]/T_m^2+a[12]/T_m^3)/(5*V_m^5)
         + a[13]/(2*T_m^3*a[14])*(2-(2-a[14]/V_m^2)*exp(-a[14]/V_m^2)))/Z;*/
  //  Modelica.Utilities.Streams.print("z="+String(z));
  //  Modelica.Utilities.Streams.print("V="+String(V));
  //  PowerPlant.Components.PipeStuff.print_msg(phi,"phi_N2=");
  end fugacity_N2_Duan2006;

  function fugacity_H2O_Duan2006 "Calculation of fugacity coefficient according to (Mao&Duan 2006 'A thermodynamic model for calculating nitrogen solubility, gasphase
composition and density of the N2?H2O?NaCl system')"
    extends partial_fugacity_pTX;

  protected
    Partial_Units.Pressure_bar P_1;
    Real[:] a = {1.86357885E-03,
                 1.17332094E-02,
                 7.82682497E-07,
                -1.15662779E-05,
                -3.13619739,
                -1.29464029E-03};
  algorithm
    phi := exp(a[1] + a[2]*p_bar + a[3]*p_bar^2 + a[4]*p_bar*T + a[5]*p_bar/T + a[6]*p_bar^2/T) "equ. 5";
  end fugacity_H2O_Duan2006;

  function degassingPressure_N2_Duan2006
    "calculates degassing pressure from concentration of dissolved gas"
    extends partial_degassingPressure_pTX;
  /*protected 
    Modelica.SIunits.MassFraction solu_soll=X[end-3];*/

  algorithm
  //  Modelica.Utilities.Streams.print("degassingPressure_N2_Duan2006("+String(p)+","+String(T)+","+String(X[end-2])+")="+String(p_gas)+" (degassingPressure_N2_Duan2006)");
  /*  p_gas := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
      function solubility_N2_Duan2006_res(p=p,T=T,X=X,MM_vec=MM_vec),
      0,
      1e8,
      1e-8);*/

    p_gas := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
        function solubility_res(
          solufun=function solubility_N2_pTX_Duan2006(),p=p,T=T,X=X,MM_vec=MM_vec,
          c_gas=X[end-2]),
        0,
        1e10,
        1e-8);
    assert(size(X,1)==9,"Wenn es nicht 9 Komponenten gibt, dann haut die Konzentrationszuordnung hier nich hin(degassingPressure_N2_Duan2006)");
  //  Modelica.Utilities.Streams.print("p_sat_N2("+String(X[end-2])+")="+String(p_gas)+" (degassingPressure_N2_Duan2006)");

  end degassingPressure_N2_Duan2006;

  function solubility_CH4_pTX_Duan2006 "Duan ZH, Mao SD. (2006) A thermodynamic model for calculating methane solubility, density and gas phase composition of methane-bearing aqueous fluids from 273 to 523 K and from 1 to 2000 bar. Geochimica et Cosmochimica Acta, 70 (13): 3369-3386. 
  http://geochem-model.org/Publications/43-GCA_2006_3369.pdf
  http://dx.doi.org/10.1016/j.gca.2006.03.018TODO Umrechnung andere Salz in NaCl"
    extends partial_solubility_pTX;
  //  output Modelica.SIunits.MassFraction c_gas "gas concentration in kg_gas/kg_H2O";

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

    Modelica.SIunits.MolarMass M_H2O = MM_vec[end];
    Partial_Units.Molality molalities[size(X,1)];
    Partial_Units.Molality molefractions[size(X,1)];
    Partial_Units.Molality m_Cl;
    Partial_Units.Molality m_Na;
    Partial_Units.Molality m_K;
    Partial_Units.Molality m_Ca;
    Partial_Units.Molality m_Mg;
    Partial_Units.Molality m_SO4;
  //  Real X_NaCl = molalities[NaCl]*M_H2O "mole fraction of NaCl in liquid phase";

  //  Modelica.SIunits.Pressure p_H2O = Modelica.Media.Water.WaterIF97_base.saturationPressure(T);
    Modelica.SIunits.Pressure p_H2O = p_sat_H2O_Duan2003(T);
    Partial_Units.Pressure_bar p_bar=Modelica.SIunits.Conversions.to_bar(p);
    Modelica.SIunits.MolarVolume v_l_H2O=M_H2O/Modelica.Media.Water.WaterIF97_base.density_pT(p,T);
    Real phi_H2O = fugacity_H2O_Duan2006b(p,T);
    Real y_H20 "mole fraction of H2O in vapor phase";
    Real y_CH4;
    Real phi_CH4 = fugacity_CH4_Duan1992(p,T);
    Real mu_l0_CH4_RT = Par_CH4_Duan2006(p,T,mu_l0_CH4_RT_c);
    Real lambda_CH4_Na = Par_CH4_Duan2006(p,T,lambda_CH4_Na_c);
    Real xi_CH4_NaCl = Par_CH4_Duan2006(p,T,xi_CH4_NaCl_c);
  algorithm

   if outOfRangeMode==1 then
     if (273>T or T>273+250) then
        Modelica.Utilities.Streams.print("T="+String(T)+" K, but  CH4 solubility  calculation is only valid for temperatures between 0 and 250°C (Partial_Gas_Data.solubility_CH4_pTX_Duan1992)");
     end if;
     if (p>1600e5) then
        Modelica.Utilities.Streams.print("p="+String(p/1e5)+" bar, but CH4 fugacity calculation only valid for pressures between 1 and 1600 bar (Partial_Gas_Data.solubility_CH4_pTX_Duan1992)");
     end if;
   elseif outOfRangeMode==2 then
     assert(273.15<=T and T<=273+250, "T="+String(T-273.15)+"°C, but CH4 solubility calculation is only valid for temperatures between 0 and 250°C (solubility_CH4_pTX_Duan1992)");
     assert(p<=1600e5, "p="+String(p/1e5)+"bar, but CH4 fugacity calculation only valid for pressures up to 1600 bar (solubility_CH4_pTX_Duan1992)");
   end if;

  //  (molefractions,molalities):=massFractionsToMoleFractions(X, MM);
    molalities:=massFractionsToMolalities(X, MM_vec);
  // Modelica.Utilities.Streams.print("molefractions[NaCl]="+String(molefractions[NaCl])+" (Partial_Gas_Data.solubility_CH4_pTX_Duan1992)");
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
    y_CH4 :=p_gas/p;

    //equ. 10
      solu := y_CH4*phi_CH4*p_bar*exp(-mu_l0_CH4_RT - 2*lambda_CH4_Na*(m_Na + m_K + 2*m_Ca + 2*m_Mg) - xi_CH4_NaCl*(m_Na + m_K + 2*m_Ca + 2*m_Mg)*(m_Cl + 2*m_SO4) - 4*0.0332*m_SO4);

  //  solu := max(0, solu) "algorithm can return negative values";
  //  solu := p_H2O;
  //  solu := 0;
  //  c_gas:=solu*M_CH4 "kg_gas / kg_H2O";
    X_gas :=solu*M_CH4*X[end];
    //    Modelica.Utilities.Streams.print("mola_CH4("+String(p_gas)+","+String(T-273.16)+")="+String(c_gas)+"->k="+String(c_gas/max(1,p_gas))+" (solubility_CH4_pTX_Duan2006)");
  end solubility_CH4_pTX_Duan2006;

  function solubility_CH4_pTX_Harting "..."
  //from Harting1981, considering only NaCl
    extends partial_solubility_pTX;

  /*  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_H2O";
  input Modelica.SIunits.MolarMass MM[:] "molar masses of components";
  output Modelica.SIunits.MassFraction c_gas "gas concentration in kg_gas/kg_H2O";
*/
  protected
    Molality molalities[size(X,1)]=massFractionsToMolalities(X,MM_vec);
    Modelica.SIunits.Temp_C T_C = Modelica.SIunits.Conversions.to_degC(T);
    Real L_0=.454 "CH4-Löslichkeit in H2O bei 25 atm 75°C";
    Real L_rel_p "pressure influence";
    Real L_rel_c "salinity influence";
    Real L_rel_T "Temperature influence";
  //  Real c = X[1]/Salt_Data.M_NaCl/X[end];
    Real c = sum(molalities[1:2])+sum(molalities[3:5])*1.8;
    Real p_atm = p_gas/101325;
  algorithm
  // Modelica.Utilities.Streams.print("mola_CH4("+String(p_gas)+","+String(T-273.16)+") (solubility_CH4_pTX_Duan2006)");

  //  assert(273<=T and T<=400, "T="+String(T)+" K, but N2 solubility calculation is only valid for temperatures between 0 and 127°C");
   if 273>T or  T>400 then
      Modelica.Utilities.Streams.print("T="+String(T)+" K, but N2 solubility calculation is only valid for temperatures between 0 and 127°C (Partial_Gas_Data.solubility_N2_pTX_Duan2006())");
    end if; /**/
    assert(p>=1e5 and p<=600e5, "p="+String(p/1e5)+" bar, but N2 solubility calculation only valid for pressures between 1 and 600 bar");
  //  assert(molalities[NaCl]<6, "mola[NaCl]="+String(molalities[NaCl])+" mol/kg, but N2 solubility calculation only valid for salinities up to 6 mol/kg");
    if molalities[NaCl]>6 then
      Modelica.Utilities.Streams.print("mola[NaCl]="+String(molalities[NaCl])+" mol/kg, but N2 solubility calculation only valid for salinities up to 6 mol/kg");
    end if;

  //page 19
    L_rel_p :=(.4009*p_atm - 7.454e-4*p_atm^2 + 5.985e-7*p_atm^3)/100*10.4537157651017 "CH4";

    L_rel_c :=exp(-.315*c +.01452 *c^2) "S. 15";

    L_rel_T := -0.0000003493*(T_C-75)^3 + 0.0001054*(T_C-75)^2 - 0.000293*(T_C-75) + 1.01
      "fitted with Excel";
    solu :=L_0*L_rel_p*L_rel_c*L_rel_T "l/kg_H2O";

    X_gas :=solu/22.4*M_CH4*X[end];
  /*  Modelica.Utilities.Streams.print("L_rel_p="+String(L_rel_p)+" (solubility_CH4_pTX_Harting)");
  Modelica.Utilities.Streams.print("L_rel_c(c="+String(c)+")="+String(L_rel_c)+" (solubility_CH4_pTX_Harting)");
  Modelica.Utilities.Streams.print("L_rel_T="+String(L_rel_T)+" (solubility_CH4_pTX_Harting)");
  Modelica.Utilities.Streams.print("mola_CH4("+String(p_gas)+","+String(T-273.16)+","+String(molalities[1])+")="+String(solu)+" (solubility_CH4_pTX_Harting)");
  */
  end solubility_CH4_pTX_Harting;

  function Par_CH4_Duan2006 "Duan,Sun(2003)"
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input Real[:] c;
    output Real Par;
  protected
    Partial_Units.Pressure_bar p_bar=Modelica.SIunits.Conversions.to_bar(p);
  algorithm
    //eq. 7
    Par :=  c[1] + c[2]*T + c[3]/T + c[4]*T^2 + c[5]/T^2 + c[6]*p_bar + c[7]*p_bar*T + c[8]*p_bar/T + c[9]*p_bar/T^2+c[10]*p_bar^2*T;
  end Par_CH4_Duan2006;

  function fugacity_H2O_Duan2006b
    "Calculation of fugacity coefficient according to Duan ZH, Mao SD. (2006)"
   /*A thermodynamic model for calculating methane solubility, density and gas phase composition of methane-bearing aqueous fluids from 273 to 523 K and from 1 to 2000 bar. Geochimica et Cosmochimica Acta, 70 (13): 3369-3386.
 */
    extends partial_fugacity_pTX;

  protected
    Partial_Units.Pressure_bar P_1;
    Real[:] a = {-1.42006707E-02,
                  1.08369910E-02,
                 -1.59213160E-06,
                 -1.10804676E-05,
                 -3.14287155E+00,
                  1.06338095E-03};
  algorithm
    phi := exp(a[1] + a[2]*p_bar + a[3]*p_bar^2 + a[4]*p_bar*T + a[5]*p_bar/T + a[6]*p_bar^2/T) "equ. 6";
  end fugacity_H2O_Duan2006b;

  function fugacity_CH4_Duan1992 "Nullstellensuche mit EOS aus Duan1992"
    extends partial_fugacity_pTX;
  protected
    Modelica.SIunits.SpecificVolume V_neu=.024 "Startwert";
    Modelica.SIunits.SpecificVolume V=0;
    Real a[:]= {8.72553928E-02,
            -7.52599476E-01,
            3.75419887E-01,
            1.07291342E-02,
            5.49626360E-03,
            -1.84772802E-02,
            3.18993183E-04,
            2.11079375E-04,
            2.01682801E-05,
            -1.65606189E-05,
            1.19614546E-04,
            -1.08087289E-04,
            4.48262295E-02,
            7.53970000E-01,
            7.71670000E-02};
    Real alpha=a[13];
    Real beta=a[14];
    Real gamma = a[15];
    Modelica.SIunits.Temperature T_c = 190.6;
    Modelica.SIunits.Pressure P_c = 46.41e5;
    Real P_r = p/P_c;
    Real T_r = T/T_c;
     Real B = a[1]+a[2]/T_r^2+a[3]/T_r^3;
     Real C = a[4]+a[5]/T_r^2+a[6]/T_r^3;
     Real D = a[7]+a[8]/T_r^2+a[9]/T_r^3;
     Real E = a[10]+a[11]/T_r^2+a[12]/T_r^3;
     Real F = alpha/T_r^3;
    Real ln_phi;
    Real Z;
    Real V_r;
    Real G;
    Integer z=0 "only a counter to avoid getting caught in the iteration loop";
    Real d=.7 " dampening factor 0=no dampening, 1=no progress";
  algorithm

  /*  V := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
      function fugacity_CH4_Duan1992_res(p=p,T=T,a=a,B=B,C=C,D=D,E=E,F=F,P_c=P_c,T_c=T_c),
      1e-6,
      1e1,
      1e-6);*/
  //  PowerPlant.Components.PipeStuff.print_msg(,"V_end=")

   while abs(V-V_neu)>1e-8 loop
  //    V:=(1-d)*V_neu+d*V;
      V:=if z<5 then V_neu else (1-d)*V_neu+d*V "dampened";
  //    d:= min(0.95,z/100);    V:= (1-d)*V_neu+d*V;
      V_r:=V/(Modelica.Constants.R*T_c/P_c);
      G := F/(2*gamma)*(beta+1-(beta+1+gamma/V_r^2)*exp(-gamma/V_r^2));
      Z:= 1 + B/V_r + C/V_r^2 + D/V_r^4 + E/V_r^5 + F/V_r^2*(beta + gamma/V_r^2)*exp(-gamma/V_r^2);
      V_neu :=Z/p*Modelica.Constants.R*T;
  //    Modelica.Utilities.Streams.print("V("+String(z)+")="+String(V_neu));
      z:=z + 1;
      assert(z<1000," Reached maximum number of iterations for CH4 fugacity calculation.(fugacity_CH4_Duan1992)");
    end while;

  /*  V_r := V/(Modelica.Constants.R*T_c/P_c);
  Z := 1 + B/V_r + C/V_r^2 + D/V_r^4 + E/V_r^5 + F/V_r^2*(beta + gamma/V_r^2)*exp(-gamma/V_r^2);*/
  /*  ln_phi := Z-1-log(Z) + B/V_r + C/(2*V_r^2) + D/(4*V_r^4) + E/(5*V_r^5) + G;
  phi := exp(ln_phi) "fugacity coefficient";*/
  //  G :=F/(2*gamma)*(beta + 1 - (beta + 1 + gamma/V_r^2)*exp(-gamma/V_r^2));
    phi := exp(Z-1 + B/V_r + C/(2*V_r^2) + D/(4*V_r^4) + E/(5*V_r^5) + G)/Z
      "fugacity coefficient";
  //  PowerPlant.Components.PipeStuff.print_msg(phi,"phi_CH4=");
  end fugacity_CH4_Duan1992;

  function degassingPressure_CH4_Duan2006
    "calculates degassing pressure from concentration of dissolved gas"
    extends partial_degassingPressure_pTX;
  /*protected 
    Modelica.SIunits.MassFraction solu_soll=X[end-3];*/

  algorithm
  //  Modelica.Utilities.Streams.print("p_sat_CH4("+String(X[end-1])+") (degassingPressure_CH4_Duan2006)");
  /*  p_gas := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
      function solubility_CH4_Duan2006_res(p=p,T=T,X=X,MM_vec=MM_vec),
      0,
      1e8,
      1e-8);
*/

    p_gas := Modelica.Math.Nonlinear.solveOneNonlinearEquation(
        function solubility_res(
          solufun=function solubility_CH4_pTX_Duan2006(),p=p,T=T,X=X,MM_vec=MM_vec,
          c_gas=X[end-1]),
        0,
        1e10,
        1e-8);

  //  Modelica.Utilities.Streams.print("p_sat_CH4("+String(X[end-1])+")="+String(p_gas)+" (degassingPressure_CH4_Duan2006)");

  end degassingPressure_CH4_Duan2006;

  partial function partial_solutionEnthalpy
    "template for calculation of solution enthalpy"
    input Modelica.SIunits.Temp_K T;
    output Modelica.SIunits.SpecificEnthalpy Delta_h_solution;
  protected
    Modelica.SIunits.Pressure p_H2O;

  end partial_solutionEnthalpy;

  function solutionEnthalpy_CO2_Duan2003
    "calculation of solution enthalpy of CO2+H2O according to Duan(2003) equ. 8"
    extends partial_solutionEnthalpy;
  /*  input Modelica.SIunits.Temp_K T;
  output Modelica.SIunits.SpecificEnthalpy Delta_h_solution;
 //TODO: Was passiert wenn weniger CO2 zur Verfügung steht als gelöst werden kann?

  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
  output Modelica.SIunits.SpecificEnthalpy h;

  output Modelica.SIunits.SpecificEnthalpy Delta_h_solution 
    "solution enthalpy per mole";
  Modelica.SIunits.SpecificEnthalpy h_l;
  Modelica.SIunits.SpecificEnthalpy h_g;
  Modelica.SIunits.SpecificEnthalpy h_brine;
  Modelica.SIunits.SpecificEnthalpy h_CO2_dissoluted;
  Modelica.SIunits.SpecificEnthalpy h_H2O_g;
  Modelica.SIunits.SpecificEnthalpy h_CO2_g;
  Modelica.SIunits.SpecificEnthalpy Delta_h_mix "mixing enthalpy in gas phase";*/
  protected
    Modelica.SIunits.Pressure p_H2O;

    Real[:] c = { 28.9447706,
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
  algorithm
  /*  h_brine := Brine_Driesner.specificEnthalpy_pTX(p,T,X);

  h_CO2_dissoluted := 0;
*/
    p_H2O := 0 "Modelica.Media.Water.WaterIF97_base.saturationPressure(T)";
    Delta_h_solution := Modelica.Constants.R*T^2 * (c[2]-c[3]/T^2+2*c[4]*T+c[5]/(630-T)^2
      + c[7]*p_H2O/T - c[8]*p_H2O/T^2 + c[9]*p_H2O/(630-T)^2
      + 2*c[10]*p_H2O^2/(630-T)^3) "eq. 8 Duan 2003";

  /*  h_H2O_g := 0 
    "Modelica.Media.Water.WaterIF97_base.dewEnthalpy(Modelica.Media.Water.WaterIF97_base.setSat_p(p))";

  h_CO2_g := 0 
    "Modelica.Media.IdealGases.SingleGases.CO2.specificEnthalpy(Modelica.Media.IdealGases.SingleGases.CO2.specificEnthalpy(p,T))";

  Delta_h_mix :=0 "ideal mixing";

  h_l := h_brine + h_CO2_dissoluted + Delta_h_solution 
    "specific enthalpy of liquid phase";
  h_g := h_H2O_g + h_CO2_g + Delta_h_mix "specific enthalpy of gas phase";

  h :=((1 - GVF)*d_l*h_l + GVF*d_g*h_g)/d "summarize according to mass ratio";*/
  end solutionEnthalpy_CO2_Duan2003;

  function solutionEnthalpy_N2 "calculation of solution enthalpy of N2+H2O"
    extends partial_solutionEnthalpy;
  algorithm
    Delta_h_solution := 0 "TODO";
  end solutionEnthalpy_N2;

  function solutionEnthalpy_CH4 "calculation of solution enthalpy of CH4+H2O"
    extends partial_solutionEnthalpy;
  algorithm
    Delta_h_solution := 0 "TODO";
  end solutionEnthalpy_CH4;

end Partial_Gas_Data;
