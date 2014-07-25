within BrineProp;
package BrineDriesner "NaCl solution using Driesner density and enthalpy function"
  extends BrineProp.PartialBrineMultiSaltOnePhase(
      redeclare package Salt_data = BrineProp.SaltDataDuan,
      final saltNames = {"sodium chloride"},
      final MM_salt = {Salt_data.M_NaCl},
      final nM_salt = {Salt_data.nM_NaCl});


  redeclare function extends density_pTX
  "density calculation according to Driesner et al: 10-1000degC; 0.1-500MPa; 0.25-5 mol/kg"
  /*  input SI.Pressure p;
  input SI.Temp_K T;
  input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
  output SI.Density d; */

  /*  constant Real M_NaCl = 0.058443 "molar mass in [kg/mol]";
  constant Real M_H2O = 0.018015 "molar mass in [kg/mol]";
*/
public
    constant Real M_NaCl=BrineProp.SaltData.M_NaCl "molar mass in [kg/mol]";
  //  constant Real M_H2O= BrineProp.M_H2O "molar mass in [kg/mol] TODO";

protected
    PartialUnits.Molality mola=X[1]/M_NaCl "molality b (mol_NaCl/kg_sol)";
    SI.Temp_C T_C = SI.Conversions.to_degC(T);
    SI.Temp_C T_Scale_V;
    PartialUnits.Pressure_bar p_bar=SI.Conversions.to_bar(p);
    Real n_21;
    Real n_22;
    Real n_20;
    Real n_23;
    Real n_2;
    Real n_300;
    Real n_301;
    Real n_302;
    Real n_30;
    Real n_310;
    Real n_311;
    Real n_312;
    Real n_31;
    Real D;
    Real n_11;
    Real n_10;
    Real n_12;
    Real n_1;
    Modelica.Media.Water.WaterIF97_base.ThermodynamicState state_H2O;
    Real x_NaCl;
    SI.MolarMass M_Solution "[kg/mol]";
  algorithm
    p_bar := SI.Conversions.to_bar(p);
    assert(T_C>=0 and T_C<=1000, "Temperature must be between 0 and 1000degC");
    assert(p_bar>=1 and p_bar<=5000, "Pressure must be between 1 and 5000 bar");
  //  assert(mola>=.25 and mola<=5, "Molality must be between 0.25 and 5 mol/kg");

  //Salinity conversion
    if X[1]==0 then
      x_NaCl := 0;
    else
      x_NaCl := 1/(M_NaCl/M_H2O*(1/X[1]-1)+1) "mol fraction";
    end if;
    M_Solution := x_NaCl*M_NaCl + (1-x_NaCl)* M_H2O;

  //CALCULATION OF EQUIVALENT TEMPERATURE_V
    //Calculation n2
    n_21 := -2.6142 - 0.000239092*p_bar;
    n_22 := 0.0356828 + 4.37235E-6*p_bar + 2.0566E-9*p_bar^2;

    n_20 := 1-n_21*sqrt(n_22) "x_NaCl = 0 results in n_2=1";
    n_23 := -n_20 - n_21*sqrt(1+n_22) -0.0370751 + 0.00237723*sqrt(p_bar) + 5.42049E-5*p_bar + 5.84709E-9*p_bar^2 - 5.99373E-13*p_bar^3
    "x_NaCl = 1 is pure NaCl";

    n_2 := n_20 + n_21*sqrt(x_NaCl+n_22) + n_23*x_NaCl;
    //End of Calculation n2

    //Calculation D
    n_300 := 7.60664E6 / (p_bar + 472.051)^2;
    n_301 := -50 - 86.144 * exp(-6.21128E-4 * p_bar);
    n_302 := 294.318 * exp(5.66735E-3 * p_bar);
    n_30 := n_300 * (exp(n_301*x_NaCl) - 1) + n_302*x_NaCl;

    n_310 := -0.0732761 * exp(-2.3772E-3 * p_bar) - 5.2948E-5 * p_bar;
    n_311 := -47.2747 + 24.3653 * exp(-1.25533E-3 * p_bar);
    n_312 := -0.278529 - 0.00081381 * p_bar;
    n_31 := n_310 * exp(n_311 * x_NaCl) + n_312 * x_NaCl;

    D := n_30 * exp(n_31 * T_C);
    //End of Calculation D

    //Calculation n1
    n_11 := -54.2958 - 45.7623*exp(-9.44785E-4*p_bar);
    n_10 := 330.47 + 0.942876*sqrt(p_bar) + 0.0817193*p_bar - 2.47556E-8*p_bar^2 + 3.45052E-10*p_bar^3
    "x_NaCl = 1 is pure NaCl";
    n_12 := -n_10 -n_11 "x_NaCl = 0 results in n_1=0";

    n_1 := n_10 + n_11* (1-x_NaCl) + n_12*(1-x_NaCl)^2;
    //End of Calculation n1

    T_Scale_V := T_C*n_2 + n_1 + D;
    //END OF CALCULATION OF EQUIVALENT TEMPERATURE

    state_H2O := Modelica.Media.Water.WaterIF97_base.setState_pTX(p, SI.Conversions.from_degC(T_Scale_V), fill(0,0));
    d := Modelica.Media.Water.WaterIF97_base.density(state_H2O)*M_Solution/M_H2O;
  //  print("T_Scale_V: "+String(T_Scale_V)+"(density_Driesner_pTX)");
  //  print("Density: "+String(d)+"(density_Driesner_pTX)");
  end density_pTX;


  redeclare function extends specificEnthalpy_pTX
  "enthalpy calculation according to Driesner 2007 et al: 0-1000degC; 0.1-500MPa (doi:10.1016/j.gca.2007.05.026)"
  //Pressure limited to 100 MPa by Modelica Water property function
  /*  input SI.Pressure p;
  input SI.Temp_K T;
  input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
  output SI.SpecificEnthalpy h;*/

  algorithm
    h := BrineProp.SpecificEnthalpies.specificEnthalpy_pTX_Driesner(p,T,X[1]);

  //  print("Brine_Driesner.specificEnthalpy_pTX: "+String(p*1e-5)+"bar."+String(T_Scale_h)+"degC->"+String(h)+" J/kg");
  end specificEnthalpy_pTX;
end BrineDriesner;
