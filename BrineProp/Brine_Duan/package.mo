within BrineProp;
package Brine_Duan "NaCl solution using Duan density"
  extends BrineProp.PartialBrine_MultiSalt_1Phase(
                                        redeclare package Salt_data =
      BrineProp.SaltData_Duan);

  constant Salt_data.SaltConstants salt = Salt_data.saltConstants[NaCl];


  redeclare function extends dynamicViscosity_pTX
   //  constant Real M_NaCl=0.058443 "molar mass in [kg/mol]";
    /*  public 
     constant SI.MolarMass M_NaCl = salt.M_salt; 
      "[kg/mol]";*/
protected
    Molality mola = X[1]/(salt.M_salt*(1-X[1])) "molality b (mol_NaCl/kg_H2O)";
    SI.Temp_C T_C = SI.Conversions.to_degC(T_K);
    Pressure_bar p_bar= SI.Conversions.to_bar(p_Pa);

    Real a_0_NaCl = -0.21319213;
    Real a_1_NaCl = +0.13651589E-2;
    Real a_2_NaCl = - 0.12191756E-5;

    Real b_0_NaCl = +0.69161945E-1;
    Real b_1_NaCl = -0.27292263E-3;
    Real b_2_NaCl = +0.20852448E-6;

    Real c_0_NaCl = -0.25988855E-2;
    Real c_1_NaCl = +0.77989227E-5;

    Real A_NaCl;
    Real B_NaCl;
    Real C_NaCl;

    Real eta_relative;

    SI.DynamicViscosity eta_H2O;
    Modelica.Media.Water.WaterIF97_base.ThermodynamicState state_H2O;
  algorithm
    assert(T_C>=0 and T_C<=300, "Temperature is "+String(SI.Conversions.to_degC(T_K))+"°C, but must be between 10 and 350°C");
    assert(p_bar>=1 and p_bar<=1000, "Pressure must be between 1 and 500 bar");
    assert(mola>=0 and mola<=6, "Molality must be between 0.25 and 5 mol/kg");
    //factors
    A_NaCl := a_0_NaCl + a_1_NaCl*T_K + a_2_NaCl*T_K^2;
    B_NaCl := b_0_NaCl + b_1_NaCl*T_K + b_2_NaCl*T_K^2;
    C_NaCl := c_0_NaCl + c_1_NaCl*T_K;
    eta_relative := exp(A_NaCl*mola + B_NaCl*mola^2 + C_NaCl*mola^3);

    //viscosity calculation
    state_H2O := Modelica.Media.Water.WaterIF97_base.setState_pTX(p_Pa, T_K, X);
    eta_H2O := Modelica.Media.Water.WaterIF97_base.dynamicViscosity(state_H2O);
    eta := eta_relative * eta_H2O;
  //  print("mola: "+String(mola));
  //  print("eta_relative: "+String(eta_relative));
  //  print("eta: "+String(eta));

  end dynamicViscosity_pTX;


  redeclare function extends specificEnthalpy_pTX
  "enthalpy calculation according to Driesner et al: 0-1000°C; 0.1-500MPa"
  algorithm
    h :=Brine_Driesner.specificEnthalpy_pTX(p,T,X);
  //  h :=T;

  end specificEnthalpy_pTX;
end Brine_Duan;
