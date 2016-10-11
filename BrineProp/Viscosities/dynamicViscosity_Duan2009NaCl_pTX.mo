within BrineProp.Viscosities;
function dynamicViscosity_Duan2009NaCl_pTX
  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X[:] "mass fraction m_NaCl/m_Sol";
  output SI.DynamicViscosity eta;
protected
  Types.Molality mola=X[1]/(BrineProp.SaltData.M_NaCl*(1 - X[1]))
  "molality b (mol_NaCl/kg_H2O)";
  SI.Temp_C T_C = SI.Conversions.to_degC(T);
  Types.Pressure_bar p_bar=SI.Conversions.to_bar(p);

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
  Modelica.Media.Water.WaterIF97_pT.ThermodynamicState state_H2O;
algorithm
  if debugmode then
    print("p="+String(p)+" Pa, T_K"+String(T)+"X[1]="+String(X[1])+" K (Brine.Viscosities.dynamicViscosity_Duan2009NaClKCl_pTX)");
  end if;
  assert(T_C>=0 and T_C<=350, "Temperature is "+String(SI.Conversions.to_degC(T))+"degC, but must be between 10 and 350degC");
  assert(p_bar>=1 and p_bar<=1000, "Pressure must be between 1 and 500 bar");
  assert(mola>=0 and mola<=6, "Molality must be between 0.25 and 5 mol/kg");
  //factors
  A_NaCl := a_0_NaCl + a_1_NaCl*T + a_2_NaCl*T^2;
  B_NaCl := b_0_NaCl + b_1_NaCl*T + b_2_NaCl*T^2;
  C_NaCl := c_0_NaCl + c_1_NaCl*T;
  eta_relative := exp(A_NaCl*mola + B_NaCl*mola^2 + C_NaCl*mola^3);

  //viscosity calculation
  state_H2O := Modelica.Media.Water.WaterIF97_pT.setState_pTX(p, T, X);
  eta_H2O := Modelica.Media.Water.WaterIF97_pT.dynamicViscosity(state_H2O);
  eta := eta_relative * eta_H2O;
//  print("mola: "+String(mola));
//  print("eta_relative: "+String(eta_relative));
//  print("eta: "+String(eta));

end dynamicViscosity_Duan2009NaCl_pTX;
