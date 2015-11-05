within BrineProp.SpecificEnthalpies;
function T_Scale_h_Driesner
  "enthalpy calculation according to Driesner 2007 et al: 0-1000degC; 0.1-500MPa (doi:10.1016/j.gca.2007.05.026)"
//Pressure limited to 100 MPa by Modelica Water property function
  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X_NaCl "mass fraction m_NaCl/m_Sol";
  output SI.Temp_K T_Scale_h;
  output Real q_2;
          //public
protected
  constant Real M_NaCl=BrineProp.SaltData.M_NaCl "molar mass in [kg/mol]";
//  constant Real M_H2O =  PartialBrine.M_H2O "molar mass in [kg/mol] TODO";

  constant Types.Pressure_bar p_min=1;
  constant Types.Pressure_bar p_max=1000;
  constant SI.Temp_C T_min=0;
  constant SI.Temp_C T_max=1000;
  SI.Temp_C T_C = SI.Conversions.to_degC(T);
//  SI.Temp_C T_Scale_h;
  Types.Pressure_bar p_bar=SI.Conversions.to_bar(p);
  Real q_21;
  Real q_22;
  Real q_20;
  Real q_23;
//  Real q_2;
  Real q_11;
  Real q_10;
  Real q_12;
  Real q_1;
//  Modelica.Media.Water.WaterIF97_pT.ThermodynamicState state_H2O;
  Real x_NaCl "mol fraction";
  //  SI.MolarMass M_Solution "[kg/mol]";
algorithm
  assert(size(ignoreLimitSalt_p,1)>0,"Wrong length of ignoreLimitSalt_T ("+String(size(ignoreLimitSalt_p,1))+")"); //needed here, because flag vector with fewer than nX_salts elements causes "out of bounds" and is not caught elsewere
  assert(size(ignoreLimitSalt_T,1)>0,"Wrong length of ignoreLimitSalt_p ("+String(size(ignoreLimitSalt_p,1))+")"); //should be in PartialFlags, but asserts can't be in packages

  if AssertLevel>0 then
    assert(ignoreLimitSalt_p[iNaCl] or (p_bar>=p_min and p_bar<=p_max),"Pressure p=" + String(p/1e5) + " bar is out of validity range ["+String(p_min/1e5)+"..."+String(p_max/1e5)+"]bar.\nTo ignore set ignoreLimitSalt_p["+String(iNaCl)+"]=true",aLevel);
    assert(ignoreLimitSalt_T[iNaCl] or (T_C>=T_min and T_C<=T_max),"Temperature T=" + String(T-273.15) + " C out of validity range ["+String(T_min)+"..."+String(T_max)+"]C.\nTo ignore set ignoreLimitSalt_T["+String(iNaCl)+"]=true",aLevel);
//    assert(ignoreLimitSalt_b[1] or (mola>=0.25 and mola<=5),"NaCl-Molality is out of validity range: m[i]=" + String(mola) + " mol/kg.\nTo ignore set ignoreLimitSalt_b["+String(iNaCl)+"]=true\n",aLevel);
  end if;

//Salinity conversion
//  if X[1]==0 then
  if X_NaCl==0 then
    x_NaCl := 0;
  else
//    x_NaCl := 1/(M_NaCl/M_H2O*(1/sum(X[1:5])-1)+1) "mol fraction";
    x_NaCl := 1/(M_NaCl/M_H2O*(1/X_NaCl-1)+1) "mol fraction";
  end if;
//  M_Solution := x_NaCl*M_NaCl + (1-x_NaCl)* M_H2O;

//CALCULATION OF EQUIVALENT TEMPERATURE_h
  q_21 := -1.69513-4.52781E-4*p_bar-6.04279E-8*p_bar^2;
  q_22 := 0.0612567+1.88082E-5*p_bar;

  q_20 := 1-q_21*sqrt(q_22) "x_NaCl = 0 results in q_2=1";

  q_23 := -q_20-q_21*sqrt(1+q_22) + 0.241022 + 3.45087E-5*p_bar - 4.28356E-9*p_bar^2
    "x_NaCl = 1 is pure NaCl";
  q_2 := q_20+q_21*sqrt(x_NaCl+q_22)+q_23*x_NaCl;

  q_10 := 47.9048 - 9.36994E-3*p_bar + 6.51059E-6*p_bar^2;
  q_11 := -32.1724+0.0621255*p_bar;

  q_12 := -q_10 - q_11;
  //x_NaCl=0 results in q_1=0 / x_NaCl = 1 is pure NaCl
  q_1 := q_10 + q_11*(1-x_NaCl) + q_12*(1-x_NaCl)^2;

  T_Scale_h := SI.Conversions.from_degC(T_C*q_2 + q_1);

//  p_check:=max(p, Modelica.Media.Water.WaterIF97_pT.saturationPressure(T_Scale_h)) "To make sure its liquid";

  //END OF CALCULATION OF EQUIVALENT TEMPERATURE

/*  state_H2O := Modelica.Media.Water.WaterIF97_pT.setState_pTX(p, SI.Conversions.from_degC(T_Scale_h), fill(0,0));
  h := Modelica.Media.Water.WaterIF97_pT.specificEnthalpy(state_H2O);*/
//  h := Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p, SI.Conversions.from_degC(T_Scale_h));

//  print("Brine_Driesner.specificEnthalpy_pTX: "+String(p*1e-5)+"bar."+String(T_Scale_h)+"degC->"+String(h)+" J/kg");
end T_Scale_h_Driesner;
