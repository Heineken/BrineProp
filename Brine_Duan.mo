within Brine;
package Brine_Duan "NaCl solution using Duan density"
  extends PartialBrine(redeclare package Salt_data =
        Salt_Data_Duan);

  constant Salt_data.SaltConstants salt = Salt_data.saltConstants[NaCl];

  function density_Duan2008_pTX "density calculation of an aqueous salt solution according to Shide Mao and Zhenhao Duan (2008) 0-300°C; 0.1-100MPa; 0-6 mol/kg
  http://dx.doi.org/10.1016/j.jct.2008.03.005
  http://www.geochem-model.org/wp-content/uploads/2009/09/55-JCT_40_1046.pdf
  Problems: Brine has the same evaporation temperature as pure water, only different  "

      input Modelica.SIunits.Pressure p;
      input Modelica.SIunits.Temp_K T;
      input Modelica.SIunits.MassFraction X[:] "mass fractions m_NaCl/m_Sol";
      input Modelica.SIunits.MolarMass MM[:] "molar masses of components";

    //  input saltCoefficients salt;
      output Modelica.SIunits.Density d;
      //    constant Modelica.SIunits.MolarMass M_NaCl=0.058443 "molar mass in [kg/mol]";
    //  constant Modelica.SIunits.MolarMass M_H2O = 0.018015 "[kg/mol]";
  public
      final constant Real b = 1.2;
      final constant Real U[:] = {
         3.4279E2,
        -5.0866E-3,
         9.4690E-7,
        -2.0525,
         3.1159E3,
        -1.8289E2,
        -8.0325E3,
         4.2142E6,
       2.1417};  //dielectric constant D of pure water according to Bradley and Pitzer (1979)
      final constant Real N_0(final unit="1/mol") = Modelica.Constants.N_A
      "Avogadro constant in [1/mol]";
    //  e := 1.60217733E-19 [C] "elementary charge in Coulomb";
      final constant Real e = 1.60217733E-19 *10*299792458
      "elementary charge in [esu]";
    //k := 1.3806505E-23 "Boltzmann constant in [J/K]";
      final constant Real k = 1.3806505E-16 "Boltzmann constant in [erg/K]";
      final constant Real R = Modelica.Constants.R "Gas constant [J/mol*K]";
    /*  constant Integer nX_salt =  size(X,1) 
    "TODO: diese Zeile und alle Verweise au nX_salt entfernen";*/
  protected
      constant Integer nX_salt=5;
      Real mola "molality b (mol_NaCl/kg_sol)";
      Modelica.SIunits.MassFraction w_salt;
    //  Modelica.SIunits.Temp_C T_C = Modelica.SIunits.Conversions.to_degC(T);
      Pressure_bar p_bar=Modelica.SIunits.Conversions.to_bar(p);
      Pressure_MPa p_MPa=p*1e-6;
      Real v;
     // Modelica.Media.Water.WaterIF97_base.ThermodynamicState state_H2O;
      Modelica.SIunits.Mass m;
      Real I;
      Real I_mr;
      Modelica.SIunits.Density rho_sol_r;
      Modelica.SIunits.Density rho_H2O;
      Modelica.SIunits.Density rho_H2O_plus;
      Modelica.SIunits.Density rho_H2O_minus;
      Real p_plus_bar;
      Real p_minus_bar;
      Real D_plus;
      Real D_minus;
      Real A_Phi_plus;
      Real A_Phi_minus;
      Real B_v;
      Real C_v;
      Real V_m_r;
      Real Bb;
      Real Cc;
      Real D_1000;
      Real D;

      Real A_Phi;
      Real dp;
      Real dA_Phi;
      Real A_v;
      Real V_o_Phi;
      Real V_Phi;
      Real h;
      Real h_mr;

      Real M_salt;
      Real m_r;
      Real z_plus;
      Real z_minus;
      Real v_plus;
      Real v_minus;
      Real[23] c;

      Modelica.SIunits.Density[nX_salt] rho;
      SaltConstants salt;
    //  constant Molality[:] molalities=Partial_Gas_Data.massFractionsToMolalities(X,MM[1:nX_salt])
      constant Molality[:] molalities=PowerPlant.Media.Brine.Partial_Gas_Data.massFractionsToMolalities_(
                                                                                 X,MM)
      "TODO: nich sauber der Querverweis";
    //  constant Molality[nXi] molalities = massFractionsToMolality(X[1:nXi],salt.MM_salt[1:nXi]);

  algorithm
    //  p_bar := Modelica.SIunits.Conversions.to_bar(p_Pa);
      //Density of pure water
    /*  state_H2O := Modelica.Media.Water.WaterIF97_base.setState_pTX(p, T, fill(0,0));
  rho_H2O := Modelica.Media.Water.WaterIF97_base.density(state_H2O) * 1e-3 "kg/m³->kg/dm³";*/
    //  rho_H2O := Modelica.Media.Water.WaterIF97_base.density(Modelica.Media.Water.WaterIF97_base.setState_pTX(p, T, fill(0,0))) * 1e-3 "kg/m³->kg/dm³";

      assert(Modelica.Media.Water.WaterIF97_base.saturationPressure(T)<p,"T="+String(T-273.15)+"°C is above evaporation temperature at p="+String(p/1e5)+" bar!");
    /*  if (Modelica.Media.Water.WaterIF97_base.saturationPressure(T)>p) then
    d:=-1 "if above evaporation temperature";
    Modelica.Utilities.Streams.print("above evaporation temperature!");
    return;
  end if;*/

      rho_H2O := Modelica.Media.Water.WaterIF97_base.density_pT(p, T) * 1e-3
      "kg/m³->kg/dm³";
    //   Modelica.Utilities.Streams.print("rho_H2O=" +String(rho_H2O)+" kg/dm³");

      //for pure water skip the whole calculation and return water density
      if max(X[1:nX_salt]) <= 1e-12 then
        d:= rho_H2O*1000;
        return;
      end if;
      for i in 1:nX_salt loop
        if X[i]>0 then
          salt := saltConstants[i];

    //      Modelica.Utilities.Streams.print(String(X[i])+"");

          assert(T>=salt.T_min_rho and T<=salt.T_max_rho, "Temperature is "+String(Modelica.SIunits.Conversions.to_degC(T))+"°C, but for "+salt.name+" must be between "+String(Modelica.SIunits.Conversions.to_degC(salt.T_min_rho))+"°C and "+String(Modelica.SIunits.Conversions.to_degC(salt.T_max_rho))+"°C");
          assert(p>=salt.p_min_rho and p<=salt.p_max_rho, "Pressure is "+String(p_bar)+" bar, but for "+salt.name+" must be between "+String(salt.p_min_rho*1e-5)+" bar and "+String(salt.p_max_rho*1e-5)+" bar");
          assert(molalities[i]>=0 and molalities[i]<=salt.mola_max_rho, "Molality of "+salt.name+" is "+String(molalities[i])+", but must be between 0 and "+String(salt.mola_max_rho)+" mol/kg");

          M_salt := salt.M_salt * 1000 "in g/mol";
          m_r :=salt.m_r;
          z_plus :=salt.z_plus;
          z_minus :=salt.z_minus;
          v_plus :=salt.v_plus;
          v_minus :=salt.v_minus;
          c :=salt.C;

      //    Modelica.Utilities.Streams.print(salt.name+": "+String(X[i]));

          v := v_plus + v_minus;

          //Conversion to mass and mol fraction
        //        w_salt := (m*M_salt*Convert('g';'kg'))/(1+m*M_salt*Convert('g';'kg'));
          w_salt :=X[i];
          m := w_salt/(M_salt*1e-3*(1-w_salt));
        //  x_salt := (m*M_H2O *Convert('g';'kg'))/(1+m*M_H2O *Convert('g';'kg'));

        //---------------------------------------------------

          //Equation 3: Ionic strength
          I    := 1/2 * (m  *v_plus*z_plus^2 + m  *v_minus*z_minus^2);
          I_mr := 1/2 * (m_r*v_plus*z_plus^2 + m_r*v_minus*z_minus^2);

          //Equation 4:
          h    := Modelica.Math.log10(1+b*I   ^(0.5))/(2*b);
          h_mr := Modelica.Math.log10(1+b*I_mr^(0.5))/(2*b);

        //---------------------------------------------------
        // equations using empirically fitted coefficients
        //---------------------------------------------------

          //Equation 10: solution volume at reference molality
          V_m_r :=  c[01] + c[02]*T + c[03]*T^2 + c[04]*T^3
                  + p_MPa*(c[05] + c[06]*T + c[07]*T^2 + c[08]*T^3);

          //Check: solution density at reference molality
          rho_sol_r := (1000 + m_r*M_salt)/V_m_r;

          //Equation 11: second virial coefficient. depends on temperature and pressure
          B_v :=     c[09]/(T-227) + c[10] + c[11]*T + c[12]*T^2 + c[13]/(647-T)
                + p_MPa*( c[14]/(T-227) + c[15] + c[16]*T + c[17]*T^2 + c[18]/(647-T));

          //Equation 12: third virial coefficient. depends on temperature
          C_v :=     c[19]/(T-227) + c[20] + c[21]*T + c[22]*T^2 + c[23]/(647-T);

        //---------------------------------------------------
        // Appendix A: Debye-Hückel limiting law slopes"
        //---------------------------------------------------

          Bb := U[07] + U[08]/T + U[09]*T;
          Cc := U[04] + U[05]/(U[06]+T);
          D_1000 := U[01]*exp(U[02]*T + U[03]*T^2);
          D := D_1000 + Cc*log((Bb + p_bar)/(Bb+1000));

          //DH-slope for osmotic coefficient according to Bradley and Pitzer (1979)
          A_Phi := 1/3 * ((2*Modelica.Constants.pi*N_0*rho_H2O)/1000)^(1/2) * (e^2/(D*k*T))^(3/2);

          //numeric differentiation per dp
          dp := 1E-3 * p_bar;
          p_plus_bar  := p_bar + dp/2;
          p_minus_bar := p_bar - dp/2;
          D_plus  := D_1000 + Cc*log((Bb + p_plus_bar) /(Bb + 1000));
          D_minus := D_1000 + Cc*log((Bb + p_minus_bar)/(Bb + 1000));

          /*state_H2O := Modelica.Media.Water.WaterIF97_base.setState_pTX(Modelica.SIunits.Conversions.from_bar(p_plus), T, fill(0,0));
      rho_H2O_plus := Modelica.Media.Water.WaterIF97_base.density(state_H2O) * 1e-3;*/
          rho_H2O_plus := Modelica.Media.Water.WaterIF97_base.density_pT(p_plus_bar*1e5, T) * 1e-3
          "kg/m³->kg/dm³";

          /*state_H2O := Modelica.Media.Water.WaterIF97_base.setState_pTX(Modelica.SIunits.Conversions.from_bar(p_minus), T, fill(0,0));
      rho_H2O_minus := Modelica.Media.Water.WaterIF97_base.density(state_H2O) * 1e-3;*/
          rho_H2O_minus := Modelica.Media.Water.WaterIF97_base.density_pT(p_minus_bar*1e5, T) * 1e-3
          "kg/m³->kg/dm³";
          A_Phi_plus  := 1/3 * (2*Modelica.Constants.pi*N_0*rho_H2O_plus /1000)^(1/2)* (e^2/(D_plus *k*T))^(3/2);
          A_Phi_minus := 1/3 * (2*Modelica.Constants.pi*N_0*rho_H2O_minus/1000)^(1/2)* (e^2/(D_minus*k*T))^(3/2);
          dA_Phi := (A_Phi_plus - A_Phi_minus);

          //DH-slope for apparent molar volume according to Rogers and Pitzer (1982)
          A_v := 23 * (-4*R*T * dA_Phi/dp) "where does the 23 come from??";

        //---------------------------------------------------
        // Solution 1: using V_o_Phi and V_Phi
        //---------------------------------------------------

                        //Equation 13: apparent molar Volume at infinite dilution in cm³/mol
                        V_o_Phi := (V_m_r/m_r - 1000/(m_r*rho_H2O) - v*abs(z_plus*z_minus)*A_v*h_mr
                        - 2*v_plus*v_minus*R*T * (B_v*m_r + v_plus*z_plus*C_v*m_r^2));

                        //Equation 2: apparent molar Volume in cm^3/mol
                        V_Phi := V_o_Phi + v*abs(z_plus*z_minus)*A_v*h + 2*v_plus*v_minus*m*R*T * (B_v + v_plus*z_plus*m*C_v);

                        //Equation 1: density of the solution
                        rho[i] := ((1000+m*M_salt)*rho_H2O) / (1000+m*V_Phi*rho_H2O)*1000;

        //---------------------------------------------------
        // Solution 2: using V_m
        //---------------------------------------------------

                        /*
                    //Equation 8: solution volume
                    //V_m = m*( V_m_r/m_r + 1000/rho_H2O *(1/m - 1/m_r) +              v*abs(z_plus*z_minus) * A_v*(h-h_mr) + 2*v_plus*v_minus*R*T_K* (B_v*(m-m_r) + v_plus*z_plus*C_v*(m^2-m_r^2)) )
                    V_m =     V_m_r*m/m_r + 1000/rho_H2O - (1000*m)/(rho_H2O*m_r) + m*(v*abs(z_plus*z_minus) * A_v*(h-h_mr) + 2*v_plus*v_minus*R*T_K* (B_v*(m-m_r) + v_plus*z_plus*C_v*(m^2-m_r^2)) )
    
                    //density of the solution
                    Density_Duan = (1000 + m*M_salt)/V_m
                    */
        end if;
      end for;
    //  d := X[1:end-1]*rho/(1-X[end]) "linear mixture (matrix multiplication)";
    //  Modelica.Utilities.Streams.print("rho: "+String(rho[2]));
      d := molalities[1:nX_salt]*rho/(sum(molalities[1:nX_salt]))
      "linear mixture (matrix multiplication)";
    //  Modelica.Utilities.Streams.print("Density: "+String(d));
  end density_Duan2008_pTX;

  redeclare function extends dynamicViscosity_pTX
   //  constant Real M_NaCl=0.058443 "molar mass in [kg/mol]";
    /*  public 
     constant Modelica.SIunits.MolarMass M_NaCl = salt.M_salt; 
      "[kg/mol]";*/
  protected
    Molality mola = X[1]/(salt.M_salt*(1-X[1])) "molality b (mol_NaCl/kg_H2O)";
    Modelica.SIunits.Temp_C T_C = Modelica.SIunits.Conversions.to_degC(T_K);
    Pressure_bar p_bar= Modelica.SIunits.Conversions.to_bar(p_Pa);

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

    Modelica.SIunits.DynamicViscosity eta_H2O;
    Modelica.Media.Water.WaterIF97_base.ThermodynamicState state_H2O;
  algorithm
    assert(T_C>=0 and T_C<=300, "Temperature is "+String(Modelica.SIunits.Conversions.to_degC(T_K))+"°C, but must be between 10 and 350°C");
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
  //  Modelica.Utilities.Streams.print("mola: "+String(mola));
  //  Modelica.Utilities.Streams.print("eta_relative: "+String(eta_relative));
  //  Modelica.Utilities.Streams.print("eta: "+String(eta));

  end dynamicViscosity_pTX;

  redeclare function extends specificEnthalpy_pTX
    "enthalpy calculation according to Driesner et al: 0-1000°C; 0.1-500MPa"
  algorithm
    h :=Brine_Driesner.specificEnthalpy_pTX(p,T,X);
  //  h :=T;

  end specificEnthalpy_pTX;
end Brine_Duan;
