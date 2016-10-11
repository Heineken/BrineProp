within BrineProp.Viscosities;
function dynamicViscosity_Duan2009NaClKCl_pTX
  "Multisalt-Version of viscosity calculation according to Duan et al 2009: Considers NaCl and KCL, with geometric mixture rule"
  //doi:10.1007/s10765-009-0646-7

  input SI.Pressure p_Pa;
  input SI.Temp_K T_K;
  input SI.MassFraction X[:] "mass fraction m_NaCl/m_Sol";
  input SI.MolarMass MM[:];
  input BrineProp.SaltDataDuan.SaltConstants[:] Salt_Constants;
  output SI.DynamicViscosity eta;
protected
  SI.Temp_C T_C = SI.Conversions.to_degC(T_K);
  Pressure_bar p_bar=SI.Conversions.to_bar(p_Pa);

  Real A;
  Real B;
  Real C;

  Real phi;
  Real b;
  Real eta_relative;

  Integer nX_salt = size(Salt_Constants,1);
  SI.DynamicViscosity eta_H2O;
  Modelica.Media.Water.WaterIF97_pT.ThermodynamicState state_H2O;

  SaltDataDuan.SaltConstants salt;
  constant Molality[:] molalities=Utilities.massFractionsToMolalities(X,MM);
 String msg;
algorithm
  if debugmode then
    print("p="+String(p_Pa)+" Pa, T_K"+String(T_K)+"X[1]="+String(X[1])+" K (Brine.Viscosities.dynamicViscosity_Duan2009NaClKCl_pTX)");
  end if;
  assert(T_C>=0 and T_C<=350, "Temperature must be between 10 and 350degC");
  assert(p_bar>=1 and p_bar<=1000, "Pressure must be between 1 and 500 bar");

  //viscosity calculation
  state_H2O := Modelica.Media.Water.WaterIF97_pT.setState_pTX(p_Pa, T_K, X);
  eta_H2O := Modelica.Media.Water.WaterIF97_pT.dynamicViscosity(state_H2O);
//  print("eta_H2O= "+String(eta_H2O)+" Pa.s");

  //for pure water skip the whole calculation and return water viscosity
  if max(X[1:nX_salt])<1e-8 then
    eta:= eta_H2O;
    return;
  end if;

  eta:= eta_H2O "^X[end]";
//  print("eta_H2O^X[end]= "+String(eta_H2O)+"^"+String(X[end]) + " -> "+String(eta)+" Pa.s");

 for i in 1:nX_salt loop
    if X[i]>0 then
       salt := Salt_Constants[i];
    if AssertLevel>0 then
     assert(ignoreLimitSalt_b[i] or (molalities[i]>=0 and molalities[i]<=salt.mola_max_eta), "Molality of " + salt.name + " is " + String(molalities[i]) + "(X="
           + String(X[i]) + "), but must be between 0 and " + String(salt.mola_max_eta)
           + " mol/kg.\nTo ignore set ignoreLimitSalt_b["+String(i)+"]=true",aLevel);
     assert(max(cat(1,salt.a,salt.b,salt.c))>0,"No coefficients for " + salt.name+".",aLevel);
    end if;

//      print(salt.name+" content = "+String(molalities[i])+" (Viscosities.dynamicViscosity_Duan_pTX)");
      //factors
      phi:=molalities[i]/sum(molalities[1:nX_salt])
        "geometric mean mixture rule weighted with mass fraction (as in Laliberte)";
      b:=molalities[i]/phi;
      A := salt.a[1] + salt.a[2]*T_K + salt.a[3]*T_K^2;
      B := salt.b[1] + salt.b[2]*T_K + salt.b[3]*T_K^2;
      C := salt.c[1] + salt.c[2]*T_K;
      eta_relative := exp(A*b + B*b^2 + C*b^3)
        "Mixture is composed of binary solutions of the same molality";
//      print("Viscosity: "+String(etas[i])+"->"+String(eta)+" Pa.s (BrineProp.Viscosities.dynamicViscosity_Duan_pTX)");
    end if;
      eta:= eta*eta_relative^phi;
  end for;
//  print("Viscosity("+String(p_Pa)+","+String(T_K)+"): "+String(eta)+" Pa.s (Partial_Viscosity_Phillips.dynamicViscosity_Duan_pTX)");
end dynamicViscosity_Duan2009NaClKCl_pTX;
