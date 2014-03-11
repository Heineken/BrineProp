within BrineProp.Viscosities;
function dynamicViscosity_DuanZhang_pTXd
  "Multisalt-Version of viscosity calculation according to Duan et al 2009 and Zhang et al 1997: Considers NaCl and KCL, with geometric mixture rule"
  //doi:10.1007/s10765-009-0646-7

  input Modelica.SIunits.Pressure p_Pa;
  input Modelica.SIunits.Temp_K T_K;
  input Modelica.SIunits.MassFraction X[:] "mass fraction m_NaCl/m_Sol";
  input Modelica.SIunits.Density d;
  input Modelica.SIunits.MolarMass MM[:];
  input BrineProp.SaltData_Duan.SaltConstants[:] Salt_Constants;
  output Modelica.SIunits.DynamicViscosity eta;
protected
  Modelica.SIunits.Temp_C T_C = Modelica.SIunits.Conversions.to_degC(T_K);
  Pressure_bar p_bar=Modelica.SIunits.Conversions.to_bar(p_Pa);

  Real A;
  Real B;
  Real C;

  Real eta_relative;

//  Integer nXi = size(X,1);
  Integer nX_salt = size(Salt_Constants,1);
//  Modelica.SIunits.DynamicViscosity[nX_salt] etas=fill(0,nX_salt);
  Modelica.SIunits.DynamicViscosity eta_H2O;
  Modelica.Media.Water.WaterIF97_base.ThermodynamicState state_H2O;

  SaltData_Duan.SaltConstants salt;
  constant Molality[:] molalities=massFractionsToMolalities(X,MM);
 // constant Partial_Units.Molality[:] molalities=X[1:nX_salt] ./ MM[1:nX_salt]/      X[end];
   Molarity_molperliter c;
  Molality b "component molality";
  Real phi "mixing weight";
protected
  constant Pressure_bar p_min=1;
  constant Pressure_bar p_max=1000;
  constant Modelica.SIunits.Temp_C T_min=0;
  constant Modelica.SIunits.Temp_C T_max=400;
algorithm
//   Modelica.Utilities.Streams.print("X[3]="+String(X[3])+" (Brine.Viscosities.dynamicViscosity_DuanZhuang_pTXd)");
  if debugmode then
    Modelica.Utilities.Streams.print("p="+String(p_Pa)+" Pa, T_K"+String(T_K)+" K (Brine.Viscosities.dynamicViscosity_DuanZhuang_pTXd)");
  end if;

   if outOfRangeMode==1 then
      if not (p_bar>=p_min and p_bar<=p_max) then
        Modelica.Utilities.Streams.print("Pressure is " + String(p_bar) + " bar, but must be between " + String(p_min) + " bar and " + String(p_max) + " bar (BrineProp.Viscosities.dynamicViscosity_DuanZhang_pTXd)");
      end if;
      if not (T_C>=T_min and T_C<=T_max) then
        Modelica.Utilities.Streams.print("Temperature is "+String(T_C) + "�C, but must be between " + String(T_min) + "�C and " + String(T_max) + "�C (BrineProp.Viscosities.dynamicViscosity_DuanZhang_pTXd)");
      end if;
   elseif outOfRangeMode==2 then
    assert(p_bar>=p_min and p_bar<=p_max, "p="+String(p_bar)+", but must be between "+String(p_min)+" and "+String(p_max)+" bar");
    assert(T_C>=T_min and T_C<=T_max, "T="+String(T_C)+", but must be between "+String(T_min)+" and "+String(T_max)+"�C");
   end if;

  //viscosity calculation
  state_H2O := Modelica.Media.Water.WaterIF97_base.setState_pTX(p_Pa, T_K, X);
  eta_H2O := Modelica.Media.Water.WaterIF97_base.dynamicViscosity(state_H2O);
//  Modelica.Utilities.Streams.print("eta_H2O= "+String(eta_H2O)+" Pa�s");

  //for pure water skip the whole calculation and return water viscosity

  eta:= eta_H2O "^X[end]";
//  Modelica.Utilities.Streams.print("eta_H2O^X[end]= "+String(eta_H2O)+"^"+String(X[end]) + " -> "+String(eta)+" Pa�s");
  if max(X[1:nX_salt]) <= 1e-8 then
    //pure water -> skip the rest
    return;
  end if;

 for i in 1:nX_salt loop
    if X[i]>0 then
      salt := Salt_Constants[i];
      if outOfRangeMode==1 then
        if molalities[i]>salt.mola_max_eta then
          Modelica.Utilities.Streams.print(salt.name+" content exceeds limit in Viscosities.dynamicViscosity_Duan_pTX");
//          molalities[i]=min(molalities[i],salt.mola_max_eta);
        end if;
      elseif outOfRangeMode==2 then
          assert(ignoreLimitSalt_visc[i] or (molalities[i]>=0 and molalities[i]<=salt.mola_max_eta), "Molality of "+salt.name+" is "+String(molalities[i])+"(X="+String(X[i])+"), but must be between 0 and "+String(salt.mola_max_eta)+" mol/kg");
      end if;

//      Modelica.Utilities.Streams.print(salt.name+" content = "+String(molalities[i])+" (Viscosities.dynamicViscosity_Duan_pTX)");
//      if salt.name <> "NaCl" then
      //factors
      //MIXING WEIGHT
    //  phi:=X[i]/sum(X[1:nX_salt]) "geometric mean mixture rule weighted with mass fraction (as in Lalibert�)";
      phi:=molalities[i]/sum(molalities[1:nX_salt])
        "geometric mean mixture rule weighted with mass fraction (as in Lalibert�)";

        if i==3 then
        //Zhang (available for NaCl, KCl and CaCl)
       c :=X[i]/MM[i]*d/1000/phi "component molarity";
         eta_relative := 1 + salt.Zh_A*c^0.5 + salt.Zh_B*c + salt.Zh_D*c^2 + 1e-4*salt.Zh_E*c^3.5 + 1e-5*salt.Zh_F*c^7;
      else
        //Duan (available for NaCl and KCl)
        b:=molalities[i]/phi;
        A := salt.a[1] + salt.a[2]*T_K + salt.a[3]*T_K^2;
        B := salt.b[1] + salt.b[2]*T_K + salt.b[3]*T_K^2;
        C := salt.c[1] + salt.c[2]*T_K;
//        Modelica.Utilities.Streams.print(salt.name+" Duan. A="+String(A)+" B="+String(B)+" C="+String(C)+" b="+String(molalities[i]));
//        eta_relative := exp(A*molalities[i] + B*molalities[i]^2 + C*molalities[i]^3);
        eta_relative := exp(A*b + B*b^2 + C*b^3)
          "Mixture is composed of binary solutions of the same molality";
      end if;
//      etas[i] := eta_relative "* eta_H2O";
//    Modelica.Utilities.Streams.print("molaMulti["+String(i)+"]: "+String(molalities[i]));
//    Modelica.Utilities.Streams.print("eta_H2O: "+String(eta_H2O));
//    Modelica.Utilities.Streams.print("eta_relative: "+String(eta_relative));
//    Modelica.Utilities.Streams.print("etas["+String(i)+"]: "+String(etas[i]));

//      eta:= eta*etas[i] ^X[i] "geometric mean mixture rule (as in Lalibert�)";
      eta:= eta*eta_relative^phi;

//      Modelica.Utilities.Streams.print("Viscosity "+salt.name+" phi="+String(phi)+": "+String(eta_relative)+"->"+String(eta)+" Pa�s (BrineProp.Viscosities.dynamicViscosity_Duan_pTX)");

    end if;
//    eta := eta + etas[i]*molalities[i];
  end for;
//  eta := eta * (1+ sum(etas)) "additive mixing rule (Zhang 1997)";
//  eta := etas*molalities[1:nX_salt]/(sum(molalities[1:nX_salt])) "linear mixing rule molality weighted (as in Duan2009)";
//  eta := prod(cat(1,{eta_H2O},etas).^X[1:nX_salt]) "geometric mean mixture rule (as in Lalibert�)";

//  Modelica.Utilities.Streams.print("Viscosity("+String(p_Pa)+","+String(T_K)+"): "+String(eta)+" Pa�s (Partial_Viscosity_Phillips.dynamicViscosity_Duan_pTX)");
end dynamicViscosity_DuanZhang_pTXd;
