within BrineProp.Viscosities;
function dynamicViscosity_Duan_pTX
  "Multisalt-Version of viscosity calculation according to Duan et al 2009: Considers NaCl and KCL, with linear mixing"
  //doi:10.1007/s10765-009-0646-7

  input Modelica.SIunits.Pressure p_Pa;
  input Modelica.SIunits.Temp_K T_K;
  input Modelica.SIunits.MassFraction X[:] "mass fraction m_NaCl/m_Sol";
  input Modelica.SIunits.MolarMass MM[:];
  input BrineProp.SaltData_Duan.SaltConstants[:]
                                        Salt_Constants;
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
  Modelica.SIunits.DynamicViscosity[nX_salt] etas=fill(0,nX_salt);
  Modelica.SIunits.DynamicViscosity eta_H2O;
  Modelica.Media.Water.WaterIF97_base.ThermodynamicState state_H2O;

  SaltData_Duan.SaltConstants salt;
  constant Molality[:] molalities=massFractionsToMolalities(X,MM);
 // constant Partial_Units.Molality[:] molalities=X[1:nX_salt] ./ MM[1:nX_salt]/      X[end];
algorithm
  if debugmode then
    Modelica.Utilities.Streams.print("p="+String(p_Pa)+" Pa, T_K"+String(T_K)+" K (Brine.Viscosities.dynamicViscosity_Duan_pTX)");
  end if;
  assert(T_C>=0 and T_C<=400, "Temperature must be between 10 and 350°C");
  assert(p_bar>=1 and p_bar<=1000, "Pressure must be between 1 and 500 bar");

  //viscosity calculation
  state_H2O := Modelica.Media.Water.WaterIF97_base.setState_pTX(p_Pa, T_K, X);
  eta_H2O := Modelica.Media.Water.WaterIF97_base.dynamicViscosity(state_H2O);
//  Modelica.Utilities.Streams.print("eta_H2O= "+String(eta_H2O)+" Pa·s");

  //for pure water skip the whole calculation and return water viscosity
  if max(X[1:nX_salt])<1e-8 then
    eta:= eta_H2O;
    return;
  end if;

  eta:= eta_H2O^X[end];
//  Modelica.Utilities.Streams.print("eta_H2O^X[end]= "+String(eta_H2O)+"^"+String(X[end]) + " -> "+String(eta)+" Pa·s");
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
      A := salt.a[1] + salt.a[2]*T_K + salt.a[3]*T_K^2;
      B := salt.b[1] + salt.b[2]*T_K + salt.b[3]*T_K^2;
      C := salt.c[1] + salt.c[2]*T_K;
      eta_relative := exp(A*molalities[i] + B*molalities[i]^2 + C*molalities[i]^3);
      etas[i] := eta_relative * eta_H2O;
//    Modelica.Utilities.Streams.print("molaMulti["+String(i)+"]: "+String(molalities[i]));
//    Modelica.Utilities.Streams.print("eta_H2O: "+String(eta_H2O));
//    Modelica.Utilities.Streams.print("eta_relative: "+String(eta_relative));
//    Modelica.Utilities.Streams.print("etas["+String(i)+"]: "+String(etas[i]));
      eta:= eta*etas[i]^X[i] "geometric mean mixture rule (as in Laliberté)";
//      Modelica.Utilities.Streams.print("Viscosity: "+String(etas[i])+"->"+String(eta)+" Pa·s (BrineProp.Viscosities.dynamicViscosity_Duan_pTX)");
    end if;
//    eta := eta + etas[i]*molalities[i];
  end for;
//  eta := etas*molalities[1:nX_salt]/(sum(molalities[1:nX_salt])) "linear mixing rule molality weighted (as in Duan2009)";
//  eta := prod(cat(1,{eta_H2O},etas).^X[1:nX_salt]) "geometric mean mixture rule (as in Laliberté)";

//  Modelica.Utilities.Streams.print("Viscosity("+String(p_Pa)+","+String(T_K)+"): "+String(eta)+" Pa·s (Partial_Viscosity_Phillips.dynamicViscosity_Duan_pTX)");
end dynamicViscosity_Duan_pTX;
