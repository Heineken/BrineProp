within BrineProp.SpecificEnthalpies;
function specificHeatCapacityCp_pTX_liq_Francke
  "calculation of liquid specific heat capacity for NaCl-KCl-CaCl2-liquid from apparent molar heat capacities"
  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X[:] "mass fractions m_i/m_Sol";
  input SI.MolarMass MM_vec[:];
  output SI.SpecificHeatCapacity cp;
//  extends BrineProp.SaltDataDuan.defineSaltOrder_;

protected
  Integer nX_salt = size(X, 1)-1;
  Types.Molality b[nX_salt+1]
    "=Utilities.massFractionsToMolalities(X_,cat(1,MM_vec,fill(-1, size(X_, 1) - size(MM_vec, 1))))";

/*  Real cp_by_cpWater[:]={0,
      SpecificEnthalpies.HeatCapacityRatio_KCl_White(T, b[KCl]),
      SpecificEnthalpies.HeatCapacityRatio_CaCl2_White(T, b[CaCl2]),
      0,0} "cp/cp_H2O of salt solutions";*/
  Types.PartialMolarHeatCapacity[nX_salt] Cp_appmol
    "Apparent molar enthalpy of salts";

  SI.SpecificHeatCapacity cp_Driesner;

  //  SI.SpecificHeatCapacity cp_H2O=Modelica.Media.Water.IF97_Utilities.cp_pT(p,T);

//  SI.MassFraction X_[:]=state.X "mass fraction m_NaCl/m_Sol";
//    SI.MassFraction X_[:]=cat(1,state.X[1:end-1],{1-sum(state.X[1:end-1])}) "Doesn't work in function in OM";
//    SI.MassFraction X_[size(state.X,1)] "OM workaround for cat";
    SI.MassFraction X_[nX_salt+1]
    "OM workaround for cat TODO: still necessary?";
algorithm
    if debugmode then
      print("Running specificHeatCapacityCp_pTX_liq_Francke("+String(p/1e5)+" bar,"+String(T-273.15)+"degC, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
    end if;
    X_[1:end-1]:=X[1:end-1] "OM workaround for cat";
    X_[end]:=1-sum(X[1:end-1]) "OM workaround for cat";
    b:=Utilities.massFractionsToMolalities(X_, cat(
    1,MM_vec,fill(-1, size(X_, 1) - size(MM_vec, 1))));
//    assert(X[end]>0, "No water in brine.");
    cp_Driesner:=specificHeatCapacity_pTX_Driesner(p,T,X_[1]/(X_[1] + X_[end]));

  Cp_appmol:=cat(1,{0,
      if b[iKCl] > 0 then appMolarHeatCapacity_KCl_White(T, b[iKCl]) else 0,
      if b[iCaCl2] > 0 then appMolarHeatCapacity_CaCl2_White(T,b[iCaCl2]) else 0},
      fill(0,nX_salt-3)) "Apparent molar enthalpy of salts";
//    Cp_appmol:={(if b[i] > 0 and cp_by_cpWater[i]>0 then ((1 .+ MM_vec[i] .* b[i]) .* cp_by_cpWater[i] .- 1)*cp_H2O ./ b[i] else 0) for i in 1:5};

    cp := (X_[iNaCl]+X_[end])*cp_Driesner + X_[end]*b[2:nX_salt]*(Cp_appmol[2:nX_salt]);
    annotation (Documentation(info="<html>
                                <p>In the two phase region this function returns the interpolated heat capacity between the
                                liquid and vapour state heat capacities.</p>
                                </html>"));
end specificHeatCapacityCp_pTX_liq_Francke;
