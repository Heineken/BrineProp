within BrineProp.SpecificEnthalpies;
function specificEnthalpy_pTX_liq_Francke_cp "enthalpy calculation DIY"
//based on Driesner enthalpy function NaCl and pressure independent 2D-fits (T,b) for cp measurement data for KCl and CaCl2 solutions
//TODO: Add Cp_appmol for (NaCl),MgCl and SrCl
  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X[:] "mass fractions m_i/m_Sol";
  input SI.MolarMass MM_vec[:];
  input Boolean ignoreTlimit=false "activated by temperature_phX";
  output SI.SpecificEnthalpy h;
//  constant Real M_NaCl=Salt_Data.M_NaCl "molar mass in [kg/mol]";
//  output Real val2=H_appmol[CaCl2];
//  extends BrineProp.SaltDataDuan.defineSaltOrder_;

//  SI.MolarMass MM_vec[:]=BrineProp.SaltData.MM_salt[1:5];
//  SI.MolarMass MM_vec[:];

protected
  SI.MolarInternalEnergy Delta_h_solution_NaCl = 0
    "HeatOfSolution_NaCl_Sanahuja1986(T) [J/mol_NaCl]";
 //  constant SI.MolarInternalEnergy Delta_h_solution_NaCl = -3460 "[J/mol_NaCl] @313.15K Sanahuja1984 http://dx.doi.org/10.1016/0021-9614(84)90192-7";
// SI.MolarInternalEnergy Delta_h_solution_KCl = HeatOfSolution_KCl_Sanahuja1986(T) "[J/mol_NaCl]";
// constant SI.MolarInternalEnergy Delta_h_solution_KCl = -17000 "[J/mol_KCl] http://webserver.dmt.upm.es/~isidoro/dat1/Heat%20of%20solution%20data.htm";
//  constant SI.MolarInternalEnergy Delta_h_solution_CaCl2 = 82900 "[J/mol_CaCl2] http://webserver.dmt.upm.es/~isidoro/dat1/Heat%20of%20solution%20data.htm";
  constant SI.MolarInternalEnergy Delta_h_solution_CaCl2 = 81850
    "[J/mol_CaCl2] Sinke1985 http://dx.doi.org/10.1016/0021-9614(85)90083-7";
  constant SI.MolarInternalEnergy Delta_h_solution_MgCl2 = 0 "[J/mol_MgCl2]";
  constant SI.MolarInternalEnergy Delta_h_solution_SrCl2 = 0 "[J/mol_SrCl]";
/*  constant SI.MolarInternalEnergy[:] Delta_h_solution = {
    Delta_h_solution_NaCl,
    Delta_h_solution_KCl,
    Delta_h_solution_CaCl2,
    Delta_h_solution_MgCl2,
    Delta_h_solution_SrCl2};*/

  Types.Molality b[size(X, 1)]=
      Utilities.massFractionsToMolalities(X,cat(1,MM_vec,fill(-1, size(X, 1) - size(MM_vec, 1))));

  SI.Temperature T2 "fixed temperature (workaround for input NaN input";
//  SI.SpecificEnthalpy h_H2O =  Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p, T);
  SI.SpecificEnthalpy h_Driesner;

  BrineProp.Types.PartialMolarEnthalpy[5] H_appmol
    "TODO: remove absolute indices";
  parameter Integer[:] otherSalts=cat(1,
  if iKCl>0 then {iKCl} else fill(0,0),
  if iCaCl2>0 then {iCaCl2} else fill(0,0),
  if iMgCl2>0 then {iMgCl2} else fill(0,0),
  if iSrCl2>0 then {iSrCl2} else fill(0,0));
algorithm
    if debugmode then
      print("\nRunning specificEnthalpy_pTX_liq_Francke_cp("+String(p/1e5)+" bar,"+String(T-273.15)+"degC, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
    end if;

  assert(size(ignoreLimitSalt_p,1)==size(X,1)-1,"Length of ignoreLimitSalt_p ("+String(size(ignoreLimitSalt_p,1))+" and X ("+String(size(X,1)-1)+") don't match."); //needed here, because flag vector with fewer than nX_salts elements causes "out of bounds" and is not caught elsewere
  assert(size(ignoreLimitSalt_T,1)==size(X,1)-1,"Length of ignoreLimitSalt_T ("+String(size(ignoreLimitSalt_T,1))+" and X ("+String(size(X,1)-1)+") don't match."); //should be in PartialFlags, but asserts can't be in packages
  assert(size(ignoreLimitSalt_b,1)==size(X,1)-1,"Length of ignoreLimitSalt_b ("+String(size(ignoreLimitSalt_b,1))+" and X ("+String(size(X,1)-1)+") don't match.");

  if String(T)== "-1.#IND" then
    assert(false, "T > 1e99, probably division by zero. Setting T=300 K.",AssertionLevel.warning);
    T2:=350;
  else
    T2 := T;
  end if;

  h_Driesner :=specificEnthalpy_pTX_Driesner(p,T2,X[1]/(X[1] + X[end]));

  H_appmol:={0,if b[iKCl] > 0 then appMolarEnthalpy_KCl_White(
    T2,
    b[iKCl],
    ignoreTlimit) else 0,if b[iCaCl2] > 0 then appMolarEnthalpy_CaCl2_White(
    T2,
    b[iCaCl2],
    ignoreTlimit) else 0,0,0} "TODO: remove absolute indices";

  h := (X[iNaCl]+X[end])*h_Driesner + X[end]*b[otherSalts]*H_appmol[otherSalts]
    "TODO: remove absolute indices";

//  print("MM_vec      :"+PowerPlant.vector2string(MM_vec));
//  print("b                :"+PowerPlant.vector2string(b));
//  print("otherSalts:"+Modelica.Math.Matrices.toString({otherSalts}));
//  print("H_appmol         :"+PowerPlant.vector2string(H_appmol));

//  print("Brine.specificEnthalpy_pTX_Francke: "+String(p*1e-5)+"bar. "+String(T)+"degC->"+String(h)+" J/kg");
end specificEnthalpy_pTX_liq_Francke_cp;
