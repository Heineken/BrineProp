within BrineProp.SpecificEnthalpies;
function specificEnthalpy_pTX_liq_Francke_cp "enthalpy calculation DIY"
//based on Driesner enthalpy function NaCl and pressure independent 2D-fits (T,b) for cp measurement data for KCl and CaCl2 solutions
//TODO: Add Cp_appmol for (NaCl),MgCl and SrCl
  input SI.Pressure p;
  input SI.Temp_K T;
  input SI.MassFraction X[:] "mass fractions m_i/m_Sol";
  output SI.SpecificEnthalpy h;/**/
//  constant Real M_NaCl=Salt_Data.M_NaCl "molar mass in [kg/mol]";
//  output Real val2=H_appmol[CaCl2];
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

  Partial_Units.Molality b[size(X,1)]=massFractionsToMolalities(X,cat(1,MM_vec_salt,fill(-1,size(X,1)-size(MM_vec_salt,1))));
//  Partial_Units.Molality b[5]=mola[1:5];

//  SI.SpecificEnthalpy h_H2O =  Modelica.Media.Water.WaterIF97_base.specificEnthalpy_pT(p, T);
  SI.SpecificEnthalpy h_Driesner = specificEnthalpy_pTX_Driesner(p, T, X[1]/(X[1]+X[end]));

  SI.MolarMass MM_vec_salt[:]=BrineProp.SaltData.MM_salt[1:5];
  /*SI.Temp_C T_C = SI.Conversions.to_degC(T);
  Pressure_bar p_bar=SI.Conversions.to_bar(p); */
/*
//  Modelica.Media.Water.WaterIF97_base.ThermodynamicState state_H2O;
//  SI.MolarMass M_Solution "[kg/mol]";
//  SI.Pressure p_check;
    Real int_cp_by_cpWater[5]={0,
    int_HeatCapacity_KCl_White_unused(T,b[KCl]),
    int_HeatCapacity_CaCl2_White_unused(T,b[CaCl2]),
    0,0} "cp ratios integrated over T";
 SI.SpecificHeatCapacity cp_H2O=Modelica.Media.Water.IF97_Utilities.cp_pT(p,T,1);
 
 BrineProp.Partial_Units.PartialMolarEnthalpy[:] H_appmol={(if b[i] > 0 then 
     ((1 .+ MM_vec_salt[i] .* b[i]) .* int_cp_by_cpWater[i] .- T)*cp_H2O ./ b[i]
       else 0) for i in 1:5} .+ Delta_h_solution 
    "Apparent molar enthalpy of salts";*/

  BrineProp.Partial_Units.PartialMolarEnthalpy[:] H_appmol={0,
   if b[KCl] >0 then appMolarEnthalpy_KCl_White(T,b[KCl]) else 0,
   if b[CaCl2] >0 then appMolarEnthalpy_CaCl2_White(T,b[CaCl2]) else 0,
   0,0};
algorithm
//  print("specificEnthalpy_pTX_Driesner: "+String(h_Driesner)+" J/kg");
//print("h_H2O: "+String(h_H2O)+" J/kg");
//print("h_salt: "+String(h_salt[1]) +" J/kg");
//print("cp_salt: "+String(size(cp_salt,1))+"");
//print("Delta_h_solution_NaCl: "+String(Delta_h_solution_NaCl)+" J/kg");
//print("b[CaCl2]: "+String(b[CaCl2])+" mol/kg");

//  p_bar := SI.Conversions.to_bar(p);
//  assert(T_C>=0 and T_C<=1000, "T="+String(T-273.15)+", but must be between 0 and 1000°C");
//  assert(p_bar>=1 and p_bar<=1000, "P="+String(p/1e5)+" bar, but must be between 1 and 1000 bar");
//  assert(mola>=.25 and mola<=5, "Molality must be between 0.25 and 5 mol/kg");

//Salinity conversion
/*
  if X[1]==0 then
    x_NaCl := 0;
  else
    x_NaCl := 1/(M_NaCl/M_H2O*(1/sum(X[1:5])-1)+1) "mol fraction";
  end if;
  M_Solution := x_NaCl*M_NaCl + (1-x_NaCl)* M_H2O;
*/

  h := (X[NaCl]+X[end])*h_Driesner + X[end]*b[2:5]*H_appmol[2:5];
//  h := h_H2O;
//  print("MM_vec_salt      :"+PowerPlant.vector2string(MM_vec_salt));
//  print("b                :"+PowerPlant.vector2string(b));
//  print("int_cp_by_cpWater:"+Modelica.Math.Matrices.toString({int_cp_by_cpWater}));
//  print("H_appmol         :"+PowerPlant.vector2string(H_appmol));

//  print("Brine.specificEnthalpy_pTX_Francke: "+String(p*1e-5)+"bar. "+String(T)+"°C->"+String(h)+" J/kg");
end specificEnthalpy_pTX_liq_Francke_cp;
