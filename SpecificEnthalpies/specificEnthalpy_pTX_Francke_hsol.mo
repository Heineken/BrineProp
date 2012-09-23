within Brine.SpecificEnthalpies;
function specificEnthalpy_pTX_Francke_hsol
  "enthalpy calculation using solution enthalpies"
//Pressure limited to 100 MPa by Modelica Water property function
  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input Modelica.SIunits.MassFraction X[:] "mass fractions m_i/m_Sol";
  output Modelica.SIunits.SpecificEnthalpy h;/**/
//  constant Real M_NaCl=Salt_Data.M_NaCl "molar mass in [kg/mol]";
protected
  Modelica.SIunits.MolarMass MM_vec_salt[:]=Salt_Data.MM_salt[1:5];

//   constant Modelica.SIunits.SpecificHeatCapacity cp_NaCl = 50.5/M_NaCl "[J/(kg·K)] http://hyperphysics.phy-astr.gsu.edu/hbase/tables/therprop.html";
//   constant Modelica.SIunits.SpecificHeatCapacity cp_NaCl = 36.79/MM_vec_salt[NaCl] "[J/(kg·K)] http://en.wikipedia.org/wiki/Sodium_chloride";
   constant Modelica.SIunits.SpecificHeatCapacity cp_NaCl = 892.2568
    "[J/(kg·K)] aus Driesner 6molar 0...100°C";
   constant Modelica.SIunits.SpecificHeatCapacity cp_KCl = 690
    "[J/(kg·K)] http://www.korth.de/index.php/material-detailansicht/items/16.html";
   constant Modelica.SIunits.SpecificHeatCapacity cp_CaCl2 = 72.59/MM_vec_salt[CaCl2]
    "[J/(kg·K)] http://hyperphysics.phy-astr.gsu.edu/hbase/tables/therprop.html";
   constant Modelica.SIunits.SpecificHeatCapacity cp_MgCl2 = 0 "[J/(kg·K)]";
   constant Modelica.SIunits.SpecificHeatCapacity cp_SrCl2 = 0 "[J/(kg·K)]";
  constant Modelica.SIunits.SpecificHeatCapacity[:] cp_salt = {
    cp_NaCl,
    cp_KCl,
    cp_CaCl2,
    cp_MgCl2,
    cp_SrCl2};

  constant Modelica.SIunits.MolarInternalEnergy Delta_h_solution_NaCl = -3880
    "[J/mol_NaCl]";
  constant Modelica.SIunits.MolarInternalEnergy Delta_h_solution_KCl = -17000
    "[J/mol_KCl] http://webserver.dmt.upm.es/~isidoro/dat1/Heat%20of%20solution%20data.htm";
  constant Modelica.SIunits.MolarInternalEnergy Delta_h_solution_CaCl2 = 82900
    "[J/mol_CaCl2] http://webserver.dmt.upm.es/~isidoro/dat1/Heat%20of%20solution%20data.htm";
  constant Modelica.SIunits.MolarInternalEnergy Delta_h_solution_MgCl2 = 0
    "[J/mol_MgCl2]";
  constant Modelica.SIunits.MolarInternalEnergy Delta_h_solution_SrCl2 = 0
    "[J/mol_SrCl]";
  constant Modelica.SIunits.MolarInternalEnergy[:] Delta_h_solution = {
    Delta_h_solution_NaCl,
    Delta_h_solution_KCl,
    Delta_h_solution_CaCl2,
    Delta_h_solution_MgCl2,
    Delta_h_solution_SrCl2};

  constant Modelica.SIunits.SpecificEnthalpy[:] h_salt_ref = {286935,0,0,0,0}
    "salt enthalpies @ T_ref";
  constant Modelica.SIunits.Temp_K[:] T_ref_salt = {373.16,0,0,0,0};

  Modelica.SIunits.SpecificEnthalpy h_H2O =  Modelica.Media.Water.WaterIF97_base.specificEnthalpy_pT(p, T);
  Modelica.SIunits.SpecificEnthalpy[5] h_salt = h_salt_ref - (T_ref_salt .- T).*cp_salt
    "J/mol_salt solid enthalpy";

  Partial_Units.Molality mola[size(X,1)]=massFractionsToMolalities(X,cat(1,MM_vec_salt,fill(-1,size(X,1)-size(MM_vec_salt,1))));
  Partial_Units.Molality mola_salt[5]=mola[1:5];
  Modelica.SIunits.Temp_C T_C = Modelica.SIunits.Conversions.to_degC(T);
  Pressure_bar p_bar=Modelica.SIunits.Conversions.to_bar(p);

//  Modelica.Media.Water.WaterIF97_base.ThermodynamicState state_H2O;
//  Modelica.SIunits.MolarMass M_Solution "[kg/mol]";
//  Modelica.SIunits.Pressure p_check;

algorithm
Modelica.Utilities.Streams.print("h_H2O: "+String(h_H2O)+" J/kg");
Modelica.Utilities.Streams.print("h_salt: "+String(h_salt[1]) +" J/kg");
//Modelica.Utilities.Streams.print("cp_salt: "+String(size(cp_salt,1))+"");
//Modelica.Utilities.Streams.print("Delta_h_solution_NaCl: "+String(Delta_h_solution_NaCl)+" J/kg");
//Modelica.Utilities.Streams.print("mola_salt[NaCl]: "+String(mola_salt[NaCl])+" J/kg");

//  p_bar := Modelica.SIunits.Conversions.to_bar(p);
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

/*  state_H2O := Modelica.Media.Water.WaterIF97_base.setState_pTX(p, Modelica.SIunits.Conversions.from_degC(T_Scale_h), fill(0,0));
  h := Modelica.Media.Water.WaterIF97_base.specificEnthalpy(state_H2O);*/
  h := X[end]*(h_H2O+Delta_h_solution*mola_salt) +X[1:5]*h_salt;
//  h := h_H2O;

  Modelica.Utilities.Streams.print("Brine.specificEnthalpy_pTX_Francke: "+String(p*1e-5)+"bar."+String(T)+"°C->"+String(h)+" J/kg");
end specificEnthalpy_pTX_Francke_hsol;
