within BrineProp.Partial_Gas_Data;
function solubility_CH4_pTX_Harting "..."
//from Harting1981, considering only NaCl
  extends partial_solubility_pTX;

/*  input Modelica.SIunits.Pressure p;
  input Modelica.SIunits.Temp_K T;
  input Modelica.SIunits.MassFraction X[:] "mass fractions m_x/m_H2O";
  input Modelica.SIunits.MolarMass MM[:] "molar masses of components";
  output Modelica.SIunits.MassFraction c_gas "gas concentration in kg_gas/kg_H2O";
*/
protected
  Molality molalities[size(X,1)]=massFractionsToMolalities(X,MM_vec);
  Modelica.SIunits.Temp_C T_C = Modelica.SIunits.Conversions.to_degC(T);
  Real L_0=.454 "CH4-Löslichkeit in H2O bei 25 atm 75°C";
  Real L_rel_p "pressure influence";
  Real L_rel_c "salinity influence";
  Real L_rel_T "Temperature influence";
//  Real c = X[1]/Salt_Data.M_NaCl/X[end];
  Real c = sum(molalities[1:2])+sum(molalities[3:5])*1.8;
  Real p_atm = p_gas/101325;
algorithm
// Modelica.Utilities.Streams.print("mola_CH4("+String(p_gas)+","+String(T-273.16)+") (solubility_CH4_pTX_Duan2006)");

//  assert(273<=T and T<=400, "T="+String(T)+" K, but N2 solubility calculation is only valid for temperatures between 0 and 127°C");
 if 273>T or  T>400 then
    Modelica.Utilities.Streams.print("T="+String(T)+" K, but N2 solubility calculation is only valid for temperatures between 0 and 127°C (Partial_Gas_Data.solubility_N2_pTX_Duan2006())");
  end if; /**/
  assert(p>=1e5 and p<=600e5, "p="+String(p/1e5)+" bar, but N2 solubility calculation only valid for pressures between 1 and 600 bar");
//  assert(molalities[NaCl]<6, "mola[NaCl]="+String(molalities[NaCl])+" mol/kg, but N2 solubility calculation only valid for salinities up to 6 mol/kg");
  if molalities[NaCl]>6 then
    Modelica.Utilities.Streams.print("mola[NaCl]="+String(molalities[NaCl])+" mol/kg, but N2 solubility calculation only valid for salinities up to 6 mol/kg");
  end if;

//page 19
  L_rel_p :=(.4009*p_atm - 7.454e-4*p_atm^2 + 5.985e-7*p_atm^3)/100*10.4537157651017 "CH4";

  L_rel_c :=exp(-.315*c +.01452 *c^2) "S. 15";

  L_rel_T := -0.0000003493*(T_C-75)^3 + 0.0001054*(T_C-75)^2 - 0.000293*(T_C-75) + 1.01
    "fitted with Excel";
  solu :=L_0*L_rel_p*L_rel_c*L_rel_T "l/kg_H2O";

  X_gas :=solu/22.4*M_CH4*X[end];
/*  Modelica.Utilities.Streams.print("L_rel_p="+String(L_rel_p)+" (solubility_CH4_pTX_Harting)");
  Modelica.Utilities.Streams.print("L_rel_c(c="+String(c)+")="+String(L_rel_c)+" (solubility_CH4_pTX_Harting)");
  Modelica.Utilities.Streams.print("L_rel_T="+String(L_rel_T)+" (solubility_CH4_pTX_Harting)");
  Modelica.Utilities.Streams.print("mola_CH4("+String(p_gas)+","+String(T-273.16)+","+String(molalities[1])+")="+String(solu)+" (solubility_CH4_pTX_Harting)");
  */
end solubility_CH4_pTX_Harting;
