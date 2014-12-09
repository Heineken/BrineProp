within BrineProp.GasData;
function solubility_CH4_pTX_Harting
  //from Harting1981 (not used)
  //made for NaCl, considering other salts via conversion factors
  extends partial_solubility_pTX;

protected
  Types.Molality molalities[size(X, 1)]=
      Utilities.massToMoleFractions(X,
      MM_vec);
  SI.Temp_C T_C = SI.Conversions.to_degC(T);
  Real L_0=0.454 "CH4 solubility in H2O at 25, 75degC";
  Real L_rel_p "pressure influence";
  Real L_rel_c "salinity influence";
  Real L_rel_T "Temperature influence";
//  Real c = X[1]/Salt_Data.M_NaCl/X[end];
  Real c = sum(molalities[1:2])+sum(molalities[3:5])*1.8;
  Real p_atm = p_gas/101325;
algorithm
// print("mola_CH4("+String(p_gas)+","+String(T-273.16)+") (solubility_CH4_pTX_Duan2006)");

  if AssertLevel>0 then
     assert(ignoreLimitCH4_T or (273<T and T<400),"\nTemperature out of validity range T=" + String(T) + " K.\nTo ignore set ignoreLimitCH4_T=true",aLevel);
     assert(ignoreLimitCH4_p or (1e5<p or p<600e5),"\nPressure out of validity range: p=" + String(p/1e5) + " bar.\nTo ignore set ignoreLimitCH4_p=true",aLevel);
     assert(ignoreLimitCH4_b or molalities[1]<6,"\nMolality out of validity range: b=" + String(molalities[1]) + " mol/kg.\nTo ignore set ignoreLimitCH4_p=true",aLevel);
  end if;
//page 19
  L_rel_p :=(0.4009
                  *p_atm - 7.454e-4*p_atm^2 + 5.985e-7*p_atm^3)/100*10.4537157651017 "CH4";

  L_rel_c :=exp(-0.315
                     *c +0.01452*c^2) "S. 15";

  L_rel_T := -0.0000003493*(T_C-75)^3 + 0.0001054*(T_C-75)^2 - 0.000293*(T_C-75) + 1.01
    "fitted with Excel";
  solu :=L_0*L_rel_p*L_rel_c*L_rel_T "l/kg_H2O";

  X_gas :=solu/22.4*M_CH4*X[end];
/*  print("L_rel_p="+String(L_rel_p)+" (solubility_CH4_pTX_Harting)");
  print("L_rel_c(c="+String(c)+")="+String(L_rel_c)+" (solubility_CH4_pTX_Harting)");
  print("L_rel_T="+String(L_rel_T)+" (solubility_CH4_pTX_Harting)");
  print("mola_CH4("+String(p_gas)+","+String(T-273.16)+","+String(molalities[1])+")="+String(solu)+" (solubility_CH4_pTX_Harting)");
  */
end solubility_CH4_pTX_Harting;
