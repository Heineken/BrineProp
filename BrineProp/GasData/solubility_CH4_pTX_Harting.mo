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
  String msg;
algorithm
// print("mola_CH4("+String(p_gas)+","+String(T-273.16)+") (solubility_CH4_pTX_Duan2006)");
  if 273>T or  T>400 then
    msg:="T=" + String(T) + " K, but CH4 solubility calculation is only valid for temperatures between 0 and 127degC (Partial_Gas_Data.solubility_CH4_pTX_Harting)";
  end if;
  if (p<1e5 or p>600e5) then
    msg:="p=" + String(p/1e5) + " bar, but CH4 solubility calculation only valid for pressures between 1 and 600 bar (Partial_Gas_Data.solubility_CH4_pTX_Harting)";
  end if;
  if molalities[1]>6 then
    msg :="mola[NaCl]=" + String(molalities[1]) + " mol/kg, but N2 solubility calculation only valid for salinities up to 6 mol/kg (Partial_Gas_Data.solubility_CH4_pTX_Harting)";
  end if;

  if msg<>"" then
    if outOfRangeMode==1 then
     print(msg);
    elseif outOfRangeMode==2 then
     assert(false,msg);
    end if;
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
