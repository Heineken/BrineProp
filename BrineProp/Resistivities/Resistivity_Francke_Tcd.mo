within BrineProp.Resistivities;
function Resistivity_Francke_Tcd "mixture ansatz"
  //returns conductivity instead of resistivity to allow easy handling of zero salinities
  input SI.Temperature T;
  input Types.Molarity[:] c_vec
    "molar concentration (molarity) in mol / litre solution";
  input SI.Density d;
  input Integer AssertLevel=2
    "when out of validity range: 0-do nothing, 1-show warnings, 2-throw error";
  output SI.Resistivity rho;
protected
  constant AssertionLevel aLevel = if AssertLevel==1 then AssertionLevel.warning else AssertionLevel.error;

  SI.Conductivity[3] gamma_vec=Conductivity_Ucok1980_Tcd(T,c_vec,d);
  Integer n_vec[:] = {1,1,2};
  SI.Conductivity gamma;
algorithm
  if AssertLevel>0 then
/*    assert(ignoreLimitSalt_p[i] or (p >= salt.p_min_rho and p <= salt.p_max_rho),"Pressure p=" + String(p/1e5) + " bar is out of validity range  ["+String(salt.p_min_rho/1e5)+"..."+String(salt.p_max_rho/1e5)+"]bar for "+salt.name + ":.\nTo ignore set ignoreLimitSalt_p["+String(i)+"]=true",aLevel);
    assert(ignoreLimitSalt_T[i] or (T >= salt.T_min_rho and T <= salt.T_max_rho),"Temperature  T=" + String(T-273.15) + " C is out of validity range ["+String(salt.T_min_rho-273.15)+"..."+String(salt.T_max_rho-273.15)+"] for "+salt.name + ".\nTo ignore set ignoreLimitSalt_T["+String(i)+"]=true",aLevel);
    assert(ignoreLimitSalt_b[i] or (m[i] >= 0 and m[i] <= salt.mola_max_rho),salt.name + "Molality is out of validity range: m[i]=" + String(m[i]) + " mol/kg.\nTo ignore set ignoreLimitSalt_b["+String(i)+"]=true",aLevel);
    */
    assert(max(c_vec)>0, "#infinite resistivity for pure water",aLevel);
  end if;

  //mixture ansatz: average weighted with molar concentration and ion valence
  gamma := c_vec.* n_vec/sum(c_vec.* n_vec) *gamma_vec;
  rho:=1/gamma;
end Resistivity_Francke_Tcd;
