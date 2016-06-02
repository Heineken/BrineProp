within BrineProp.Resistivities;
function Resistivity_Francke_TXd "mixture ansatz"
  //Ucok1980 function + mixture ansatz

  input SI.Temperature T;
  input SI.MassFraction X[:] "mass fraction m_NaCl/m_brine";
//  input Types.Molarity[:] c_vec "molar concentration (molarity) in mol / litre solution";
  input SI.Density d;
  input SI.MolarMass MM_vec[:];
  input Integer AssertLevel=2
    "when out of validity range: 0-do nothing, 1-show warnings, 2-throw error";
  output SI.Resistivity rho;
protected
  SI.Conductivity[3] gamma_vec=Conductivity_Ucok1980_Tcd(T,c_vec_total,d,AssertLevel);
  Integer n_vec[:] = {1,1,2};
  SI.Conductivity gamma;
  Types.Molarity c_vec[:]=X[1:end-1] ./ MM_vec[1:end-1]*d/1000;
  Types.Molarity c_vec_total[:]=fill(sum(c_vec), size(c_vec,1))
    "Total molality to be passed to conductivity function";
  constant AssertionLevel aLevel = if AssertLevel==1 then AssertionLevel.warning else AssertionLevel.error;
/* constant c_min = 0.03 "minimum molar concentration";
 constant c_max = 0.26 "maximum molar concentration";*/
  SI.MassFraction X_total "total salinity m_XCl/m_brine";
algorithm
  print("c_vec="+Modelica.Math.Matrices.toString(transpose([c_vec])));
  print("c_vec_total="+Modelica.Math.Matrices.toString(transpose([c_vec_total])));
//  print("MM_vec="+Modelica.Math.Matrices.toString(transpose([MM_vec])));
  print("gamma_vec="+Modelica.Math.Matrices.toString(transpose([gamma_vec])));

  if AssertLevel>0 then
    assert(X[end]<1, "#infinite resistivity for pure water",aLevel);
    assert(T-273.15 >= 22 and T-273.15 <= 375,"Temperature  T=" + String(T-273.15) + " C is out of validity range.\nTo ignore validity range call with AssertLevel=1 or 0",aLevel);
    for i in 1:3 loop
      X_total :=c_vec_total[i]*MM_vec[i];
      assert(X_total==0 or (X_total >= 0.03 and X[i] <= 0.26),"total salinity X="+String(X_total) + " is out of validity range.\nTo ignore validity range call with AssertLevel=1 or 0",aLevel);
    end for;
  end if;

  //mixture ansatz: average weighted with molar concentration and ion valence
  gamma := c_vec.* n_vec/sum(c_vec.* n_vec) *gamma_vec;
  rho:=1/gamma;
end Resistivity_Francke_TXd;
