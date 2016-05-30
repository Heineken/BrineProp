within BrineProp.Resistivities;
function Resistivity_Francke_Tcd "mixture ansatz"
  //returns conductivity instead of resistivity to allow easy handling of zero salinities
  input SI.Temperature T;
  input Types.Molarity[:] c_vec
    "molar concentration (molarity) in mol / litre solution";
  input SI.Density d;
  output SI.Resistivity rho;
protected
  SI.Conductivity[3] gamma_vec=Conductivity_Ucok1980_Tcd(T,c_vec,d);
  Integer n_vec[:] = {1,1,2};
  SI.Conductivity gamma;
algorithm
  //mixture ansatz: average weighted with molar concentration and ion valence
  gamma := c_vec.* n_vec/sum(c_vec.* n_vec) *gamma_vec;
  rho:=1/gamma;
end Resistivity_Francke_Tcd;
