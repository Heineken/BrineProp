within BrineProp.Examples;
model SalinityFromDensity
  "Invert density function to determine NaCl content"
  //needs "Advanced.PedanticModelica:=false" to run

  package Medium = Brine3salts "specify medium";
  Medium.BaseProperties props;
  Real X_NaCl; //(min=0,max=0.25,start=0);
equation
//specify thermodynamic state
  props.p = 1e5;
  props.T = 45+273.15;
  props.d = 1048;
//SPECIFY MEDIUM COMPOSITION {NaCl, KCl, CaCl2}
//  props.Xi = {min(0.25,X_NaCl),0,0} this limitation is ignored by the solver
  props.Xi = {X_NaCl,0,0};
end SalinityFromDensity;
