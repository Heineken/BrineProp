within BrineProp.Examples;
model BrineProps1PhasePureNaCl "Pure salt example - only h and cp"
  //needs "Advanced.PedanticModelica:=false" to run
  // Driesner EOS support also pure NaCl density, but as it is used here only for h and cp only those values can be calculatesd

  package Medium = Brine3salts(ignoreLimitSalt_p={false,true,true})
    "specify medium";
//  Medium.BaseProperties props;
  SI.SpecificEnthalpy h=Medium.specificEnthalpy_pTX(1013250,20+273.15,{1,0,0,0});
  SI.SpecificHeatCapacity c_p_brine= Medium.specificHeatCapacityCp(Medium.ThermodynamicState(p=1013250,T=20+273.15,X={1,0,0,0},d=-1))
    "specific heat capacity";
equation
//specify thermodynamic state
/*  props.p = 100e5;
  props.T = 245+273.15;
  
//SPECIFY MEDIUM COMPOSITION {NaCl, KCl, CaCl2, MgCl2, SrCl2}
  props.Xi = {1,0,0} "pure NaCl";
*/
end BrineProps1PhasePureNaCl;
