within BrineProp.Examples;
model BrineProps1PhaseMinimal
  "Minimal example for 1-phase brine property model"
  package Medium = Brine5salts(ignoreLimitSalt_p={false,true,true,false,false})
    "specify medium";
  Medium.BaseProperties props;
equation
//specify thermodynamic state
  props.p = 100e5;
  props.T = 245+273.15;

//SPECIFY MEDIUM COMPOSITION {NaCl, KCl, CaCl2, MgCl2, SrCl2}
  props.Xi = {0.0839077010751,0.00253365118988,0.122786737978,0,0}
    "Feldbusch 2-2013 1.1775g/ml V2";
end BrineProps1PhaseMinimal;
