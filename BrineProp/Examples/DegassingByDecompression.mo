within BrineProp.Examples;
model DegassingByDecompression
  "Degassing example for 2-phase brine property model"
  //needs "Advanced.PedanticModelica:=false" to run
  //should be run with 500 steps on equidistant time grid

//SPECIFY MEDIUM and COMPOSITION
  package Medium = BrineProp.Brine3salts3gas(AssertLevel=2);
//DEFINE BRINE COMPOSITION (NaCl, KCl, CaCl2, CO2, N2, CH4)
  Real[Medium.nXi] Xi = {0.0839077010751,0.00253365118988,0.122786737978,7.2426359111e-05,0.000689505657647,6.14906384726e-05}
    "GrSk brine (Feldbusch 2-2013 1.1775g/ml V2)";

  Medium.BaseProperties props;

  Real GVF = props.GVF "<- PLOT ME!";
  Real p_bar = props.p/1e5 "pressure in bar (For Plotting)";
equation
  //SPECIFY THERMODYNAMIC STATE
  //degassing by decompression starting at reservoir conditions
  props.p = (100-98.5*time)*1e5;
  props.T = 125+273.15;

  //specify brine composition
  props.Xi = {0.0839077010751,0.00253365118988,0.122786737978,7.2426359111e-05,0.000689505657647,6.14906384726e-05}
    "GrSk brine (Feldbusch 2-2013 1.1775g/ml V2)";
  annotation (__Dymola_Commands(file="Resources/Scripts/DegassingByDecompression.mos"
        "Plot degassing"), experiment(__Dymola_NumberOfIntervals=100));
end DegassingByDecompression;
