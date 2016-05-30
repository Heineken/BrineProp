within BrineProp.Resistivities;
package Examples
  extends Modelica.Icons.ExamplesPackage;

  model ResistivityUcok1980
    //needs "Advanced.PedanticModelica:=false" to run

    parameter SI.Temp_C T_C=25;
    SI.Temperature T=T_C+273.15;
    parameter SI.MassFraction x_NaCl = 0;
    parameter SI.MassFraction x_KCl = 0;
    parameter SI.MassFraction x_CaCl2 = 0.03;

    package Medium = Brine3salts(ignoreLimitSalt_p={false,true,true},ignoreLimitSalt_T={false,true,true})
      "specify medium";
    Medium.BaseProperties props;
    Types.Molarity c[:]=props.X[1:end-1] ./ Medium.MM_vec[1:end-1]*props.d/1000;

  //  SI.Conductivity gamma[:]=Conductivity_Ucok1980_Tbd(props.T,c,props.d);
    SI.Resistivity rho=Resistivity_Francke_Tcd(props.T,c,props.d);
  //  SI.Resistivity Resistivity_Ucok1980_TdX(props.T,props.d, Utilities.massFractionsToMolalities(props.X,Medium(MM_vec)));

  equation
  //specify thermodynamic state
    props.p =1e5;
    props.T = T;

  //SPECIFY MEDIUM COMPOSITION {NaCl, KCl, CaCl2}
    props.Xi = {x_NaCl,x_KCl,x_CaCl2};
  end ResistivityUcok1980;
end Examples;
