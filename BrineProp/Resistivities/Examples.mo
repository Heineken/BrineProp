within BrineProp.Resistivities;
package Examples
  extends Modelica.Icons.ExamplesPackage;

  model ResistivityUcok1980
    //needs "Advanced.PedanticModelica:=false" to run

    parameter SI.Temp_C T_C=225;
    SI.Temperature T=T_C+273.15;
    parameter SI.MassFraction x_NaCl = 0.03;
    parameter SI.MassFraction x_KCl = 0.03;
    parameter SI.MassFraction x_CaCl2 = 0;

    package Medium = Brine3salts(ignoreLimitSalt_p={false,true,true},ignoreLimitSalt_T={false,true,true},ignoreLimitSalt_b={true,true,true})
      "specify medium";
    Medium.BaseProperties props;
    SI.Resistivity rho=Resistivity_Francke_TXd(props.T,props.X,props.d,Medium.MM_vec,Medium.AssertLevel);

  equation
  //specify thermodynamic state
    props.p =1e5;
    props.T = T;

  //SPECIFY MEDIUM COMPOSITION {NaCl, KCl, CaCl2}
    props.Xi = {x_NaCl,x_KCl,x_CaCl2};
  end ResistivityUcok1980;
end Examples;
