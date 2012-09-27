within BrineProp;
partial package PartialBrine_MultiSalt_2Phase_unused "Template package for two-phase brines without property functions "

constant Modelica.Media.Interfaces.PartialMedium.FluidConstants[1]
  simpleWaterConstants(
     each chemicalFormula = "H2O+NaCl",
     each structureFormula="H2O+NaCl",
     each casRegistryNumber="007",
     each iupacName="Aqueous NaCl-Solution",
     each molarMass=0);


 extends Modelica.Media.Interfaces.PartialMixtureMedium(
   final mediumName="Geothermal Brine",
   final substanceNames={"sodium chloride","water"},
   final reducedX =  true,
   final singleState=false,
   reference_X={0.01,0.99},
   fluidConstants = waterConstants);


  extends Partial_Units;

/*
,
   AbsolutePressure(
     min=1,
     max=500),
   Temperature(
     min=283.15,
     max=623.15,
     start=323)
   final fixedX = false
     */

//   extends Partial_Salt_Data;


  replaceable package Salt_data = BrineProp.SaltData;

                              //definition of molar masses
/*
   
   constant Modelica.SIunits.MolarMass M_H2O = 0.018015 "[kg/mol]";
   constant Modelica.SIunits.MolarMass M_NaCl = 0.058443 "[kg/mol]";
   constant Modelica.SIunits.MolarMass M_KCl = 0.074551 "[kg/mol]";
   constant Modelica.SIunits.MolarMass M_CaCl = 0.1109840 "[kg/mol]";
   constant Modelica.SIunits.MolarMass M_MgCl2 = 0.095236 "[kg/mol]";
   constant Modelica.SIunits.MolarMass M_SrCl2 = 0.158536 "[kg/mol]";
*/


 redeclare model extends BaseProperties "Base properties of medium"

    //    PowerPlant.Media.TableLookup Table;
    //  protected
    /*     constant Modelica.SIunits.MolarMass M_H2O = PartialBrine.M_H2O "[kg/mol]";
     constant Modelica.SIunits.MolarMass M_NaCl = PartialBrine.M_NaCl 
        "[kg/mol]";*/
 equation
   d = density_pTX(p,T,X);
   h = specificEnthalpy_pTX(p,T,X);
 //  T = temperature_phX(p,h,X);
   u = 1 "h - p/d";
    MM = X[1]*PartialBrine_MultiSalt_2Phase_unused.Salt_data.M_NaCl + X[2]*
      BrineProp.PartialBrine_MultiSalt_2Phase_unused.M_H2O;
   R  = 8.3144/MM;

   state.p = p;
   state.T = T;
   state.X = X;

   state.s = 0 "specificEntropy_phX(p,h,X)";
   state.h = h;
   state.d = d;

   annotation (Documentation(revisions="<html>

</html>"));
 end BaseProperties;

 /* Provide implementations of the following optional properties.
     If not available, delete the corresponding function.
     The record "ThermodynamicState" contains the input arguments
     of all the function and is defined together with the used
     type definitions in PartialMedium. The record most often contains two of the
     variables "p, T, d, h" (e.g. medium.T)
  */

// redeclare replaceable record ThermodynamicState


redeclare record extends ThermodynamicState
  "a selection of variables that uniquely defines the thermodynamic state"
/*  AbsolutePressure p "Absolute pressure of medium";
  Temperature T(unit="K") "Temperature of medium";
  MassFraction X[nX] "Mass fraction of NaCl in kg/kg";*/
  SpecificEnthalpy h "Specific enthalpy";
  SpecificEntropy s "Specific entropy";
  Density d(start=300) "density";

   annotation (Documentation(info="<html>

</html>"));
end ThermodynamicState;


redeclare function extends setState_pTX
  "Return thermodynamic state of water as function of p and T"
algorithm
  state := ThermodynamicState(p= p,
                              h= specificEnthalpy_pTX(p,T,X),
                              X= X,
                              T= T,
                              s= 0,
                              d= density_pTX(p,T,X));
end setState_pTX;


redeclare function extends setState_phX
  "Return thermodynamic state of water as function of p and T"
algorithm
  state := ThermodynamicState(p= p,
                              h= h,
                              X= X,
                              T= 0,
                              s= 0,
                              d= 0);
end setState_phX;


redeclare function extends setState_psX
  "Return thermodynamic state of water as function of p and T"
algorithm
  state := ThermodynamicState(p= p,
                              h= specificEnthalpy_ps(p,s,X),
                              X= X,
                              T= temperature_psX(p,s,X),
                              s= s,
                              d= density_psX(p,s,X));
end setState_psX;


  redeclare replaceable partial function density_pTX
  "Return density from p, T, and X or Xi"
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[:] "Mass fractions";
    input Modelica.SIunits.MolarMass MM[:]={1} "molar masses of components";
    output Density d "Density";
    annotation(Documentation(info="<html></html>"));
  end density_pTX;


  redeclare function extends temperature "return temperature of ideal gas"
  algorithm
    T := state.T;
  end temperature;


  redeclare function extends specificEnthalpy "Return specific enthalpy"
    extends Modelica.Icons.Function;
  algorithm
    h := state.h;
  end specificEnthalpy;


  redeclare function extends dynamicViscosity
  algorithm
    eta :=dynamicViscosity_pTX(
        state.p,
        state.T,
        state.X);
  end dynamicViscosity;


  redeclare replaceable function specificEnthalpy_pTX
     input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
    output Modelica.SIunits.SpecificEnthalpy h;

  /*algorithm 
  h := 4*T;
*/
  end specificEnthalpy_pTX;


  replaceable function dynamicViscosity_pTX "viscosity calculation"
    input Modelica.SIunits.Pressure p;
    input Modelica.SIunits.Temp_K T;
    input MassFraction X[:] "mass fraction m_NaCl/m_Sol";
    output Modelica.SIunits.DynamicViscosity eta;
  //  constant Real M_NaCl=0.058443 "molar mass in [kg/mol]";
  end dynamicViscosity_pTX;


annotation (Documentation(info="<html><h1>PartialBrine</h1>
A base package for the Modelica.Media model for aqueous NaCl-Solutions.<br>
Gives density and viscosity.</html>"));
end PartialBrine_MultiSalt_2Phase_unused;
