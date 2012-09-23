within ;
package Brine 
 import SI = Modelica.SIunits;
 constant Boolean debugmode = false "print messages in functions";
 constant Boolean ignoreLimitN2_T=true;
 constant Boolean[5] ignoreLimitSalt_p={false,true,false,false,false};
 constant Boolean[5] ignoreLimitSalt_visc={false,false,true,false,false};
 constant Integer outOfRangeMode=2
  "when out of validity range: 0-do nothing, 1-show warnings, 2-throw error";

  constant Integer NaCl=1 "reference number";
  constant Integer KCl=2 "reference number";
  constant Integer CaCl2=3 "reference number";
  constant Integer MgCl2=4 "reference number";
  constant Integer SrCl2=5 "reference number";

  constant Modelica.SIunits.MolarMass M_H2O = 0.018015 "[kg/mol] TODO";


  replaceable function massFractionsToMolalities
  "Calculate molalities (mole_i per kg H2O) from mass fractions X"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.MassFraction X[:] "Mass fractions of mixture";
    input Modelica.SIunits.MolarMass MM[:] "molar masses of components";
    output Partial_Units.Molality molalities[size(X, 1)]
    "Molalities moles/m_H2O";
  //  Real invMMX[size(X, 1)] "inverses of molar weights";
  //  Modelica.SIunits.MolarMass Mmix "molar mass of mixture";
  algorithm
   assert(size(X, 1)==size(MM, 1), "Inconsistent vectors for mass fraction("+String(size(X, 1))+") and molar masses("+String(size(MM, 1))+")");
  // Modelica.Utilities.Streams.print(String(size(X,1))+" "+String(X[end]));

    if X[end]>0 then
      for i in 1:size(X, 1) loop
    // Modelica.Utilities.Streams.print("MM["+String(i)+"]="+String(MM[i]));
    //   Modelica.Utilities.Streams.print("X["+String(i)+"]="+String(X[i]));
    //    molalities[i] := if X[end]>0 then X[i]/(MM[i]*X[end]) else -1;
        molalities[i] := X[i]/(MM[i]*X[end]);
      end for;
    else
       molalities:=fill(-1, size(X, 1));
    end if;
    annotation(smoothOrder=5);
  end massFractionsToMolalities;



  annotation (Documentation(info="<html>
<p>
<b>Brine</b> is a package that provides properties of a specified brine, i.e. an aqueous solution of salts and gases, with a potential gas phase, therefore
 including de/gassing and evaporation/condensation.
It is based on an extension to and therefore largely compatible to the Modelica.Media library. This necessary extension is PartialMixtureTwoPhaseMedium 
(not included in this package).
This package has been developed and tested in Dymola up to 2012 FD01.
</p>
<p>
All files in this library, including the C source files are released under the Modelica License 2.
</p>

<h2>Installation</h2>
<p>
This package needs the package MediaTwoPhaseMixture.PartialMixtureTwoPhaseMedium which is not included in this package. It is found in <a href=\"https://github.com/Heineken/REFPROP2Modelica\">
REFPROP2Modelica</a>.
This provided, the examples under Examples should work right away.
</p>

<h2>Usage</h2>
As it is based on Modelica.Media, the usage is little different from the usage of the two-phase water model:<br/>
Create an Instance of the Medium:
<pre>
  package Medium = Brine_Duan_Multi_TwoPhase_ngas_3;
</pre>
Create an Instance of Medium.Baseproperties:
<pre>
  Medium.BaseProperties props;
</pre>
You can then use the BaseProperties model to define the actual brine composition(Xi or X), to define the thermodynamic state and calculate the corresponding properties.
<pre>
  props.p = 1e5;
  props.T = 300;
  props.Xi = {0.08, 0.004, 0.12, 0.001, 0.002, 1-4, 7e-4, 6e-005} \"NaCl, KCl, CaCl2, MgCl2, SrCl2, CO2, N2, CH4\";
  d = props.d;
</pre>
<p>Pressure and temperature as well pressure and specific enthalpy can be used to define a thermodynamic state. 
</p>
<p>All calculated values are returned in SI-Units and are mass based.
</p>


<h2>Details</h2>
  Brine is a mixture of components and, if gases are involved, potentially has two phases. As this is not possible with the components provided in Modelica.Media
  a new Medium template had to be created by merging Modelica.Media.Interfaces.PartialMixtureMedium and Modelica.Media.Interfaces.PartialTwoPhaseMedium of the 
  Modelica Standard Library 3.1.
  <!-- Alternatively, there is a version of this package limited to single phase fluid (Brine_Duan_Multi) which uses the  
  template PartialBrine_Multi based on the standard template Modelica.Media.Interfaces.PartialMixtureMedium. -->
  <br>
  The model is explicit for p and T, but for h(p,T) the inverse function T(p,h) is defined. T(p,h) is inverts h(p,T) numerically by bisection, stopping at a given tolerance.<br>
  In order to calculate h(p,T), the vapour-liquid-equilibrium (VLE) is determined, i.e. the gas mass fraction q and the compositions of the liquid phases X_l. 
  Only h is returned, due to the limitation of DYMOLA/Modelica not allowing inverse functions of functions that are returning an array. As q (gas mass fraction) and X_l (composition of liquid phase) are of interest themselves and
  required to calculate density and viscosity, the VLE calculation is conducted one more time, this time with T known. This additional calculation doubles the workload when 
  p,h are given. When p,T are given, however, it adds only one more calculation to the multiple iterations of the bisection algorithm.
<p>

<p>
<h2>TODO:</h2>
<ul>
<li></li>
</ul>

</p>


<h3> Created by</h3>
Henning Francke<br/>
Helmholtz Centre Potsdam<br/>
GFZ German Research Centre for Geosciences<br/>
Telegrafenberg, D-14473 Potsdam<br/>
Germany
<p>
<a href=mailto:francke@gfz-potsdam.de>francke@gfz-potsdam.de</a>
</html>
",
 revisions="<html>

</html>"), version="0.1", versionDate="2012-08-01", uses(Modelica(version="3.2"), MediaTwoPhaseMixture(version="0.2")));
end Brine;
