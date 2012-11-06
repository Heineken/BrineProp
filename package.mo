within ;
package BrineProp "Media models for p-h-T-rho-eta properties of aqueous solutions of multiple salts and gases"
import SI = Modelica.SIunits;
 constant Boolean debugmode = false "print messages in functions";
 constant Boolean ignoreLimitN2_T=true;
 constant Boolean ignoreLimitN2_p=true;
 constant Boolean[5] ignoreLimitSalt_p={false,true,true,false,false};
 constant Boolean[5] ignoreLimitSalt_visc={false,false,true,false,false};
 constant Integer outOfRangeMode=2
  "when out of validity range: 0-do nothing, 1-show warnings, 2-throw error";

  constant Integer NaCl=1 "reference number";
  constant Integer KCl=2 "reference number";
  constant Integer CaCl2=3 "reference number";
  constant Integer MgCl2=4 "reference number";
  constant Integer SrCl2=5 "reference number";

  constant Modelica.SIunits.MolarMass M_H2O = Modelica.Media.Water.waterConstants[1].molarMass
  "0.018015 [kg/mol]";


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


  replaceable function massFractionsToMoleFractions
  "Return mole_i/sum(mole_i) from mass fractions X"
    extends Modelica.Icons.Function;
    input Modelica.SIunits.MassFraction X[:] "Mass fractions of mixture";
    input Modelica.SIunits.MolarMass MM[:] "molar masses of components";
    output Partial_Units.Molality molefractions[size(X, 1)] "Molalities";
    output Partial_Units.Molality molalities[size(X, 1)]
    "Molalities moles/m_H2O";
  //  Real invMMX[size(X, 1)] "inverses of molar weights";
  //  Modelica.SIunits.MolarMass Mmix "molar mass of mixture";
  //  Real n[size(X, 1)] "number of Moles";
protected
    Real n_total;
  algorithm
   assert(size(X, 1)==size(MM, 1), "Inconsistent vectors for mass fraction("+String(size(X, 1))+") and molar masses("+String(size(MM, 1))+")");
  // Modelica.Utilities.Streams.print(String(size(X,1))+" "+String(X[end]));
  //  printVector(MM);
    for i in 1:size(X, 1) loop
  // Modelica.Utilities.Streams.print("MMX["+String(i)+"]="+String(MMX[i]));
      molalities[i] := if X[end]>0 then X[i]/(MM[i]*X[end]) else -1;
  //    n[i] := X[i]/MMX[i];
    end for;
    n_total :=sum(molalities);
    for i in 1:size(X, 1) loop
      molefractions[i] := molalities[i]/n_total;
    end for;
    annotation(smoothOrder=5);
  end massFractionsToMoleFractions;


  annotation (Documentation(info="<html>
<p><b>BrineProp</b> is a package that provides properties of a specified brine, i.e. an aqueous solution of salts and gases, with a potential gas phase, therefore including de/gassing and evaporation/condensation. It is based on an extension to and therefore largely compatible to the Modelica.Media library. This necessary extension is PartialMixtureTwoPhaseMedium (not included in this package). This package has been developed and tested in Dymola up to 2012 FD01. </p>
<p>All files in this library, including the C source files are released under the Modelica License 2. </p>
<p><b></font><font style=\"font-size: 12pt; \">Installation</b></p>
<p>The sub-package <code>PartialBrine_ngas_Newton</code> / <code>Brine_Duan_Multi_TwoPhase_ngas_3</code> needs the package <code>MediaTwoPhaseMixture.PartialMixtureTwoPhaseMedium</code> which is not included in this package. It is found in <a href=\"https://github.com/Heineken/REFPROP2Modelica\">REFPROP2Modelica</a>. This provided, the examples under Examples should work right away.</p>
<!-- <p><h4>Packages</h4></p>
<dl>
 <dt><pre>Examples</pre></dt><dd>Usage examples</dd>
 <dt><pre>PartialUnits</pre></dt><dd>Definition of additional units used in BrineProp</dd>
 <dt><pre>SaltData</pre></dt><dd>Molar masses and mole numbers of the contained salts</dd>
 <dt><pre>Densities</pre></dt><dd>Density functions</dd>
 <dt><pre>SpecificEnthalpies</pre></dt><dd>Enthalpy functions</dd>
 <dt><pre>Viscosities</pre></dt><dd>Viscosity functions</dd>
 <dt><pre>PartialBrine_MultiSalt_1Phase</pre></dt><dd>Template for one-phase (liquid) brine based on PartialMediaMixtureMedium</dd>
 <dt><pre>SaltData_Duan</pre></dt><dd>Coefficients used in Duan density calculation</dd>
 <dt><pre>Brine_Duan</pre></dt><dd>One-phase (liquid) aqueous NaCl solution using functions by Duan.</dd>
 <dt><pre>Brine_Driesner</pre></dt><dd>One-phase (liquid) aqueous NaCl solution using functions by Driesner</dd>
 <dt><pre>Brine_5salts_nogas</pre></dt><dd>One-phase (liquid) multisalt brine solution</dd>
 <dt><pre>PartialGasData</pre></dt><dd>Molar masses and ion numbers of the contained gases</dd>
 <dt><pre>PartialBrine_ngas_Newton</pre></dt><dd>Template for multisalt-multigas brines with vapor-liquid-equilibrium calculation</dd>
 <dt><pre>Brine_5salts_TwoPhase_3gas</pre></dt><dd>Brine with 5 salts and 3 gases</dd>
</dl>
-->
<p><b>Usage</b></p>
Check the (non-partial) Brine packages or <pre>Examples</pre> for instructions.

<p>All calculated values are returned in SI-Units and are mass based. </p>
<p><h4>TODO:</h4></p>
<p><ul>
<li>Add apparent molar heat capacity for (NaCl,) MgCl2 and SrCl2</li>
<li>Add multi-salt viscosity </li>
</ul></p>
<font style=\"font-size: 10pt; \">
<p><i>Created by</i><br/>
Henning Francke<br/>
Helmholtz Centre Potsdam<br/>
GFZ German Research Centre for Geosciences<br/>
Telegrafenberg, D-14473 Potsdam<br/>
Germany</p></font>
<p><a href=\"mailto:francke@gfz-potsdam.de\">francke@gfz-potsdam.de</a> </p>
</html>",
 revisions="<html>

</html>"), version="0.1", versionDate="2012-08-01", uses(Modelica(version="3.2"), MediaTwoPhaseMixture(version="0.2")));
end BrineProp;
