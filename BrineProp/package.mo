within ;
package BrineProp "Media models for p-h-T-rho-eta properties of aqueous solutions of multiple salts and gases"

 import SI = Modelica.SIunits;
 import Modelica.Utilities.Streams.print;

 constant Boolean debugmode = false "print messages in functions";
 constant Boolean ignoreLimitN2_T=true;
 constant Boolean ignoreLimitN2_p=true;
 constant Boolean ignoreLimit_h_KCl_Tmin=true
  "ignore Tmin in appMolarEnthalpy_KCl_White and appMolarHeatCapacity_KCl_White";
 constant Boolean ignoreLimit_h_CaCl2_Tmin=true
  "ignore Tmin in appMolarEnthalpy_CaCl2_White and appMolarHeatCapacity_CaCl2_White";
 constant Boolean[5] ignoreLimitSalt_p={false,true,true,false,false};
 constant Boolean[5] ignoreLimitSalt_T={false,false,false,false,false};
 constant Boolean[5] ignoreLimitSalt_visc={false,false,true,false,false};
 constant Integer outOfRangeMode=2
  "when out of validity range: 0-do nothing, 1-show warnings, 2-throw error";

  constant Integer NaCl=1 "reference number TODO: does not belong here";
  constant Integer KCl=2 "reference number";
  constant Integer CaCl2=3 "reference number";
  constant Integer MgCl2=4 "reference number";
  constant Integer SrCl2=5 "reference number";

  constant SI.MolarMass M_H2O = Modelica.Media.Water.waterConstants[1].molarMass
  "0.018015 [kg/mol]";


  function massFractionsToMolalities
  "Calculate molalities (mole_i per kg H2O) from mass fractions X"
    extends Modelica.Icons.Function;
    input SI.MassFraction X[:] "Mass fractions of mixture";
    input SI.MolarMass MM_vec[:] "molar masses of components";
    output Partial_Units.Molality molalities[size(X, 1)]
    "Molalities moles/m_H2O";

  algorithm
   assert(size(X, 1)==size(MM_vec, 1), "Inconsistent vectors for mass fraction("+String(size(X, 1))+") and molar masses("+String(size(MM_vec, 1))+")");

    if X[end]>0 then
      for i in 1:size(X, 1) loop
        molalities[i] := if X[i]<1e-6 then 0 else X[i]/(MM_vec[i]*X[end])
        "numerical errors my create X[i]>0, this prevents it";
      end for;
    else
       molalities:=fill(-1, size(X, 1));
    end if;
    annotation(smoothOrder=5);
  end massFractionsToMolalities;


  function massFractionsToMoleFractions
  "Return mole_i/sum(mole_i) from mass fractions X"
    extends Modelica.Icons.Function;
    input SI.MassFraction X[:] "Mass fractions of mixture";
    input SI.MolarMass MM_vec[:] "molar masses of components";
    output Partial_Units.Molality molefractions[size(X, 1)] "Molalities";
    output Partial_Units.Molality molalities[size(X, 1)]
    "Molalities moles/m_H2O";
protected
    Real n_total;
    Integer n=size(X, 1);
  algorithm
   assert(n==size(MM_vec, 1), "Inconsistent vectors for mass fraction("+String(n)+") and molar masses("+String(size(MM_vec, 1))+")");
  // print(String(size(X,1))+" "+String(X[end]));
  //  printVector(MM);
    for i in 1:n loop
  // print("MMX["+String(i)+"]="+String(MMX[i]));
      molalities[i] := if X[end]>0 then X[i]/(MM_vec[i]*X[end]) else -1;
  //    n[i] := X[i]/MMX[i];
    end for;
    n_total :=sum(molalities);
    for i in 1:n loop
      molefractions[i] := molalities[i]/n_total;
    end for;
    annotation(smoothOrder=5);
  end massFractionsToMoleFractions;


  replaceable function Xi2X "calculates the full mass vector X from Xi"
    extends Modelica.Icons.Function;
    input SI.MassFraction Xi[:] "Mass fractions of mixture";
    output SI.MassFraction X[size(Xi,1)+1] "Mass fractions of mixture";
  algorithm
    X:=cat(
      1,
      Xi,
      {1 - sum(Xi)});

  end Xi2X;


  annotation (Documentation(info="<html>
<p><b>BrineProp</b> is a modelica package that calculates the thermodynamic properties of a specified brine, i.e. an aqueous solution of salts and gases, with a potential gas phase, including degassing/evaporation and solution/condensation.</p>
<p>It was developed as a part of a PhD projected, documented in the thesis &QUOT;<a href=\"http://nbn-resolving.de/urn:nbn:de:kobv:83-opus4-47126\">Thermo-hydraulic model of the two-phase flow in the brine circuit of a geothermal power plant</a>&QUOT;. </p>
<p>This package has been developed and tested in Dymola up to 2014 FD01.</p>
<p><b>Licensed by the </b>Henning Francke<b> under the Modelica License 2</b></p>
<p>Copyright &copy; 2009-2014 Helmholtz Centre Potsdam, GFZ German Research Centre for Geosciences.</p>
<p><br><i>This Modelica package is <u>free</u> software and the use is completely at <u>your own risk</u>; it can be redistributed and/or modified under the terms of the Modelica License 2. For license conditions (including the disclaimer of warranty) see <a href=\"modelica://Modelica.UsersGuide.ModelicaLicense2\">Modelica.UsersGuide.ModelicaLicense2</a> or visit <a href=\"http://www.modelica.org/licenses/ModelicaLicense2\">http://www.modelica.org/licenses/ModelicaLicense2</a>.</i> </p>
<h4><span style=\"color:#008000\">Installation</span></h4>
<p>Make sure the package directory is named BrineProp. This provided, the examples under <code>Examples</code> should work right away. </p>
<h4>Usage</h4>
<p>Check the (non-partial) Brine packages or <code>BrineProp/Examples </code>for instructions. </p>
<p>All calculated values are returned in SI-Units and are mass based. </p>
<h4>TODO:</h4>
<ul>
<li>Add differentials</li>
<li>make 1phase-transient calculation work</li>
<li>Add apparent molar heat capacity/enthalpy for (NaCl,) MgCl2 and SrCl2</li>
<li>Make it work in JModelica/OpenModelica </li>
</ul>
<h5>Created by</h5>
<p>Henning Francke</p><p>Helmholtz Centre Potsdam GFZ German Research Centre for Geosciences</p><p>Telegrafenberg, D-14473 Potsdam</p><p>Germany</p>
<p><a href=\"mailto:info@xrg-simulation.de\">francke@gfz-potsdam.de</a></p>
</html>", revisions="<html>

</html>"),
    version="0.2.0",
    versionDate="2014-04-02",
    uses(Modelica(version="3.2"), MediaTwoPhaseMixture(version="0.2")));
end BrineProp;
