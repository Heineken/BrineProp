within BrineProp;
package Vacuum "Medium with zero properties"
  extends PartialBrineGas(
    final substanceNames={"carbondioxide","nitrogen","methane","water"},
    final MM_vec = {M_CO2,M_N2,M_CH4, M_H2O},
    final nM_vec = {nM_CO2,nM_N2,nM_CH4, nM_CH4});


 redeclare model extends BaseProperties
 //Dummy for OM
 end BaseProperties;

/* redeclare record extends ThermodynamicState
 //Dummy for OM
 end ThermodynamicState;
*/
  constant Boolean waterSaturated=false "activates water saturation";


  redeclare function extends density "water-saturated density from state"

  algorithm
    d := density_pTX(
      p=state.p,
      T=state.T,
      X= state.X);
  end density;


  redeclare function extends density_pTX
  "Density of an ideal mixture of ideal gases"
  algorithm
    if debugmode then
      print("Running density_pTX("+String(p/1e5)+" bar,"+String(T-273.15)+" degC, X="+Modelica.Math.Matrices.toString(transpose([X]))+")");
    end if;

    d :=0;
  end density_pTX;




  redeclare function extends dynamicViscosity
  "water-saturated  thermal conductivity of water"
  //very little influence of salinity
  algorithm
    eta := dynamicViscosity_pTX(
          p=state.p,
          T=state.T,
          X= state.X);
  //  else state.X[end - nX + 1:end]);
  //  assert(lambda>0,"lambda="+String(lambda));
  end dynamicViscosity;


  redeclare function extends dynamicViscosity_pTX
  "calculation of gas dynamic Viscosity"
  /*  import NG = Modelica.Media.IdealGases.Common.SingleGasNasa;
  input SI.Pressure p;
  input SI.Temperature T;
  input SI.MassFraction[nX] X "Mass fractions of mixture";
  output SI.DynamicViscosity eta;*/
  algorithm
    eta:=Modelica.Media.Air.MoistAir.dynamicViscosity(
      Modelica.Media.Air.MoistAir.ThermodynamicState(
      p=0,
      T=T,
      X={0,0}));
  end dynamicViscosity_pTX;


  redeclare function extends thermalConductivity
  algorithm
    lambda := thermalConductivity_pTX(
          p=state.p,
          T=state.T,
          X= state.X);

  end thermalConductivity;


  redeclare function extends thermalConductivity_pTX
  "calculation of gas thermal conductivity"
  algorithm
    lambda:=0;
  end thermalConductivity_pTX;


  redeclare function specificEnthalpy_pTX
  "calculation of specific enthalpy of gas mixture"
  //  import Modelica.Media.IdealGases.Common.SingleGasNasa;
    import Modelica.Media.IdealGases.SingleGases;
    extends Modelica.Icons.Function;
    input AbsolutePressure p "Pressure";
    input Temperature T "Temperature";
    input MassFraction X[:]=reference_X "Mass fractions";
    output SpecificEnthalpy h "Specific enthalpy";
  algorithm
    h := 0;
  end specificEnthalpy_pTX;


  redeclare function extends specificEnthalpy
  "water-saturated specific enthalpy of gas phase"
  algorithm
       h := specificEnthalpy_pTX(
          p=state.p,
          T=state.T,
          X=state.X);
  //  else state.X[end - nX + 1:end]);

  end specificEnthalpy;


  annotation (Documentation(info="<html>
<p><b>BrineGas_3Gas</b> is a medium package that, based on Brine.PartialBrineGas, defines a brine with 3 gases (CO<sub>2</sub>, N<sub>2</sub>, CH<sub>4</sub>), which are the main gases in the geofluid in Gross Schoenebeck, Germany.</p>
<h4>Usage</h4>
<p>It is based on Modelica.Media, the usage is accordingly:</p>
<p>Create an instance of the Medium: </p>
<pre>  package Medium = BrineGas_3Gas;</pre>
<p>Create an instance of Medium.Baseproperties: </p>
<pre>  Medium.BaseProperties props;</pre>
<p>Use the BaseProperties model to define the actual brine composition(Xi or X), to define the thermodynamic state and calculate the corresponding properties. </p>
<pre>  props.p = 1e5;
  props.T = 300;
  props.Xi = {1-4, 7e-4, 6e-005} \"CO2, N2, CH4\"
  d = props.d;
</pre>

<p>See <code><a href=\"Modelica://BrineProp.Examples.BrineGas\">BrineProp.Examples.BrineGas</a></code> for more usage examples.</p>
<p>Returns properties for given composition when _pTX functions are called directly.
  Returns properties for given gas composition + saturated water when called via state functions (e.g. density)
</p>  
<p>All calculated values are returned in SI units and are mass based.</p>
<h4>Potential speedup:</h4>
<p>Calculate water saturated composition externally once (instead of separately in each property function) and pass on.</p>
</html>"));
end Vacuum;
