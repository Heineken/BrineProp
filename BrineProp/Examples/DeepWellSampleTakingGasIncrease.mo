within BrineProp.Examples;
model DeepWellSampleTakingGasIncrease
  "What happens when gas content is increased?"
//SPECIFY MEDIUM and COMPOSITION

//   package Medium = BrineProp.Brine5salts3gas(ignoreLimitN2_T=true);
   package Medium = BrineProp.Brine5salts3gas(AssertLevel=2);

//DEFINE BRINE COMPOSITION (NaCl, KCl, CaCl2, MgCl2, SrCl2, CO2, N2, CH4)
  Real[Medium.nXi] Xi = {f*0.0839077010751,f*0.00253365118988,f*0.122786737978,0,0,g*7.2426359111e-05,0*g*0.000689505657647,0*g*6.14906384726e-05}
    "GrSk brine (Feldbusch 2-2013 1.1775g/ml V2)";

/*  package Medium = BrineProp.Water_MixtureTwoPhase_pT;
    Real[Medium.nXi] Xi= fill(0,0);*/

  Medium.BaseProperties props;
  Real GVF = props.GVF;

  Real ratioGasLiquid_STP = V_l*props.d_l/1000
    "<- PLOT ME! gas-liquid volume ratio at standard conditions, fully degassed";

  Real g=time;
  Real f=1;
protected
  Modelica.SIunits.MolarVolume Vm0 = Modelica.Constants.R*273.15/101325
    "molar volume at STP";
  Modelica.SIunits.MolarVolume V_l = sum(props.X_l[6:8]./Medium.MM_gas)*Vm0*1000/props.X_l[end]
    "Liter of dissolved gas per kg_brine would have after complete degassing at standard conditions";

equation
  //SPECIFY THERMODYNAMIC STATE
/*  //degassing by heating starting at STP
  props.p = 1.01325e5;
  props.T = 20+273.15+time;
*/
  //degassing by decompression starting at reservoir conditions
  props.p = 200*1e5;
  props.T = 120+273.15;

  //specify brine composition
  props.Xi = Xi;
  annotation (experiment(
      StartTime=1,
      StopTime=300,
      __Dymola_NumberOfIntervals=100),
      __Dymola_experimentSetupOutput,
    Documentation(info="<html>
<p>This calculation aims at reproducing taking a fluid sample in a deep geothermal well.</p>
<p>The sample is taken in an stationary (idle for long time) well at 200 bar, 120 C (only liquid phase). Then the content in dissolved gases is determined by completely degassing (using ultrasonic) at standard conditions (STP) and measuring the gas volume (neglecting water vapour).</p>
<p>The fluid is assumed to be in-situ saturated with gas and in thermal equilibrium with the formation, this situation having been created by stopping the production of reservoir fluid of a given composition, which was saturated or even contained free gas, at the respective depth in the fluid column, where it degassed until saturation. Hence in the calculation fluid with an initial mass composition is subjected to in-situ conditions, the VLE is calculated and the liquid part is extracted. Its gas mass is converted to a gas volume and to a gas-liquid volume ratio (GLVR) at STP.</p>
<p>If this is done over the whole depth of the well, a profile of the saturated GLVR can be plotted. Comparing this to the measured values, however, shows apparently over- and undersaturated, which cannot be explained by the process above.</p>
<h4><span style=\"color:#008000\">Observation</span></h4>
<p>[Plot <code>ratioGasLiquid_STP</code> vs. <code>g</code>]</p>
<p>In order to show the effect of varying volumes of free gas in this model the original gas content is multiplied with an increasing factor (1...300) <i><b>g</b></i>. The resulting GLVR are plotted and show an increase far beyond the initial value, reaching a maximum and decreasing.</p>
<p>[Plot <code>props.p_gas[1...4]</code> vs. <code>g</code>]</p>
<p>Plotting the partial pressures shows that CO2 and CH4 replace N2 in the gas phase in-situ. p_CO2 reaches a maximum and decreases again. As partial pressures correspond to dissolved gas concentrations the same behaviour can be seen there.</p>
<p>The partial pressure of water is nearly constant because it is tied to the saturation pressure which only depends on temperature and salination.</p>
<p>[Plot <code>GVF</code> vs. <code>g</code>]</p>
<p>The in-situ gas volume fraction obviously increases when more gas is present.</p>
<h4><span style=\"color:#008000\">Explanation</span></h4>
<p>For low initial gas contents most of the CO2 dissolves in the liquid phase due to its low content and good solubility (50 times the one of N2). What is dissolved is removed from the gas phase reducing partially pressure and consequentially dissolved concentration. When more total gas mass is present, the CO2 partial pressure remains higher inspite of increased dissolved mass. The same holds for CH4 to a lower extent. The effect on N2 is inverse, because it dissolves only to a small amount and is amply present in the gas phase already initially.</p>
<p>As the partial pressure of gas is constant and the in-situ gas fraction increases more water mass is in the gas phase. This increases salinity in the liquid phase, in turn reducing gas solubilities.</p>
<p>Consequentially:</p>
<ul>
<li>Without salt, the GLVR keeps rising for increased gas multiplier <i><b>g</b></i>.</li>
<li>If only one gas is present the GLVR decreases by salting-out once a gas phase is present for increased gas multiplier <i><b>g</b></i>. <br>The maximum values for 200 bar, 120 &deg;C are: CO2 12:1, N2 0.75:1, CH4 1,5:1</li>
</ul>
<p><br>This means that all for saturated GLVR al values between 0.75 and 12:1 are possible for a free combination of CO2 and N2. (Lower if undersaturated.)</p>
</html>"));
end DeepWellSampleTakingGasIncrease;
