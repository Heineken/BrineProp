within BrineProp;
package Brine5salts3gas "Two-phase aqueous solution of NaCl, KCl, CaCl2, MgCl2, SrCl2, N2, CO2, CH4"

//TODO: use Fluid limits


  extends PartialBrineMultiSaltMultiGasTwoPhase(
    redeclare package Salt_data = BrineProp.SaltDataDuan,
    final gasNames = {"carbondioxide","nitrogen","methane"},
    final saltNames = {"sodium chloride","potassium chloride","calcium chloride","magnesium chloride","strontium chloride"},
    final MM_gas = {M_CO2,M_N2,M_CH4},
    final nM_gas = {nM_CO2,nM_N2,nM_CH4},
    final MM_salt = Salt_data.MM_salt,
    final nM_salt = Salt_data.nM_salt);



  redeclare function extends setState_pTX "to avoid check error"
  end setState_pTX;


  redeclare function extends setState_phX "to avoid check error"
  end setState_phX;


  redeclare function extends solubilities_pTX
  "solubility calculation of CO2 in seawater Duan, Sun(2003), returns gas concentration in kg/kg H2O"
  algorithm
  //  print("p="+String(p)+" bar, T=("+String(T)+") (solubilities_pTX)");
  //  if gasname =="carbondioxide" then
      solu[1] := if X[nX_salt+1]>0 then solubility_CO2_pTX_Duan2006(p,T,X_l,MM_vec,p_gas[1]) else -1
    "aus Partial_Gas_Data, mol/kg_H2O -> kg_CO2/kg_H2O";
  //  elseif gasname =="nitrogen" then
      solu[2] := if X[nX_salt+2]>0 then solubility_N2_pTX_Duan2006(p,T,X_l,MM_vec,p_gas[2]) else -1
    "aus Partial_Gas_Data, mol/kg_H2O -> kg_N2/kg_H2O";
  //    solu[2] := if X[nX_salt+2]>0 then solubility_N2_pTX_Harting(p,T,X_l,MM_vec,p_gas[2]) else -1
  //  elseif gasname =="methane" then
      solu[3] := if X[nX_salt+3]>0 then solubility_CH4_pTX_Duan2006(p,T,X_l,MM_vec,p_gas[3]) else -1
    "aus Partial_Gas_Data, mol/kg_H2O -> kg_CH4/kg_H2O";
  //    solu[3] := if X[nX_salt+3]>0 then solubility_CH4_pTX_Harting(p,T,X_l,MM_vec,p_gas[3]) else -1
  //  end if;

  //  print("p_gas={"+String(p_gas[1])+", "+String(p_gas[2])+", "+String(p_gas[3])+"} (solubilities_pTX)");
  //  print("k={"+String(solu[1]/p_gas[1])+", "+String(solu[2]/p_gas[2])+", "+String(solu[3]/p_gas[3])+"}(solubilities_pTX)");
  //  print("solu={"+String(solu[1])+", "+String(solu[2])+", "+String(solu[3])+"}(solubilities_pTX)");
  //  print(Modelica.Math.Matrices.toString({MM_vec}));
  end solubilities_pTX;


  redeclare function extends density_liq_pTX
  //  PowerPlant.Media.Brine.Salt_Data_Duan.density_Duan2008_pTX;
protected
    parameter Integer[:] liqIndex=cat(1,1:nX_salt,{nX});
  algorithm
    d := Densities.density_Duan2008_pTX(p,T,X[liqIndex],MM[liqIndex]);
  //   print("density_liquid_pTX: "+String(p*1e-5)+" bar,"+String(T)+" K->"+String(d)+"kg/m^3");
  end density_liq_pTX;


 redeclare function extends specificEnthalpy_liq_pTX
 // Partial_Units.Molality molalities = massFractionsToMoleFractions(X, MM_vec);
 //  SI.SpecificEnthalpy h_H2O := Modelica.Media.Water.WaterIF97_pT.specificEnthalpy_pT(p, T) "H2O";
 algorithm
 //    h := SpecificEnthalpies.specificEnthalpy_pTX_Driesner(p,T,X);
     h := SpecificEnthalpies.specificEnthalpy_pTX_liq_Francke_cp(p,T,X);
 //  print(String(p*1e-5)+" bar,"+String(T)+" K->"+String(h)+" J/kg (Brine_Duan_Multi_TwoPhase_ngas_3.specificEnthalpy_liq_pTX)");
 end specificEnthalpy_liq_pTX;


 redeclare function extends specificEnthalpy_gas_pTX

 algorithm
     h :=BrineGas3Gas.specificEnthalpy_pTX(
         p,
         T,
         X);
 end specificEnthalpy_gas_pTX;


 redeclare function extends dynamicViscosity_liq
protected
   SI.Temperature T_corr;
 algorithm
  if state.T<273.16 then
     print("T="+String(state.T)+" too low (<0 degC), setting to 0 degC in BrineProp.Brine5salts3gas.dynamicViscosity_liq");
  end if;
  T_corr:= max(273.16,state.T);

  /*eta := Viscosities.dynamicViscosity_Duan_pTX(
    state.p,
    T_corr,
    state.X_l,
    MM_vec,
    Salt_data.saltConstants);*/
    eta := Viscosities.dynamicViscosity_DuanZhang_pTXd(
       state.p,
       T_corr,
       state.X,
       state.d,
       MM_vec,
       Salt_data.saltConstants);
       assert(eta>0,"Error in liquid viscosity calculation.");
 end dynamicViscosity_liq;


 redeclare function extends dynamicViscosity_gas
 algorithm
   eta  :=BrineGas3Gas.dynamicViscosity(BrineGas3Gas.ThermodynamicState(
         state.p,
         state.T,
         state.X_g));
   assert(eta>0,"Error in gas viscosity calculation.");
 end dynamicViscosity_gas;


  redeclare function extends saturationPressures
  algorithm

  //  if gasname =="carbondioxide" then
      p_sat[1] := if X[nX_salt+1]>0 then degassingPressure_CO2_Duan2006(p,T,X,MM_vec) else 0
    "aus Partial_Gas_Data TODO: use numeral";
  //  elseif gasname =="nitrogen" then
      p_sat[2] := if X[nX_salt+2]>0 then degassingPressure_N2_Duan2006(p,T,X,MM_vec) else 0
    "aus Partial_Gas_Data";
  //  elseif gasname =="methane" then
      p_sat[3] := if X[nX_salt+3]>0 then degassingPressure_CH4_Duan2006(p,T,X,MM_vec) else 0
    "aus Partial_Gas_Data";
  //  end if;
    if debugmode then
      print("saturationPressures("+String(p)+","+String(T)+")={"+Modelica.Math.Matrices.toString({p_sat})+"}");
    end if;
  end saturationPressures;


  redeclare function extends thermalConductivity
  "Thermal conductivity of water"
  algorithm
    lambda := Modelica.Media.Water.IF97_Utilities.thermalConductivity(
        state.d,
        state.T,
        state.p,
        state.phase);
  assert(lambda>0,"lambda = " + String(lambda) + "W/(m.K)");
  /*if lambda<0 then
    print("lambda = " + String(lambda) + "W/(m.K)");
  end if;*/
  end thermalConductivity;


  redeclare function extends surfaceTension
  algorithm
     sigma:=Modelica.Media.Water.WaterIF97_pT.surfaceTension(sat)
    "TODO http://www.if.ufrgs.br/~levin/Pdfs.dir/6756.pdf";
  end surfaceTension;


  redeclare function extends specificHeatCapacityCp_liq
  "calculation of liquid specific heat capacity from apparent molar heat capacities"
    extends BrineProp.SaltDataDuan.defineSaltOrder;

protected
    SI.MolarMass MM_vec_salt[:]=BrineProp.SaltData.MM_salt[1:5];
    SI.Pressure p=state.p;
    SI.Temperature T=state.T;
    Types.Molality b[size(X, 1)]=
        Utilities.massToMoleFractions(X,
        cat(1,
            MM_vec_salt,
            fill(-1, size(X, 1) - size(MM_vec_salt, 1))));

  /*  Real cp_by_cpWater[:]={0,
      SpecificEnthalpies.HeatCapacityRatio_KCl_White(T, b[KCl]),
      SpecificEnthalpies.HeatCapacityRatio_CaCl2_White(T, b[CaCl2]),
      0,0} "cp/cp_H2O of salt solutions";*/
    Types.PartialMolarHeatCapacity[5] Cp_appmol
    "Apparent molar enthalpy of salts";

    SI.SpecificHeatCapacity cp_Driesner
    "=SpecificEnthalpies.specificHeatCapacity_pTX_Driesner(p,T,X[1]/(X[1]+X[end]))";

    //  SI.SpecificHeatCapacity cp_H2O=Modelica.Media.Water.IF97_Utilities.cp_pT(p,T);

  //  SI.MassFraction X[:]=state.X "mass fraction m_NaCl/m_Sol";
  //    SI.MassFraction X[:]=cat(1,state.X[1:end-1],{1-sum(state.X[1:end-1])}) "Doesn't work in function in OM";
      SI.MassFraction X[size(state.X,1)] "OM workaround for cat";
  algorithm
      if debugmode then
        print("Running specificHeatCapacityCp_liq("+String(p/1e5)+" bar,"+String(T-273.15)+"degC, X="+Modelica.Math.Matrices.toString(transpose([state.X]))+")");
      end if;
      X[1:end-1]:=state.X[1:end-1] "OM workaround for cat";
      X[end]:=1-sum(state.X[1:end-1]) "OM workaround for cat";
  //    assert(state.X[end]>0, "No water in brine.");
      cp_Driesner:=SpecificEnthalpies.specificHeatCapacity_pTX_Driesner(p,T,X[1]/(X[1] + X[end]));

    Cp_appmol:={0,if b[KCl] > 0 then
      SpecificEnthalpies.appMolarHeatCapacity_KCl_White(T, b[KCl]) else 0,if b[
      CaCl2] > 0 then SpecificEnthalpies.appMolarHeatCapacity_CaCl2_White(T, b[
      CaCl2]) else 0,0,0} "Apparent molar enthalpy of salts";
  //    Cp_appmol:={(if b[i] > 0 and cp_by_cpWater[i]>0 then ((1 .+ MM_vec_salt[i] .* b[i]) .* cp_by_cpWater[i] .- 1)*cp_H2O ./ b[i] else 0) for i in 1:5};

      cp := (X[NaCl]+X[end])*cp_Driesner + X[end]*b[2:5]*(Cp_appmol[2:5]);

  //  cp:=(specificEnthalpy_pTX(state.p,state.T+.1,state.X)-state.h)/.1;
  //  cp := Modelica.Media.Water.IF97_Utilities.cp_pT(state.p, state.T)+mola[1:size(MM_vec_salt,1)];
  //  print("Cp_appmol: "+PowerPlant.vector2string(Cp_appmol)+" J/kg/K");
  //  print("cp_Driesner("+String(cp_Driesner)+")= J/(kg.K)");

      annotation (Documentation(info="<html>
                                <p>In the two phase region this function returns the interpolated heat capacity between the
                                liquid and vapour state heat capacities.</p>
                                </html>"));
  end specificHeatCapacityCp_liq;


  redeclare function extends specificHeatCapacityCp_gas
  "calculation of gas specific heat capacity"
  import SG = Modelica.Media.IdealGases.SingleGases;
  algorithm
    if state.x>0 then
  /*    cp_vec:= {SG.CO2.specificHeatCapacityCp(state),
             SG.N2.specificHeatCapacityCp(state),
             SG.CH4.specificHeatCapacityCp(state),
             SG.H2O.specificHeatCapacityCp(state)};
    cp:=X_g[end - nX_gas:end]*cp_vec;*/

  //    cp :=specificHeatCapacityCp_gas_TX(
      cp :=BrineGas3Gas.specificHeatCapacityCp_pTX(
            p=state.p,
            T=state.T,
            X=X_g[end - nX_gas:end]);
    else
      cp:=-1;
    end if;

      annotation (Documentation(info="<html>
                                <p>In the two phase region this function returns the interpolated heat capacity between the
                                liquid and vapour state heat capacities.</p>
                                </html>"));
  end specificHeatCapacityCp_gas;


  annotation (Documentation(info="<html>
<p><b>Brine_Duan_Multi_TwoPhase_ngas_3</b> is a medium package that, based on Brine.PartialBrine_ngas_Newton, defines a brine property model with five salts (NaCl, KCl, CaCl<sub>2</sub>, MgCl<sub>2</sub>, SrCl<sub>2</sub>) and 3 gases (CO<sub>2</sub>, N<sub>2</sub>, CH<sub>4</sub>), which are the main constituents of the geofluid in Gross Schoenebeck, Germany.</p>
<p>It was used for the calculations documented in this <a href=\"http://nbn-resolving.de/urn:nbn:de:kobv:83-opus4-47126\">PhD thesis</a>.</p>
<h4>Usage</h4>
<p>As it is based on Modelica.Media, the usage differs little from the usage of the two-phase water model:</p>
<p>Create an Instance of the Medium: </p>
<pre>  package Medium = Brine_Duan_Multi_TwoPhase_ngas_3;</pre>
<p>Create an Instance of Medium.Baseproperties: </p>
<pre>  Medium.BaseProperties props;</pre>
<p>Use the BaseProperties model to define the actual brine composition(Xi or X), to define the thermodynamic state and calculate the corresponding properties. </p>
<pre>  props.p = 1e5;
  props.T = 300;
  props.Xi = {0.08, 0.004, 0.12, 0.001, 0.002, 1-4, 7e-4, 6e-005} \"NaCl, KCl, CaCl2, MgCl2, SrCl2, CO2, N2, CH4\"
  d = props.d;
</pre>

<p>See <code><a href=\"Modelica://BrineProp.Examples.BrineProps2phase\">BrineProp.Examples.BrineProps2phase</a></code> for more usage examples.</p>

<p>All calculated values are returned in SI units and are mass based.</p>


<h4>Specific Enthalpy:</h4>
<p>The enthalpy is calculated as a mass fraction weighted average of the enthalpies of the two phases.</p>
<pre align=\"center\">h = x&middot;h_g + (1-x)&middot;h_l</pre>
<h5>Specific enthalpy of gas phase:</h5>
<p>Enthalpy of the gas phase is modelled as the enthalpy of an ideal mixture of ideal gases, i.e. it is calculated as the mass weighted average of the individual gas enthalpies including water. </p>
<p align=\"center\"><code>h_g = sum(h\"<sub>i</sub>&middot;X\"<sub>i</sub>)</code></p>
<p>The individual gas enthalpies are calculated using the functions for ideal gases in Modelica.Media.IdealGases.SingleGases.</p>
<h5>Specific enthalpy of liquid phase: </h5>
<p>Enthalpy of the liquid phase is assembled from the enthalpy of a NaCl-solution (Driesner) and the apparent molar enthalpies of the salts multiplied by their respective molalities</p>
<p align=\"center\"><code>h_l = h<sub>Driesner</sub>+sum(H<sub>i</sub><sup>app</sup>&middot;b<sub>i</sub>)</code></p>
<p>The apparent molar enthalpies of KCl and CaCl2 are calculated from apparent molar heat capacities which are calculated from a 2D fit for data from literature. The contributions of MgCl2 and SrCl2 are neglected, due to their small influence in the GrSbk fluid. </p>
<h4>Density</h4>
<p>The total density <code>d</code> of the fluid is calculated by combining the densities of both phases (<code>d_g </code>and <code>d_l</code>) according to their volume fractions. The gas phase is assumed to be an Density of the gas phase is assumed to be an ideal mixture of ideal gases. </p>
<pre align=\"center\">d = 1/(x/d_g + (1 - x)/d_l)</pre>
<h5>Density of the gas phase:</h5>
<p>Gas density is calculated using the ideal gas law.</p>
<h5>Density of liquid phase:</h5>
<p>Density of the liquid phase is calculated by combining the densities of solutions of single salts. The density model by Duan for single salt solutions is adapted for multi-salt solutions, resulting in an approach with apparent molar volumes analogous to the mixing rule for enthalpy. </p>
<h5>Created by</h5>
<div>Henning Francke<br/>
Helmholtz Centre Potsdam GFZ German Research Centre for Geosciences<br/>
Telegrafenberg, D-14473 Potsdam<br/>
Germany</div>
<p><a href=\"mailto:info@xrg-simulation.de\">francke@gfz-potsdam.de</a></p>
</html>",
 revisions="<html>

</html>"));
end Brine5salts3gas;
