within BrineProp;
package Brine3salts "One-phase (liquid) multisalt brine solution"
  extends SaltDataDuan;// "for the molar masses below"
  extends BrineProp.PartialBrineMultiSaltOnePhase(
    redeclare package Salt_data = BrineProp.SaltDataDuan,
    saltNames = {"sodium chloride","potassium chloride","calcium chloride"},
    saltConstants = {
      saltConstants_NaCl,
      saltConstants_KCl,
      saltConstants_CaCl2},
    MM_salt = {M_NaCl,M_KCl,M_CaCl2},
    nM_salt = {nM_NaCl,nM_KCl,nM_CaCl2},
      final iNaCl=1,
      final iKCl=2,
      final iCaCl2=3);

    /*    redeclare package Salt_data = BrineProp.SaltDataDuan,
    final MM_salt = Salt_data.MM_salt,
    final nM_salt = Salt_data.nM_salt*/

  redeclare function extends density_pTX
  //  extends density_Duan2008_pTX(MM_vec=cat(1,MM_salt, {M_H2O}));
     //TODO should take MM_vec;
  //  extends Densities.density_Duan2008_pTX(MM_vec=MM_vec);
  algorithm
  //   print("density_liquid_pTX: "+String(p*1e-5)+" bar,"+String(T)+" K->"+String(d)+"kg/m^3");
    d := density_Duan2008_pTX(p,T,X,MM_vec,saltConstants);
  //  d := Brine_Driesner.density_pTX(p,T,X[1:nX_salt],MM_salt,saltConstants);
  //  d := Modelica.Media.Water.WaterIF97_pT.density_pT(p,T)  "*(1+sum(X[1:nX_salt]))/X[end]";
  end density_pTX;

 redeclare function specificEnthalpy_pTX
   extends specificEnthalpy_pTX_liq_Francke_cp(MM_vec=MM_salt);
 end specificEnthalpy_pTX;

 redeclare function extends dynamicViscosity_pTXd
  protected
   SI.Temperature T_corr;
 algorithm
 //print("p="+String(p)+" Pa, T="+String(T)+" K (BrineProp.Brine_5salts_noGas.dynamicViscosity_pTX)");

  if T<273.16 then
     print("T="+String(T)+" too low (<0 degC), setting to 0 degC in PartialBrine_ngas_Newton.quality_pTX()");
  end if;
     T_corr:= max(273.16,T);

    eta := dynamicViscosity_DuanZhang_pTXd(
       p,
       T_corr,
       X,
       d,
       MM_vec,
       saltConstants);
 end dynamicViscosity_pTXd;

  redeclare function extends specificHeatCapacityCp
    "calculation of liquid specific heat capacity from apparent molar heat capacities"
  algorithm
      cp:=specificHeatCapacityCp_pTX_liq_Francke(p=state.p,T=state.T,X=state.X,
          MM_vec=MM_salt);
  end specificHeatCapacityCp;

  redeclare function extends surfaceTension_T
  algorithm
     sigma:=Modelica.Media.Water.IF97_Utilities.surfaceTension(T)
      "TODO http://www.if.ufrgs.br/~levin/Pdfs.dir/6756.pdf";
  end surfaceTension_T;

  redeclare function extends thermalConductivity
    "Thermal conductivity of water TODO"
  algorithm
    lambda := Modelica.Media.Water.IF97_Utilities.thermalConductivity(
        Modelica.Media.Water.IF97_Utilities.rho_props_pT(state.p, state.T, Modelica.Media.Water.IF97_Utilities.waterBaseProp_pT(state.p, state.T, 1)),
        state.T,
        state.p,
        1);
  end thermalConductivity;

 redeclare function extends dynamicViscosity_pTX
  protected
   SI.Temperature T_corr;
 algorithm
 //print("p="+String(p)+" Pa, T="+String(T)+" K (BrineProp.Brine_5salts_noGas.dynamicViscosity_pTX)");

  if T<273.16 then
     print("T="+String(T)+" too low (<0 degC), setting to 0 degC in PartialBrine_ngas_Newton.quality_pTX()");
  end if;
     T_corr:= max(273.16,T);
    eta := dynamicViscosity_Duan_pTX(
       p,
       T_corr,
       X,
       MM_vec,
       saltConstants);
 end dynamicViscosity_pTX;

  annotation (Documentation(info="<html>
<p><b>BrineProp.Brine_5salts</b> is a medium package that provides properties of one-phase solution of five salts (NaCl, KCl, CaCl<sub>2</sub>, MgCl<sub>2</sub>, SrCl<sub>2</sub>).</p>
<h4>Usage</h4>
<p>It is based on Modelica.Media, the usage is accordingly:</p>
<p>Create an instance of the Medium (optionally deactivating range checks, for all options see .PartialFlags): </p>
<pre>  package Medium = Brine_5salts(AssertLevel=1,ignoreLimitSalt_T={false,true,true,false,false});</pre>
<p>Create an instance of Medium.Baseproperties: </p>
<pre>  Medium.BaseProperties props;</pre>
<p>Use the BaseProperties model to define the actual brine composition(Xi or X), to define the thermodynamic state and calculate the corresponding properties. </p>
<pre>  props.p = 1e5;
  props.T = 300;
  props.Xi = {0.08, 0.004, 0.12, 0.001, 0.002} \"NaCl, KCl, CaCl2, MgCl2, SrCl2\"
  d = props.d;</pre>
<p>Pressure and temperature as well as pressure and specific enthalpy can be used to define a thermodynamic state.</p>
<p>All calculated values are returned in SI-Units and are mass based. </p>
<p><h4>Details</h4></p>
<p>The model is explicit for p and T, but for h(p,T) the inverse function T(p,h) is defined. T(p,h) is inverts h(p,T) numerically by bisection, stopping at a given tolerance.</p>
<p>Density and enthalpy are calculated like the liquid phase properties in <code>BrineProp.Brine_5salts_TwoPhase_3gas</code> .</p>
</html>"));
end Brine3salts;
