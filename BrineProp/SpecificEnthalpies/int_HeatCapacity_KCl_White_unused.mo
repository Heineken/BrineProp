within BrineProp.SpecificEnthalpies;
function int_HeatCapacity_KCl_White_unused
  "cp/cp_H2O integrated over Temperature"
//2D-fit Reproduction of measurements of heat capacity of KCl solution
  extends BrineProp.SpecificEnthalpies.PartialCpRatio_KCl_White;
  output Real int_cp_by_cpWater=(a + b*bn + c*Tn/2 + d*bn^2 + e*bn*Tn/2 + f*Tn^2/3 + g*bn^2*Tn/2 + h*bn*Tn^2/3 + i*Tn^3/4)*T_std*Tn;
//  SI.SpecificHeatCapacity cp_Water =  Modelica.Media.Water.IF97_Utilities.cp_pT(p, T);
/*protected 
  Real cp_by_cpWater;*/
algorithm
//cp_by_cpWater :=p00 + p10*T + p01*b + p20*T^2 + p11*T*b + p02*b^2;
//int_cp_by_cpWater:=(P00 + P01*b + P10/2*T + P20/3*T^2 + P11/2*T*b + P02*b^2)*T "cp ratio integrated over T";
//cp_by_cpWater :=a + b*bn + c*Tn + d*bn^2 + e*bn*Tn + f*Tn^2 + g*bn^2*Tn + h*bn*Tn^2 + i*Tn^3;
//print("KCl   bn: "+String(bn)+" Tn: "+String(Tn)+"");
//print("KCl   cp_by_cpWater("+String(mola)+" mol/kg): "+String(cp_by_cpWater)+"");
//int_cp_by_cpWater :=(a + b*bn + c*Tn/2 + d*bn^2 + e*bn*Tn/2 + f*Tn^2/3 + g*bn^2*Tn/2 + h*bn*Tn^2/3 + i*Tn^3/4)*T_std*Tn;

//cp = cp_by_cpWater*cp_Water;
end int_HeatCapacity_KCl_White_unused;
