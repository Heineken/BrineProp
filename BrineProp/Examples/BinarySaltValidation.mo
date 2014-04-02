within BrineProp.Examples;
model BinarySaltValidation "Calculation of Values given by Zhang1997"
package Medium = Brine_Duan_Multi_TwoPhase_ngas_3;

//NaCl
/*  Real b_t[:]={0.0805,0.1884,0.1959,0.262,0.4753,0.4912,0.7652,0.8386,1.2958,1.625,1.6953,1.7648,2.0102,3.6734,2.1492,2.2029,2.4938,2.7129,3.3438,3.5449,3.8629,4.0519,4.0913,4.2702,4.7142};
  Real f]=3.0218 "mole ratio of salts";
*/
  /*  
  Real b_t[:]={0.0505,0.1247,0.2137,0.3606,0.463,0.761,0.8437,1.1412,1.2591,1.8499,1.8933,2.664,2.0628,2.2922,2.3136,2.3985,2.6436,2.8055,3.3025,3.4498,3.6082,3.8786};
  Real f=1.0001 "mole ratio of salts";
   

  Real b_t[:]={0.0373,0.044 ,0.1157,0.152 ,0.2659,0.3694,0.401 ,0.551 ,0.6233,1.1907,1.4846,1.5343,1.646 ,3.2797,1.8152,2.0704,2.3294,2.407 ,2.4175,3.0112,3.4989,3.5791,3.6924,3.8449,3.9103,4.0365};
  Real f=1.0029 "mole ratio of salts";

  Real b_t[:]={0.0145,0.0362,0.0544,0.0783,0.1273,0.1514,0.2453,0.406,0.4845,0.5209,0.6562,0.9876,1.349,1.8159,2.2903,2.6288,2.8655,3.3164,3.5427,3.8457,4.046,4.1158,4.3216,4.4712};
  Real f=0.333400013 "mole ratio of salts";
*/

//KCl
/*
      Real b_t[:]={0.0043,0.0053,0.0131,0.0397,0.0463,0.0521,0.0887,0.1064,0.1656,0.1793,0.2799,0.2851,0.4813,2.2659,0.5846,0.7406,1.06,1.1123,1.7861,2.0733,3.2759,3.6716,3.9096,4.1326,4.3093};
  Real f=3.0138 "mole ratio of salts";

  Real b_t[:]={0.0075,0.056,0.2,0.3673,0.733,1.3793,2.28,2.9884,3.9055,4.0287,4.2295};
  Real f=3.0155 "mole ratio of salts";
   

  Real b_t[:]={0.013,0.0333,0.0534,0.0809,0.1048,0.1371,0.176,0.2251,0.2817,0.3029,0.4183,0.6443,0.9416,3.3067,1.077,1.2289,1.3517,1.5162,2.4113,3.0097,3.5053,3.7728,3.9678,4.1343,4.2724};
  Real f=1 "mole ratio of salts";

  Real b_t[:]={0.0328,0.051,0.0772,0.1255,0.2093,0.2295,0.4151,0.5273,0.5405,0.6981,0.8709,1.8217,0.9542,1.0324,1.0489,1.1445,1.4745,2.0021,2.377,2.5219,2.7118,2.8382};
  Real f=0.9996 "mole ratio of salts";

  Real b_t[:]={0.0071,0.0278,0.0452,0.2042,0.3003,0.4389,0.4991,0.7645,1.0255,1.3779,1.6425,2.9588,1.8372,1.9389,2.2741,2.3654,2.6535,3.9322,4.2837,4.5742,4.7571};
  Real f=0.33411293 "mole ratio of salts";
*/
//  Real b_t[:]={0.0177,0.0491,0.0778,0.15,0.1587,0.252,0.2905,1.9623,0.6285,1.2731,1.6893,2.6506,4.6633,4.7419};
  Real b_t[:]={0.0177};
  Real f=0.334896182 "mol ratio of salts";
/**/

  Medium.BaseProperties[size(b_t,1)] props;
  Real n_CaCl2[size(b_t,1)];
  SI.Density d[size(b_t,1)];
equation

//  for j in 1:size(f,1) loop
    for i in 1:size(b_t,1) loop
      props[i].p = 101325;
      props[i].T = 273.16+25;
      n_CaCl2[i] = b_t[i]/(1+f) "mol/l solvent";
      //NaCl
//      props[i].Xi = {f*n_CaCl2[i]*Medium.Salt_data.M_NaCl,0,n_CaCl2[i]*Medium.Salt_data.M_CaCl2,0,0,0,0,1e-3}/(1+f*n_CaCl2[i]*Medium.Salt_data.M_NaCl+n_CaCl2[i]*Medium.Salt_data.M_CaCl2);
      //KCl
      props[i].Xi = {0,f*n_CaCl2[i]*Medium.Salt_data.M_KCl,n_CaCl2[i]*Medium.Salt_data.M_CaCl2,0,0,0,0,1e-3}/(1+f*n_CaCl2[i]*Medium.Salt_data.M_KCl+n_CaCl2[i]*Medium.Salt_data.M_CaCl2);
      d[i] = props[i].d;
    end for;

algorithm
//  print("rho="+String(d)+" kg/m³, TDS = " + String(TDS) + " g/l -> "+ String(f*265/TDS));
end BinarySaltValidation;
