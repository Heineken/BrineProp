within Brine.Brine_Duan_Multi_TwoPhase_ngas_3;
function fugacity_H2O
  "Calculation of fugacity coefficient according to (Duan 2003)"
  extends fugacity_pTX;

protected
  Pressure_bar P_1;
  Real[:] a = {1.86357885E-03,
               1.17332094E-02,
               7.82682497E-07,
              -1.15662779E-05,
              -3.13619739,
              -1.29464029E-03};
algorithm
  phi := exp(a[1] + a[2]*p_bar + a[3]*p_bar^2 + a[4]*p_bar*T_K + a[5]*p_bar/T_K + a[6]*p_bar^2/T_K);
end fugacity_H2O;
