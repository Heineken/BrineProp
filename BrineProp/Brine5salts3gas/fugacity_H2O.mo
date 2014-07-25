within BrineProp.Brine5salts3gas;
function fugacity_H2O
  "Calculation of fugacity coefficient according to (Duan 2003)"
//  extends fugacity_pTX;
  input SI.Pressure p;
  input SI.Temp_K T;
  output Real phi;

protected
  Types.Pressure_bar p_bar=SI.Conversions.to_bar(p);
  Real[:] a = {1.86357885E-03,
               1.17332094E-02,
               7.82682497E-07,
              -1.15662779E-05,
              -3.13619739,
              -1.29464029E-03};
algorithm
  phi := exp(a[1] + a[2]*p_bar + a[3]*p_bar^2 + a[4]*p_bar*T + a[5]*p_bar/T + a[6]*p_bar^2/T);
end fugacity_H2O;
