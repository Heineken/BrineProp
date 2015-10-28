Attribute VB_Name = "GasData"
' Ideal gas coefficients from Modelica.Media.IdealGases.GasData
' Solubility functions
' TODO: put limits into the record
' TODO: improve limits management (change on function call)


' by Henning Francke francke@gfz-potsdam.de
' 2014 GFZ Potsdam

Option Explicit
Option Base 1
Public Const M_CO2 = 0.0440095 '[kg/mol]
Public Const M_N2 = 0.0280134 '[kg/mol]
Public Const M_CH4 = 0.01604246 '[kg/mol]
Public Const M_H2O = 0.018015 '[kg/mol]

Type DataRecord 'Coefficient data record for properties of ideal gases based on NASA source
  name As String 'Name of ideal gas
  MM As Double 'Molar mass
  Hf As Double 'Enthalpy of formation at 298.15K
  h0 As Double 'H0(298.15K) - H0(0K)
  Tlimit As Double 'Temperature limit between low and high data sets
  alow As Variant '(7) 'Low temperature coefficients a
  blow As Variant '(2) 'Low temperature constants b
  ahigh As Variant '(7) 'High temperature coefficients a
  bhigh As Variant '(2) 'High temperature constants b
  R As Double 'Gas constant
  nM As Integer 'number of ions per molecule
End Type


Function MM_CO2() As Double
    MM_CO2 = M_CO2
End Function
Function MM_N2() As Double
    MM_N2 = M_N2
End Function
Function MM_CH4() As Double
    MM_CH4 = M_CH4
End Function
Function MM_H2O() As Double
    MM_H2O = M_H2O
End Function

Function CO2() As DataRecord
    With CO2
        .name = "CO2"
        .MM = M_CO2
        .Hf = -8941478.54440518
        .h0 = 212805.621513537
        .Tlimit = 1000
        .alow = Array(49436.5054, -626.411601, 5.30172524, 0.002503813816, -2.127308728E-07, -7.68998878E-10, 2.849677801E-13)
        .blow = Array(-45281.9846, -7.04827944)
        .ahigh = Array(117696.2419, -1788.791477, 8.29152319, -0.0000922315678, 4.86367688E-09, -1.891053312E-12, 6.33003659E-16)
        .bhigh = Array(-39083.5059, -26.52669281)
        .R = 188.924482214067
        .nM = 1
    End With
End Function

Function N2() As DataRecord
    With N2
        .name = "N2"
        .MM = M_N2
        .Hf = 0
        .h0 = 309498.454311151
        .Tlimit = 1000
        .alow = Array(22103.71497, -381.846182, 6.08273836, -0.00853091441, 0.00001384646189, -9.62579362E-09, 2.519705809E-12)
        .blow = Array(710.846086, -10.76003744)
        .ahigh = Array(587712.406, -2239.249073, 6.06694922, -0.00061396855, 1.491806679E-07, -1.923105485E-11, 1.061954386E-15)
        .bhigh = Array(12832.10415, -15.86640027)
        .R = 296.803386950531
        .nM = 1
    End With
End Function
    
Function CH4() As DataRecord
    With CH4
        .name = "CH4"
        .MM = M_N2
        .Hf = -4650159.63885838
        .h0 = 624355.740952447
        .Tlimit = 1000
        .alow = Array(-176685.0998, 2786.18102, -12.0257785, 0.0391761929, -0.0000361905443, 2.026853043E-08, -4.97670549E-12)
        .blow = Array(-23313.1436, 89.0432275)
        .ahigh = Array(3730042.76, -13835.01485, 20.49107091, -0.001961974759, 0.000000472731304, -3.72881469E-11, 1.623737207E-15)
        .bhigh = Array(75320.6691, -121.9124889)
        .R = 518.279116793809
        .nM = 1
    End With
End Function

Function H2O() As DataRecord
    With H2O
        .name = "H2O"
        .MM = M_H2O
        .Hf = -13423382.8172529
        .h0 = 549760.647628014
        .Tlimit = 1000
        .alow = Array(-39479.6083, 575.573102, 0.931782653, 0.00722271286, -0.00000734255737, 4.95504349E-09, -1.336933246E-12)
        .blow = Array(-33039.7431, 17.24205775)
        .ahigh = Array(1034972.096, -2412.698562, 4.64611078, 0.002291998307, -0.000000683683048, 9.42646893E-11, -4.82238053E-15)
        .bhigh = Array(-13842.86509, -7.97814851)
        .R = 461.523329085088
        .nM = 1
    End With
End Function

Function solubility_CO2_pTX_Duan2006(p As Double, T As Double, X, p_gas) 'CO2 solubility in aqueous saltsolutions
'' Zhenhao Duan et al. (2006) An improved model for the calculation of CO2 solubility in aqueous
'' solutions containing Na+,K+,Ca2+,Mg2+,Cl-, and SO4_2-. Marine Chemistry 98131-139.
'' fugacity from doi10.1016/j.marchem.2005.09.001

 Dim solu As Double 'CO2 solubility in mol_CO2/kg H2O
 Dim mu_l0_CO2_RT_c, lambda_CO2_Na_c, zeta_CO2_NaCl_c
 mu_l0_CO2_RT_c = Array(28.9447706, -0.0354581768, -4770.67077, 0.0000102782768, 33.8126098, 0.0090403714, -0.00114934031, -0.307405726, -0.0907301486, 0.000932713393, 0)

 lambda_CO2_Na_c = Array(-0.411370585, 0.000607632013, 97.5347708, 0, 0, 0, 0, -0.0237622469, 0.0170656236, 0, 0.0000141335834)

 zeta_CO2_NaCl_c = Array(0.000336389723, -0.000019829898, 0, 0, 0, 0, 0, 0.0021222083, -0.00524873303, 0, 0)

 Dim p_H2O As Double
 p_H2O = IAPWS.Waterpsat_T(T)
 Dim phi  As Double
 Dim mu_l0_CO2_RT  As Double
 Dim lambda_CO2_Na  As Double
 Dim zeta_CO2_NaCl As Double

 'constant
 Dim molalities '() As Double ReDim molalities(nX)
 molalities = massFractionsToMolalities(X, Brine.MM_vec)
 Dim m_Cl As Double, m_Na As Double, m_K As Double, m_Ca As Double, m_Mg As Double, m_SO4 As Double
 m_Cl = molalities(i_NaCl) + molalities(i_KCl) + 2 * molalities(i_CaCl2)
 m_Na = molalities(i_NaCl)
 m_K = molalities(i_KCl)
 m_Ca = molalities(i_CaCl2)
 
 If Not p_gas > 0 Then
    solubility_CO2_pTX_Duan2006 = 0
 Else
    Dim msg As String
    If T < 273 Or T > 533 Then
        msg = "T=" & T - 273.15 & "°C, CO2 solubility only valid for 0<T<260°C (GasData.solubility_CO2_pTX_Duan2003)"
    End If
    If (p < 0 Or p > 2000 * 10 ^ 5) Then
         msg = "p=" & p / 10 ^ 5 & " bar, CO2 fugacity only valid for 0<p<2000 bar (GasData.solubility_CO2_pTX_Duan2003)"
    End If
    If Len(msg) > 0 Then
        If outOfRangeMode = 1 Then
            Debug.Print msg
        ElseIf outOfRangeMode = 2 Then
            solubility_CO2_pTX_Duan2006 = msg
            Exit Function
        End If
    End If

     'equ. 9
     phi = fugacity_CO2_Duan2006(p_gas + p_H2O, T)
     mu_l0_CO2_RT = Par_CO2_Duan2003(p_gas + p_H2O, T, mu_l0_CO2_RT_c)
     lambda_CO2_Na = Par_CO2_Duan2003(p_gas + p_H2O, T, lambda_CO2_Na_c)
     zeta_CO2_NaCl = Par_CO2_Duan2003(p_gas + p_H2O, T, zeta_CO2_NaCl_c)
    
     solu = phi * p_gas / 10 ^ 5 * Exp(-mu_l0_CO2_RT - 2 * lambda_CO2_Na * (m_Na + m_K + 2 * m_Ca + 2 * m_Mg) - zeta_CO2_NaCl * m_Cl * (m_Na + m_K + m_Mg + m_Ca) + 0.07 * m_SO4 * 0)
     solubility_CO2_pTX_Duan2006 = solu * M_CO2 * X(Brine.nX) 'molality->mass fraction
 End If
End Function



Function fugacity_CO2_Duan2006(p As Double, T As Double) 'Calculation of fugacity coefficient according to (Duan 2006)
  'doi:10.1016/j.marchem.2005.09.001

  Dim p_bar As Double, P_1 As Double
  Dim c
  p_bar = p / 10 ^ 5
  
  If outOfRangeMode = 1 Then
    If T < 273 Or T > 573 Then
      Debug.Print "T=" & T & "K, but CO2 solubility calculation is only valid for temperatures between 0 and 260°C (GasData.fugacity_CO2_Duan2006)"
    End If
   If (p < 0 Or p > 2000 * 10 ^ 5) Then
       Debug.Print "p=" & p / 10 ^ 5 & " bar, but CO2 fugacity calculation only valid for pressures between 0 and 2000 bar (Partial_Gas_Data.fugacity_CO2_Duan2006)"
   End If
  ElseIf outOfRangeMode = 2 Then
    If Not (T > 273 And T < 573) Then
        fugacity_CO2_Duan2006 = "#T=" & T - 273.15 & "°C out of range(0...300°C) for CO2 fugacity calculation (fugacity_CO2_Duan2006)"
    End If
    If Not p < 2000 * 10 ^ 5 Then
        fugacity_CO2_Duan2006 = "#p=" & p / 10 ^ 5 & " bar out of range for CO2 fugacity calculation (fugacity_CO2_Duan2006)"
      End If
  End If

    If T < 305 Then
      P_1 = p_sat_CO2(T) / 10 ^ 5
    ElseIf T < 405 Then
      P_1 = 75 + (T - 305) * 1.25
    Else
      P_1 = 200
    End If

  If p_bar < P_1 Then
    '1 273<T<573 and p_bar<P_1
      c = Array(1#, 0.0047586835, -0.0000033569963, 0, -1.3179396, -0.0000038389101, 0, 0.0022815104, 0, 0, 0, 0, 0, 0, 0)
  ElseIf T > 435 Then
          '6 T>435 and p_bar>P_1
          c = Array(-0.1569349, 0.00044621407, -0.00000091080591, 0, 0, 0.00000010647399, 2.4273357E-10, 0, 0.35874255, 0.00006331971, -249.89661, 0, 0, 888.768, -0.00000066348003)

  ElseIf p_bar < 1000 Then
    'P_1<p_bar<1000 and 273<T<435
    If T < 340 Then
     '2 273<T<340 and P_1<p_bar<1000
           c = Array(-0.71734882, 0.00015985379, -0.00000049286471, 0, 0, -0.00000027855285, 1.1877015E-09, 0, 0, 0, 0, -96.539512, 0.44774938, 101.81078, 0.0000053783879)
    Else
      'T>340
          '4 340<T<435 and P_1<p_bar<1000
          c = Array(5.0383896, -0.0044257744, 0, 1.9572733, 0, 0.0000024223436, 0, -0.00093796135, -1.502603, 0.003027224, -31.377342, -12.847063, 0, 0, -0.000015056648)
    End If

      Else
    'p_bar>1000 bar and 273<T<435
        If T < 340 Then
          '3 273<T<340 and p_bar>1000
          c = Array(-0.065129019, -0.00021429977, -0.000001144493, 0#, 0#, -0.00000011558081, 0.000000001195237, 0#, 0#, 0#, 0#, -221.34306, 0#, 71.820393, 0.0000066089246)
        Else 'T>340
          '5 340<T<435 and p_bar>1000
          c = Array(-16.063152, -0.002705799, 0, 0.14119239, 0, 0.00000081132965, 0, -0.00011453082, 2.3895671, 0.00050527457, -17.76346, 985.92232, 0, 0, -0.00000054965256)
        End If
        'Region 6 omitted
  End If
    c = ToDouble(c)
  fugacity_CO2_Duan2006 = c(1) + (c(2) + c(3) * T + c(4) / T + c(5) / (T - 150)) * p_bar + (c(6) + c(7) * T + c(8) / T) * p_bar ^ 2 + (c(9) + c(10) * T + c(11) / T) * Log(p_bar) + (c(12) + c(13) * T) / p_bar + c(14) / T + c(15) * T ^ 2
End Function

  Function Par_CO2_Duan2003(p As Double, T As Double, c) 'Duan,Sun(2003)
    Dim p_bar As Double
    p_bar = p / 10 ^ 5
   
    'eq. 7
    c = ToDouble(c)
    Par_CO2_Duan2003 = c(1) + c(2) * T + c(3) / T + c(4) * T ^ 2 + c(5) / (630 - T) + c(6) * p_bar + c(7) * p_bar * Log(T) + c(8) * p_bar / T + c(9) * p_bar / (630 - T) + c(10) * p_bar ^ 2 / (630 - T) ^ 2 + c(11) * T * Log(p_bar)
  End Function

  Function p_sat_CO2(T As Double) 'calculates saturation pressure, polynom derived from EES calculations
    If Not T < 305 Then
        p_sat_CO2 = "#Temperature above critical Temperature (304.1 K)"
    Else
        p_sat_CO2 = 1178.4 * T ^ 2 - 555378 * T + 7 * 10 ^ 7
    End If
  End Function

 Function solubility_N2_pTX_Duan2006(p As Double, T As Double, X, p_gas) 'solubility calculation of N2 in seawater Mao&Duan(2006)
    ' Shide Mao and Zhenhao Duan (2006) A thermodynamic model for calculating nitrogen solubility, gas phase composition and density of the H2O-N2-NaCl system. Fluid Phase Equilibria, 248 (2): 103-114
    ' 273-400 K, 1-600 bar and 0-6 mol/kg
    ' http://dx.doi.org/10.1016/j.fluid.2006.07.020
    ' http://www.geochem-model.org/wp-content/uploads/2009/09/46-FPE_248_103.pdf
    
    'Dim M_H2O As Double:M_H2O = H2O.MM
    Dim molalities
    molalities = ToDouble(massFractionsToMolalities(X, Brine.MM_vec))
    If VarType(molalities) = vbString Then
        solubility_N2_pTX_Duan2006 = molalities
        Exit Function
    End If
    Dim m_Cl As Double, m_Na As Double, m_K As Double, m_Ca As Double, m_Mg As Double, m_SO4 As Double
    m_Cl = molalities(i_NaCl) + molalities(i_KCl) + 2 * molalities(i_CaCl2) ' + 2 * molalities(i_MgCl2)
    m_Na = molalities(i_NaCl)
    m_K = molalities(i_KCl)
    m_Ca = molalities(i_CaCl2)
    m_Mg = 0 ' molalities(i_MgCl2)
    m_SO4 = 0 ' molalities(i_MgCl2)

    Dim p_H2O As Double, x_NaCl As Double
    p_H2O = IAPWS.Waterpsat_T(T)
    x_NaCl = molalities(i_NaCl) * M_H2O 'mole fraction of NaCl in liquid phase
    Dim v_l_H2O As Double 'MolarVolume
    v_l_H2O = M_H2O / IAPWS.Density_pT(p, T)
    Dim phi_H2O As Double
    phi_H2O = fugacity_H2O_Duan2006N2(p, T)
    Const R = 83.14472 'bar.cm3/(mol.K) Molar gas constant"
    Dim phi_N2
    Dim mu_l0_N2_RT As Double
    Dim lambda_N2_Na As Double
    Dim xi_N2_NaCl As Double
   
    If Not p_gas > 0 Then
      solubility_N2_pTX_Duan2006 = 0
    Else
        Dim msg As String
     If outOfRangeMode > 0 Then
       If Not ignoreLimitN2_T And (273 > T Or T > 400) Then
          msg = "#T=" & (T - 273.15) & " °C, N2 solubility only valid for 0<T<127°C (GasData.solubility_N2_pTX_Duan2006)"
       End If
       If (p < 10 ^ 5 Or p > 600 * 10 ^ 5) Then
          msg = "#p=" & (p / 10 ^ 5) & " bar, N2 solubility only valid for 1<p<600 bar (GasData.solubility_N2_pTX_Duan2006)"
       End If
       If molalities(i_NaCl) > 6 Then
         msg = "#mola(i_NaCl)=" & (molalities(i_NaCl)) & " mol/kg, but N2 solubility only valid up to 6 mol/kg (GasData.solubility_N2_pTX_Duan2006)"
       End If
         If Len(msg) > 0 Then
            If outOfRangeMode = 1 Then
                Debug.Print msg
            ElseIf outOfRangeMode = 2 Then
                solubility_N2_pTX_Duan2006 = msg
                Exit Function
            End If
        End If
     End If

    Dim mu_l0_N2_RT_c, lambda_N2_Na_c, xi_N2_NaCl_c
    mu_l0_N2_RT_c = Array(-23.093813, 0.056048525, 9880.8898, -0.000051091621, -1322029.8, -0.00049542866, 0.0000012698747, 0.51411144, -0.000064733978)
    lambda_N2_Na_c = Array(-2.4434074, 0.0036351795, 447.47364, 0, 0, -0.000013711527, 0, 0, 0.0000071037217)
    xi_N2_NaCl_c = Array(-0.0058071053, 0, 0, 0, 0, 0, 0, 0, 0)
      
      phi_N2 = fugacity_N2_Duan2006(p_gas + p_H2O, T)
      If VarType(phi_N2) = vbString Then
        solubility_N2_pTX_Duan2006 = phi_N2
        Exit Function
      End If
      
      mu_l0_N2_RT = Par_N2_Duan2006(p_gas + p_H2O, T, mu_l0_N2_RT_c)
      If VarType(mu_l0_N2_RT) = vbString Then
        solubility_N2_pTX_Duan2006 = mu_l0_N2_RT
        Exit Function
      End If
      
      lambda_N2_Na = Par_N2_Duan2006(p_gas + p_H2O, T, lambda_N2_Na_c)
      If VarType(lambda_N2_Na) = vbString Then
        solubility_N2_pTX_Duan2006 = lambda_N2_Na
        Exit Function
      End If
      
      xi_N2_NaCl = Par_N2_Duan2006(p_gas + p_H2O, T, xi_N2_NaCl_c)
      If VarType(xi_N2_NaCl) = vbString Then
        solubility_N2_pTX_Duan2006 = xi_N2_NaCl
        Exit Function
      End If

    'equ. 9
      Dim solu As Double
      solu = p_gas / 10 ^ 5 * phi_N2 * Exp(-mu_l0_N2_RT - 2 * lambda_N2_Na * (m_Na + m_K + 2 * m_Ca + 2 * m_Mg) - xi_N2_NaCl * (m_Cl + 2 * m_SO4) * (m_Na + m_K + 2 * m_Ca + 2 * m_Mg) - 4 * 0.0371 * m_SO4)
      solubility_N2_pTX_Duan2006 = solu * M_N2 * X(Brine.nX) 'molality->mass fraction
    End If
End Function

Function fugacity_N2_Duan2006(p As Double, T As Double)  'Zero search with EOS from Duan2006
'doi:10.1016/j.?uid.2006.07.020
'Shide Mao, Zhenhao Duan:A thermodynamic model for calculating nitrogen solubility, gas phase composition and density of the N2?H2O?NaCl system
    Dim p_bar As Double
    p_bar = p / 10 ^ 5
    
    Dim V_neu As Double, V As Double  'SpecificVolume
    V_neu = 0.024 'start value SpecificVolume
    Dim a
    a = Array(0.0375504388, -10873.0273, 1109648.61, 0.000541589372, 112.094559, -5921.91393, 0.00000437200027, 0.495790731, -164.902948, -7.07442825E-08, 0.00965727297, 0.487945175, 16225.7402, 0.00899)
    Dim sigma As Double
    sigma = 3.63
    Dim epsilon As Double 'Temperature
    epsilon = 101
    Dim T_m As Double, b As Double, c As Double, d As Double, E As Double, f As Double
    T_m = 154 * T / epsilon
    b = a(1) + a(2) / T_m ^ 2 + a(3) / T_m ^ 3
    c = a(4) + a(5) / T_m ^ 2 + a(6) / T_m ^ 3
    d = a(7) + a(8) / T_m ^ 2 + a(9) / T_m ^ 3
    E = a(10) + a(11) / T_m ^ 2 + a(12) / T_m ^ 3
    f = a(13) / T_m ^ 3
    
    Dim P_m As Double, G As Double, ln_phi As Double, V_m As Double, Z As Double
    P_m = 3.0626 * sigma ^ 3 * p_bar / epsilon
    Dim z_ As Integer
    Dim d_ As Double 'only a counter to avoid getting caught in the iteration loop
    d_ = 0.7 'dampening factor 0=no dampening, 1=no progress
    
    'iterative solution
    While Abs(V - V_neu) > 10 ^ -8
        V = IIf(z_ < 5, V_neu, (1 - d_) * V_neu + d_ * V)
        V_m = V * 10 ^ 6 / (1000# * (sigma / 3.691) ^ 3)
        Z = 1 + b / V_m + c / V_m ^ 2 + d / V_m ^ 4 + E / V_m ^ 5 + f / V_m ^ 2 * (1 + a(14) / V_m ^ 2) * Exp(-a(14) / V_m ^ 2)
        V_neu = Z / p * Constants.R * T
        '    print("V("&(z)&")="&(V_neu))
        z_ = z_ + 1
        
        If z_ >= 1000 Then
            fugacity_N2_Duan2006 = "#Reached maximum number of iterations for fugacity calculation.(fugacity_N2_Duan2006)"
            Exit Function
        End If
    Wend
    
    V_m = 1000# * V / (sigma / 3.691) ^ 3 'm³/mol -> dm³/mol"
    Z = 1 + (a(1) + a(2) / T_m ^ 2 + a(3) / T_m ^ 3) / V_m _
        + (a(4) + a(5) / T_m ^ 2 + a(6) / T_m ^ 3) / V_m ^ 2 _
        + (a(7) + a(8) / T_m ^ 2 + a(9) / T_m ^ 3) / V_m ^ 4 _
        + (a(10) + a(11) / T_m ^ 2 + a(12) / T_m ^ 3) / V_m ^ 5 _
        + a(13) / T_m ^ 3 / V_m ^ 2 * (1 + a(14) / V_m ^ 2) * Exp(-a(14) / V_m ^ 2)
    G = a(13) / T_m ^ 3 / (2 * a(14)) * (2 - (2 + a(14) / V_m ^ 2) * Exp(-a(14) / V_m ^ 2))
    fugacity_N2_Duan2006 = Exp(Z - 1 + b / V_m + c / (2 * V_m ^ 2) + d / (4 * V_m ^ 4) + E / (5 * V_m ^ 5) + G) / Z 'fugacity coefficient
End Function
  
Function Par_N2_Duan2006(p As Double, T As Double, c) 'Mao,Duan(2006)
    Dim p_bar As Double
    p_bar = p / 10 ^ 5
    'eq. 7
    Par_N2_Duan2006 = c(1) + c(2) * T + c(3) / T + c(4) * T ^ 2 + c(5) / T ^ 2 + c(6) * p_bar + c(7) * p_bar * T + c(8) * p_bar / T + c(9) * p_bar ^ 2 / T
End Function

 Function fugacity_H2O_Duan2006N2(p As Double, T As Double) 'Calculation of fugacity coefficient
' according to (Mao&Duan 2006 'A thermodynamic model for calculating nitrogen solubility, gasphase composition and density of the N2?H2O?NaCl system')"
    Dim p_bar As Double, a, phi As Double
    p_bar = p / 10 ^ 5
    a = Array(0.00186357885, 0.0117332094, 0.000000782682497, -0.0000115662779, -3.13619739, -0.00129464029)
    phi = Exp(a(1) + a(2) * p_bar + a(3) * p_bar ^ 2 + a(4) * p_bar * T + a(5) * p_bar / T + a(6) * p_bar ^ 2 / T) 'equ. 5
 End Function


Function solubility_CH4_pTX_Duan2006(p As Double, T As Double, X, p_gas) 'Duan ZH, Mao SD. (2006) A thermodynamic model for calculating methane solubility, density and gas phase composition of methane-bearing aqueous fluids from 273 to 523 K and from 1 to 2000 bar. Geochimica et Cosmochimica Acta, 70 (13): 3369-3386.
' http://geochem-model.org/Publications/43-GCA_2006_3369.pdf
' http://dx.doi.org/10.1016/j.gca.2006.03.018TODO Umrechnung andere Salz in NaCl"
'  output SI.MassFraction c_gas "gas concentration in kg_gas/kg_H2O"

    Dim mu_l0_CH4_RT_c, lambda_CH4_Na_c, xi_CH4_NaCl_c
    mu_l0_CH4_RT_c = Array(8.3143711, -0.00072772168, 2148.9858, -0.000014019672, -667434.49, 0.007698589, -0.0000050253331, -3.0092013, 484.68502, 0)
    lambda_CH4_Na_c = Array(-0.81222036, 0.0010635172, 188.94036, 0, 0, 0.000044105635, 0, 0, 0, -4.6797718E-11)
    xi_CH4_NaCl_c = Array(-0.0029903571, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    
    'Dim M_H2O As Double: M_H2O = H2O.MM
    Dim p_H2O As Double: p_H2O = p_sat_H2O_Duan2003(T)
    Dim v_l_H2O As Double: v_l_H2O = M_H2O / IAPWS.Density_pT(p, T) 'MolarVolume
    Dim phi_H2O As Double: phi_H2O = fugacity_H2O_Duan2006CH4(p, T)
    Dim phi_CH4 As Double, mu_l0_CH4_RT As Double, lambda_CH4_Na As Double, xi_CH4_NaCl As Double
    
    If Not p_gas > 0 Then
        solubility_CH4_pTX_Duan2006 = 0
    Else
        Dim msg As String
        If outOfRangeMode > 0 Then
            If 273 > T Or T > 273 + 250 Then
                msg = "#T=" & (T) & " K, CH4 solubility only valid for 0<T<250°C (GasData.solubility_CH4_pTX_Duan2006)"
            End If
            If p < 10 ^ 5 Or p > 2000 * 10 ^ 5 Then
                msg = "#p=" & (p / 10 ^ 5) & " bar, but CH4 fugacity only valid for 1<p<1600 bar (GasData.solubility_CH4_pTX_Duan2006)"
            End If
        
            If Len(msg) > 0 Then
                If outOfRangeMode = 1 Then
                   Debug.Print msg
                ElseIf outOfRangeMode = 2 Then
                    solubility_CH4_pTX_Duan2006 = msg
                    Exit Function
                End If
            End If
        End If

        Dim molalities '(nX)
        molalities = massFractionsToMolalities(X, Brine.MM_vec)
        Dim m_Cl As Double, m_Na As Double, m_K As Double, m_Ca As Double, m_Mg As Double, m_SO4 As Double                 'Molality
        m_Cl = molalities(i_NaCl) + molalities(i_KCl) + 2 * molalities(i_CaCl2)  '+ 2 * molalities(i_MgCl2)
        m_Na = molalities(i_NaCl)
        m_K = molalities(i_KCl)
        m_Ca = molalities(i_CaCl2)
        m_Mg = 0 ' molalities(i_MgCl2)
        m_SO4 = 0 ' molalities(i_MgCl2)
        
        phi_CH4 = fugacity_CH4_Duan1992(p_gas + p_H2O, T)
        mu_l0_CH4_RT = Par_CH4_Duan2006(p_gas + p_H2O, T, mu_l0_CH4_RT_c)
        lambda_CH4_Na = Par_CH4_Duan2006(p_gas + p_H2O, T, lambda_CH4_Na_c)
        xi_CH4_NaCl = Par_CH4_Duan2006(p_gas + p_H2O, T, xi_CH4_NaCl_c)
        
        'equ. 10
        Dim solu As Double
        solu = p_gas / 10 ^ 5 * phi_CH4 * Exp(-mu_l0_CH4_RT - 2 * lambda_CH4_Na * (m_Na + m_K + 2 * m_Ca + 2 * m_Mg) - xi_CH4_NaCl * (m_Na + m_K + 2 * m_Ca + 2 * m_Mg) * (m_Cl + 2 * m_SO4) - 4 * 0.0332 * m_SO4)
        
        solubility_CH4_pTX_Duan2006 = solu * M_CH4 * X(Brine.nX) 'molality->mass fraction
    End If
End Function

  
Function fugacity_CH4_Duan1992(p As Double, T As Double)    'Zero search with EOS from Duan1992
    Dim V_neu As Double, V As Double  'SpecificVolume
    V_neu = 0.024 'start value SpecificVolume
    Dim a
    a = Array(0.0872553928, -0.752599476, 0.375419887, 0.0107291342, 0.0054962636, -0.0184772802, 0.000318993183, 0.000211079375, 0.0000201682801, -0.0000165606189, 0.000119614546, -0.000108087289, 0.0448262295, 0.75397, 0.077167)
    Dim alpha As Double, beta As Double, gamma As Double, T_C As Double
    alpha = a(13)
    beta = a(14)
    gamma = a(15)
    T_C = 190.6
    Dim P_c As Double, P_r As Double, T_r As Double
    P_c = 46.41 * 10 ^ 5
    P_r = p / P_c
    T_r = T / T_C
    Dim b As Double, c As Double, d As Double, E As Double, f As Double
    b = a(1) + a(2) / T_r ^ 2 + a(3) / T_r ^ 3
    c = a(4) + a(5) / T_r ^ 2 + a(6) / T_r ^ 3
    d = a(7) + a(8) / T_r ^ 2 + a(9) / T_r ^ 3
    E = a(10) + a(11) / T_r ^ 2 + a(12) / T_r ^ 3
    f = alpha / T_r ^ 3
    Dim ln_phi As Double, d_ As Double, G As Double, Z As Double, V_r As Double
    Dim z_ As Integer 'counter to avoid getting caught in the iteration loop
    d_ = 0.7 'dampening factor 0=no dampening, 1=no progress
    
    While Abs(V - V_neu) > 10 ^ -8
        V = IIf(z_ < 5, V_neu, (1 - d_) * V_neu + d_ * V) 'dampened
        V_r = V / (Constants.R * T_C / P_c)
        G = f / (2 * gamma) * (beta + 1 - (beta + 1 + gamma / V_r ^ 2) * Exp(-gamma / V_r ^ 2))
        Z = 1 + b / V_r + c / V_r ^ 2 + d / V_r ^ 4 + E / V_r ^ 5 + f / V_r ^ 2 * (beta + gamma / V_r ^ 2) * Exp(-gamma / V_r ^ 2)
        V_neu = Z / p * Constants.R * T
        z_ = z_ + 1
        If z_ >= 1000 Then
            fugacity_CH4_Duan1992 = "#Reached maximum number of iterations for fugacity calculation.(fugacity_CH4_Duan1992)"
            Exit Function
        End If
    Wend
    
    fugacity_CH4_Duan1992 = Exp(Z - 1 + b / V_r + c / (2 * V_r ^ 2) + d / (4 * V_r ^ 4) + E / (5 * V_r ^ 5) + G) / Z 'fugacity coefficient
End Function
  
Function Par_CH4_Duan2006(p As Double, T As Double, c) 'Duan,Sun(2003)
    Dim p_bar
    p_bar = p / 10 ^ 5
    'eq. 7
    Par_CH4_Duan2006 = c(1) + c(2) * T + c(3) / T + c(4) * T ^ 2 + c(5) / T ^ 2 + c(6) * p_bar + c(7) * p_bar * T + c(8) * p_bar / T + c(9) * p_bar / T ^ 2 + c(10) * p_bar ^ 2 * T
End Function

Function p_sat_H2O_Duan2003(T) 'calculates saturation pressure of water, with the equation given in Duan2003
'An improved model calculating CO2 solubility in pure water and aqueous NaCl solutions from 273 to 533 K and from 0 to 2000 bar
    Dim c, P_c As Double, T_C As Double, t_ As Double
    c = Array(-38.640844, 5.894842, 59.876516, 26.654627, 10.637097)
    P_c = 220.85 * 10 ^ 5
    T_C = 647.29
    t_ = (T - T_C) / T_C
    
    p_sat_H2O_Duan2003 = (P_c * T / T_C) * (1 + c(1) * (-t_) ^ 1.9 + c(2) * t_ + c(3) * t_ ^ 2 + c(4) * t_ ^ 3 + c(5) * t_ ^ 4)
End Function

Function fugacity_H2O_Duan2006CH4(p As Double, T As Double) 'Calculation of fugacity coefficient according to Duan ZH, Mao SD. (2006)
'A thermodynamic model for calculating methane solubility, density and gas phase composition of methane-bearing aqueous fluids from 273 to 523 K and from 1 to 2000 bar. Geochimica et Cosmochimica Acta, 70 (13): 3369-3386.
    Dim p_bar As Double, a, phi As Double
    p_bar = p / 10 ^ 5
    a = Array(-0.0142006707, 0.010836991, -0.0000015921316, -0.0000110804676, -3.14287155, 0.00106338095)
    'equ. 6
    phi = Exp(a(1) + a(2) * p_bar + a(3) * p_bar ^ 2 + a(4) * p_bar * T + a(5) * p_bar / T + a(6) * p_bar ^ 2 / T)
End Function

Function MoistAirDynamicViscosity(T)
    'From Modelica.Media.Air.MoistAir.dynamicViscosity()
    Dim msg As String
    If T < 240 Or T > 400 Then
        msg = "T out of range (240..400 K) for MoistAir viscosity calculation (GasData.MoistAirDynamicViscosity)"
        If outOfRangeMode = 1 Then
            Debug.Print msg
        ElseIf outOfRangeMode = 2 Then
            MoistAirDynamicViscosity = msg
            Exit Function
        End If
    End If

    Dim T_C As Double: T_C = T - 273.15
    MoistAirDynamicViscosity = -4.96717436974791 * 10 ^ -11 * T_C ^ 2 + 5.06626785714286 * 10 ^ -8 * T_C + 1.72937731092437 * 10 ^ -5
End Function
