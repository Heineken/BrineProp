Attribute VB_Name = "IAPWS"
' IAPWS functions transcoded from MATLAB code of X-Steam
' http://www.mathworks.com/matlabcentral/fileexchange/9817-x-steam-thermodynamic-properties-of-water-and-steam

' Added non-private functions with Prefix "Water", that take an return values in SI units (function names as in Modelica.Media.Water.IF97_Utilities.BaseIF97.Basic.psat)
Option Explicit
Option Base 1

Private Const r = 0.461526 'kJ/(kg K)

' COEFFICIENTS FOR REGION 2 ?
Private Function J0()
    J0 = Array(0, 1, -5, -4, -3, -2, -1, 2, 3)
End Function
Private Function n0()
    n0 = Array(-9.6927686500217, 10.086655968018, -0.005608791128302, 0.071452738081455, -0.40710498223928, 1.4240819171444, -4.383951131945, -0.28408632460772, 0.021268463753307)
End Function
Private Function Ir()
    Ir = Array(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24)
End Function
Private Function Jr()
    Jr = Array(0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58)
End Function
Private Function nr()
    nr = Array(-1.7731742473213E-03, -0.017834862292358, -0.045996013696365, -0.057581259083432, -0.05032527872793, -3.3032641670203E-05, -1.8948987516315E-04, -3.9392777243355E-03, -0.043797295650573, -2.6674547914087E-05, 2.0481737692309E-08, 4.3870667284435E-07, -3.227767723857E-05, -1.5033924542148E-03, -0.040668253562649, -7.8847309559367E-10, 1.2790717852285E-08, 4.8225372718507E-07, 2.2922076337661E-06, -1.6714766451061E-11, -2.1171472321355E-03, -23.895741934104, -5.905956432427E-18, -1.2621808899101E-06, -0.038946842435739, 1.1256211360459E-11, -8.2311340897998, 1.9809712802088E-08, 1.0406965210174E-19, -1.0234747095929E-13, -1.0018179379511E-09, -8.0882908646985E-11, 0.10693031879409, -0.33662250574171, 8.9185845355421E-25, 3.0629316876232E-13, -4.2002467698208E-06, -5.9056029685639E-26, 3.7826947613457E-06, -1.2768608934681E-15, 7.3087610595061E-29, 5.5414715350778E-17, -9.436970724121E-07)
End Function

' COEFFICIENTS FOR REGION 1
Private Function J1()
    J1 = Array(-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41)
End Function
Private Function I1()
    I1 = Array(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32)
End Function
Private Function n1()
    n1 = Array(0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26)
End Function

' COEFFICIENTS FOR REGION 5
Private Function Ji0()
    Ji0 = Array(0, 1, -3, -2, -1, 2)
End Function
Private Function ni0()
    ni0 = Array(-13.179983674201, 6.8540841634434, -0.024805148933466, 0.36901534980333, -3.1161318213925, -0.32961626538917)
End Function
Private Function Iir()
    Iir = Array(1, 1, 1, 2, 3)
End Function
Private Function Jir()
    Jir = Array(0, 1, 3, 9, 3)
End Function
Private Function nir()
    nir = Array(-1.2563183589592E-04, 2.1774678714571E-03, -0.004594282089991, -3.9724828359569E-06, 1.2919228289784E-07)
End Function


Function region_pT(p_MPa As Double, ByVal T As Double) 'From MSL
  'return the current region (valid values: 1,2,3,5) in IF97, given pressure and temperature
    Dim p As Double: p = p_MPa * 10 ^ 6
'  output Integer region    "region (valid values: 1,2,3,5) in IF97, region 4 is impossible!";
    If p < 16529200# Then
      'test for regions 1,2,5
      If T > 1073.15 Then
        region_pT = 5
      ElseIf T > WaterTsat_p(p) Then
        region_pT = 2
      Else
        region_pT = 1
      End If
    Else
      'test for regions 1,2,3
      If T < 623.15 Then
        region_pT = 1
      ElseIf T < boundary23ofp(p) Then
        region_pT = 3
      Else
        region_pT = 2
      End If
    End If
    'Debug.Print "Region:" & region_pT
End Function
Function boundary23ofp(p As Double)
  'boundary function for region boundary between regions 2 and 3 (input pressure)
  Dim n: n = Array(348.05185628969, -1.1671859879975, 1.0192970039326E-03, 572.54459862746, 13.91883977887) 'polynomial coefficients for boundary between regions 2 and 3
  Const ptriple = 611.657
  Dim Pi As Double: Pi = p / 10 ^ 6 'dimensionless pressure
  If p < ptriple Then boundary23ofp = "IF97 medium function boundary23ofp called with too low pressure\n" _
                                    & "p = " & p & " Pa <= " & ptriple & " Pa (triple point pressure)"
  T = n(4) + ((Pi - n(5)) / n(3)) ^ 0.5
End Function


'Function region_pT(p, T) ' p in MPa, T in K
'Dim ps As Double
'    If T > 1073.15 And p < 10 And T < 2273.15 And p > 0.000611 Then
'        region_pT = 5
'    ElseIf T <= 1073.15 And T > 273.15 And p <= 100 And p > 0.000611 Then
'        If T > 623.15 Then
'            If p > B23p_T(T) Then
'                region_pT = 3
'                If T < 647.096 Then
'                    ps = p4_T(T)
'                    If Abs(p - ps) < 0.00001 Then
'                        region_pT = 4
'                    End If
'                End If
'            Else
'                region_pT = 2
'            End If
'        Else
'            ps = p4_T(T)
'            If Abs(p - ps) < 0.00001 Then
'                region_pT = 4
'            ElseIf p > ps Then
'                region_pT = 1
'            Else
'                region_pT = 2
'            End If
'        End If
'    Else
'        region_pT = "#IAPWS Region Error"
'    End If
'    If DebugMode Then
'        Debug.Print "Region: " & region_pT
'    End If
'End Function


Function Density_pT(p_Pa, ByVal T) ' p in Pa, T in K
    'Density_pT = 1 / v1_pT(p_Pa / 10 ^ 6, T) 'Pa->MPa
    Dim p As Double: p = p_Pa / 10 ^ 6 'Pa->MPa
    
    'Dim r: r = region_pT(p, T)
    Select Case region_pT(p, T)
        Case 1
            Density_pT = 1 / v1_pT(p, T)
        Case 2
            Density_pT = 1 / v2_pT(p, T)
        Case 3
            Density_pT = 1 / v3_ph(p, h3_pT(p, T))
        Case 4
            Density_pT = CVErr(xlErrNA)
        Case 5
            Density_pT = 1 / v5_pT(p, T)
        Case Else
            Density_pT = CVErr(xlErrNA)
    End Select
End Function


Private Function v1_pT(p, T)
' p in MPa, T in K
' Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
' 5 Equations for Region 1, Section. 5.1 Basic Equation
' Eqution 7, Table 3, Page 6
    'Dim I1, J1, n1
    ' I1 = Array(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32)
    ' J1 = Array(-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41)
    ' n1 = Array(0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26)
    Dim Pi As Double, tau As Double, gamma_der_pi As Double, i As Integer
    'Const R = 0.461526 ' kJ/(kg K)
    Pi = p / 16.53
    tau = 1386 / T
    gamma_der_pi = 0
    For i = 1 To 34
        gamma_der_pi = gamma_der_pi - n1(i) * I1(i) * (7.1 - Pi) ^ (I1(i) - 1) * (tau - 1.222) ^ J1(i)
    Next i
    v1_pT = r * T / p * Pi * gamma_der_pi / 1000
End Function

Function v2_pT(p, T)
    'Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
    '6 Equations for Region 2, Section. 6.1 Basic Equation
    'Table 11 and 12, Page 14 and 15
    Dim Pi As Double: Pi = p
    Dim tau  As Double: tau = 540 / T
    Dim g0_pi As Double: g0_pi = 1 / Pi
    Dim gr_pi As Double: gr_pi = 0
    Dim i As Integer
    For i = 1 To 43
        gr_pi = gr_pi + nr(i) * Ir(i) * Pi ^ (Ir(i) - 1) * (tau - 0.5) ^ Jr(i)
    Next i
    v2_pT = r * T / p * Pi * (g0_pi + gr_pi) / 1000
End Function

Function v3_ph(p, h)
    'Revised Supplementary Release on Backward Equations for the functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam
    '2004
    'Section 3.3 Backward Equations T(p,h) and v(p,h) for Subregions 3a and 3b
    'Boundary equation, Eq 1 Page 5
    Dim h3ab As Double, ps  As Double, hs As Double, vs As Double
    h3ab = 2014.64004206875 + 3.74696550136983 * p - 2.19921901054187E-02 * p ^ 2 + 8.7513168600995E-05 * p ^ 3
    Dim Ii, Ji, ni
    If h < h3ab Then
        ' Subregion 3a
        ' Eq 4, Table 6, Page 9
        Ii = Array(-12, -12, -12, -12, -10, -10, -10, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1, 0, 0, 1, 1, 1, 2, 2, 3, 4, 5, 8)
        Ji = Array(6, 8, 12, 18, 4, 7, 10, 5, 12, 3, 4, 22, 2, 3, 7, 3, 16, 0, 1, 2, 3, 0, 1, 0, 1, 2, 0, 2, 0, 2, 2, 2)
        ni = Array(5.29944062966028E-03, -0.170099690234461, 11.1323814312927, -2178.98123145125, -5.06061827980875E-04, 0.556495239685324, -9.43672726094016, -0.297856807561527, 93.9353943717186, 1.92944939465981E-02, 0.421740664704763, -3689141.2628233, -7.37566847600639E-03, -0.354753242424366, -1.99768169338727, 1.15456297059049, 5683.6687581596, 8.08169540124668E-03, 0.172416341519307, 1.04270175292927, -0.297691372792847, 0.560394465163593, 0.275234661176914, -0.148347894866012, -6.51142513478515E-02, -2.92468715386302, 6.64876096952665E-02, 3.52335014263844, -1.46340792313332E-02, -2.24503486668184, 1.10533464706142, -4.08757344495612E-02)
        ps = p / 100
        hs = h / 2100
        vs = 0
        Dim i As Integer
        For i = 1 To 32
            vs = vs + ni(i) * (ps + 0.128) ^ Ii(i) * (hs - 0.727) ^ Ji(i)
        Next i
        v3_ph = vs * 0.0028
    Else
        ' Subregion 3b
        ' Eq 5, Table 7, Page 9
        Ii = Array(-12, -12, -8, -8, -8, -8, -8, -8, -6, -6, -6, -6, -6, -6, -4, -4, -4, -3, -3, -2, -2, -1, -1, -1, -1, 0, 1, 1, 2, 2)
        Ji = Array(0, 1, 0, 1, 3, 6, 7, 8, 0, 1, 2, 5, 6, 10, 3, 6, 10, 0, 2, 1, 2, 0, 1, 4, 5, 0, 0, 1, 2, 6)
        ni = Array(-2.25196934336318E-09, 1.40674363313486E-08, 2.3378408528056E-06, -3.31833715229001E-05, 1.07956778514318E-03, -0.271382067378863, 1.07202262490333, -0.853821329075382, -2.15214194340526E-05, 7.6965608822273E-04, -4.31136580433864E-03, 0.453342167309331, -0.507749535873652, -100.475154528389, -0.219201924648793, -3.21087965668917, 607.567815637771, 5.57686450685932E-04, 0.18749904002955, 9.05368030448107E-03, 0.285417173048685, 3.29924030996098E-02, 0.239897419685483, 4.82754995951394, -11.8035753702231, 0.169490044091791, -1.79967222507787E-02, 3.71810116332674E-02, -5.36288335065096E-02, 1.6069710109252)
        ps = p / 100
        hs = h / 2800
        vs = 0
        For i = 1 To 30
            vs = vs + ni(i) * (ps + 0.0661) ^ Ii(i) * (hs - 0.72) ^ Ji(i)
        Next i
        v3_ph = vs * 0.0088
    End If
End Function

Function v5_pT(p, T)
    'Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
    'Basic Equation for Region 5
    'Eq 32,33, Page 36, Tables 37-41
    Dim tau As Double: tau = 1000 / T
    Dim Pi  As Double: Pi = p
    Dim gamma0_pi As Double: gamma0_pi = 1 / Pi
    Dim gammar_pi As Double: gammar_pi = 0
    Dim i As Integer
    For i = 1 To 5
        gammar_pi = gammar_pi + nir(i) * Iir(i) * Pi ^ (Iir(i) - 1) * tau ^ Jir(i)
    Next i
    v5_pT = r * T / p * Pi * (gamma0_pi + gammar_pi) / 1000
End Function

Public Function dynamicViscosity_pT(p_Pa, T) ' p in MPa, T in K
    Dim p As Double: p = p_Pa / 10 ^ 6 'Pa->MPa

    Dim h0: h0 = Array(0.5132047, 0.3205656, 0, 0, -0.7782567, 0.1885447)
    Dim h1: h1 = Array(0.2151778, 0.7317883, 1.241044, 1.476783, 0, 0)
    Dim h2: h2 = Array(-0.2818107, -1.070786, -1.263184, 0, 0, 0)
    Dim h3: h3 = Array(0.1778064, 0.460504, 0.2340379, -0.4924179, 0, 0)
    Dim h4: h4 = Array(-0.0417661, 0, 0, 0.1600435, 0, 0)
    Dim h5: h5 = Array(0, -0.01578386, 0, 0, 0, 0)
    Dim h6: h6 = Array(0, 0, 0, -0.003629481, 0, 0)
    
    'Calculate density.
    Dim rho As Double: rho = Density_pT(p_Pa, T)
    
    Dim rhos As Double: rhos = rho / 317.763
    Dim Ts As Double: Ts = T / 647.226
    Dim ps As Double: ps = p / 22.115
    
    'Check valid area
    If T > 900 + 273.15 Or (T > 600 + 273.15 And p > 300) Or (T > 150 + 273.15 And p > 350) Or p > 500 Then
        dynamicViscosity_pT = CVErr(xlErrNA)
        Return
    End If
    Dim my0 As Double
    my0 = Ts ^ 0.5 / (1 + 0.978197 / Ts + 0.579829 / (Ts ^ 2) - 0.202354 / (Ts ^ 3))
    Dim Sum As Double: Sum = 0
    Dim i As Integer
    For i = 0 To 5
        Sum = Sum + h0(i + 1) * (1 / Ts - 1) ^ i _
        + h1(i + 1) * (1 / Ts - 1) ^ i * (rhos - 1) ^ 1 _
        + h2(i + 1) * (1 / Ts - 1) ^ i * (rhos - 1) ^ 2 _
        + h3(i + 1) * (1 / Ts - 1) ^ i * (rhos - 1) ^ 3 _
        + h4(i + 1) * (1 / Ts - 1) ^ i * (rhos - 1) ^ 4 _
        + h5(i + 1) * (1 / Ts - 1) ^ i * (rhos - 1) ^ 5 _
        + h6(i + 1) * (1 / Ts - 1) ^ i * (rhos - 1) ^ 6
    Next i
    Dim my1 As Double: my1 = Exp(rhos * Sum)
    Dim mys As Double: mys = my0 * my1
    dynamicViscosity_pT = mys * 0.000055071
End Function


Function SpecificHeatCapacityCp_pT(p_Pa As Double, T_K As Double)
    Dim p As Double
    p = p_Pa / 10 ^ 6 'Pa->MPa
    Dim T As Double: T = CDbl(T_K) 'toSIunit_T(In2);
    Dim Region As Integer
    Region = region_pT(p, T)
        Select Case Region
        Case 1
            SpecificHeatCapacityCp_pT = cp1_pT(p, T) * 10 ^ 3 'kJ->J
        Case 2
            SpecificHeatCapacityCp_pT = cp2_pT(p, T) * 10 ^ 3 'kJ->J
        Case 3
            Dim hs As Double, rhos As Double
            hs = h3_pT(p, T)
            rhos = 1 / v3_ph(p, hs)
            SpecificHeatCapacityCp_pT = cp3_rhoT(rhos, T)
        Case 4
            SpecificHeatCapacityCp_pT = "#No cp for Region 4" 'CVErr(xlErrNA)
        Case 5
            SpecificHeatCapacityCp_pT = cp5_pT(p, T) * 10 ^ 3 'kJ->J
        Case Else
            SpecificHeatCapacityCp_pT = "#Unknown Region in (IAPWS.SpecificHeatCapacityCp_pT)" 'CVErr(xlErrNA)
    End Select
End Function

'Function SpecificHeatCapacityCp_pT(p_Pa As Double, T As Double) ' p in MPa, T in K
'   ' Dim p As Double: p = p_Pa / 10 ^ 6 'Pa->MPa
'    SpecificHeatCapacityCp_pT = h1_pT(p_Pa / 10 ^ 6, T) 'Pa->MPa
'End Function
Private Function cp1_pT(p As Double, T As Double)
    ' Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
    ' 5 Equations for Region 1, Section. 5.1 Basic Equation
    ' Eqution 7, Table 3, Page 6
     
'    Dim I1, J1, n1
'    I1 = Array(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32)
'    J1 = Array(-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41)
'    n1 = Array(0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26)
    Dim Pi As Double, tau As Double, gamma_der_tautau As Double, i As Integer
    Const r = 0.461526 'kJ/(kg K)
    Pi = p / 16.53
    tau = 1386 / T
    gamma_der_tautau = 0
    For i = 1 To 34
        gamma_der_tautau = gamma_der_tautau + (n1(i) * (7.1 - Pi) ^ I1(i) * J1(i) * (J1(i) - 1) * (tau - 1.222) ^ (J1(i) - 2))
    Next i
    cp1_pT = -r * tau ^ 2 * gamma_der_tautau
End Function

Function cp2_pT(p, T)
' Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
' 6 Equations for Region 2, Section. 6.1 Basic Equation
' Table 11 and 12, Page 14 and 15
Dim Pi As Double, tau As Double, g0_tautau As Double, gr_tautau As Double
    Pi = p
    tau = 540 / T
    g0_tautau = 0
    Dim i As Integer
    For i = 1 To 9
        g0_tautau = g0_tautau + n0(i) * J0(i) * (J0(i) - 1) * tau ^ (J0(i) - 2)
    Next i
    gr_tautau = 0
    For i = 1 To 43
        gr_tautau = gr_tautau + nr(i) * Pi ^ Ir(i) * Jr(i) * (Jr(i) - 1) * (tau - 0.5) ^ (Jr(i) - 2)
    Next i
    cp2_pT = -r * tau ^ 2 * (g0_tautau + gr_tautau)
End Function

Function cp3_rhoT(rho, T)
' Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
' 7 Basic Equation for Region 3, Section. 6.1 Basic Equation
' Table 30 and 31, Page 30 and 31
    Dim Ii, Ji, ni
    Ii = Array(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 7, 8, 9, 9, 10, 10, 11)
    Ji = Array(0, 0, 1, 2, 7, 10, 12, 23, 2, 6, 15, 17, 0, 2, 6, 7, 22, 26, 0, 2, 4, 16, 26, 0, 2, 4, 26, 1, 3, 26, 0, 2, 26, 2, 26, 2, 26, 0, 1, 26)
    ni = Array(1.0658070028513, -15.732845290239, 20.944396974307, -7.6867707878716, 2.6185947787954, -2.808078114862, 1.2053369696517, -8.4566812812502E-03, -1.2654315477714, -1.1524407806681, 0.88521043984318, -0.64207765181607, 0.38493460186671, -0.85214708824206, 4.8972281541877, -3.0502617256965, 0.039420536879154, 0.12558408424308, -0.2799932969871, 1.389979956946, -2.018991502357, -8.2147637173963E-03, -0.47596035734923, 0.0439840744735, -0.44476435428739, 0.90572070719733, 0.70522450087967, 0.10770512626332, -0.32913623258954, -0.50871062041158, -0.022175400873096, 0.094260751665092, 0.16436278447961, -0.013503372241348, -0.014834345352472, 5.7922953628084E-04, 3.2308904703711E-03, 8.0964802996215E-05, -1.6557679795037E-04, -4.4923899061815E-05)
   
    Dim tc As Double, pc As Double, rhoc As Double, Delta As Double, tau As Double, fitautau As Double, fidelta As Double, fideltadelta As Double
    tc = 647.096 ' K
    pc = 22.064 ' MPa
    rhoc = 322 ' kg/m3
    Delta = rho / rhoc
    tau = tc / T
    fitautau = 0
    fidelta = 0
    Dim fideltatau As Double '= 0
    fideltadelta = 0
    Dim i As Integer
    For i = 2 To 40
        fitautau = fitautau + ni(i) * Delta ^ Ii(i) * Ji(i) * (Ji(i) - 1) * tau ^ (Ji(i) - 2)
        fidelta = fidelta + ni(i) * Ii(i) * Delta ^ (Ii(i) - 1) * tau ^ Ji(i)
        fideltatau = fideltatau + ni(i) * Ii(i) * Delta ^ (Ii(i) - 1) * Ji(i) * tau ^ (Ji(i) - 1)
        fideltadelta = fideltadelta + ni(i) * Ii(i) * (Ii(i) - 1) * Delta ^ (Ii(i) - 2) * tau ^ Ji(i)
    Next i
    fidelta = fidelta + ni(1) / Delta
    fideltadelta = fideltadelta - ni(1) / (Delta ^ 2)
    cp3_rhoT = r * (-tau ^ 2 * fitautau + (Delta * fidelta - Delta * tau * fideltatau) ^ 2 / (2 * Delta * fidelta + Delta ^ 2 * fideltadelta))
End Function

Function cp5_pT(p As Double, T As Double)
' Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
' Basic Equation for Region 5
' Eq 32,33, Page 36, Tables 37-41
    Dim tau As Double, Pi As Double, gamma0_tautau As Double, gamma0_tau As Double, gammar_tau As Double, gammar_tautau As Double
    tau = 1000 / T
    Pi = p
    'gamma0_tautau = 0
    Dim i As Integer
    For i = 1 To 6
        gamma0_tautau = gamma0_tautau + ni0(i) * Ji0(i) * (Ji0(i) - 1) * tau ^ (Ji0(i) - 2)
    Next i
    'gammar_tautau = 0
    For i = 1 To 5
        gammar_tautau = gammar_tautau + nir(i) * Pi ^ Iir(i) * Jir(i) * (Jir(i) - 1) * tau ^ (Jir(i) - 2)
    Next i
    cp5_pT = -r * tau ^ 2 * (gamma0_tautau + gammar_tautau)
End Function

Function SpecificEnthalpy_pT(p_Pa, T_K)
    Dim p As Double
    p = p_Pa / 10 ^ 6 'Pa->MPa
    Dim T As Double: T = CDbl(T_K) ' toSIunit_T(In2)
    Select Case region_pT(p, T)
    Case 1
        SpecificEnthalpy_pT = h1_pT(p, T) * 10 ^ 3 'kJ->J
    Case 2
        SpecificEnthalpy_pT = h2_pT(p, T) * 10 ^ 3 'kJ->J
    Case 3
        SpecificEnthalpy_pT = h3_pT(p, T) * 10 ^ 3 'kJ->J
    Case 4
        SpecificEnthalpy_pT = CVErr(xlErrNA)
    Case 5
        SpecificEnthalpy_pT = h5_pT(p, T) * 10 ^ 3 'kJ->J
    Case Else
        SpecificEnthalpy_pT = CVErr(xlErrNA)
    End Select
End Function

Private Function h1_pT(p, T) As Double ' p in MPa, T in K
    ' Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
    ' 5 Equations for Region 1, Section. 5.1 Basic Equation
    ' Eqution 7, Table 3, Page 6
'    Dim I1, J1, n1
'    I1 = Array(0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 5, 8, 8, 21, 23, 29, 30, 31, 32)
'    J1 = Array(-2, -1, 0, 1, 2, 3, 4, 5, -9, -7, -1, 0, 1, 3, -3, 0, 1, 3, 17, -4, 0, 6, -5, -2, 10, -8, -11, -6, -29, -31, -38, -39, -40, -41)
'    n1 = Array(0.14632971213167, -0.84548187169114, -3.756360367204, 3.3855169168385, -0.95791963387872, 0.15772038513228, -0.016616417199501, 8.1214629983568E-04, 2.8319080123804E-04, -6.0706301565874E-04, -0.018990068218419, -0.032529748770505, -0.021841717175414, -5.283835796993E-05, -4.7184321073267E-04, -3.0001780793026E-04, 4.7661393906987E-05, -4.4141845330846E-06, -7.2694996297594E-16, -3.1679644845054E-05, -2.8270797985312E-06, -8.5205128120103E-10, -2.2425281908E-06, -6.5171222895601E-07, -1.4341729937924E-13, -4.0516996860117E-07, -1.2734301741641E-09, -1.7424871230634E-10, -6.8762131295531E-19, 1.4478307828521E-20, 2.6335781662795E-23, -1.1947622640071E-23, 1.8228094581404E-24, -9.3537087292458E-26)
    Dim Pi As Double, tau As Double, gamma_der_tau As Double, i As Integer
    
    'Const R = 0.461526 'kJ/(kg K)
    Pi = p / 16.53
    tau = 1386 / T
    gamma_der_tau = 0
    For i = 1 To 34
        gamma_der_tau = gamma_der_tau + (n1(i) * (7.1 - Pi) ^ I1(i) * J1(i) * (tau - 1.222) ^ (J1(i) - 1))
        'Debug.Print gamma_der_tau
    Next i
    h1_pT = r * T * tau * gamma_der_tau
End Function

Private Function h2_pT(p, T) As Double
    ' Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
    ' 6 Equations for Region 2, Section. 6.1 Basic Equation
    ' Table 11 and 12, Page 14 and 15
    Dim tau As Double, g0_tau As Double, Pi As Double, gr_tau As Double
'    Dim J0, n0, Ir, Jr, nr, R As Double
'    J0 = Array(0, 1, -5, -4, -3, -2, -1, 2, 3)
'    n0 = Array(-9.6927686500217, 10.086655968018, -0.005608791128302, 0.071452738081455, -0.40710498223928, 1.4240819171444, -4.383951131945, -0.28408632460772, 0.021268463753307)
'    Ir = Array(1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 5, 6, 6, 6, 7, 7, 7, 8, 8, 9, 10, 10, 10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24)
'    Jr = Array(0, 1, 2, 3, 6, 1, 2, 4, 7, 36, 0, 1, 3, 6, 35, 1, 2, 3, 7, 3, 16, 35, 0, 11, 25, 8, 36, 13, 4, 10, 14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58)
'    nr = Array(-1.7731742473213E-03, -0.017834862292358, -0.045996013696365, -0.057581259083432, -0.05032527872793, -3.3032641670203E-05, -1.8948987516315E-04, -3.9392777243355E-03, -0.043797295650573, -2.6674547914087E-05, 2.0481737692309E-08, 4.3870667284435E-07, -3.227767723857E-05, -1.5033924542148E-03, -0.040668253562649, -7.8847309559367E-10, 1.2790717852285E-08, 4.8225372718507E-07, 2.2922076337661E-06, -1.6714766451061E-11, -2.1171472321355E-03, -23.895741934104, -5.905956432427E-18, -1.2621808899101E-06, -0.038946842435739, 1.1256211360459E-11, -8.2311340897998, 1.9809712802088E-08, 1.0406965210174E-19, -1.0234747095929E-13, -1.0018179379511E-09, -8.0882908646985E-11, 0.10693031879409, -0.33662250574171, 8.9185845355421E-25, 3.0629316876232E-13, -4.2002467698208E-06, -5.9056029685639E-26, 3.7826947613457E-06, -1.2768608934681E-15, 7.3087610595061E-29, 5.5414715350778E-17, -9.436970724121E-07)
'    R = 0.461526 'kJ/(kg K)
    
    
    Pi = p
    tau = 540 / T
    g0_tau = 0
    Dim i As Integer
    For i = 1 To 9
        g0_tau = g0_tau + n0(i) * J0(i) * tau ^ (J0(i) - 1)
    Next i
    gr_tau = 0
    For i = 1 To 43 'CHANGED HERE: reduced by 1, because VBA-Arrays start at 0
        gr_tau = gr_tau + nr(i) * Pi ^ Ir(i) * Jr(i) * (tau - 0.5) ^ (Jr(i) - 1)
    Next i
    h2_pT = r * T * tau * (g0_tau + gr_tau)
End Function

Private Function h3_pT(p As Double, ByVal T As Double) As Double
    'Not avalible with if 97
    'Solve function T3_ph-T=0 with half interval method.
    'ver2.6 Start corrected bug
    Dim High_Bound As Double, Low_Bound As Double, Ts As Double, hs As Double
    If p < 22.06395 Then  'Bellow tripple point
      Ts = T4_p(p)    'Saturation temperature
      If T <= Ts Then     'Liquid side
        High_Bound = h4L_p(p) 'Max h är liauid h.
        Low_Bound = h1_pT(p, 623.15)
      Else
        Low_Bound = h4V_p(p)  'Min h är Vapour h.
        High_Bound = h2_pT(p, B23T_p(p))
      End If
    Else                  'Above tripple point. R3 from R2 till R3
      Low_Bound = h1_pT(p, 623.15)
      High_Bound = h2_pT(p, B23T_p(p))
    End If
    'ver2.6 End corrected bug
    Ts = T + 1
    While Abs(T - Ts) > 0.00001
        hs = (Low_Bound + High_Bound) / 2
        Ts = T3_ph(p, hs)
        If Ts > T Then
            High_Bound = hs
        Else
            Low_Bound = hs
        End If
    Wend
    h3_pT = hs
End Function

Function T3_ph(p As Double, h As Double)
    'Revised Supplementary Release on Backward Equations for the functions T(p,h), v(p,h) and T(p,s), v(p,s) for Region 3 of the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam
    '2004
    'Section 3.3 Backward Equations T(p,h) and v(p,h) for Subregions 3a and 3b
    'Boundary equation, Eq 1 Page 5
    h3ab = 2014.64004206875 + 3.74696550136983 * p - 2.19921901054187E-02 * p ^ 2 + 8.7513168600995E-05 * p ^ 3
    If h < h3ab Then
        'Subregion 3a
        'Eq 2, Table 3, Page 7
        Ii = Array(-12, -12, -12, -12, -12, -12, -12, -12, -10, -10, -10, -8, -8, -8, -8, -5, -3, -2, -2, -2, -1, -1, 0, 0, 1, 3, 3, 4, 4, 10, 12)
        Ji = Array(0, 1, 2, 6, 14, 16, 20, 22, 1, 5, 12, 0, 2, 4, 10, 2, 0, 1, 3, 4, 0, 2, 0, 1, 1, 0, 1, 0, 3, 4, 5)
        ni = Array(-1.33645667811215E-07, 4.55912656802978E-06, -1.46294640700979E-05, 6.3934131297008E-03, 372.783927268847, -7186.54377460447, 573494.7521034, -2675693.29111439, -3.34066283302614E-05, -2.45479214069597E-02, 47.8087847764996, 7.64664131818904E-06, 1.28350627676972E-03, 1.71219081377331E-02, -8.51007304583213, -1.36513461629781E-02, -3.84460997596657E-06, 3.37423807911655E-03, -0.551624873066791, 0.72920227710747, -9.92522757376041E-03, -0.119308831407288, 0.793929190615421, 0.454270731799386, 0.20999859125991, -6.42109823904738E-03, -0.023515586860454, 2.52233108341612E-03, -7.64885133368119E-03, 1.36176427574291E-02, -1.33027883575669E-02)
        ps = p / 100
        hs = h / 2300
        Ts = 0
        For i = 1 To 31
            Ts = Ts + ni(i) * (ps + 0.24) ^ Ii(i) * (hs - 0.615) ^ Ji(i)
        Next i
        T3_ph = Ts * 760
    Else
        'Subregion 3b
        'Eq 3, Table 4, Page 7,8
        Ii = Array(-12, -12, -10, -10, -10, -10, -10, -8, -8, -8, -8, -8, -6, -6, -6, -4, -4, -3, -2, -2, -1, -1, -1, -1, -1, -1, 0, 0, 1, 3, 5, 6, 8)
        Ji = Array(0, 1, 0, 1, 5, 10, 12, 0, 1, 2, 4, 10, 0, 1, 2, 0, 1, 5, 0, 4, 2, 4, 6, 10, 14, 16, 0, 2, 1, 1, 1, 1, 1)
        ni = Array(3.2325457364492E-05, -1.27575556587181E-04, -4.75851877356068E-04, 1.56183014181602E-03, 0.105724860113781, -85.8514221132534, 724.140095480911, 2.96475810273257E-03, -5.92721983365988E-03, -1.26305422818666E-02, -0.115716196364853, 84.9000969739595, -1.08602260086615E-02, 1.54304475328851E-02, 7.50455441524466E-02, 2.52520973612982E-02, -6.02507901232996E-02, -3.07622221350501, -5.74011959864879E-02, 5.03471360939849, -0.925081888584834, 3.91733882917546, -77.314600713019, 9493.08762098587, -1410437.19679409, 8491662.30819026, 0.861095729446704, 0.32334644281172, 0.873281936020439, -0.436653048526683, 0.286596714529479, -0.131778331276228, 6.76682064330275E-03)
        hs = h / 2800
        ps = p / 100
        Ts = 0
        For i = 1 To 33
            Ts = Ts + ni(i) * (ps + 0.298) ^ Ii(i) * (hs - 0.72) ^ Ji(i)
        Next i
        T3_ph = Ts * 860
    End If
End Function

Function h4L_p(p As Double)
    If (p > 0.000611657 & p < 22.06395) = 1 Then
        Ts = T4_p(p)
        If p < 16.529 Then
            h4L_p = h1_pT(p, Ts)
        Else
            'Iterate to find the the backward solution of p3sat_h
            Low_Bound = 1670.858218
            High_Bound = 2087.23500164864
            ps = -1000
            While Abs(p - ps) > 0.00001
                hs = (Low_Bound + High_Bound) / 2
                ps = p3sat_h(hs)
                If ps > p Then
                    High_Bound = hs
                Else
                    Low_Bound = hs
                End
            Wend
                
            h4L_p = hs
        End
    Else
        h4L_p = -99999
    End
End Function

Function B23T_p(p As Double)
'Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam 1997
'Section 4 Auxiliary Equation for the Boundary between Regions 2 and 3
'Eq 6, Page 6
    B23T_p = 572.54459862746 + ((p - 13.91883977887) / 1.0192970039326E-03) ^ 0.5
End Function


Private Function h5_pT(p, T) As Double
    ' Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
    ' Basic Equation for Region 5
    ' Eq 32,33, Page 36, Tables 37-41
    'R = 0.461526 'kJ/(kg K)
    Dim tau As Double, Pi As Double, gamma0_tau As Double, gammar_tau As Double
    tau = 1000 / T
    Pi = p
    gamma0_tau = 0
    For i = 1 To 6
        gamma0_tau = gamma0_tau + ni0(i) * Ji0(i) * tau ^ (Ji0(i) - 1)
    Next i
    gammar_tau = 0
    For i = 1 To 5
        gammar_tau = gammar_tau + nir(i) * Pi ^ Iir(i) * Jir(i) * tau ^ (Jir(i) - 1)
    Next
    h5_pT = r * T * tau * (gamma0_tau + gammar_tau)
End Function


Function Waterpsat_T(T) As Double
    Waterpsat_T = p4_T(T) * 10 ^ 6 'MPa->Pa
End Function
Private Function p4_T(T) As Double
    ' Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
    ' Section 8.1 The Saturation-Pressure Equation
    ' Eq 30, Page 33
    Dim teta As Double, a As Double, b As Double, c As Double
    teta = T - 0.23855557567849 / (T - 650.17534844798)
    a = teta ^ 2 + 1167.0521452767 * teta - 724213.16703206
    b = -17.073846940092 * teta ^ 2 + 12020.82470247 * teta - 3232555.0322333
    c = 14.91510861353 * teta ^ 2 - 4823.2657361591 * teta + 405113.40542057
    p4_T = (2 * c / (-b + (b ^ 2 - 4 * a * c) ^ 0.5)) ^ 4 'MPa
End Function

Function WaterTsat_p(p_Pa) 'Saturation temperature
    WaterTsat_p = T4_p(p_Pa / 10 ^ 6) 'Pa->MPa
End Function
Private Function T4_p(p) As Double
    ' Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam, September 1997
    ' Section 8.2 The Saturation-Temperature Equation
    ' Eq 31, Page 34
    Dim beta As Double, E As Double, f As Double, G As Double, d As Double
    beta = p ^ 0.25
    E = beta ^ 2 - 17.073846940092 * beta + 14.91510861353
    f = 1167.0521452767 * beta ^ 2 + 12020.82470247 * beta - 4823.2657361591
    G = -724213.16703206 * beta ^ 2 - 3232555.0322333 * beta + 405113.40542057
    d = 2 * G / (-f - (f ^ 2 - 4 * E * G) ^ 0.5)
    T4_p = (650.17534844798 + d - ((650.17534844798 + d) ^ 2 - 4 * (-0.23855557567849 + 650.17534844798 * d)) ^ 0.5) / 2
End Function

Function Waterh_satv_p(p_Pa) 'Saturated vapour enthalpy
    Waterh_satv_p = h4V_p(p_Pa / 10 ^ 6) * 10 ^ 3 'Pa->MPa, kJ->J
End Function
Private Function h4V_p(p)
    Dim Ts As Double, Low_Bound As Double, High_Bound As Double, ps As Double, hs As Double
    If p > 0.000611657 And p < 22.06395 Then
        Ts = T4_p(p)
        If p < 16.529 Then
            h4V_p = h2_pT(p, Ts)
        Else
            'Iterate to find the the backward solution of p3sat_h
            Low_Bound = 2087.23500164864
            High_Bound = 2563.592004 + 5
            ps = -1000
            Do While Abs(p - ps) > 0.000001
                hs = (Low_Bound + High_Bound) / 2
                ps = p3sat_h(hs)
                If ps < p Then
                    High_Bound = hs
                Else
                    Low_Bound = hs
                End If
            Loop
            h4V_p = hs
        End If
    Else
        h4V_p = -99999
    End If
End Function

Private Function p3sat_h(h)
    ' Revised Supplementary Release on Backward Equations for the functions T(p,h), v(p,h)  & T(p,s), v(p,s) for Region 3 of the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water  & Steam 2004
    ' Section 4 Boundary Equations psat(h)  & psat(s) for the Saturation Lines of Region 3
    ' Se pictures Page 17, Eq 10, Table 17, Page 18
    Dim Ii: Ii = Array(0, 1, 1, 1, 1, 5, 7, 8, 14, 20, 22, 24, 28, 36)
    Dim Ji: Ji = Array(0, 1, 3, 4, 36, 3, 0, 24, 16, 16, 3, 18, 8, 24)
    Dim ni: ni = Array(0.600073641753024, -9.36203654849857, 24.6590798594147, -107.014222858224, -91582131580576.8, -8623.32011700662, -23.5837344740032, 2.52304969384128E+17, -3.89718771997719E+18, -3.33775713645296E+22, 35649946963.6328, -1.48547544720641E+26, 3.30611514838798E+18, 8.13641294467829E+37)
    Dim hs As Double: hs = h / 2600
    Dim ps As Double ': ps = 0
    Dim i As Interior
    For i = 0 To 13 'CHANGED HERE: reduced by 1, because VBA-Arrays start at 0
        ps = ps + ni(i) * (hs - 1.02) ^ Ii(i) * (hs - 0.608) ^ Ji(i)
    End
    p3sat_h = ps * 22
End Function

Private Function B23p_T(T)
    ' Release on the IAPWS Industrial formulation 1997 for the Thermodynamic Properties of Water and Steam
    ' 1997
    ' Section 4 Auxiliary Equation for the Boundary between Regions 2 and 3
    'Eq 5, Page 5
    B23p_T = 348.05185628969 - 1.1671859879975 * T + 1.0192970039326E-03 * T ^ 2
End Function

