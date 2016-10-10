Attribute VB_Name = "Brine_liq"
' Properties of liquid phase (H2O + NaCl + KCl + CaCl) properties
' (density, viscosity, specific heat capacity, specific enthalpy)
' uses IAPWS pure water properties
'
' by Henning Francke francke@gfz-potsdam.de
' 2014 GFZ Potsdam

Option Explicit
Option Base 1
Public Const nX_salt = 3
Public Const nX = nX_salt + 1


Private ignoreLimitSalt_p(1 To nX_salt) As Boolean 'TODO Limithandling vereinheitlichen
Private ignoreLimitSalt_T(1 To nX_salt) As Boolean
Private ignoreLimitSalt_visc(1 To nX_salt) As Boolean
Private Const ignoreLimit_h_KCl_Tmin = True 'ignore Tmin in appMolarEnthalpy_KCl_White and appMolarHeatCapacity_KCl_White
Private Const ignoreLimit_h_CaCl2_Tmin = True 'ignore Tmin in appMolarEnthalpy_CaCl2_White and appMolarHeatCapacity_CaCl2_White

Public Salts(1 To 3) As SaltProps

Private Sub init()
    'DefineSalts
    Salts(1) = NaCl
    Salts(2) = KCl
    Salts(3) = CaCl2
    DefineLimits
End Sub

Private Sub DefineLimits()
    ignoreLimitSalt_p(1) = False
    ignoreLimitSalt_p(2) = True
    ignoreLimitSalt_p(3) = True
    ignoreLimitSalt_T(1) = False
    ignoreLimitSalt_T(2) = False
    ignoreLimitSalt_T(3) = False
    ignoreLimitSalt_visc(1) = False
    ignoreLimitSalt_visc(2) = False
    ignoreLimitSalt_visc(3) = True
End Sub

Function specificHeatCapacityCp(p As Double, T As Double, Xin)  'calculation of liquid specific heat capacity from apparent molar heat capacities
    If DebugMode Then
        Debug.Print "Running specificHeatCapacityCp(" & p / 10 ^ 5 & " bar," & T - 273.15 & " °C, X={" & Xin(1) & Xin(2) & Xin(3); "})"
    End If
  
    If Salts(1).MM = 0 Then 'if Salt data not defined
        init
    End If
    
    Dim X
    X = CheckMassVector(Xin, nX)
    If VarType(X) = vbString Then
        specificHeatCapacityCp = X & " (Brine_liq.specificHeatCapacityCp)"
        Exit Function
    End If
    
    Dim cp_Driesner
    cp_Driesner = specificHeatCapacity_pTX_Driesner(p, T, X(1) / (X(1) + X(nX_salt + 1)))
    If VarType(cp_Driesner) = vbString Then
          specificHeatCapacityCp = cp_Driesner
          Exit Function
    End If

    If X(2) = 0 And X(3) = 0 Then 'Only NaCl
        specificHeatCapacityCp = cp_Driesner
        Exit Function
    End If
    
    Dim b: b = massFractionsToMolalities(X, MM_vec)
'    If b(1) = -1 Then
'      specificHeatCapacityCp = "#Error in massFractionsToMolalities (no water nor pure NaCl)"
'      Exit Function
'    End If

    Dim Cp_appmol: ReDim Cp_appmol(2 To nX_salt) 'Apparent molar heat capacity of salts
    'Cp_appmol(1) = 0
    If b(2) > 0 Then
        Cp_appmol(2) = appMolarHeatCapacity_KCl_White(T, b(2))
        If VarType(Cp_appmol(2)) = vbString Then
            specificHeatCapacityCp = Cp_appmol(2)
            Exit Function
        End If

    End If
    If b(3) > 0 Then
        Cp_appmol(3) = appMolarHeatCapacity_CaCl2_White(T, b(3))
        If VarType(Cp_appmol(3)) = vbString Then
            specificHeatCapacityCp = Cp_appmol(3)
            Exit Function
        End If
    End If
    
    specificHeatCapacityCp = (X(i_NaCl) + X(nX)) * cp_Driesner + X(nX) * (b(2) * Cp_appmol(2) + b(3) * Cp_appmol(3)) 'for length(H_appmol)=2
End Function


Private Function appMolarHeatCapacity_KCl_White(T As Double, ByVal mola As Double) '2D-fit Reproduction of measurements of heat capacity of KCl solution
    appMolarHeatCapacity_KCl_White = appMolarX_KCl_White(T, mola, "cp")
End Function
Private Function appMolarEnthalpy_KCl_White(T As Double, ByVal mola As Double) '2D-fit Reproduction of measurements of heat capacity of KCl solution
    appMolarEnthalpy_KCl_White = appMolarX_KCl_White(T, mola, "h")
End Function
Private Function appMolarX_KCl_White(T As Double, ByVal mola As Double, what As String) '2D-fit Reproduction of measurements of heat capacity of KCl solution
'White, D., Ryan, M., Armstrong, M.C., Gates, J. and Wood, R. (1987b) 'Heat capacities of aqueous KCl
'from 325 to 600 K at 17.9 MPa', The Journal of Chemical Thermodynamics, vol. 19, no. 10, oct, pp. 1023-1030,
'DOI: 10.1016/0021-9614(87)90012-7.
    Dim b As Double: b = 0.09818
    Dim c As Double: c = -1.244
    Dim k As Double: k = -327.9
    Dim l As Double: l = -1.31 * 10 ^ 5
    Dim M As Double: M = 628.8
    
    If outOfRangeMode > 0 Then
        Dim T_min As Double: T_min = 306
        Dim T_max As Double: T_max = 603

        Dim MsgTxt As String
        If Not ((ignoreLimit_h_KCl_Tmin Or T >= T_min) And T <= T_max) Then
            MsgTxt = "#T=" & T - 273.15 & "°C, must be within " & T_min - 273.15 & "..." & T_max - 273.15 & "°C  (appMolarX_KCl_White)"
            Debug.Print MsgTxt
            If outOfRangeMode = 2 Then
                appMolarX_KCl_White = MsgTxt
                Exit Function
            End If
        End If
    End If
    
    If what = "cp" Then
        appMolarX_KCl_White = (mola ^ b + c) * (k - l * (M - T) ^ (-1))
    ElseIf what = "h" Then
        Const T0 = 293.16 'Temperature at which HeatOfSolution is taken
        appMolarX_KCl_White = HeatOfSolution_KCl_Sanahuja1986(T0) + (mola ^ b + c) * (k * (T - T0) + l * Log((M - T) / (M - T0)))
    End If
End Function

Private Function appMolarHeatCapacity_CaCl2_White(T As Double, ByVal mola As Double) '2D-fit Reproduction of measurements of heat capacity of KCl solution
    appMolarHeatCapacity_CaCl2_White = appMolarX_CaCl2_White(T, mola, "cp")
End Function
Private Function appMolarEnthalpy_CaCl2_White(T As Double, ByVal mola As Double) '2D-fit Reproduction of measurements of heat capacity of KCl solution
    appMolarEnthalpy_CaCl2_White = appMolarX_CaCl2_White(T, mola, "h")
End Function
Private Function appMolarX_CaCl2_White(T As Double, ByVal mola As Double, what As String) '2D-fit Reproduction of measurements of heat capacity of KCl solution
'White, D., Doberstein, A., Gates, J., Tillett, D. and Wood, R. (1987a) 'Heat capacity of aqueous CaCl2
'from 306 to 603 K at 17.5 MPa', The Journal of Chemical Thermodynamics, vol. 19, no. 3, mar, pp. 251-259,
'DOI: 10.1016/0021-9614(87)90132-7.
    Dim b  As Double: b = -0.001977
    Dim c  As Double: c = -0.9958
    Dim k  As Double: k = 1373
    Dim l  As Double: l = 6736000#
    Dim M  As Double: M = 628
    
    If outOfRangeMode > 0 Then
        Dim T_min As Double: T_min = 306
        Dim T_max As Double: T_max = 603

        Dim MsgTxt As String
        If Not ((ignoreLimit_h_CaCl2_Tmin Or T >= T_min) And T <= T_max) Then
            MsgTxt = "#T=" & T - 273.15 & "°C, must be within " & T_min - 273.15 & "..." & T_max - 273.15 & "°C  (appMolarX_CaCl2_White)"
            Debug.Print MsgTxt
            If outOfRangeMode = 2 Then
                appMolarX_CaCl2_White = MsgTxt
                Exit Function
            End If
        End If
    End If
    
    If what = "cp" Then
        appMolarX_CaCl2_White = (mola ^ b + c) * (k - l * (M - T) ^ (-1))
    ElseIf what = "h" Then
        Const T0 = 298.15 'Temperature at which HeatOfSolution is take
        Const Delta_h_solution_CaCl2 = 81850 '[J/mol_CaCl2] @ 298.15K Sinke1985 http://dx.doi.org/10.1016/0021-9614(85)90083-7
        appMolarX_CaCl2_White = Delta_h_solution_CaCl2 + (mola ^ b + c) * (k * (T - T0) + l * Log((M - T) / (M - T0)))
    End If
End Function


Function specificHeatCapacity_pTX_Driesner(p As Double, T As Double, x_NaCl As Double) 'cp calculation according to Driesner 2007 et al: 0-1000°C; 0.1-500MPa (doi:10.1016/j.gca.2007.05.026)"
  Dim T_Scale_h
  Dim q_2 As Double
  T_Scale_h = T_Scale_h_Driesner(p, T, x_NaCl, q_2)
  If VarType(T_Scale_h) = vbString Then
    specificHeatCapacity_pTX_Driesner = T_Scale_h
  Else
      specificHeatCapacity_pTX_Driesner = q_2 * IAPWS.SpecificHeatCapacityCp_pT(p, CDbl(T_Scale_h)) 'J/(kg·K)
  End If
End Function

Private Function T_Scale_h_Driesner(p As Double, T As Double, X_NaCl_ As Double, Optional ByRef q_2 As Double) 'enthalpy calculation according to Driesner 2007 et al: 0-1000°C; 0.1-500MPa (doi:10.1016/j.gca.2007.05.026)"
  If Not X_NaCl_ > 0 Then 'pure water -> trivial
    T_Scale_h_Driesner = T
    Exit Function
  End If
  
  Const p_min = 1 'bar
  Const p_max = 1000 'bar
  Const T_min = 0 '°C
  Const T_max = 1000 '°C

  Dim T_C As Double: T_C = T - 273.15
  Dim p_bar As Double: p_bar = p / 10 ^ 5
  Dim q_21 As Double
  Dim q_22 As Double
  Dim q_20 As Double
  Dim q_23 As Double
  Dim q_11 As Double
  Dim q_10 As Double
  Dim q_12 As Double
  Dim q_1 As Double
  Dim x_NaCl  As Double 'mol fraction

   If outOfRangeMode > 0 Then
    Dim MsgTxt As String
      If Not (p_bar >= p_min And p_bar <= p_max) Then
        MsgTxt = "#Pressure is " & p_bar & " bar, but must be between " & p_min & " bar and " & p_max & " bar  (T_Scale_h_Driesner)"
        Debug.Print MsgTxt
        If outOfRangeMode = 2 Then
            T_Scale_h_Driesner = MsgTxt
            Exit Function
        End If
      End If
      If Not (T_C >= T_min And T_C <= T_max) Then
        MsgTxt = "#Temperature is " & T_C & "°C, but must be between " & T_min & "°C and " & T_max & "°C  (T_Scale_h_Driesner)"
        Debug.Print MsgTxt
        If outOfRangeMode = 2 Then
            T_Scale_h_Driesner = MsgTxt
            Exit Function
        End If
      End If
   End If


  If X_NaCl_ = 0 Then
    x_NaCl = 0
  Else
    x_NaCl = 1 / (M_NaCl / M_H2O * (1 / X_NaCl_ - 1) + 1) 'mass fraction -> mol fraction
  End If

'CALCULATION OF EQUIVALENT TEMPERATURE_h
  q_21 = -1.69513 - 0.000452781 * p_bar - 0.0000000604279 * p_bar ^ 2
  q_22 = 0.0612567 + 0.0000188082 * p_bar

  q_20 = 1 - q_21 * Sqr(q_22) 'x_NaCl = 0 results in q_2=1

  q_23 = -q_20 - q_21 * Sqr(1 + q_22) + 0.241022 + 0.0000345087 * p_bar - 0.00000000428356 * p_bar ^ 2 'x_NaCl = 1 is pure NaCl
  q_2 = q_20 + q_21 * Sqr(x_NaCl + q_22) + q_23 * x_NaCl

  q_10 = 47.9048 - 0.00936994 * p_bar + 0.00000651059 * p_bar ^ 2
  q_11 = -32.1724 + 0.0621255 * p_bar

  q_12 = -q_10 - q_11
  
  q_1 = q_10 + q_11 * (1 - x_NaCl) + q_12 * (1 - x_NaCl) ^ 2 'x_NaCl=0 results in q_1=0 / x_NaCl = 1 is pure NaCl

  T_Scale_h_Driesner = T_C * q_2 + q_1 + 273.15 ' K
End Function


Function specificEnthalpy(p As Double, T As Double, Xin) 'enthalpy calculation
'based on Driesner enthalpy function NaCl and pressure independent 2D-fits (T,b) for cp measurement data for KCl and CaCl2 solutions
    If DebugMode Then
        Debug.Print "Running SpecificEnthalpy_vec(" & p / 100000# & " bar," & T - 273.15 & " °C, X={" & Xin(1) & Xin(2) & Xin(3); "})"
    End If

    If Salts(1).MM = 0 Then 'if Salt data not defined
      init
    End If
    Dim X
    X = CheckMassVector(Xin, nX)
    If VarType(X) = vbString Then
        specificEnthalpy = X & " (Brine_liq.specificEnthalpy)"
        Exit Function
    End If
    
    Dim b '() As Double:
    b = massFractionsToMolalities(X, MM_vec)
    
    Dim h_Driesner: h_Driesner = SpecificEnthalpy_Driesner(p, T, X(1) / (X(1) + X(nX)))
    
    If VarType(h_Driesner) = vbString Then 'Error
        specificEnthalpy = h_Driesner
        Exit Function
    End If
    
    If X(2) = 0 And X(3) = 0 Then 'Only NaCl
        specificEnthalpy = h_Driesner
        Exit Function
    End If

'    If b = -1 Then
'      specificEnthalpy = "#Error in massFractionsToMolalities (no water nor pure NaCl)"
'      Exit Function
'    End If
    
    Dim H_appmol: ReDim H_appmol(2 To nX - 1) 'Apparent molar heat capacity of salts
    
    'H_appmol(1) = 0 'included in Driesner enthalpy
    If X(2) > 0 Then
        H_appmol(2) = appMolarEnthalpy_KCl_White(T, b(2))
        If VarType(H_appmol(2)) = vbString Then
            specificEnthalpy = H_appmol(2)
            Exit Function
        End If
    End If
    If X(3) > 0 Then
        H_appmol(3) = appMolarEnthalpy_CaCl2_White(T, b(3))
        If VarType(H_appmol(3)) = vbString Then
            specificEnthalpy = H_appmol(3)
            Exit Function
        End If
    End If

'    SpecificEnthalpy = (X(i_NaCl) + X(nX)) * h_Driesner + X(nX) * ScalProd(SubArray(b, 2, 3), SubArray(H_appmol, 2, 3))
    specificEnthalpy = (X(i_NaCl) + X(nX)) * h_Driesner + X(nX) * (b(2) * H_appmol(2) + b(3) * H_appmol(3)) 'for length(H_appmol)=2

End Function



Private Function SpecificEnthalpy_Driesner(p As Double, T As Double, x_NaCl As Double) 'enthalpy calculation according to Driesner 2007 et al: 0-1000°C; 0.1-500MPa (doi:10.1016/j.gca.2007.05.026)"
  Dim T_Scale_h:  T_Scale_h = T_Scale_h_Driesner(p, T, x_NaCl)
    If VarType(T_Scale_h) = vbString Then
        SpecificEnthalpy_Driesner = T_Scale_h
    Else
        SpecificEnthalpy_Driesner = IAPWS.SpecificEnthalpy_pT(p, CDbl(T_Scale_h)) 'J/(kg)
    End If
End Function





Private Function HeatOfSolution_KCl_Sanahuja1986(T As Double)
'2nd degree polynomial fit from Sanahuja1986 http://dx.doi.org/10.1016/0021-9614(86)90063-7
  HeatOfSolution_KCl_Sanahuja1986 = -0.0102 * T ^ 2 + 6.12926 * T - 904.206
End Function


Function MM_vec()
 'generates double vector of molar masses
    If Salts(1).MM = 0 Then 'if Salt data not defined
        init
    End If

    Dim MM(1 To nX_salt + 1) As Double
    Dim i As Integer
    For i = 1 To nX_salt
        MM(i) = Salts(i).MM
    Next i
    MM(nX_salt + 1) = M_H2O
    MM_vec = MM
End Function

Function nM_vec() As Double()
 'generates double vector of molar masses
    Dim nM(1 To nX_salt + 1) As Double
    Dim i As Integer
    For i = 1 To nX_salt
        nM(i) = Salts(i).nM
    Next i
    nM(nX_salt + 1) = H2O.nM
    nM_vec = nM
End Function

'ABOVE COMPOSITION SPECIFIC
'BELOW GENERIC=============================================================================

Function dynamicViscosity(p_Pa As Double, T As Double, Xin) 'brine density
  'Multisalt-Version of viscosity calculation according to Duan et al 2009 and Zhang et al 1997: Considers NaCl and KCL, with geometric mixture rule"
  'doi:10.1007/s10765-009-0646-7
    If DebugMode Then
      Debug.Print "Running Viscosity(" & p_Pa / 10 ^ 5 & " bar, " & T - 273.15 & " °C, X={" & Xin(1) & Xin(2) & Xin(3); "})"
    End If
    
    If Salts(1).MM = 0 Then
      init
    End If
    
    Dim T_C As Double: T_C = T - 273.15
    Dim p_bar As Double: p_bar = p_Pa / 10 ^ 5
    
    Dim X
    X = CheckMassVector(Xin, nX)
    If VarType(X) = vbString Then
        dynamicViscosity = X & " (Brine_liq.dynamicViscosity)"
        Exit Function
    End If
    
    Dim molalities: molalities = massFractionsToMolalities(X, MM_vec)
'    If molalities(1) = -1 Then
'      dynamicViscosity = "#Error in massFractionsToMolalities"
'      Exit Function
'    End If
    
     'viscosity calculation
    Dim eta_H2O As Double: eta_H2O = IAPWS.dynamicViscosity_pT(p_Pa, T)
    
     'for pure water skip the whole calculation and return water viscosity
     Dim eta As Double: eta = eta_H2O
     If Application.Max(SubArray(X, 1, nX_salt)) <= 0.1 ^ 8 Then
       dynamicViscosity = eta_H2O
       Exit Function 'pure water -> skip the rest
     End If
    
    Dim p_min As Double: p_min = 1
    Dim p_max As Double: p_max = 1000
    Dim T_min As Double: T_min = 0
    Dim T_max As Double: T_max = 400
    
    Dim MsgTxt As String
     
      If outOfRangeMode > 0 Then
         If Not (p_bar >= p_min And p_bar <= p_max) Then
           MsgTxt = "#Pressure " & p_bar & " bar out of limits {" & p_min & "..." & p_max & "} bar (dynamicViscosity_DuanZhang_pTXd)"
           Debug.Print MsgTxt
           If outOfRangeMode = 2 Then
               dynamicViscosity = MsgTxt
               Exit Function
           End If
         End If
         If Not (T_C > T_min And T_C <= T_max) Then
           MsgTxt = "#Temperature " & T_C & "°C, out of limits {" & T_min & "..." & T_max & "}°C (dynamicViscosity_DuanZhang_pTXd)"
           Debug.Print MsgTxt
           If outOfRangeMode = 2 Then
               dynamicViscosity = MsgTxt
               Exit Function
           End If
         End If
      End If
    

    
    Dim i As Integer
    For i = 1 To nX_salt
       If X(i) > 0 Then
         If outOfRangeMode > 0 Then
           If molalities(i) < 0 Or molalities(i) > Salts(i).mola_max_eta Then
               MsgTxt = (Salts(i).name + " molality " & molalities(i) & " out of limits {0..." & Salts(i).mola_max_eta & "} in (Viscosities.dynamicViscosity_Duan_pTX)")
               Debug.Print MsgTxt
               If outOfRangeMode = 2 Then
                   dynamicViscosity = MsgTxt
                   Exit Function
               End If
           End If
         End If
    
    Dim c As Double 'Molarity_molperliter
    Dim b As Double 'Molality
    Dim phi As Double:
    Dim A_ As Double
    Dim B_ As Double
    Dim C_ As Double
    Dim rho: rho = density(p_Pa, T, X)
    If VarType(rho) = vbString Then
       dynamicViscosity = rho
       Exit Function
    End If
    
    Dim eta_relative As Double
    
    'MIXING WEIGHT
         phi = molalities(i) / (Application.Sum(molalities) - molalities(UBound(molalities))) 'geometric mean mixture rule weighted with mass fraction (as in Laliberté)
           If i = 3 Then
           'Zhang (available for NaCl, KCl and CaCl)
          c = X(i) / Salts(i).MM * CDbl(rho) / 1000 / phi 'component molarity (molperliter)
            eta_relative = 1 + Salts(i).Zh_A * c ^ 0.5 + Salts(i).Zh_B * c + Salts(i).Zh_D * c ^ 2 + 0.0001 * Salts(i).Zh_E * c ^ 3.5 + 0.00001 * Salts(i).Zh_F * c ^ 7
         Else
           'Duan (available for NaCl and KCl)
           b = molalities(i) / phi
           A_ = Salts(i).a(1) + Salts(i).a(2) * T + Salts(i).a(3) * T ^ 2
           B_ = Salts(i).b(1) + Salts(i).b(2) * T + Salts(i).b(3) * T ^ 2
           C_ = Salts(i).c(1) + Salts(i).c(2) * T
           eta_relative = Exp(A_ * b + B_ * b ^ 2 + C_ * b ^ 3) 'Mixture is composed of binary solutions of the same molality
         End If
         eta = eta * eta_relative ^ phi
    
       End If
     Next i
     dynamicViscosity = eta
End Function



Function density(p As Double, T As Double, Xin) ', p_sat_MPa As Double) 'Brine density'
' density calculation of an aqueous salt solution according to Shide Mao and Zhenhao Duan (2008) 0-300degC; 0.1-100MPa; 0-6 mol/kg"
'  Mixing rule acc. to Laliberte&Cooper2004
' http://dx.doi.org/10.1016/j.jct.2008.03.005
' http://www.geochem-model.org/wp-content/uploads/2009/09/55-JCT_40_1046.pdf
' Problems: Brine has the same evaporation temperature as pure water"
    If DebugMode Then
        Debug.Print "Running Density_vec(" & p / 10 ^ 5 & " bar," & T - 273.15 & " °C" ', X={" & Xin(1) & Xin(2) & Xin(3); "})"
    End If
  
    If Salts(1).MM = 0 Then 'if Salt data not defined
        init
    End If
  
    If p = 0 Then
      density = "#p is 0"
      Exit Function
    End If
    
    If T = 0 Then
      density = "#T is 0"
      Exit Function
    End If
    
    Dim X
    X = CheckMassVector(Xin, nX)
    If VarType(X) = vbString Then
        density = X & " (Brine_liq.density)"
        Exit Function
    End If
    
    Dim M: M = massFractionsToMolalities(X, MM_vec) 'Molality
'    If M = -1 Then
'      density = "#Error in massFractionsToMolalities (no water)"
'      Exit Function
'    End If
    
    Dim v_ As Double
    Dim b As Double: b = 1.2
    Dim U(1 To 9) As Double 'dielectric constant D of pure water according to Bradley and Pitzer (1979)
    U(1) = 342.79
    U(2) = -0.0050866
    U(3) = 0.0000009469
    U(4) = -2.0525
    U(5) = 3115.9
    U(6) = -182.89
    U(7) = -8032.5
    U(8) = 4214200#
    U(9) = 2.1417
              
    Dim N_0 As Double: N_0 = 6.0221415E+23    'Avogadro constant in [1/mol]
    Dim E As Double: E = 1.60217733E-19 * 10 * 299792458 'elementary charge in [esu]
    Dim k As Double: k = 1.3806505E-16 'Boltzmann constant in [erg/K]
    Dim R As Double: R = 8.314472 'Gas constant [J/mol*K]
    Dim p_bar As Double: p_bar = p / 10 ^ 5
    Dim p_MPa As Double: p_MPa = p / 10 ^ 6
    Dim V As Double
    Dim I_ As Double
    Dim I_mr As Double
    Dim rho_sol_r As Double
    Dim rho_H2O As Double: rho_H2O = Density_pT(p, T) / 1000
    
    Dim rho_H2O_plus As Double
    Dim rho_H2O_minus As Double
    Dim p_plus_bar As Double
    Dim p_minus_bar As Double
    Dim D_plus As Double
    Dim D_minus As Double
    Dim A_Phi_plus As Double
    Dim A_Phi_minus As Double
    Dim B_v As Double
    Dim C_v As Double
    Dim V_m_r As Double
    Dim BB As Double
    Dim Cc As Double
    Dim D_1000 As Double
    Dim d As Double
    
    Dim A_Phi As Double
    Dim dp As Double
    Dim dA_Phi As Double
    Dim A_v As Double
    Dim V_o_Phi As Double
    Dim V_Phi() As Double
    ReDim V_Phi(1 To nX_salt)
    Dim h As Double
    Dim h_mr As Double
    
    Dim M_salt() As Double
    ReDim M_salt(1 To nX_salt)
    
    Dim m_r As Double
    Dim z_plus As Double
    Dim z_minus As Double
    Dim v_plus As Double
    Dim v_minus As Double
    
    Dim c() As Double
    ReDim c(1 To 23)
    
    
    If Application.Max(X) <= 0.1 ^ 12 Then 'for pure water skip the whole calculation and return water density
      density = rho_H2O * 1000
      Exit Function
    End If
    
    Dim j As Integer
    For j = 1 To nX_salt
'        If Len(Salts(j).name) = 0 Then
'            Density = "#Salt name not defined."
'            Exit Function
'        End If
        If Not X(j) > 0 Then
            M_salt(j) = 1 'To avoid division by zero at the end
        Else
            If outOfRangeMode > 0 Then
              Dim MsgTxt As String
                If Not (M(j) >= 0 And M(j) <= Salts(j).mola_max_rho) Then
                  MsgTxt = Salts(j).name & " molality " & M(j) & " out of limits {0..." & Salts(j).mola_max_rho & "} mol/kg (Density)"
                  Debug.Print MsgTxt
                  If outOfRangeMode = 1 Then
                    Debug.Print MsgTxt
                  ElseIf outOfRangeMode = 2 Then
                    density = MsgTxt
                    Exit Function
                  End If
                End If
                If Not (ignoreLimitSalt_p(j) Or (p >= Salts(j).p_min_rho And p <= Salts(j).p_max_rho)) Then
                  MsgTxt = "#p=" & p_bar & " bar, but for " & Salts(j).name & " must be within " & Salts(j).p_min_rho * 0.00001 & "..." & Salts(j).p_max_rho * 0.00001 & " bar (Density)"
                  Debug.Print MsgTxt
                  If outOfRangeMode = 1 Then
                    Debug.Print MsgTxt
                  ElseIf outOfRangeMode = 2 Then
                    density = MsgTxt
                    Exit Function
                  End If
                End If
                If Not (ignoreLimitSalt_T(j) Or (T >= Salts(j).T_min_rho And T <= Salts(j).T_max_rho)) Then
                  MsgTxt = "#T=" & T - 273.15 & "°C but for " & Salts(j).name & " must be within " & Salts(j).T_min_rho - 273.15 & "..." & Salts(j).T_max_rho - 273.15 & "°C (Density)"
                  If outOfRangeMode = 1 Then
                    Debug.Print MsgTxt
                  ElseIf outOfRangeMode = 2 Then
                    density = MsgTxt
                    Exit Function
                  End If
                End If
            End If
            
            M_salt(j) = Salts(j).MM * 1000 'in g/mol
            m_r = Salts(j).m_r
            z_plus = Salts(j).z_plus
            z_minus = Salts(j).z_minus
            v_plus = Salts(j).v_plus
            v_minus = Salts(j).v_minus
            c = Salts(j).Coeff
            
            V = v_plus + v_minus
            
            '---------------------------------------------------
            
            'Equation 3: Ionic strength
            I_ = 1 / 2 * (M(j) * v_plus * z_plus ^ 2 + M(j) * v_minus * z_minus ^ 2)
            I_mr = 1 / 2 * (m_r * v_plus * z_plus ^ 2 + m_r * v_minus * z_minus ^ 2)
            
            'Equation 4:
            h = Application.Log10(1 + b * I_ ^ (0.5)) / (2 * b)
            h_mr = Application.Log10(1 + b * I_mr ^ (0.5)) / (2 * b)
            
            '---------------------------------------------------
            ' equations using empirically fitted coefficients
            '---------------------------------------------------
            
            'Equation 10: solution volume at reference molality
            V_m_r = c(1) + c(2) * T + c(3) * T ^ 2 + c(4) * T ^ 3 + p_MPa * (c(5) + c(6) * T + c(7) * T ^ 2 + c(8) * T ^ 3)
            
            'Check: solution density at reference molality
            rho_sol_r = (1000 + m_r * M_salt(j)) / V_m_r
            
            'Equation 11: second virial coefficient. depends on temperature and pressure
            B_v = c(9) / (T - 227) + c(10) + c(11) * T + c(12) * T ^ 2 + c(13) / (647 - T) + p_MPa * (c(14) / (T - 227) + c(15) + c(16) * T + c(17) * T ^ 2 + c(18) / (647 - T))
            
            'Equation 12: third virial coefficient. depends on temperature
            C_v = c(19) / (T - 227) + c(20) + c(21) * T + c(22) * T ^ 2 + c(23) / (647 - T)
            
            '---------------------------------------------------
            ' Appendix A: Debye-Hückel limiting law slopes'
            '---------------------------------------------------
            
            BB = U(7) + U(8) / T + U(9) * T
            Cc = U(4) + U(5) / (U(6) + T)
            D_1000 = U(1) * Exp(U(2) * T + U(3) * T ^ 2)
            d = D_1000 + Cc * Log((BB + p_bar) / (BB + 1000))
            
            'DH-slope for osmotic coefficient according to Bradley and Pitzer (1979)
            A_Phi = 1 / 3 * ((2 * Application.Pi() * N_0 * rho_H2O) / 1000) ^ (1 / 2) * (E ^ 2 / (d * k * T)) ^ (3 / 2)
            
            'numeric differentiation per dp
            dp = 0.001 * p_bar
            p_plus_bar = p_bar + dp '/2
            p_minus_bar = p_bar '- dp/2
            D_plus = D_1000 + Cc * Log((BB + p_plus_bar) / (BB + 1000))
            D_minus = D_1000 + Cc * Log((BB + p_minus_bar) / (BB + 1000))
            
            rho_H2O_plus = Density_pT(p_plus_bar * 10 ^ 5, T) / 1000 'kg/m³->kg/dm³
            
            rho_H2O_minus = rho_H2O
              'Modelica.Media.Water.WaterIF97_base.density_pT(p_minus_bar*1e5, T) * 1e-3 kg/m³->kg/dm³
            A_Phi_plus = 1 / 3 * (2 * Application.Pi() * N_0 * rho_H2O_plus / 1000) ^ (1 / 2) * (E ^ 2 / (D_plus * k * T)) ^ (3 / 2)
            A_Phi_minus = 1 / 3 * (2 * Application.Pi() * N_0 * rho_H2O_minus / 1000) ^ (1 / 2) * (E ^ 2 / (D_minus * k * T)) ^ (3 / 2)
            dA_Phi = (A_Phi_plus - A_Phi_minus)
            
            'DH-slope for apparent molar volume according to Rogers and Pitzer (1982)
            A_v = 23 * (-4 * R * T * dA_Phi / dp) 'where does the 23 come from??
            
            '---------------------------------------------------
            ' Solution 1: using V_o_Phi and V_Phi
            '---------------------------------------------------
            
            'Equation 13: apparent molar Volume at infinite dilution in cm³/mol
            V_o_Phi = (V_m_r / m_r - 1000 / (m_r * rho_H2O) - V * Abs(z_plus * z_minus) * A_v * h_mr - 2 * v_plus * v_minus * R * T * (B_v * m_r + v_plus * z_plus * C_v * m_r ^ 2))
            
                          'Equation 2: apparent molar Volume in cm^3/mol
            V_Phi(j) = V_o_Phi + V * Abs(z_plus * z_minus) * A_v * h + 2 * v_plus * v_minus * M(j) * R * T * (B_v + v_plus * z_plus * M(j) * C_v)
            
        End If
  Next j
  
        v_ = X(UBound(X)) / (rho_H2O * 1000)
    For j = 1 To nX - 1
        v_ = v_ + X(j) * (V_Phi(j) / 10 ^ 6 / (M_salt(j) / 1000)) 'Mixing rule Laliberte&Cooper2004 equ. 5&6
    Next j
    density = 1 / v_
End Function



Function resistivity(p As Double, T As Double, Xin) 'electrical density'
    Dim d: d = density(p, T, Xin)
    If VarType(d) = vbString Then
        resistivity = d & " (Brine_liq.resistivity)"
        Exit Function
    End If
    
    Dim X
    X = CheckMassVector(Xin, nX) ' TODO get from density, checked there already
    If Not X(UBound(X)) < 1 Then
        resistivity = "#infinite resistivity for pure water"
        Exit Function
    Else
        Dim i As Integer
        Dim gamma_vec As Variant
        Dim n_vec As Variant: n_vec = Array(1, 1, 2) 'ion valence TODO get from elsewhere
        Dim c_vec As Variant: c_vec = Array(0, 0, 0)
        For i = 1 To 3
            c_vec(i) = X(i) / MM_vec(i) * CDbl(d) / 1000
        Next i

        Dim c_total As Double: c_total = c_vec(1) + c_vec(2) + c_vec(3) ' WorksheetFunction.Sum(SubArray(X,1,3))
        Dim j As Integer
        Dim X_total As Double
        For j = 1 To nX_salt
            If outOfRangeMode > 0 Then
                X_total = c_total * MM_vec(j)
                Dim MsgTxt As String
                If Not (X_total = 0 Or (X_total >= 0.03 And X_total <= 0.26)) Then
                  MsgTxt = "Total salinity X=" & X_total & " out of " & Salts(j).name & " limits {0.03...0.26} kg/kg (resistivity)"
                  Debug.Print MsgTxt
                  If outOfRangeMode = 1 Then
                    Debug.Print MsgTxt
                  ElseIf outOfRangeMode = 2 Then
                    resistivity = MsgTxt
                    Exit Function
                  End If
                End If
                If Not (ignoreLimitSalt_T(j) Or (T - 273.15 >= 22 And T - 273.15 <= 375)) Then
                  MsgTxt = "#T=" & T - 273.15 & "°C but for resistivity must be within 22...375 °C"
                  If outOfRangeMode = 1 Then
                    Debug.Print MsgTxt
                  ElseIf outOfRangeMode = 2 Then
                    resistivity = MsgTxt
                    Exit Function
                  End If
                End If
            End If
        Next j
    End If
        
    gamma_vec = Conductivity_Ucok1980_Tcd(T, fill(c_total, nX_salt), CDbl(d))
    
    Dim gamma As Double
    With WorksheetFunction
        gamma = .SumProduct(c_vec, n_vec, gamma_vec) / .SumProduct(c_vec, n_vec)
    End With
    resistivity = 1 / gamma
End Function

Function Conductivity_Ucok1980_Tcd(T As Double, c_vec, d As Double) As Variant 'electrical density'
  Dim B_NaCl As Variant: B_NaCl = Array( _
    Array(3.47, -59.21, 0.4551, -0.0000935, -0.00000177), _
    Array(-6.65, 198.1, -0.2058, 0.0000737, 0.000000877), _
    Array(2.633, -64.8, 0.005799, 0.0000674, -0.000000214))
  Dim B_KCl As Variant: B_KCl = Array( _
    Array(5.783, -59.23, 0.2051, 0.000182, -0.00000109), _
    Array(-6.607, 149.7, 0.1064, -0.000704, 0.00000108), _
    Array(1.665, -31.21, -0.03418, 0.000154, -0.000000195))
  Dim B_CaCl2 As Variant: B_CaCl2 = Array( _
    Array(-34.62, 780.3, 1.05, -0.002459, 0.000000999), _
    Array(24.64, -492.3, -0.5922, 0.001461, -0.000000711), _
    Array(-3.907, 64.59, 0.06735, -0.000122, -0.00000000473))
  Dim BB As Variant: BB = Array(B_NaCl, B_KCl, B_CaCl2)
    
  Dim T_C As Double: T_C = T - 273.15
  Dim T_vec As Variant: T_vec = Array(1, 1 / T_C, T_C, T_C ^ 2, T_C ^ 3)
  Dim c As Double
  Dim C_ As Variant
  Dim D_ As Variant
  Dim i As Integer
  Dim gamma(3) As Double
  Dim offset As Integer: offset = LBound(c_vec) - 1
  
    For i = 1 To 3
    c = c_vec(i + offset)
    If Not c > 0 Then
      gamma(i) = 0 ' gamma
    Else
      C_ = Array(c, c * Sqr(c), c * c * Log(c))
      D_ = WorksheetFunction.MMult(C_, BB(i))
      gamma(i) = WorksheetFunction.SumProduct(D_, T_vec)  ' gamma
    End If
    Next i
    Conductivity_Ucok1980_Tcd = gamma
End Function
