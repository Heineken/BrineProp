Attribute VB_Name = "SaltData"
Option Explicit
Option Base 1

Type SaltProps
    name As String
    MM As Double 'molar mass in [kg/mol]
    nM As Integer 'ions per salt molecule

    Coeff(1 To 23) As Double '
    m_r As Double 'reference molality
    z_plus As Double 'cation charge +
    z_minus As Double 'anion charge -
    v_plus As Double 'moles of cations per mol salt
    v_minus As Double 'moles of anions per mol salt
    mola_max_rho As Double 'maximum molality for which density function is valid
    T_min_rho As Double
    T_max_rho As Double
    p_min_rho As Double
    p_max_rho As Double
 'VISCOSITY
     a(1 To 3)  As Double '
     b(1 To 3)  As Double '
     c(1 To 2)  As Double '
     Zh_A As Double 'Zhang viscosity coefficients
     Zh_B As Double
     Zh_D As Double
     Zh_E As Double
     Zh_F As Double

     mola_max_eta As Double 'maximum molality for which viscosity function is valid
End Type
Public Const M_NaCl = 0.058443
Public Const M_KCl = 0.074551
Public Const M_CaCl2 = 0.110984

Function MM_NaCl() As Double
    MM_NaCl = M_NaCl
End Function
Function MM_KCl() As Double
    MM_KCl = M_KCl
End Function
Function MM_CaCl2() As Double
    MM_CaCl2 = M_CaCl2
End Function

Function NaCl() As SaltProps
    With NaCl
         .name = "NaCl"
         .MM = M_NaCl
         .nM = 2
         .m_r = 6
         .z_plus = 1
         .z_minus = -1
         .v_plus = 1
         .v_minus = 1
         'ReDim .Coeff(1 To 23)
         .Coeff(1) = 1066.07098
         .Coeff(2) = -0.00839622456
         .Coeff(3) = 0.000535429127
         .Coeff(4) = 0.000000755373789
         .Coeff(5) = -0.419512335
         .Coeff(6) = 0.00145082899
         .Coeff(7) = -0.00000347807732
         .Coeff(8) = 0#
         .Coeff(9) = 0.0110913788
         .Coeff(10) = 0.00114498252
         .Coeff(11) = -0.0000055118127
         .Coeff(12) = 7.05483955E-09
         .Coeff(13) = -0.0505734723
         .Coeff(14) = -0.000132747828
         .Coeff(15) = 0.00000477261581
         .Coeff(16) = -1.76888377E-08
         .Coeff(17) = 0#
         .Coeff(18) = 0.000640541237
         .Coeff(19) = 0.000307698827
         .Coeff(20) = -0.000164042763
         .Coeff(21) = 0.000000706784935
         .Coeff(22) = -6.50338372E-10
         .Coeff(23) = -0.000450906014
        .mola_max_rho = 6
        .T_min_rho = 273
        .T_max_rho = 573
        .p_min_rho = 100000#
        .p_max_rho = 100000000#
        .a(1) = -0.21319213
        .a(2) = 0.0013651589
        .a(3) = -0.0000012191756
        .b(1) = 0.069161945
        .b(2) = -0.00027292263
        .b(3) = 0.00000020852448
        .c(1) = -0.0025988855
        .c(2) = 0.0000077989227
        .Zh_A = 0.0061
        .Zh_B = 0.0799
        .Zh_D = 0.0104
        .Zh_E = 7.56
        .mola_max_eta = 6
    End With
End Function

Function KCl() As SaltProps
    With KCl
       .name = "KCl"
       .MM = M_KCl
       .nM = 2
       .m_r = 6
       .z_plus = 1
       .z_minus = -1
       .v_plus = 1
       .v_minus = 1
       .Coeff(1) = 290.812061
       .Coeff(2) = 6.54111195
       .Coeff(3) = -0.0161831978
       .Coeff(4) = 0.0000146290384
       .Coeff(5) = 14.1397987
       .Coeff(6) = -0.10726623
       .Coeff(7) = 0.000264506021
       .Coeff(8) = -0.000000219189708
       .Coeff(9) = 0.0302182158
       .Coeff(10) = -0.00215621394
       .Coeff(11) = 0.00000924163206
       .Coeff(12) = -1.10089434E-08
       .Coeff(13) = 0.0287018859
       .Coeff(14) = -0.000673119697
       .Coeff(15) = 0.000168332473
       .Coeff(16) = -0.00000079964564
       .Coeff(17) = 1.1188156E-09
       .Coeff(18) = -0.00659292385
       .Coeff(19) = -0.00202369103
       .Coeff(20) = -0.000170609099
       .Coeff(21) = 0.00000100510108
       .Coeff(22) = -1.86624642E-09
       .Coeff(23) = 0.0191919166
      .mola_max_rho = 4.5
      .T_min_rho = 273
      .T_max_rho = 543
      .p_min_rho = 100000#
      .p_max_rho = 50000000#
        .a(1) = -0.42122934
        .a(2) = 0.0018286059
        .a(3) = -0.0000013603098
        .b(1) = 0.011380205
        .b(2) = 0.0000047541391
        .b(3) = -0.000000099280575
        .c(1) = 0
        .c(2) = 0
    .Zh_A = 0.0051
    .Zh_B = -0.0152
    .Zh_D = 0.00725
    .Zh_E = 0.8
    .mola_max_eta = 4.5
  End With
End Function

Function CaCl2() As SaltProps
  With CaCl2
     .name = "CaCl2"
     .MM = M_CaCl2
     .nM = 3
     .m_r = 5
     .z_plus = 2
     .z_minus = -1
     .v_plus = 1
     .v_minus = 2
    .Coeff(1) = 1120.80057
    .Coeff(2) = -0.261669538
    .Coeff(3) = 0.0015204296
    .Coeff(4) = -0.000000689131095
    .Coeff(5) = -0.511802652
    .Coeff(6) = 0.00222234857
    .Coeff(7) = -0.00000566464544
    .Coeff(8) = 2.92950266E-09
    .Coeff(9) = 0.0243934633
    .Coeff(10) = -0.00142746873
    .Coeff(11) = 0.00000735840529
    .Coeff(12) = -9.4361548E-09
    .Coeff(13) = -0.0518606814
    .Coeff(14) = -0.0000616536928
    .Coeff(15) = -0.0000104523561
    .Coeff(16) = 4.52637296E-08
    .Coeff(17) = -1.05076158E-10
    .Coeff(18) = 0.00231544709
    .Coeff(19) = -0.00109663211
    .Coeff(20) = 0.000190836111
    .Coeff(21) = -0.000000925997994
    .Coeff(22) = 1.54388261E-09
    .Coeff(23) = -0.0129354832
    .mola_max_rho = 6
    .T_min_rho = 273
    .T_max_rho = 523
    .p_min_rho = 100000#
    .p_max_rho = 60000000#
    .Zh_A = 0.0157
    .Zh_B = 0.271
    .Zh_D = 0.04712
    .Zh_E = 94
    .Zh_F = 3
    .mola_max_eta = 3
    End With
End Function

' Public NaCl As SaltProps, CaCl2 As SaltProps, KCL As SaltProps

'Sub DefineSalts()
''Dim SaltData As New Collection
''Set SaltData = New SaltProps
'    'With SaltData
'    With NaCl
'         .name = "NaCl"
'         .MM = 0.058443
'         .m_r = 6
'         .z_plus = 1
'         .z_minus = -1
'         .v_plus = 1
'         .v_minus = 1
'         'ReDim .Coeff(1 To 23)
'         .Coeff(1) = 1066.07098
'         .Coeff(2) = -0.00839622456
'         .Coeff(3) = 0.000535429127
'         .Coeff(4) = 0.000000755373789
'         .Coeff(5) = -0.419512335
'         .Coeff(6) = 0.00145082899
'         .Coeff(7) = -0.00000347807732
'         .Coeff(8) = 0#
'         .Coeff(9) = 0.0110913788
'         .Coeff(10) = 0.00114498252
'         .Coeff(11) = -0.0000055118127
'         .Coeff(12) = 7.05483955E-09
'         .Coeff(13) = -0.0505734723
'         .Coeff(14) = -0.000132747828
'         .Coeff(15) = 0.00000477261581
'         .Coeff(16) = -1.76888377E-08
'         .Coeff(17) = 0#
'         .Coeff(18) = 0.000640541237
'         .Coeff(19) = 0.000307698827
'         .Coeff(20) = -0.000164042763
'         .Coeff(21) = 0.000000706784935
'         .Coeff(22) = -6.50338372E-10
'         .Coeff(23) = -0.000450906014
'        .mola_max_rho = 6
'        .T_min_rho = 273
'        .T_max_rho = 573
'        .p_min_rho = 100000#
'        .p_max_rho = 100000000#
'        .a(1) = -0.21319213
'        .a(2) = 0.0013651589
'        .a(3) = -0.0000012191756
'        .b(1) = 0.069161945
'        .b(2) = -0.00027292263
'        .b(3) = 0.00000020852448
'        .c(1) = -0.0025988855
'        .c(2) = 0.0000077989227
'        .Zh_A = 0.0061
'        .Zh_B = 0.0799
'        .Zh_D = 0.0104
'        .Zh_E = 7.56
'        .mola_max_eta = 6
'    End With
'
'    With KCl
'       .name = "KCl"
'       .MM = 0.074551
'       .m_r = 6
'       .z_plus = 1
'       .z_minus = -1
'       .v_plus = 1
'       .v_minus = 1
'       .Coeff(1) = 290.812061
'       .Coeff(2) = 6.54111195
'       .Coeff(3) = -0.0161831978
'       .Coeff(4) = 0.0000146290384
'       .Coeff(5) = 14.1397987
'       .Coeff(6) = -0.10726623
'       .Coeff(7) = 0.000264506021
'       .Coeff(8) = -0.000000219189708
'       .Coeff(9) = 0.0302182158
'       .Coeff(10) = -0.00215621394
'       .Coeff(11) = 0.00000924163206
'       .Coeff(12) = -1.10089434E-08
'       .Coeff(13) = 0.0287018859
'       .Coeff(14) = -0.000673119697
'       .Coeff(15) = 0.000168332473
'       .Coeff(16) = -0.00000079964564
'       .Coeff(17) = 1.1188156E-09
'       .Coeff(18) = -0.00659292385
'       .Coeff(19) = -0.00202369103
'       .Coeff(20) = -0.000170609099
'       .Coeff(21) = 0.00000100510108
'       .Coeff(22) = -1.86624642E-09
'       .Coeff(23) = 0.0191919166
'      .mola_max_rho = 4.5
'      .T_min_rho = 273
'      .T_max_rho = 543
'      .p_min_rho = 100000#
'      .p_max_rho = 50000000#
'        .a(1) = -0.42122934
'        .a(2) = 0.0018286059
'        .a(3) = -0.0000013603098
'        .b(1) = 0.011380205
'        .b(2) = 0.0000047541391
'        .b(3) = -0.000000099280575
'        .c(1) = 0
'        .c(2) = 0
'    .Zh_A = 0.0051
'    .Zh_B = -0.0152
'    .Zh_D = 0.00725
'    .Zh_E = 0.8
'    .mola_max_eta = 4.5
'  End With
'
'  With CaCl2
'     .name = "CaCl2"
'     .MM = 0.110984
'     .m_r = 5
'     .z_plus = 2
'     .z_minus = -1
'     .v_plus = 1
'     .v_minus = 2
'    .Coeff(1) = 1120.80057
'    .Coeff(2) = -0.261669538
'    .Coeff(3) = 0.0015204296
'    .Coeff(4) = -0.000000689131095
'    .Coeff(5) = -0.511802652
'    .Coeff(6) = 0.00222234857
'    .Coeff(7) = -0.00000566464544
'    .Coeff(8) = 2.92950266E-09
'    .Coeff(9) = 0.0243934633
'    .Coeff(10) = -0.00142746873
'    .Coeff(11) = 0.00000735840529
'    .Coeff(12) = -9.4361548E-09
'    .Coeff(13) = -0.0518606814
'    .Coeff(14) = -0.0000616536928
'    .Coeff(15) = -0.0000104523561
'    .Coeff(16) = 4.52637296E-08
'    .Coeff(17) = -1.05076158E-10
'    .Coeff(18) = 0.00231544709
'    .Coeff(19) = -0.00109663211
'    .Coeff(20) = 0.000190836111
'    .Coeff(21) = -0.000000925997994
'    .Coeff(22) = 1.54388261E-09
'    .Coeff(23) = -0.0129354832
'    .mola_max_rho = 6
'    .T_min_rho = 273
'    .T_max_rho = 523
'    .p_min_rho = 100000#
'    .p_max_rho = 60000000#
'    .Zh_A = 0.0157
'    .Zh_B = 0.271
'    .Zh_D = 0.04712
'    .Zh_E = 94
'    .Zh_F = 3
'    .mola_max_eta = 3
'    End With
'
'    'SaltData.Add NaCl
'    Salts(1) = NaCl
'    Salts(2) = KCl
'    Salts(3) = CaCl2
'
''    MM(1) = Salts(1).MM
''    MM(2) = Salts(2).MM
''    MM(3) = Salts(3).MM
''    MM(4) = M_H2O
'    If DebugMode Then
'        Debug.Print "Salt properties defined"
'    End If
'End Sub
