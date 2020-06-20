C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
C     A 15-step reduced mechanism based on GRI3.0
C
C     by Tianfeng Lu
C     Email: tlu@engr.uconn.edu
C
C     Reference: 
C       Tianfeng Lu and Chung K. Law, 
C       "A criterion based on computational singular perturbation
C       for the identification of quasi steady state species: 
C       A reduced mechanism for methane oxidation with NO chemistry," 
C       Combustion and Flame, Vol.154 No.4 pp.761�774, 2008.
C
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE CKWYP  (P, T, Y, ICKWRK, RCKWRK, WDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION ICKWRK(*), RCKWRK(*), WDOT(*), Y(*)
      DIMENSION RF(184), RB(184), RKLOW(22), XQ(11), C(19)
C
      CALL YTCP(P, T, Y, C)
      CALL RATT(T, RF, RB, RKLOW)
      CALL RATX(T, C, RF, RB, RKLOW)
      CALL QSSA(RF, RB, XQ)
      CALL RDOT(RF, RB, WDOT)
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE YTCP (P, T, Y, C)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION Y(*), C(*)
      DATA SMALL/1.D-50/
C
      C(1) = Y(1)*4.96046521D-1
      C(2) = Y(2)*9.92093043D-1
      C(3) = Y(3)*6.25023433D-2
      C(4) = Y(4)*3.12511716D-2
      C(5) = Y(5)*5.87980383D-2
      C(6) = Y(6)*5.55082499D-2
      C(7) = Y(7)*3.02968146D-2
      C(8) = Y(8)*2.93990192D-2
      C(9) = Y(9)*6.65112065D-2
      C(10) = Y(10)*6.23323639D-2
      C(11) = Y(11)*3.57008335D-2
      C(12) = Y(12)*2.27221341D-2
      C(13) = Y(13)*3.33039255D-2
      C(14) = Y(14)*3.12086189D-2
      C(15) = Y(15)*3.84050525D-2
      C(16) = Y(16)*3.56453112D-2
      C(17) = Y(17)*3.32556033D-2
      C(18) = Y(18)*2.37882046D-2
      C(19) = Y(19)*3.56972032D-2
C
      SUM = 0.0
      DO K = 1, 19
         SUM = SUM + C(K)
      ENDDO
      SUM = P/(SUM*T*8.314510D7)
C
      DO K = 1, 19
         C(K) = C(K) * SUM
      ENDDO
      END
C
C----------------------------------------------------------------------C
C
      SUBROUTINE RATT (T, RF, RB, RKLOW)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (RU=8.314510D7, PATM=1.01325D6, SMALL = 1.D-200)
      DIMENSION RF(*), RB(*), RKLOW(*), EQK(184), SMH(30), EG(30)
C
      ALOGT = LOG(T)
      TI = 1.0D0/T
      TI2 = TI*TI
C
      RF(1) = 1.2D17*TI
      RF(2) = 5.D17*TI
      RF(3) = EXP(1.05635949D1 +2.7D0*ALOGT -3.15013634D3*TI)
      RF(4) = 2.D13
      RF(5) = EXP(1.60803938D1 +2.D0*ALOGT -2.01286667D3*TI)
      RF(6) = 5.7D13
      RF(7) = 8.D13
      RF(8) = 1.5D13
      RF(9) = 1.5D13
      RF(10) = 5.06D13
      RF(11) = EXP(2.07430685D1 +1.5D0*ALOGT -4.32766334D3*TI)
      RF(12) = EXP(2.36136376D1 -1.20017175D3*TI)
      RF(13) = 3.D13
      RF(14) = 3.D13
      RF(15) = EXP(3.12945828D1 -1.781387D3*TI)
      RF(16) = 1.D13
      RF(17) = 1.D13
      RF(18) = EXP(1.28687606D1 +2.5D0*ALOGT -1.55997167D3*TI)
      RF(19) = EXP(1.17752897D1 +2.5D0*ALOGT -2.51608334D3*TI)
      TMP = EXP(2.D0*ALOGT -9.56111669D2*TI )
      RF(20) = 1.35D7 * TMP
      RF(21) = 6.94D6 * TMP
      RF(22) = 3.D13
      TMP = EXP(1.83D0*ALOGT -1.10707667D2*TI )
      RF(23) = 1.25D7 * TMP
      RF(167) = 6.7D6 * TMP
      RF(24) = 2.24D13
      RF(25) = EXP(1.83130955D1 +1.92D0*ALOGT -2.86330284D3*TI)
      RF(26) = 1.D14
      TMP = EXP(-4.02573334D3*TI)
      RF(27) = 1.D13 * TMP
      RF(76) = 5.D13 * TMP
      RF(28) = EXP(2.81906369D1 -6.79342501D2*TI)
      RF(29) = EXP(2.85473118D1 -2.40537567D4*TI)
      RF(30) = EXP(3.22361913D1 -2.01286667D4*TI)
      RF(31) = EXP(4.24761511D1 -8.6D-1*ALOGT)
      TMP = EXP(-1.24D0*ALOGT)
      RF(32) = 2.08D19 * TMP
      RF(34) = 2.6D19 * TMP
      RF(33) = EXP(4.38677883D1 -7.6D-1*ALOGT)
      RF(35) = EXP(3.78159211D1 -6.707D-1*ALOGT -8.57531523D3*TI)
      RF(36) = 1.D18*TI
      RF(37) = EXP(3.90385861D1 -6.D-1*ALOGT)
      RF(38) = EXP(4.55408762D1 -1.25D0*ALOGT)
      RF(39) = 5.5D20*TI2
      RF(40) = 2.2D22*TI2
      RF(41) = EXP(2.90097872D1 -3.37658384D2*TI)
      RF(42) = EXP(3.14332293D1 -5.37435401D2*TI)
      RF(43) = EXP(3.20618379D1 -3.19542584D2*TI)
      RF(44) = EXP(1.6308716D1 +2.D0*ALOGT -2.61672667D3*TI)
      RF(45) = EXP(2.99336062D1 -1.81158D3*TI)
      RF(46) = 1.65D14
      RF(47) = 6.D14
      RF(48) = 3.D13
      RF(49) = EXP(3.71706652D1 -5.34D-1*ALOGT -2.69724134D2*TI)
      RF(50) = EXP(2.03077504D1 +1.62D0*ALOGT -5.45486868D3*TI)
      RF(51) = EXP(2.77171988D1 +4.8D-1*ALOGT +1.30836334D2*TI)
      RF(52) = 7.34D13
      RF(53) = EXP(2.7014835D1 +4.54D-1*ALOGT -1.81158D3*TI)
      RF(54) = EXP(2.7014835D1 +4.54D-1*ALOGT -1.30836334D3*TI)
      RF(55) = EXP(1.78655549D1 +1.9D0*ALOGT -1.3798201D3*TI)
      RF(56) = EXP(2.76845619D1 +5.D-1*ALOGT -4.32766334D1*TI)
      RF(57) = 2.D13
      RF(58) = EXP(2.58292113D1 +6.5D-1*ALOGT +1.42913534D2*TI)
      RF(59) = EXP(3.11214496D1 -9.D-2*ALOGT -3.06962167D2*TI)
      RF(60) = EXP(2.85189124D1 +5.15D-1*ALOGT -2.51608334D1*TI)
      RF(61) = EXP(1.7541204D1 +1.63D0*ALOGT -9.68188869D2*TI)
      RF(62) = 2.D13
      RF(63) = EXP(2.80364862D1 +5.D-1*ALOGT +5.53538334D1*TI)
      RF(64) = EXP(3.31993656D1 -2.3D-1*ALOGT -5.38441834D2*TI)
      TMP = EXP(2.1D0*ALOGT -2.45066517D3*TI )
      RF(65) = 1.7D7 * TMP
      RF(66) = 4.2D6 * TMP
      RF(67) = EXP(2.93537877D1 -1.20772D3*TI)
      RF(68) = EXP(2.94360258D1 +2.7D-1*ALOGT -1.40900667D2*TI)
      RF(69) = 3.D13
      RF(70) = EXP(2.7014835D1 +4.54D-1*ALOGT -9.15854335D2*TI)
      RF(71) = EXP(1.4096923D1 +2.53D0*ALOGT -6.15937201D3*TI)
      RF(72) = EXP(4.07945264D1 -9.9D-1*ALOGT -7.95082335D2*TI)
      RF(73) = 2.D12
      RF(74) = EXP(1.85604427D1 +1.9D0*ALOGT -3.78922151D3*TI)
      RF(75) = 1.D14
      RF(77) = EXP(3.00558238D1 -1.72502674D3*TI)
      RF(78) = EXP(1.75767107D1 +1.5D0*ALOGT -4.00560467D4*TI)
      RF(79) = EXP(1.9190789D1 +1.51D0*ALOGT -1.72603317D3*TI)
      RF(80) = EXP(3.19350862D1 -3.7D-1*ALOGT)
      RF(81) = EXP(1.0482906D1 +2.4D0*ALOGT +1.06178717D3*TI)
      RF(82) = EXP(3.03051698D1 +2.51608334D2*TI)
      RF(83) = EXP(2.83241683D1 -2.14873517D2*TI)
      RF(84) = EXP(4.19771599D1 -1.47996022D4*TI)
      RF(85) = 5.D13
      RF(86) = 3.D13
      RF(87) = 2.D13
      RF(88) = EXP(1.62403133D1 +2.D0*ALOGT -1.50965D3*TI)
      RF(89) = 3.D13
      RF(90) = EXP(4.24725733D1 -1.43D0*ALOGT -6.69278168D2*TI)
      RF(91) = EXP(1.78408622D1 +1.6D0*ALOGT -2.72743434D3*TI)
      RF(92) = EXP(4.10064751D1 -1.34D0*ALOGT -7.13058018D2*TI)
      RF(93) = EXP(1.84206807D1 +1.6D0*ALOGT -1.570036D3*TI)
      RF(94) = EXP(1.76783433D1 +1.228D0*ALOGT -3.52251667D1*TI)
      RF(95) = 5.D13
      RF(96) = EXP(2.19558261D1 +1.18D0*ALOGT +2.2493785D2*TI)
      RF(97) = 5.D12
      RF(98) = 5.D12
      RF(99) = EXP(1.41801537D1 +2.D0*ALOGT +4.22702001D2*TI)
      RF(100) = EXP(1.56560602D1 +2.D0*ALOGT -7.54825001D2*TI)
      RF(101) = EXP(-8.4310155D0 +4.5D0*ALOGT +5.03216668D2*TI)
      RF(102) = EXP(-7.6354939D0 +4.D0*ALOGT +1.00643334D3*TI)
      RF(103) = 5.D12
      RF(104) = EXP(1.50964444D1 +2.D0*ALOGT -1.25804167D3*TI)
      RF(105) = EXP(1.50796373D1 +2.12D0*ALOGT -4.37798501D2*TI)
      RF(106) = EXP(2.96459241D1 -1.00643334D3*TI)
      RF(107) = EXP(2.55908003D1 +8.20243168D2*TI)
      RF(108) = EXP(3.36712758D1 -6.03860001D3*TI)
      RF(109) = 2.D13
      RF(110) = 1.D12
      RF(111) = 3.78D13
      RF(112) = EXP(3.26416564D1 -1.18759134D4*TI)
      RF(113) = EXP(1.55382772D1 +2.D0*ALOGT -6.03860001D3*TI)
      RF(114) = EXP(3.16914641D1 -2.89852801D2*TI)
      RF(115) = 5.D13
      RF(116) = 6.71D13
      RF(117) = EXP(3.23131523D1 -1.56500384D3*TI)
      RF(118) = EXP(2.93732401D1 +3.79928584D2*TI)
      RF(119) = 4.D13
      RF(120) = 3.D13
      RF(121) = 6.D13
      RF(122) = 5.D13
      RF(123) = EXP(3.28780452D1 -7.94679762D3*TI)
      RF(124) = EXP(3.21806786D1 +2.59156584D2*TI)
      RF(125) = 5.D13
      TMP = EXP(-7.54825001D2*TI)
      RF(126) = 5.D12 * TMP
      RF(171) = 5.8D12 * TMP
      RF(172) = 2.4D12 * TMP
      RF(127) = EXP(1.31223634D1 +2.D0*ALOGT -3.63825651D3*TI)
      RF(128) = EXP(3.500878D1 -6.01041988D3*TI)
      RF(129) = 4.D13
      RF(130) = EXP(1.47156719D1 +2.D0*ALOGT -4.16160184D3*TI)
      RF(131) = EXP(2.74203001D1 +5.D-1*ALOGT -2.26950717D3*TI)
      RF(132) = 3.D13
      RF(133) = EXP(3.03390713D1 -3.01930001D2*TI)
      RF(134) = 2.8D13
      RF(135) = 1.2D13
      RF(136) = 7.D13
      RF(137) = EXP(4.07167205D1 -1.16D0*ALOGT -5.76183084D2*TI)
      RF(138) = 3.D13
      TMP = EXP(2.86833501D2*TI)
      RF(139) = 1.2D13 * TMP
      RF(140) = 1.6D13 * TMP
      RF(141) = 9.D12
      RF(142) = 7.D12
      RF(143) = 1.4D13
      RF(144) = EXP(3.13199006D1 +2.76769167D2*TI)
      RF(145) = EXP(3.12033668D1 -1.5338044D4*TI)
      RF(146) = EXP(2.84682686D1 -1.02228466D4*TI)
      RF(147) = EXP(1.01064284D1 +2.47D0*ALOGT -2.60666234D3*TI)
      RF(148) = EXP(3.87538626D1 -1.18D0*ALOGT -3.29103701D2*TI)
      RF(149) = EXP(2.95538088D1 +1.D-1*ALOGT -5.33409668D3*TI)
      RF(150) = 2.648D13
      RF(151) = EXP(8.10772006D0 +2.81D0*ALOGT -2.94884967D3*TI)
      TMP = EXP(1.5D0*ALOGT -5.00197368D3*TI )
      RF(152) = 3.D7 * TMP
      RF(153) = 1.D7 * TMP
      RF(154) = EXP(1.23327053D1 +2.D0*ALOGT -4.62959334D3*TI)
      RF(155) = EXP(1.56303353D1 +1.74D0*ALOGT -5.25861418D3*TI)
      TMP = EXP(-1.D0*ALOGT -8.55468335D3*TI )
      RF(156) = 1.5D18 * TMP
      RF(157) = 1.87D17 * TMP
      RF(158) = EXP(3.02300002D1 -2.01286667D2*TI)
      RF(159) = EXP(3.05213929D1 -4.52895001D2*TI)
      RF(160) = EXP(-2.84796532D1 +7.6D0*ALOGT +1.77635484D3*TI)
      RF(161) = EXP(3.83630605D1 -1.39D0*ALOGT -5.10764918D2*TI)
      RF(162) = EXP(2.97104627D1 +4.4D-1*ALOGT -4.36641103D4*TI)
      RF(163) = EXP(2.74566677D1 -1.94996459D3*TI)
      RF(164) = EXP(2.87941719D1 -4.29747034D2*TI)
      RF(165) = 1.D13
      RF(166) = 3.37D13
      RF(168) = EXP(3.61482143D1 -8.72074485D3*TI)
      RF(169) = EXP(2.28027074D1 +5.D-1*ALOGT +8.83145252D2*TI)
      RF(170) = EXP(2.83090547D1 +4.3D-1*ALOGT +1.86190167D2*TI)
      RF(173) = EXP(3.29293385D1 -5.52984796D3*TI)
      RF(174) = EXP(2.49457104D1 +2.5D-1*ALOGT +4.70507584D2*TI)
      RF(175) = EXP(2.64369986D1 +2.9D-1*ALOGT -5.53538334D0*TI)
      RF(176) = EXP(1.41059389D1 +1.61D0*ALOGT +1.932352D2*TI)
      RF(177) = EXP(2.69105027D1 +4.22D-1*ALOGT +8.83145252D2*TI)
      RF(178) = 1.5D14
      RF(179) = 1.81D10
      RF(180) = 2.35D10
      RF(181) = 2.2D13
      RF(182) = 1.1D13
      RF(183) = 1.2D13
      RF(184) = 3.01D13
C
      CALL RDSMH (T, SMH)
      DO N = 1, 29
          EG(N) = EXP(SMH(N))
      ENDDO
C
      PFAC = PATM / (RU*T)
      PFAC2 = PFAC*PFAC
      PFAC3 = PFAC2*PFAC
C
      EQK(1)=EG(4)/EG(3)/EG(3)/PFAC
      EQK(2)=EG(5)/EG(2)/EG(3)/PFAC
      EQK(3)=EG(2)*EG(5)/EG(1)/EG(3)
      EQK(4)=EG(4)*EG(5)/EG(3)/EG(7)
      EQK(5)=EG(5)*EG(7)/EG(3)/EG(8)
      EQK(6)=EG(2)*EG(15)/EG(3)/EG(10)
      EQK(7)=EG(2)*EG(17)/EG(3)/EG(11)
      EQK(8)=EG(1)*EG(15)/EG(3)/EG(12)
      EQK(9)=EG(2)*EG(17)/EG(3)/EG(12)
      EQK(10)=EG(2)*EG(18)/EG(3)/EG(13)
      EQK(11)=EG(5)*EG(13)/EG(3)/EG(14)
      EQK(12)=EG(16)/EG(3)/EG(15)/PFAC
      EQK(13)=EG(5)*EG(15)/EG(3)/EG(17)
      EQK(14)=EG(2)*EG(16)/EG(3)/EG(17)
      EQK(15)=EG(5)*EG(17)/EG(3)/EG(18)
      EQK(16)=EG(5)*EG(18)/EG(3)/EG(19)
      EQK(17)=EG(5)*EG(18)/EG(3)/EG(20)
      EQK(18)=EG(5)*EG(19)/EG(3)/EG(21)
      EQK(19)=EG(5)*EG(20)/EG(3)/EG(21)
      EQK(20)=EG(2)*EG(27)/EG(3)/EG(22)
      EQK(21)=EG(11)*EG(15)/EG(3)/EG(22)
      EQK(22)=EG(2)*EG(28)/EG(3)/EG(23)
      EQK(23)=EG(13)*EG(17)/EG(3)/EG(24)
      EQK(24)=EG(13)*EG(18)/EG(3)/EG(25)
      EQK(25)=EG(5)*EG(25)/EG(3)/EG(26)
      EQK(26)=EG(2)*EG(15)*EG(15)/EG(3)/EG(27)*PFAC
      EQK(27)=EG(5)*EG(27)/EG(3)/EG(28)
      EQK(28)=EG(11)*EG(16)/EG(3)/EG(28)
      EQK(29)=EG(3)*EG(16)/EG(4)/EG(15)
      EQK(30)=EG(7)*EG(17)/EG(4)/EG(18)
      EQK(31)=EG(7)/EG(2)/EG(4)/PFAC
      EQK(32)=EQK(31)
      EQK(33)=EQK(31)
      EQK(34)=EQK(31)
      EQK(35)=EG(3)*EG(5)/EG(2)/EG(4)
      EQK(36)=EG(1)/EG(2)/EG(2)/PFAC
      EQK(37)=EQK(36)
      EQK(38)=EQK(36)
      EQK(39)=EQK(36)
      EQK(40)=EG(6)/EG(2)/EG(5)/PFAC
      EQK(41)=EG(3)*EG(6)/EG(2)/EG(7)
      EQK(42)=EG(1)*EG(4)/EG(2)/EG(7)
      EQK(43)=EG(5)*EG(5)/EG(2)/EG(7)
      EQK(44)=EG(1)*EG(7)/EG(2)/EG(8)
      EQK(45)=EG(5)*EG(6)/EG(2)/EG(8)
      EQK(46)=EG(1)*EG(9)/EG(2)/EG(10)
      EQK(47)=EG(13)/EG(2)/EG(11)/PFAC
      EQK(48)=EG(1)*EG(10)/EG(2)/EG(12)
      EQK(49)=EG(14)/EG(2)/EG(13)/PFAC
      EQK(50)=EG(1)*EG(13)/EG(2)/EG(14)
      EQK(51)=EG(18)/EG(2)/EG(17)/PFAC
      EQK(52)=EG(1)*EG(15)/EG(2)/EG(17)
      EQK(53)=EG(19)/EG(2)/EG(18)/PFAC
      EQK(54)=EG(20)/EG(2)/EG(18)/PFAC
      EQK(55)=EG(1)*EG(17)/EG(2)/EG(18)
      EQK(56)=EG(21)/EG(2)/EG(19)/PFAC
      EQK(57)=EG(1)*EG(18)/EG(2)/EG(19)
      EQK(58)=EG(5)*EG(13)/EG(2)/EG(19)
      EQK(59)=EG(6)*EG(12)/EG(2)/EG(19)
      EQK(60)=EG(21)/EG(2)/EG(20)/PFAC
      EQK(61)=EG(19)/EG(20)
      EQK(62)=EG(1)*EG(18)/EG(2)/EG(20)
      EQK(63)=EG(5)*EG(13)/EG(2)/EG(20)
      EQK(64)=EG(6)*EG(12)/EG(2)/EG(20)
      EQK(65)=EG(1)*EG(19)/EG(2)/EG(21)
      EQK(66)=EG(1)*EG(20)/EG(2)/EG(21)
      EQK(67)=EG(23)/EG(2)/EG(22)/PFAC
      EQK(68)=EG(24)/EG(2)/EG(23)/PFAC
      EQK(69)=EG(1)*EG(22)/EG(2)/EG(23)
      EQK(70)=EG(25)/EG(2)/EG(24)/PFAC
      EQK(71)=EG(1)*EG(23)/EG(2)/EG(24)
      EQK(72)=EG(26)/EG(2)/EG(25)/PFAC
      EQK(73)=EG(1)*EG(24)/EG(2)/EG(25)
      EQK(74)=EG(1)*EG(25)/EG(2)/EG(26)
      EQK(75)=EG(12)*EG(15)/EG(2)/EG(27)
      EQK(76)=EG(1)*EG(27)/EG(2)/EG(28)
      EQK(77)=EG(13)*EG(15)/EG(2)/EG(28)
      EQK(78)=EG(18)/EG(1)/EG(15)/PFAC
      EQK(79)=EG(2)*EG(6)/EG(1)/EG(5)
      EQK(80)=EG(8)/EG(5)/EG(5)/PFAC
      EQK(81)=EG(3)*EG(6)/EG(5)/EG(5)
      EQK(82)=EG(4)*EG(6)/EG(5)/EG(7)
      EQK(168)=EQK(82)
      EQK(83)=EG(6)*EG(7)/EG(5)/EG(8)
      EQK(84)=EQK(83)
      EQK(85)=EG(2)*EG(15)/EG(5)/EG(9)
      EQK(86)=EG(2)*EG(17)/EG(5)/EG(10)
      EQK(87)=EG(2)*EG(18)/EG(5)/EG(11)
      EQK(88)=EG(6)*EG(10)/EG(5)/EG(11)
      EQK(89)=EG(2)*EG(18)/EG(5)/EG(12)
      EQK(90)=EG(21)/EG(5)/EG(13)/PFAC
      EQK(91)=EG(6)*EG(11)/EG(5)/EG(13)
      EQK(92)=EG(6)*EG(12)/EG(5)/EG(13)
      EQK(93)=EG(6)*EG(13)/EG(5)/EG(14)
      EQK(94)=EG(2)*EG(16)/EG(5)/EG(15)
      EQK(95)=EG(6)*EG(15)/EG(5)/EG(17)
      EQK(96)=EG(6)*EG(17)/EG(5)/EG(18)
      EQK(97)=EG(6)*EG(18)/EG(5)/EG(19)
      EQK(98)=EG(6)*EG(18)/EG(5)/EG(20)
      EQK(99)=EG(6)*EG(19)/EG(5)/EG(21)
      EQK(100)=EG(6)*EG(20)/EG(5)/EG(21)
      EQK(101)=EG(2)*EG(28)/EG(5)/EG(22)
      EQK(102)=EG(13)*EG(15)/EG(5)/EG(22)
      EQK(103)=EG(6)*EG(22)/EG(5)/EG(23)
      EQK(104)=EG(6)*EG(23)/EG(5)/EG(24)
      EQK(105)=EG(6)*EG(25)/EG(5)/EG(26)
      EQK(106)=EG(6)*EG(27)/EG(5)/EG(28)
      EQK(107)=EG(4)*EG(8)/EG(7)/EG(7)
      EQK(108)=EQK(107)
      EQK(109)=EG(5)*EG(18)/EG(7)/EG(11)
      EQK(110)=EG(4)*EG(14)/EG(7)/EG(13)
      EQK(111)=EG(5)*EG(20)/EG(7)/EG(13)
      EQK(112)=EG(5)*EG(16)/EG(7)/EG(15)
      EQK(113)=EG(8)*EG(17)/EG(7)/EG(18)
      EQK(114)=EG(3)*EG(15)/EG(4)/EG(9)
      EQK(115)=EG(2)*EG(22)/EG(9)/EG(13)
      EQK(116)=EG(3)*EG(17)/EG(4)/EG(10)
      EQK(117)=EG(2)*EG(11)/EG(1)/EG(10)
      EQK(118)=EG(2)*EG(18)/EG(6)/EG(10)
      EQK(119)=EG(2)*EG(22)/EG(10)/EG(11)
      EQK(120)=EG(2)*EG(23)/EG(10)/EG(13)
      EQK(121)=EG(2)*EG(24)/EG(10)/EG(14)
      EQK(122)=EG(27)/EG(10)/EG(15)/PFAC
      EQK(123)=EG(15)*EG(17)/EG(10)/EG(16)
      EQK(124)=EG(2)*EG(28)/EG(10)/EG(18)
      EQK(125)=EG(15)*EG(22)/EG(10)/EG(27)
      EQK(127)=EG(2)*EG(13)/EG(1)/EG(11)
      EQK(128)=EG(1)*EG(22)/EG(11)/EG(11)
      EQK(129)=EG(2)*EG(24)/EG(11)/EG(13)
      EQK(130)=EG(13)*EG(13)/EG(11)/EG(14)
      EQK(131)=EG(28)/EG(11)/EG(15)/PFAC
      EQK(132)=EG(15)*EG(23)/EG(11)/EG(27)
      EQK(133)=EG(11)/EG(12)
      EQK(138)=EQK(133)
      EQK(141)=EQK(133)
      EQK(142)=EQK(133)
      EQK(134)=EG(2)*EG(5)*EG(15)/EG(4)/EG(12)*PFAC
      EQK(135)=EG(6)*EG(15)/EG(4)/EG(12)
      EQK(136)=EG(2)*EG(13)/EG(1)/EG(12)
      EQK(137)=EG(21)/EG(6)/EG(12)/PFAC
      EQK(139)=EG(2)*EG(24)/EG(12)/EG(13)
      EQK(140)=EG(13)*EG(13)/EG(12)/EG(14)
      EQK(143)=EG(15)*EG(18)/EG(12)/EG(16)
      EQK(144)=EG(13)*EG(25)/EG(12)/EG(26)
      EQK(145)=EG(3)*EG(20)/EG(4)/EG(13)
      EQK(146)=EG(5)*EG(18)/EG(4)/EG(13)
      EQK(147)=EG(7)*EG(14)/EG(8)/EG(13)
      EQK(148)=EG(26)/EG(13)/EG(13)/PFAC
      EQK(149)=EG(2)*EG(25)/EG(13)/EG(13)
      EQK(150)=EG(14)*EG(15)/EG(13)/EG(17)
      EQK(151)=EG(14)*EG(17)/EG(13)/EG(18)
      EQK(152)=EG(14)*EG(19)/EG(13)/EG(21)
      EQK(153)=EG(14)*EG(20)/EG(13)/EG(21)
      EQK(154)=EG(14)*EG(23)/EG(13)/EG(24)
      EQK(155)=EG(14)*EG(25)/EG(13)/EG(26)
      EQK(156)=EG(2)*EG(15)/EG(17)*PFAC
      EQK(157)=EQK(156)
      EQK(158)=EG(7)*EG(15)/EG(4)/EG(17)
      EQK(159)=EG(7)*EG(18)/EG(4)/EG(19)
      EQK(160)=EG(7)*EG(18)/EG(4)/EG(20)
      EQK(161)=EG(17)*EG(18)/EG(4)/EG(23)
      EQK(162)=EG(1)*EG(22)/EG(24)*PFAC
      EQK(163)=EG(7)*EG(24)/EG(4)/EG(25)
      EQK(164)=EG(5)*EG(15)*EG(15)/EG(4)/EG(27)*PFAC
      EQK(165)=EG(15)*EG(15)*EG(22)/EG(27)/EG(27)*PFAC
      EQK(167)=EG(2)*EG(29)/EG(3)/EG(24)
      EQK(170)=EG(13)/EG(1)/EG(10)/PFAC
      EQK(172)=EG(3)*EG(18)/EG(4)/EG(11)
      EQK(175)=EG(3)*EG(29)/EG(4)/EG(23)
      EQK(176)=EG(7)*EG(22)/EG(4)/EG(23)
      EQK(177)=EG(29)/EG(2)/EG(28)/PFAC
      EQK(181)=EG(13)*EG(17)/EG(2)/EG(29)
      EQK(182)=EG(1)*EG(28)/EG(2)/EG(29)
      EQK(183)=EG(6)*EG(28)/EG(5)/EG(29)
      EQK(184)=EG(17)*EG(19)/EG(5)/EG(29)
C
      RB(1) = RF(1) / MAX(EQK(1),SMALL)
      RB(2) = RF(2) / MAX(EQK(2),SMALL)
      RB(3) = RF(3) / MAX(EQK(3),SMALL)
      RB(4) = RF(4) / MAX(EQK(4),SMALL)
      RB(5) = RF(5) / MAX(EQK(5),SMALL)
      RB(6) = RF(6) / MAX(EQK(6),SMALL)
      RB(7) = RF(7) / MAX(EQK(7),SMALL)
      RB(8) = RF(8) / MAX(EQK(8),SMALL)
      RB(9) = RF(9) / MAX(EQK(9),SMALL)
      RB(10) = RF(10) / MAX(EQK(10),SMALL)
      RB(11) = RF(11) / MAX(EQK(11),SMALL)
      RB(12) = RF(12) / MAX(EQK(12),SMALL)
      RB(13) = RF(13) / MAX(EQK(13),SMALL)
      RB(14) = RF(14) / MAX(EQK(14),SMALL)
      RB(15) = RF(15) / MAX(EQK(15),SMALL)
      RB(16) = RF(16) / MAX(EQK(16),SMALL)
      RB(17) = RF(17) / MAX(EQK(17),SMALL)
      RB(18) = RF(18) / MAX(EQK(18),SMALL)
      RB(19) = RF(19) / MAX(EQK(19),SMALL)
      RB(20) = RF(20) / MAX(EQK(20),SMALL)
      RB(21) = RF(21) / MAX(EQK(21),SMALL)
      RB(22) = RF(22) / MAX(EQK(22),SMALL)
      RB(23) = RF(23) / MAX(EQK(23),SMALL)
      RB(24) = RF(24) / MAX(EQK(24),SMALL)
      RB(25) = RF(25) / MAX(EQK(25),SMALL)
      RB(26) = RF(26) / MAX(EQK(26),SMALL)
      RB(27) = RF(27) / MAX(EQK(27),SMALL)
      RB(28) = RF(28) / MAX(EQK(28),SMALL)
      RB(29) = RF(29) / MAX(EQK(29),SMALL)
      RB(30) = RF(30) / MAX(EQK(30),SMALL)
      RB(31) = RF(31) / MAX(EQK(31),SMALL)
      RB(32) = RF(32) / MAX(EQK(32),SMALL)
      RB(33) = RF(33) / MAX(EQK(33),SMALL)
      RB(34) = RF(34) / MAX(EQK(34),SMALL)
      RB(35) = RF(35) / MAX(EQK(35),SMALL)
      RB(36) = RF(36) / MAX(EQK(36),SMALL)
      RB(37) = RF(37) / MAX(EQK(37),SMALL)
      RB(38) = RF(38) / MAX(EQK(38),SMALL)
      RB(39) = RF(39) / MAX(EQK(39),SMALL)
      RB(40) = RF(40) / MAX(EQK(40),SMALL)
      RB(41) = RF(41) / MAX(EQK(41),SMALL)
      RB(42) = RF(42) / MAX(EQK(42),SMALL)
      RB(43) = RF(43) / MAX(EQK(43),SMALL)
      RB(44) = RF(44) / MAX(EQK(44),SMALL)
      RB(45) = RF(45) / MAX(EQK(45),SMALL)
      RB(46) = RF(46) / MAX(EQK(46),SMALL)
      RB(47) = RF(47) / MAX(EQK(47),SMALL)
      RB(48) = RF(48) / MAX(EQK(48),SMALL)
      RB(49) = RF(49) / MAX(EQK(49),SMALL)
      RB(50) = RF(50) / MAX(EQK(50),SMALL)
      RB(51) = RF(51) / MAX(EQK(51),SMALL)
      RB(52) = RF(52) / MAX(EQK(52),SMALL)
      RB(53) = RF(53) / MAX(EQK(53),SMALL)
      RB(54) = RF(54) / MAX(EQK(54),SMALL)
      RB(55) = RF(55) / MAX(EQK(55),SMALL)
      RB(56) = RF(56) / MAX(EQK(56),SMALL)
      RB(57) = RF(57) / MAX(EQK(57),SMALL)
      RB(58) = RF(58) / MAX(EQK(58),SMALL)
      RB(59) = RF(59) / MAX(EQK(59),SMALL)
      RB(60) = RF(60) / MAX(EQK(60),SMALL)
      RB(61) = RF(61) / MAX(EQK(61),SMALL)
      RB(62) = RF(62) / MAX(EQK(62),SMALL)
      RB(63) = RF(63) / MAX(EQK(63),SMALL)
      RB(64) = RF(64) / MAX(EQK(64),SMALL)
      RB(65) = RF(65) / MAX(EQK(65),SMALL)
      RB(66) = RF(66) / MAX(EQK(66),SMALL)
      RB(67) = RF(67) / MAX(EQK(67),SMALL)
      RB(68) = RF(68) / MAX(EQK(68),SMALL)
      RB(69) = RF(69) / MAX(EQK(69),SMALL)
      RB(70) = RF(70) / MAX(EQK(70),SMALL)
      RB(71) = RF(71) / MAX(EQK(71),SMALL)
      RB(72) = RF(72) / MAX(EQK(72),SMALL)
      RB(73) = RF(73) / MAX(EQK(73),SMALL)
      RB(74) = RF(74) / MAX(EQK(74),SMALL)
      RB(75) = RF(75) / MAX(EQK(75),SMALL)
      RB(76) = RF(76) / MAX(EQK(76),SMALL)
      RB(77) = RF(77) / MAX(EQK(77),SMALL)
      RB(78) = RF(78) / MAX(EQK(78),SMALL)
      RB(79) = RF(79) / MAX(EQK(79),SMALL)
      RB(80) = RF(80) / MAX(EQK(80),SMALL)
      RB(81) = RF(81) / MAX(EQK(81),SMALL)
      RB(82) = RF(82) / MAX(EQK(82),SMALL)
      RB(83) = RF(83) / MAX(EQK(83),SMALL)
      RB(84) = RF(84) / MAX(EQK(84),SMALL)
      RB(85) = RF(85) / MAX(EQK(85),SMALL)
      RB(86) = RF(86) / MAX(EQK(86),SMALL)
      RB(87) = RF(87) / MAX(EQK(87),SMALL)
      RB(88) = RF(88) / MAX(EQK(88),SMALL)
      RB(89) = RF(89) / MAX(EQK(89),SMALL)
      RB(90) = RF(90) / MAX(EQK(90),SMALL)
      RB(91) = RF(91) / MAX(EQK(91),SMALL)
      RB(92) = RF(92) / MAX(EQK(92),SMALL)
      RB(93) = RF(93) / MAX(EQK(93),SMALL)
      RB(94) = RF(94) / MAX(EQK(94),SMALL)
      RB(95) = RF(95) / MAX(EQK(95),SMALL)
      RB(96) = RF(96) / MAX(EQK(96),SMALL)
      RB(97) = RF(97) / MAX(EQK(97),SMALL)
      RB(98) = RF(98) / MAX(EQK(98),SMALL)
      RB(99) = RF(99) / MAX(EQK(99),SMALL)
      RB(100) = RF(100) / MAX(EQK(100),SMALL)
      RB(101) = RF(101) / MAX(EQK(101),SMALL)
      RB(102) = RF(102) / MAX(EQK(102),SMALL)
      RB(103) = RF(103) / MAX(EQK(103),SMALL)
      RB(104) = RF(104) / MAX(EQK(104),SMALL)
      RB(105) = RF(105) / MAX(EQK(105),SMALL)
      RB(106) = RF(106) / MAX(EQK(106),SMALL)
      RB(107) = RF(107) / MAX(EQK(107),SMALL)
      RB(108) = RF(108) / MAX(EQK(108),SMALL)
      RB(109) = RF(109) / MAX(EQK(109),SMALL)
      RB(110) = RF(110) / MAX(EQK(110),SMALL)
      RB(111) = RF(111) / MAX(EQK(111),SMALL)
      RB(112) = RF(112) / MAX(EQK(112),SMALL)
      RB(113) = RF(113) / MAX(EQK(113),SMALL)
      RB(114) = RF(114) / MAX(EQK(114),SMALL)
      RB(115) = RF(115) / MAX(EQK(115),SMALL)
      RB(116) = RF(116) / MAX(EQK(116),SMALL)
      RB(117) = RF(117) / MAX(EQK(117),SMALL)
      RB(118) = RF(118) / MAX(EQK(118),SMALL)
      RB(119) = RF(119) / MAX(EQK(119),SMALL)
      RB(120) = RF(120) / MAX(EQK(120),SMALL)
      RB(121) = RF(121) / MAX(EQK(121),SMALL)
      RB(122) = RF(122) / MAX(EQK(122),SMALL)
      RB(123) = RF(123) / MAX(EQK(123),SMALL)
      RB(124) = RF(124) / MAX(EQK(124),SMALL)
      RB(125) = RF(125) / MAX(EQK(125),SMALL)
      RB(126) = 0.0
      RB(127) = RF(127) / MAX(EQK(127),SMALL)
      RB(128) = RF(128) / MAX(EQK(128),SMALL)
      RB(129) = RF(129) / MAX(EQK(129),SMALL)
      RB(130) = RF(130) / MAX(EQK(130),SMALL)
      RB(131) = RF(131) / MAX(EQK(131),SMALL)
      RB(132) = RF(132) / MAX(EQK(132),SMALL)
      RB(133) = RF(133) / MAX(EQK(133),SMALL)
      RB(134) = RF(134) / MAX(EQK(134),SMALL)
      RB(135) = RF(135) / MAX(EQK(135),SMALL)
      RB(136) = RF(136) / MAX(EQK(136),SMALL)
      RB(137) = RF(137) / MAX(EQK(137),SMALL)
      RB(138) = RF(138) / MAX(EQK(138),SMALL)
      RB(139) = RF(139) / MAX(EQK(139),SMALL)
      RB(140) = RF(140) / MAX(EQK(140),SMALL)
      RB(141) = RF(141) / MAX(EQK(141),SMALL)
      RB(142) = RF(142) / MAX(EQK(142),SMALL)
      RB(143) = RF(143) / MAX(EQK(143),SMALL)
      RB(144) = RF(144) / MAX(EQK(144),SMALL)
      RB(145) = RF(145) / MAX(EQK(145),SMALL)
      RB(146) = RF(146) / MAX(EQK(146),SMALL)
      RB(147) = RF(147) / MAX(EQK(147),SMALL)
      RB(148) = RF(148) / MAX(EQK(148),SMALL)
      RB(149) = RF(149) / MAX(EQK(149),SMALL)
      RB(150) = RF(150) / MAX(EQK(150),SMALL)
      RB(151) = RF(151) / MAX(EQK(151),SMALL)
      RB(152) = RF(152) / MAX(EQK(152),SMALL)
      RB(153) = RF(153) / MAX(EQK(153),SMALL)
      RB(154) = RF(154) / MAX(EQK(154),SMALL)
      RB(155) = RF(155) / MAX(EQK(155),SMALL)
      RB(156) = RF(156) / MAX(EQK(156),SMALL)
      RB(157) = RF(157) / MAX(EQK(157),SMALL)
      RB(158) = RF(158) / MAX(EQK(158),SMALL)
      RB(159) = RF(159) / MAX(EQK(159),SMALL)
      RB(160) = RF(160) / MAX(EQK(160),SMALL)
      RB(161) = RF(161) / MAX(EQK(161),SMALL)
      RB(162) = RF(162) / MAX(EQK(162),SMALL)
      RB(163) = RF(163) / MAX(EQK(163),SMALL)
      RB(164) = RF(164) / MAX(EQK(164),SMALL)
      RB(165) = RF(165) / MAX(EQK(165),SMALL)
      RB(166) = 0.0
      RB(167) = RF(167) / MAX(EQK(167),SMALL)
      RB(168) = RF(168) / MAX(EQK(168),SMALL)
      RB(169) = 0.0
      RB(170) = RF(170) / MAX(EQK(170),SMALL)
      RB(171) = 0.0
      RB(172) = RF(172) / MAX(EQK(172),SMALL)
      RB(173) = 0.0
      RB(174) = 0.0
      RB(175) = RF(175) / MAX(EQK(175),SMALL)
      RB(176) = RF(176) / MAX(EQK(176),SMALL)
      RB(177) = RF(177) / MAX(EQK(177),SMALL)
      RB(178) = 0.0
      RB(179) = 0.0
      RB(180) = 0.0
      RB(181) = RF(181) / MAX(EQK(181),SMALL)
      RB(182) = RF(182) / MAX(EQK(182),SMALL)
      RB(183) = RF(183) / MAX(EQK(183),SMALL)
      RB(184) = RF(184) / MAX(EQK(184),SMALL)
C
      RKLOW(1) = EXP(3.40312786D1 -1.50965D3/T)
      RKLOW(2) = EXP(5.99064331D1 -2.76D0*ALOGT -8.05146668D2/T)
      RKLOW(3) = EXP(7.69484824D1 -4.76D0*ALOGT -1.22784867D3/T)
      RKLOW(4) = EXP(5.61662604D1 -2.57D0*ALOGT -2.13867084D2/T)
      RKLOW(5) = EXP(7.39217399D1 -4.82D0*ALOGT -3.28600484D3/T)
      RKLOW(6) = EXP(6.98660102D1 -4.8D0*ALOGT -2.79788467D3/T)
      RKLOW(7) = EXP(7.28526099D1 -4.65D0*ALOGT -2.55634067D3/T)
      RKLOW(8) = EXP(9.59450043D1 -7.44D0*ALOGT -7.08529068D3/T)
      RKLOW(9) = EXP(9.34384048D1 -7.27D0*ALOGT -3.63322434D3/T)
      RKLOW(10) = EXP(6.9414025D1 -3.86D0*ALOGT -1.67067934D3/T)
      RKLOW(11) = EXP(9.61977483D1 -7.62D0*ALOGT -3.50742017D3/T)
      RKLOW(12) = EXP(9.50941235D1 -7.08D0*ALOGT -3.36400342D3/T)
      RKLOW(13) = EXP(6.37931383D1 -3.42D0*ALOGT -4.24463259D4/T)
      RKLOW(14) = EXP(4.22794408D1 -9.D-1*ALOGT +8.55468335D2/T)
      RKLOW(15) = EXP(8.42793577D1 -5.92D0*ALOGT -1.58010034D3/T)
      RKLOW(16) = EXP(6.54619238D1 -3.74D0*ALOGT -9.74227469D2/T)
      RKLOW(17) = EXP(7.69748493D1 -5.11D0*ALOGT -3.57032226D3/T)
      RKLOW(18) = EXP(8.81295053D1 -6.36D0*ALOGT -2.536212D3/T)
      RKLOW(19) = EXP(9.56297642D1 -7.03D0*ALOGT -1.38988444D3/T)
      RKLOW(20) = EXP(1.17889265D2 -9.3D0*ALOGT -4.92145901D4/T)
      RKLOW(21) = EXP(5.91374013D1 -2.8D0*ALOGT -2.96897834D2/T)
      RKLOW(22) = EXP(9.67205025D1 -7.63D0*ALOGT -1.93939704D3/T)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RDSMH  (T, SMH)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION SMH(*), TN(5)
C
      TLOG = LOG(T)
      TI = 1.0D0/T
      TN(1) = TLOG - 1.0
      TN(2) = T
      TN(3) = TN(2)*T
      TN(4) = TN(3)*T
      TN(5) = TN(4)*T
C
      IF (T .GT. 1.D3) THEN
C
      SMH(1) = -3.20502331D+00 + 9.50158922D+02*TI 
     *         + 3.33727920D+00*TN(1) - 2.47012365D-05*TN(2) 
     *         + 8.32427963D-08*TN(3) - 1.49638662D-11*TN(4) 
     *         + 1.00127688D-15*TN(5) 
      SMH(2) = -4.46682914D-01 - 2.54736599D+04*TI 
     *         + 2.50000001D+00*TN(1) - 1.15421486D-11*TN(2) 
     *         + 2.69269913D-15*TN(3) - 3.94596029D-19*TN(4) 
     *         + 2.49098679D-23*TN(5) 
      SMH(3) = 4.78433864D+00 - 2.92175791D+04*TI 
     *         + 2.56942078D+00*TN(1) - 4.29870569D-05*TN(2) 
     *         + 6.99140982D-09*TN(3) - 8.34814992D-13*TN(4) 
     *         + 6.14168455D-17*TN(5) 
      SMH(4) = 5.45323129D+00 + 1.08845772D+03*TI 
     *         + 3.28253784D+00*TN(1) + 7.41543770D-04*TN(2) 
     *         - 1.26327778D-07*TN(3) + 1.74558796D-11*TN(4) 
     *         - 1.08358897D-15*TN(5) 
      SMH(5) = 4.47669610D+00 - 3.85865700D+03*TI 
     *         + 3.09288767D+00*TN(1) + 2.74214858D-04*TN(2) 
     *         + 2.10842047D-08*TN(3) - 7.32884630D-12*TN(4) 
     *         + 5.87061880D-16*TN(5) 
      SMH(6) = 4.96677010D+00 + 3.00042971D+04*TI 
     *         + 3.03399249D+00*TN(1) + 1.08845902D-03*TN(2) 
     *         - 2.73454197D-08*TN(3) - 8.08683225D-12*TN(4) 
     *         + 8.41004960D-16*TN(5) 
      SMH(7) = 3.78510215D+00 - 1.11856713D+02*TI 
     *         + 4.01721090D+00*TN(1) + 1.11991007D-03*TN(2) 
     *         - 1.05609692D-07*TN(3) + 9.52053083D-12*TN(4) 
     *         - 5.39542675D-16*TN(5) 
      SMH(8) = 2.91615662D+00 + 1.78617877D+04*TI 
     *         + 4.16500285D+00*TN(1) + 2.45415847D-03*TN(2) 
     *         - 3.16898708D-07*TN(3) + 3.09321655D-11*TN(4) 
     *         - 1.43954153D-15*TN(5) 
      SMH(9) = 4.80150373D+00 - 8.54512953D+04*TI 
     *         + 2.49266888D+00*TN(1) + 2.39944642D-05*TN(2) 
     *         - 1.20722503D-08*TN(3) + 3.11909191D-12*TN(4) 
     *         - 2.43638946D-16*TN(5) 
      SMH(10) = 5.48497999D+00 - 7.10124364D+04*TI 
     *         + 2.87846473D+00*TN(1) + 4.85456841D-04*TN(2) 
     *         + 2.40742758D-08*TN(3) - 1.08906541D-11*TN(4) 
     *         + 8.80396915D-16*TN(5) 
      SMH(11) = 6.17119324D+00 - 4.62636040D+04*TI 
     *         + 2.87410113D+00*TN(1) + 1.82819646D-03*TN(2) 
     *         - 2.34824328D-07*TN(3) + 2.16816291D-11*TN(4) 
     *         - 9.38637835D-16*TN(5) 
      SMH(12) = 8.62650169D+00 - 5.09259997D+04*TI 
     *         + 2.29203842D+00*TN(1) + 2.32794319D-03*TN(2) 
     *         - 3.35319912D-07*TN(3) + 3.48255000D-11*TN(4) 
     *         - 1.69858183D-15*TN(5) 
      SMH(13) = 8.48007179D+00 - 1.67755843D+04*TI 
     *         + 2.28571772D+00*TN(1) + 3.61995018D-03*TN(2) 
     *         - 4.97857247D-07*TN(3) + 4.96403870D-11*TN(4) 
     *         - 2.33577197D-15*TN(5) 
      SMH(14) = 1.84373180D+01 + 9.46834459D+03*TI 
     *         + 7.48514950D-02*TN(1) + 6.69547335D-03*TN(2) 
     *         - 9.55476348D-07*TN(3) + 1.01910446D-10*TN(4) 
     *         - 5.09076150D-15*TN(5) 
      SMH(15) = 7.81868772D+00 + 1.41518724D+04*TI 
     *         + 2.71518561D+00*TN(1) + 1.03126372D-03*TN(2) 
     *         - 1.66470962D-07*TN(3) + 1.91710840D-11*TN(4) 
     *         - 1.01823858D-15*TN(5) 
      SMH(16) = 2.27163806D+00 + 4.87591660D+04*TI 
     *         + 3.85746029D+00*TN(1) + 2.20718513D-03*TN(2) 
     *         - 3.69135673D-07*TN(3) + 4.36241823D-11*TN(4) 
     *         - 2.36042082D-15*TN(5) 
      SMH(17) = 9.79834492D+00 - 4.01191815D+03*TI 
     *         + 2.77217438D+00*TN(1) + 2.47847763D-03*TN(2) 
     *         - 4.14076022D-07*TN(3) + 4.90968148D-11*TN(4) 
     *         - 2.66754356D-15*TN(5) 
      SMH(18) = 1.36563230D+01 + 1.39958323D+04*TI 
     *         + 1.76069008D+00*TN(1) + 4.60000041D-03*TN(2) 
     *         - 7.37098022D-07*TN(3) + 8.38676767D-11*TN(4) 
     *         - 4.41927820D-15*TN(5) 
      SMH(19) = 5.81043215D+00 + 3.24250627D+03*TI 
     *         + 3.69266569D+00*TN(1) + 4.32288399D-03*TN(2) 
     *         - 6.25168533D-07*TN(3) + 6.56028863D-11*TN(4) 
     *         - 3.24277101D-15*TN(5) 
      SMH(20) = 2.92957500D+00 - 1.27832520D+02*TI 
     *         + 3.77079900D+00*TN(1) + 3.93574850D-03*TN(2) 
     *         - 4.42730667D-07*TN(3) + 3.28702583D-11*TN(4) 
     *         - 1.05630800D-15*TN(5) 
      SMH(21) = 1.45023623D+01 + 2.53748747D+04*TI 
     *         + 1.78970791D+00*TN(1) + 7.04691460D-03*TN(2) 
     *         - 1.06083472D-06*TN(3) + 1.15142571D-10*TN(4) 
     *         - 5.85301100D-15*TN(5) 
      SMH(22) = -1.23028121D+00 - 2.59359992D+04*TI 
     *         + 4.14756964D+00*TN(1) + 2.98083332D-03*TN(2) 
     *         - 3.95491420D-07*TN(3) + 3.89510143D-11*TN(4) 
     *         - 1.80617607D-15*TN(5) 
      SMH(23) = 7.78732378D+00 - 3.46128739D+04*TI 
     *         + 3.01672400D+00*TN(1) + 5.16511460D-03*TN(2) 
     *         - 7.80137248D-07*TN(3) + 8.48027400D-11*TN(4) 
     *         - 4.31303520D-15*TN(5) 
      SMH(24) = 1.03053693D+01 - 4.93988614D+03*TI 
     *         + 2.03611116D+00*TN(1) + 7.32270755D-03*TN(2) 
     *         - 1.11846319D-06*TN(3) + 1.22685769D-10*TN(4) 
     *         - 6.28530305D-15*TN(5) 
      SMH(25) = 1.34624343D+01 - 1.28575200D+04*TI 
     *         + 1.95465642D+00*TN(1) + 8.69863610D-03*TN(2) 
     *         - 1.33034445D-06*TN(3) + 1.46014741D-10*TN(4) 
     *         - 7.48207880D-15*TN(5) 
      SMH(26) = 1.51156107D+01 + 1.14263932D+04*TI 
     *         + 1.07188150D+00*TN(1) + 1.08426339D-02*TN(2) 
     *         - 1.67093445D-06*TN(3) + 1.84510001D-10*TN(4) 
     *         - 9.50014450D-15*TN(5) 
      SMH(27) = -3.93025950D+00 - 1.93272150D+04*TI 
     *         + 5.62820580D+00*TN(1) + 2.04267005D-03*TN(2) 
     *         - 2.65575783D-07*TN(3) + 2.38550433D-11*TN(4) 
     *         - 9.70391600D-16*TN(5) 
      SMH(28) = 6.32247205D-01 + 7.55105311D+03*TI 
     *         + 4.51129732D+00*TN(1) + 4.50179872D-03*TN(2) 
     *         - 6.94899392D-07*TN(3) + 7.69454902D-11*TN(4) 
     *         - 3.97419100D-15*TN(5) 
      SMH(29) = -5.04525100D+00 - 4.90321800D+02*TI 
     *         + 5.97567000D+00*TN(1) + 4.06529550D-03*TN(2) 
     *         - 4.57270667D-07*TN(3) + 3.39192000D-11*TN(4) 
     *         - 1.08800850D-15*TN(5) 
C
      ELSE
C
      SMH(1) = 6.83010238D-01 + 9.17935173D+02*TI 
     *         + 2.34433112D+00*TN(1) + 3.99026037D-03*TN(2) 
     *         - 3.24635850D-06*TN(3) + 1.67976745D-09*TN(4) 
     *         - 3.68805881D-13*TN(5) 
      SMH(2) = -4.46682853D-01 - 2.54736599D+04*TI 
     *         + 2.50000000D+00*TN(1) + 3.52666409D-13*TN(2) 
     *         - 3.32653273D-16*TN(3) + 1.91734693D-19*TN(4) 
     *         - 4.63866166D-23*TN(5) 
      SMH(3) = 2.05193346D+00 - 2.91222592D+04*TI 
     *         + 3.16826710D+00*TN(1) - 1.63965942D-03*TN(2) 
     *         + 1.10717733D-06*TN(3) - 5.10672187D-10*TN(4) 
     *         + 1.05632986D-13*TN(5) 
      SMH(4) = 3.65767573D+00 + 1.06394356D+03*TI 
     *         + 3.78245636D+00*TN(1) - 1.49836708D-03*TN(2) 
     *         + 1.64121700D-06*TN(3) - 8.06774591D-10*TN(4) 
     *         + 1.62186419D-13*TN(5) 
      SMH(5) = -1.03925458D-01 - 3.61508056D+03*TI 
     *         + 3.99201543D+00*TN(1) - 1.20065876D-03*TN(2) 
     *         + 7.69656402D-07*TN(3) - 3.23427778D-10*TN(4) 
     *         + 6.82057350D-14*TN(5) 
      SMH(6) = -8.49032208D-01 + 3.02937267D+04*TI 
     *         + 4.19864056D+00*TN(1) - 1.01821705D-03*TN(2) 
     *         + 1.08673369D-06*TN(3) - 4.57330885D-10*TN(4) 
     *         + 8.85989085D-14*TN(5) 
      SMH(7) = 3.71666245D+00 - 2.94808040D+02*TI 
     *         + 4.30179801D+00*TN(1) - 2.37456025D-03*TN(2) 
     *         + 3.52638152D-06*TN(3) - 2.02303245D-09*TN(4) 
     *         + 4.64612562D-13*TN(5) 
      SMH(8) = 3.43505074D+00 + 1.77025821D+04*TI 
     *         + 4.27611269D+00*TN(1) - 2.71411208D-04*TN(2) 
     *         + 2.78892835D-06*TN(3) - 1.79809011D-09*TN(4) 
     *         + 4.31227182D-13*TN(5) 
      SMH(9) = 4.53130848D+00 - 8.54438832D+04*TI 
     *         + 2.55423955D+00*TN(1) - 1.60768862D-04*TN(2) 
     *         + 1.22298708D-07*TN(3) - 6.10195741D-11*TN(4) 
     *         + 1.33260723D-14*TN(5) 
      SMH(10) = 2.08401108D+00 - 7.07972934D+04*TI 
     *         + 3.48981665D+00*TN(1) + 1.61917771D-04*TN(2) 
     *         - 2.81498442D-07*TN(3) + 2.63514439D-10*TN(4) 
     *         - 7.03045335D-14*TN(5) 
      SMH(11) = 1.56253185D+00 - 4.60040401D+04*TI 
     *         + 3.76267867D+00*TN(1) + 4.84436072D-04*TN(2) 
     *         + 4.65816402D-07*TN(3) - 3.20909294D-10*TN(4) 
     *         + 8.43708595D-14*TN(5) 
      SMH(12) = -7.69118967D-01 - 5.04968163D+04*TI 
     *         + 4.19860411D+00*TN(1) - 1.18330710D-03*TN(2) 
     *         + 1.37216037D-06*TN(3) - 5.57346651D-10*TN(4) 
     *         + 9.71573685D-14*TN(5) 
      SMH(13) = 1.60456433D+00 - 1.64449988D+04*TI 
     *         + 3.67359040D+00*TN(1) + 1.00547588D-03*TN(2) 
     *         + 9.55036427D-07*TN(3) - 5.72597854D-10*TN(4) 
     *         + 1.27192867D-13*TN(5) 
      SMH(14) = -4.64130376D+00 + 1.02466476D+04*TI 
     *         + 5.14987613D+00*TN(1) - 6.83548940D-03*TN(2) 
     *         + 8.19667665D-06*TN(3) - 4.03952522D-09*TN(4) 
     *         + 8.33469780D-13*TN(5) 
      SMH(15) = 3.50840928D+00 + 1.43440860D+04*TI 
     *         + 3.57953347D+00*TN(1) - 3.05176840D-04*TN(2) 
     *         + 1.69469055D-07*TN(3) + 7.55838237D-11*TN(4) 
     *         - 4.52212249D-14*TN(5) 
      SMH(16) = 9.90105222D+00 + 4.83719697D+04*TI 
     *         + 2.35677352D+00*TN(1) + 4.49229839D-03*TN(2) 
     *         - 1.18726045D-06*TN(3) + 2.04932518D-10*TN(4) 
     *         - 7.18497740D-15*TN(5) 
      SMH(17) = 3.39437243D+00 - 3.83956496D+03*TI 
     *         + 4.22118584D+00*TN(1) - 1.62196266D-03*TN(2) 
     *         + 2.29665743D-06*TN(3) - 1.10953411D-09*TN(4) 
     *         + 2.16884433D-13*TN(5) 
      SMH(18) = 6.02812900D-01 + 1.43089567D+04*TI 
     *         + 4.79372315D+00*TN(1) - 4.95416685D-03*TN(2) 
     *         + 6.22033347D-06*TN(3) - 3.16071051D-09*TN(4) 
     *         + 6.58863260D-13*TN(5) 
      SMH(19) = 5.47302243D+00 + 3.19391367D+03*TI 
     *         + 3.86388918D+00*TN(1) + 2.79836152D-03*TN(2) 
     *         + 9.88786318D-07*TN(3) - 8.71100100D-10*TN(4) 
     *         + 2.18483639D-13*TN(5) 
      SMH(20) = 1.31521770D+01 - 9.78601100D+02*TI 
     *         + 2.10620400D+00*TN(1) + 3.60829750D-03*TN(2) 
     *         + 8.89745333D-07*TN(3) - 6.14803000D-10*TN(4) 
     *         + 1.03780500D-13*TN(5) 
      SMH(21) = -1.50409823D+00 + 2.56427656D+04*TI 
     *         + 5.71539582D+00*TN(1) - 7.61545645D-03*TN(2) 
     *         + 1.08740193D-05*TN(3) - 5.92339074D-09*TN(4) 
     *         + 1.30676349D-12*TN(5) 
      SMH(22) = 1.39397051D+01 - 2.64289807D+04*TI 
     *         + 8.08681094D-01*TN(1) + 1.16807815D-02*TN(2) 
     *         - 5.91953025D-06*TN(3) + 2.33460364D-09*TN(4) 
     *         - 4.25036487D-13*TN(5) 
      SMH(23) = 8.51054025D+00 - 3.48598468D+04*TI 
     *         + 3.21246645D+00*TN(1) + 7.57395810D-04*TN(2) 
     *         + 4.32015687D-06*TN(3) - 2.98048206D-09*TN(4) 
     *         + 7.35754365D-13*TN(5) 
      SMH(24) = 4.09733096D+00 - 5.08977593D+03*TI 
     *         + 3.95920148D+00*TN(1) - 3.78526124D-03*TN(2) 
     *         + 9.51650487D-06*TN(3) - 5.76323961D-09*TN(4) 
     *         + 1.34942187D-12*TN(5) 
      SMH(25) = 4.70720924D+00 - 1.28416265D+04*TI 
     *         + 4.30646568D+00*TN(1) - 2.09329446D-03*TN(2) 
     *         + 8.28571345D-06*TN(3) - 4.99272172D-09*TN(4) 
     *         + 1.15254502D-12*TN(5) 
      SMH(26) = 2.66682316D+00 + 1.15222055D+04*TI 
     *         + 4.29142492D+00*TN(1) - 2.75077135D-03*TN(2) 
     *         + 9.99063813D-06*TN(3) - 5.90388571D-09*TN(4) 
     *         + 1.34342886D-12*TN(5) 
      SMH(27) = 1.24904170D+01 - 2.00594490D+04*TI 
     *         + 2.25172140D+00*TN(1) + 8.82751050D-03*TN(2) 
     *         - 3.95485017D-06*TN(3) + 1.43964658D-09*TN(4) 
     *         - 2.53324055D-13*TN(5) 
      SMH(28) = 1.22156480D+01 + 7.04291804D+03*TI 
     *         + 2.13583630D+00*TN(1) + 9.05943605D-03*TN(2) 
     *         - 2.89912457D-06*TN(3) + 7.78664640D-10*TN(4) 
     *         - 1.00728807D-13*TN(5) 
      SMH(29) = 9.55829000D+00 - 1.52147660D+03*TI 
     *         + 3.40906200D+00*TN(1) + 5.36928700D-03*TN(2) 
     *         + 3.15248667D-07*TN(3) - 5.96548583D-10*TN(4) 
     *         + 1.43369250D-13*TN(5) 
      ENDIF
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RATX (T, C, RF, RB, RKLOW)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (SMALL = 1.D-200)
      DIMENSION C(*), RF(*), RB(*), RKLOW(*)
      DIMENSION CTB(184)
C
      ALOGT = LOG(T)
      CTOT = 0.0
      DO K = 1, 19
         CTOT = CTOT + C(K)
      ENDDO
C
      CTB(1)   = CTOT + 1.4D0*C(1) + 1.44D1*C(6) + C(10) + 7.5D-1*C(11) 
     *                + 2.6D0*C(12) + 2.D0*C(17) 
      CTB(2)   = CTOT + C(1) + 5.D0*C(6) + C(10) + 5.D-1*C(11) + C(12) 
     *                + 2.D0*C(17) 
      CTB(47)  = CTB(2)
      CTB(51)  = CTB(2)
      CTB(53)  = CTB(2)
      CTB(54)  = CTB(2)
      CTB(56)  = CTB(2)
      CTB(60)  = CTB(2)
      CTB(67)  = CTB(2)
      CTB(68)  = CTB(2)
      CTB(70)  = CTB(2)
      CTB(72)  = CTB(2)
      CTB(78)  = CTB(2)
      CTB(80)  = CTB(2)
      CTB(90)  = CTB(2)
      CTB(122) = CTB(2)
      CTB(131) = CTB(2)
      CTB(137) = CTB(2)
      CTB(148) = CTB(2)
      CTB(162) = CTB(2)
      CTB(170) = CTB(2)
      CTB(177) = CTB(2)
      CTB(12)  = CTOT + C(1) + 5.D0*C(4) + 5.D0*C(6) + C(10) 
     *                + 5.D-1*C(11) + 2.5D0*C(12) + 2.D0*C(17) 
      CTB(31)  = CTOT - C(4) - C(6) - 2.5D-1*C(11) + 5.D-1*C(12) 
     *                + 5.D-1*C(17) - C(19) 
      CTB(36)  = CTOT - C(1) - C(6) + C(10) - C(12) + 2.D0*C(17) 
      CTB(40)  = CTOT - 2.7D-1*C(1) + 2.65D0*C(6) + C(10) + 2.D0*C(17) 
      CTB(49)  = CTOT + C(1) + 5.D0*C(6) + 2.D0*C(10) + 5.D-1*C(11) 
     *                + C(12) + 2.D0*C(17) 
      CTB(157) = CTOT + C(1) - C(6) + C(10) + 5.D-1*C(11) + C(12) 
     *                + 2.D0*C(17) 
C
      PR = RKLOW(1) * CTB(12) / RF(12)
      PCOR = PR / (1.0 + PR)
      RF(12) = RF(12) * PCOR
      RB(12) = RB(12) * PCOR
C
      PR = RKLOW(2) * CTB(47) / RF(47)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 4.38D-1*EXP(-T/9.1D1) + 5.62D-1*EXP(-T/5.836D3)
     *     + EXP(-8.552D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(47) = RF(47) * PCOR
      RB(47) = RB(47) * PCOR
C
      PR = RKLOW(3) * CTB(49) / RF(49)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.17D-1*EXP(-T/7.4D1) + 7.83D-1*EXP(-T/2.941D3)
     *     + EXP(-6.964D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(49) = RF(49) * PCOR
      RB(49) = RB(49) * PCOR
C
      PR = RKLOW(4) * CTB(51) / RF(51)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.176D-1*EXP(-T/2.71D2) + 7.824D-1*EXP(-T/2.755D3)
     *     + EXP(-6.57D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(51) = RF(51) * PCOR
      RB(51) = RB(51) * PCOR
C
      PR = RKLOW(5) * CTB(53) / RF(53)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.813D-1*EXP(-T/1.03D2) + 7.187D-1*EXP(-T/1.291D3)
     *     + EXP(-4.16D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(53) = RF(53) * PCOR
      RB(53) = RB(53) * PCOR
C
      PR = RKLOW(6) * CTB(54) / RF(54)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.42D-1*EXP(-T/9.4D1) + 7.58D-1*EXP(-T/1.555D3)
     *     + EXP(-4.2D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(54) = RF(54) * PCOR
      RB(54) = RB(54) * PCOR
C
      PR = RKLOW(7) * CTB(56) / RF(56)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 4.D-1*EXP(-T/1.D2) + 6.D-1*EXP(-T/9.D4)
     *     + EXP(-1.D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(56) = RF(56) * PCOR
      RB(56) = RB(56) * PCOR
C
      PR = RKLOW(8) * CTB(60) / RF(60)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 3.D-1*EXP(-T/1.D2) + 7.D-1*EXP(-T/9.D4)
     *     + EXP(-1.D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(60) = RF(60) * PCOR
      RB(60) = RB(60) * PCOR
C
      PR = RKLOW(9) * CTB(67) / RF(67)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.493D-1*EXP(-T/9.85D1) + 7.507D-1*EXP(-T/1.302D3)
     *     + EXP(-4.167D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(67) = RF(67) * PCOR
      RB(67) = RB(67) * PCOR
C
      PR = RKLOW(10) * CTB(68) / RF(68)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.18D-1*EXP(-T/2.075D2) + 7.82D-1*EXP(-T/2.663D3)
     *     + EXP(-6.095D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(68) = RF(68) * PCOR
      RB(68) = RB(68) * PCOR
C
      PR = RKLOW(11) * CTB(70) / RF(70)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.47D-2*EXP(-T/2.1D2) + 9.753D-1*EXP(-T/9.84D2)
     *     + EXP(-4.374D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(70) = RF(70) * PCOR
      RB(70) = RB(70) * PCOR
C
      PR = RKLOW(12) * CTB(72) / RF(72)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 1.578D-1*EXP(-T/1.25D2) + 8.422D-1*EXP(-T/2.219D3)
     *     + EXP(-6.882D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(72) = RF(72) * PCOR
      RB(72) = RB(72) * PCOR
C
      PR = RKLOW(13) * CTB(78) / RF(78)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 6.8D-2*EXP(-T/1.97D2) + 9.32D-1*EXP(-T/1.54D3)
     *     + EXP(-1.03D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(78) = RF(78) * PCOR
      RB(78) = RB(78) * PCOR
C
      PR = RKLOW(14) * CTB(80) / RF(80)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.654D-1*EXP(-T/9.4D1) + 7.346D-1*EXP(-T/1.756D3)
     *     + EXP(-5.182D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(80) = RF(80) * PCOR
      RB(80) = RB(80) * PCOR
C
      PR = RKLOW(15) * CTB(90) / RF(90)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 5.88D-1*EXP(-T/1.95D2) + 4.12D-1*EXP(-T/5.9D3)
     *     + EXP(-6.394D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(90) = RF(90) * PCOR
      RB(90) = RB(90) * PCOR
C
      PR = RKLOW(16) * CTB(122) / RF(122)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 4.243D-1*EXP(-T/2.37D2) + 5.757D-1*EXP(-T/1.652D3)
     *     + EXP(-5.069D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(122) = RF(122) * PCOR
      RB(122) = RB(122) * PCOR
C
      PR = RKLOW(17) * CTB(131) / RF(131)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 4.093D-1*EXP(-T/2.75D2) + 5.907D-1*EXP(-T/1.226D3)
     *     + EXP(-5.185D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(131) = RF(131) * PCOR
      RB(131) = RB(131) * PCOR
C
      PR = RKLOW(18) * CTB(137) / RF(137)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 3.973D-1*EXP(-T/2.08D2) + 6.027D-1*EXP(-T/3.922D3)
     *     + EXP(-1.018D4/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(137) = RF(137) * PCOR
      RB(137) = RB(137) * PCOR
C
      PR = RKLOW(19) * CTB(148) / RF(148)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 3.81D-1*EXP(-T/7.32D1) + 6.19D-1*EXP(-T/1.18D3)
     *     + EXP(-9.999D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(148) = RF(148) * PCOR
      RB(148) = RB(148) * PCOR
C
      PR = RKLOW(20) * CTB(162) / RF(162)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 2.655D-1*EXP(-T/1.8D2) + 7.345D-1*EXP(-T/1.035D3)
     *     + EXP(-5.417D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(162) = RF(162) * PCOR
      RB(162) = RB(162) * PCOR
C
      PR = RKLOW(21) * CTB(170) / RF(170)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 4.22D-1*EXP(-T/1.22D2) + 5.78D-1*EXP(-T/2.535D3)
     *     + EXP(-9.365D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(170) = RF(170) * PCOR
      RB(170) = RB(170) * PCOR
C
      PR = RKLOW(22) * CTB(177) / RF(177)
      PCOR = PR / (1.0 + PR)
      PRLOG = LOG10(MAX(PR,SMALL))
      FCENT = 5.35D-1*EXP(-T/2.01D2) + 4.65D-1*EXP(-T/1.773D3)
     *     + EXP(-5.333D3/T)
      FCLOG = LOG10(MAX(FCENT,SMALL))
      XN    = 0.75 - 1.27*FCLOG
      CPRLOG= PRLOG - (0.4 + 0.67*FCLOG)
      FLOG = FCLOG/(1.0 + (CPRLOG/(XN-0.14*CPRLOG))**2)
      FC = 10.0**FLOG
      PCOR = FC * PCOR
      RF(177) = RF(177) * PCOR
      RB(177) = RB(177) * PCOR
C
      RF(1) = RF(1)*CTB(1)*C(3)*C(3)
      RF(2) = RF(2)*CTB(2)*C(3)*C(2)
      RF(3) = RF(3)*C(3)*C(1)
      RF(4) = RF(4)*C(3)*C(7)
      RF(5) = RF(5)*C(3)*C(8)
      RF(6) = RF(6)*C(3)
      RF(7) = RF(7)*C(3)
      RF(8) = RF(8)*C(3)
      RF(9) = RF(9)*C(3)
      RF(10) = RF(10)*C(3)*C(9)
      RF(11) = RF(11)*C(3)*C(10)
      RF(12) = RF(12)*C(3)*C(11)
      RF(13) = RF(13)*C(3)
      RF(14) = RF(14)*C(3)
      RF(15) = RF(15)*C(3)*C(13)
      RF(16) = RF(16)*C(3)
      RF(17) = RF(17)*C(3)
      RF(18) = RF(18)*C(3)*C(14)
      RF(19) = RF(19)*C(3)*C(14)
      RF(20) = RF(20)*C(3)*C(15)
      RF(21) = RF(21)*C(3)*C(15)
      RF(22) = RF(22)*C(3)
      RF(23) = RF(23)*C(3)*C(16)
      RF(24) = RF(24)*C(3)
      RF(25) = RF(25)*C(3)*C(17)
      RF(26) = RF(26)*C(3)
      RF(27) = RF(27)*C(3)*C(18)
      RF(28) = RF(28)*C(3)*C(18)
      RF(29) = RF(29)*C(4)*C(11)
      RF(30) = RF(30)*C(4)*C(13)
      RF(31) = RF(31)*CTB(31)*C(2)*C(4)
      RF(32) = RF(32)*C(2)*C(4)*C(4)
      RF(33) = RF(33)*C(2)*C(4)*C(6)
      RF(34) = RF(34)*C(2)*C(4)*C(19)
      RF(35) = RF(35)*C(2)*C(4)
      RF(36) = RF(36)*CTB(36)*C(2)*C(2)
      RF(37) = RF(37)*C(2)*C(2)*C(1)
      RF(38) = RF(38)*C(2)*C(2)*C(6)
      RF(39) = RF(39)*C(2)*C(2)*C(12)
      RF(40) = RF(40)*CTB(40)*C(2)*C(5)
      RF(41) = RF(41)*C(2)*C(7)
      RF(42) = RF(42)*C(2)*C(7)
      RF(43) = RF(43)*C(2)*C(7)
      RF(44) = RF(44)*C(2)*C(8)
      RF(45) = RF(45)*C(2)*C(8)
      RF(46) = RF(46)*C(2)
      RF(47) = RF(47)*C(2)
      RF(48) = RF(48)*C(2)
      RF(49) = RF(49)*C(2)*C(9)
      RF(50) = RF(50)*C(2)*C(10)
      RF(51) = RF(51)*C(2)
      RF(52) = RF(52)*C(2)
      RF(53) = RF(53)*C(2)*C(13)
      RF(54) = RF(54)*C(2)*C(13)
      RF(55) = RF(55)*C(2)*C(13)
      RF(56) = RF(56)*C(2)
      RF(57) = RF(57)*C(2)
      RF(58) = RF(58)*C(2)
      RF(59) = RF(59)*C(2)
      RF(60) = RF(60)*C(2)
      RF(61) = RF(61)*C(2)
      RF(62) = RF(62)*C(2)
      RF(63) = RF(63)*C(2)
      RF(64) = RF(64)*C(2)
      RF(65) = RF(65)*C(2)*C(14)
      RF(66) = RF(66)*C(2)*C(14)
      RF(67) = RF(67)*C(2)*C(15)
      RF(68) = RF(68)*C(2)
      RF(69) = RF(69)*C(2)
      RF(70) = RF(70)*C(2)*C(16)
      RF(71) = RF(71)*C(2)*C(16)
      RF(72) = RF(72)*C(2)
      RF(73) = RF(73)*C(2)
      RF(74) = RF(74)*C(2)*C(17)
      RF(75) = RF(75)*C(2)
      RF(76) = RF(76)*C(2)*C(18)
      RF(77) = RF(77)*C(2)*C(18)
      RF(78) = RF(78)*C(1)*C(11)
      RF(79) = RF(79)*C(5)*C(1)
      RF(80) = RF(80)*C(5)*C(5)
      RF(81) = RF(81)*C(5)*C(5)
      RF(82) = RF(82)*C(5)*C(7)
      RF(83) = RF(83)*C(5)*C(8)
      RF(84) = RF(84)*C(5)*C(8)
      RF(85) = RF(85)*C(5)
      RF(86) = RF(86)*C(5)
      RF(87) = RF(87)*C(5)
      RF(88) = RF(88)*C(5)
      RF(89) = RF(89)*C(5)
      RF(90) = RF(90)*C(5)*C(9)
      RF(91) = RF(91)*C(5)*C(9)
      RF(92) = RF(92)*C(5)*C(9)
      RF(93) = RF(93)*C(5)*C(10)
      RF(94) = RF(94)*C(5)*C(11)
      RF(95) = RF(95)*C(5)
      RF(96) = RF(96)*C(5)*C(13)
      RF(97) = RF(97)*C(5)
      RF(98) = RF(98)*C(5)
      RF(99) = RF(99)*C(5)*C(14)
      RF(100) = RF(100)*C(5)*C(14)
      RF(101) = RF(101)*C(5)*C(15)
      RF(102) = RF(102)*C(5)*C(15)
      RF(103) = RF(103)*C(5)
      RF(104) = RF(104)*C(5)*C(16)
      RF(105) = RF(105)*C(5)*C(17)
      RF(106) = RF(106)*C(5)*C(18)
      RF(107) = RF(107)*C(7)*C(7)
      RF(108) = RF(108)*C(7)*C(7)
      RF(109) = RF(109)*C(7)
      RF(110) = RF(110)*C(7)*C(9)
      RF(111) = RF(111)*C(7)*C(9)
      RF(112) = RF(112)*C(7)*C(11)
      RF(113) = RF(113)*C(7)*C(13)
      RF(114) = RF(114)*C(4)
      RF(115) = RF(115)*C(9)
      RF(116) = RF(116)*C(4)
      RF(117) = RF(117)*C(1)
      RF(118) = RF(118)*C(6)
      RF(120) = RF(120)*C(9)
      RF(121) = RF(121)*C(10)
      RF(122) = RF(122)*C(11)
      RF(123) = RF(123)*C(12)
      RF(124) = RF(124)*C(13)
      RF(126) = RF(126)*C(4)
      RF(127) = RF(127)*C(1)
      RF(129) = RF(129)*C(9)
      RF(130) = RF(130)*C(10)
      RF(131) = RF(131)*C(11)
      RF(133) = RF(133)*C(19)
      RF(134) = RF(134)*C(4)
      RF(135) = RF(135)*C(4)
      RF(136) = RF(136)*C(1)
      RF(137) = RF(137)*C(6)
      RF(138) = RF(138)*C(6)
      RF(139) = RF(139)*C(9)
      RF(140) = RF(140)*C(10)
      RF(141) = RF(141)*C(11)
      RF(142) = RF(142)*C(12)
      RF(143) = RF(143)*C(12)
      RF(144) = RF(144)*C(17)
      RF(145) = RF(145)*C(9)*C(4)
      RF(146) = RF(146)*C(9)*C(4)
      RF(147) = RF(147)*C(9)*C(8)
      RF(148) = RF(148)*C(9)*C(9)
      RF(149) = RF(149)*C(9)*C(9)
      RF(150) = RF(150)*C(9)
      RF(151) = RF(151)*C(9)*C(13)
      RF(152) = RF(152)*C(9)*C(14)
      RF(153) = RF(153)*C(9)*C(14)
      RF(154) = RF(154)*C(9)*C(16)
      RF(155) = RF(155)*C(9)*C(17)
      RF(156) = RF(156)*C(6)
      RF(157) = RF(157)*CTB(157)
      RF(158) = RF(158)*C(4)
      RF(159) = RF(159)*C(4)
      RF(160) = RF(160)*C(4)
      RF(161) = RF(161)*C(4)
      RF(162) = RF(162)*C(16)
      RF(163) = RF(163)*C(4)
      RF(164) = RF(164)*C(4)
      RF(166) = RF(166)*C(3)*C(9)
      RF(167) = RF(167)*C(3)*C(16)
      RF(168) = RF(168)*C(5)*C(7)
      RF(169) = RF(169)*C(5)*C(9)
      RF(170) = RF(170)*C(1)
      RF(171) = RF(171)*C(4)
      RF(172) = RF(172)*C(4)
      RF(174) = RF(174)*C(6)
      RF(175) = RF(175)*C(4)
      RF(176) = RF(176)*C(4)
      RF(177) = RF(177)*C(2)*C(18)
      RF(178) = RF(178)*C(3)
      RF(179) = RF(179)*C(4)
      RF(180) = RF(180)*C(4)
      RF(181) = RF(181)*C(2)
      RF(182) = RF(182)*C(2)
      RF(183) = RF(183)*C(5)
      RF(184) = RF(184)*C(5)
      RB(1) = RB(1)*CTB(1)*C(4)
      RB(2) = RB(2)*CTB(2)*C(5)
      RB(3) = RB(3)*C(2)*C(5)
      RB(4) = RB(4)*C(5)*C(4)
      RB(5) = RB(5)*C(5)*C(7)
      RB(6) = RB(6)*C(2)*C(11)
      RB(7) = RB(7)*C(2)
      RB(8) = RB(8)*C(1)*C(11)
      RB(9) = RB(9)*C(2)
      RB(10) = RB(10)*C(2)*C(13)
      RB(11) = RB(11)*C(5)*C(9)
      RB(12) = RB(12)*C(12)
      RB(13) = RB(13)*C(5)*C(11)
      RB(14) = RB(14)*C(2)*C(12)
      RB(15) = RB(15)*C(5)
      RB(16) = RB(16)*C(5)*C(13)
      RB(17) = RB(17)*C(5)*C(13)
      RB(18) = RB(18)*C(5)
      RB(19) = RB(19)*C(5)
      RB(20) = RB(20)*C(2)
      RB(21) = RB(21)*C(11)
      RB(22) = RB(22)*C(2)*C(18)
      RB(23) = RB(23)*C(9)
      RB(24) = RB(24)*C(9)*C(13)
      RB(25) = RB(25)*C(5)
      RB(26) = RB(26)*C(2)*C(11)*C(11)
      RB(27) = RB(27)*C(5)
      RB(28) = RB(28)*C(12)
      RB(29) = RB(29)*C(3)*C(12)
      RB(30) = RB(30)*C(7)
      RB(31) = RB(31)*CTB(31)*C(7)
      RB(32) = RB(32)*C(7)*C(4)
      RB(33) = RB(33)*C(7)*C(6)
      RB(34) = RB(34)*C(7)*C(19)
      RB(35) = RB(35)*C(3)*C(5)
      RB(36) = RB(36)*CTB(36)*C(1)
      RB(37) = RB(37)*C(1)*C(1)
      RB(38) = RB(38)*C(1)*C(6)
      RB(39) = RB(39)*C(1)*C(12)
      RB(40) = RB(40)*CTB(40)*C(6)
      RB(41) = RB(41)*C(3)*C(6)
      RB(42) = RB(42)*C(4)*C(1)
      RB(43) = RB(43)*C(5)*C(5)
      RB(44) = RB(44)*C(7)*C(1)
      RB(45) = RB(45)*C(5)*C(6)
      RB(46) = RB(46)*C(1)
      RB(47) = RB(47)*C(9)
      RB(48) = RB(48)*C(1)
      RB(49) = RB(49)*C(10)
      RB(50) = RB(50)*C(9)*C(1)
      RB(51) = RB(51)*C(13)
      RB(52) = RB(52)*C(1)*C(11)
      RB(55) = RB(55)*C(1)
      RB(56) = RB(56)*C(14)
      RB(57) = RB(57)*C(1)*C(13)
      RB(58) = RB(58)*C(5)*C(9)
      RB(59) = RB(59)*C(6)
      RB(60) = RB(60)*C(14)
      RB(61) = RB(61)*C(2)
      RB(62) = RB(62)*C(1)*C(13)
      RB(63) = RB(63)*C(5)*C(9)
      RB(64) = RB(64)*C(6)
      RB(65) = RB(65)*C(1)
      RB(66) = RB(66)*C(1)
      RB(68) = RB(68)*C(16)
      RB(69) = RB(69)*C(1)*C(15)
      RB(71) = RB(71)*C(1)
      RB(72) = RB(72)*C(17)
      RB(73) = RB(73)*C(1)*C(16)
      RB(74) = RB(74)*C(1)
      RB(75) = RB(75)*C(11)
      RB(76) = RB(76)*C(1)
      RB(77) = RB(77)*C(9)*C(11)
      RB(78) = RB(78)*C(13)
      RB(79) = RB(79)*C(2)*C(6)
      RB(80) = RB(80)*C(8)
      RB(81) = RB(81)*C(3)*C(6)
      RB(82) = RB(82)*C(4)*C(6)
      RB(83) = RB(83)*C(7)*C(6)
      RB(84) = RB(84)*C(7)*C(6)
      RB(85) = RB(85)*C(2)*C(11)
      RB(86) = RB(86)*C(2)
      RB(87) = RB(87)*C(2)*C(13)
      RB(88) = RB(88)*C(6)
      RB(89) = RB(89)*C(2)*C(13)
      RB(90) = RB(90)*C(14)
      RB(91) = RB(91)*C(6)
      RB(92) = RB(92)*C(6)
      RB(93) = RB(93)*C(9)*C(6)
      RB(94) = RB(94)*C(2)*C(12)
      RB(95) = RB(95)*C(6)*C(11)
      RB(96) = RB(96)*C(6)
      RB(97) = RB(97)*C(6)*C(13)
      RB(98) = RB(98)*C(6)*C(13)
      RB(99) = RB(99)*C(6)
      RB(100) = RB(100)*C(6)
      RB(101) = RB(101)*C(2)*C(18)
      RB(102) = RB(102)*C(9)*C(11)
      RB(103) = RB(103)*C(6)*C(15)
      RB(104) = RB(104)*C(6)
      RB(105) = RB(105)*C(6)
      RB(106) = RB(106)*C(6)
      RB(107) = RB(107)*C(4)*C(8)
      RB(108) = RB(108)*C(4)*C(8)
      RB(109) = RB(109)*C(5)*C(13)
      RB(110) = RB(110)*C(4)*C(10)
      RB(111) = RB(111)*C(5)
      RB(112) = RB(112)*C(5)*C(12)
      RB(113) = RB(113)*C(8)
      RB(114) = RB(114)*C(3)*C(11)
      RB(115) = RB(115)*C(2)*C(15)
      RB(116) = RB(116)*C(3)
      RB(117) = RB(117)*C(2)
      RB(118) = RB(118)*C(2)*C(13)
      RB(119) = RB(119)*C(2)*C(15)
      RB(120) = RB(120)*C(2)
      RB(121) = RB(121)*C(2)*C(16)
      RB(123) = RB(123)*C(11)
      RB(124) = RB(124)*C(2)*C(18)
      RB(125) = RB(125)*C(11)*C(15)
      RB(127) = RB(127)*C(2)*C(9)
      RB(128) = RB(128)*C(1)*C(15)
      RB(129) = RB(129)*C(2)*C(16)
      RB(130) = RB(130)*C(9)*C(9)
      RB(131) = RB(131)*C(18)
      RB(132) = RB(132)*C(11)
      RB(133) = RB(133)*C(19)
      RB(134) = RB(134)*C(2)*C(5)*C(11)
      RB(135) = RB(135)*C(11)*C(6)
      RB(136) = RB(136)*C(9)*C(2)
      RB(137) = RB(137)*C(14)
      RB(138) = RB(138)*C(6)
      RB(139) = RB(139)*C(2)*C(16)
      RB(140) = RB(140)*C(9)*C(9)
      RB(141) = RB(141)*C(11)
      RB(142) = RB(142)*C(12)
      RB(143) = RB(143)*C(11)*C(13)
      RB(144) = RB(144)*C(9)
      RB(145) = RB(145)*C(3)
      RB(146) = RB(146)*C(5)*C(13)
      RB(147) = RB(147)*C(7)*C(10)
      RB(148) = RB(148)*C(17)
      RB(149) = RB(149)*C(2)
      RB(150) = RB(150)*C(10)*C(11)
      RB(151) = RB(151)*C(10)
      RB(152) = RB(152)*C(10)
      RB(153) = RB(153)*C(10)
      RB(154) = RB(154)*C(10)
      RB(155) = RB(155)*C(10)
      RB(156) = RB(156)*C(2)*C(11)*C(6)
      RB(157) = RB(157)*CTB(157)*C(2)*C(11)
      RB(158) = RB(158)*C(7)*C(11)
      RB(159) = RB(159)*C(7)*C(13)
      RB(160) = RB(160)*C(7)*C(13)
      RB(161) = RB(161)*C(13)
      RB(162) = RB(162)*C(1)*C(15)
      RB(163) = RB(163)*C(7)*C(16)
      RB(164) = RB(164)*C(5)*C(11)*C(11)
      RB(165) = RB(165)*C(11)*C(11)*C(15)
      RB(167) = RB(167)*C(2)
      RB(168) = RB(168)*C(4)*C(6)
      RB(170) = RB(170)*C(9)
      RB(172) = RB(172)*C(3)*C(13)
      RB(175) = RB(175)*C(3)
      RB(176) = RB(176)*C(7)*C(15)
      RB(181) = RB(181)*C(9)
      RB(182) = RB(182)*C(18)*C(1)
      RB(183) = RB(183)*C(6)*C(18)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE QSSA(RF, RB, XQ)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (SMALL = 1.D-200)
      DIMENSION RF(*), RB(*), XQ(*)
C
      RF(119) = 0.D0
      RF(125) = 0.D0
      RF(128) = 0.D0
      RF(132) = 0.D0
      RF(165) = 0.D0
      RF(173) = 0.D0
      RB(180) = 0.D0
      RB(184) = 0.D0
C
C     C
      DEN = +RF( 85) +RF(114) +RF(115) +RB( 46) 
      A1_0 = ( +RB( 85) +RB(114) +RB(115) )/MAX(DEN, SMALL)
      A1_2 = ( +RF( 46) )/MAX(DEN, SMALL)
C     CH
      DEN = +RF(  6) +RF( 46) +RF( 86) +RF(116) +RF(117) 
     *  +RF(118) +RF(120) +RF(121) +RF(122) +RF(123) +RF(124) 
     *  +RF(170) +RB( 48) +RB( 88) 
      A2_0 = ( +RB(  6) +RB(118) +RB(119) +RB(121) +RB(124) 
     *  +RB(125) +RB(170) )/MAX(DEN, SMALL)
      A2_1 = ( +RB( 46) )/MAX(DEN, SMALL)
      A2_3 = ( +RF( 88) +RB(117) )/MAX(DEN, SMALL)
      A2_4 = ( +RF( 48) )/MAX(DEN, SMALL)
      A2_5 = ( +RB( 86) +RB(116) +RB(123) )/MAX(DEN, SMALL)
      A2_8 = ( +RB(120) )/MAX(DEN, SMALL)
      A2_10 = ( +RB(122) )/MAX(DEN, SMALL)
C     CH2
      DEN = +RF(  7) +RF( 47) +RF( 87) +RF( 88) +RF(109) 
     *  +RF(126) +RF(127) +RF(129) +RF(130) +RF(131) +RF(171) 
     *  +RF(172) +RB( 21) +RB( 28) +RB( 91) +RB(117) +RB(133) 
     *  +RB(138) +RB(141) +RB(142) +RB(178) 
      A3_0 = ( +RF( 21) +RF( 28) +RB( 47) +RB( 87) +RF( 91) 
     *  +RB(109) +RB(119) +RB(126) +RB(127) +RB(128) +RB(128) 
     *  +RB(129) +RB(130) +RB(131) +RB(171) +RB(172) +RB(173) 
     *  +RB(173) )/MAX(DEN, SMALL)
      A3_2 = ( +RB( 88) +RF(117) )/MAX(DEN, SMALL)
      A3_4 = ( +RF(133) +RF(138) +RF(141) +RF(142) )/MAX(DEN, SMALL)
      A3_5 = ( +RB(  7) )/MAX(DEN, SMALL)
      A3_8 = ( +RB(132) )/MAX(DEN, SMALL)
      A3_11 = ( +RF(178) )/MAX(DEN, SMALL)
C     CH2(S)
      DEN = +RF(  8) +RF(  9) +RF( 48) +RF( 89) +RF(133) 
     *  +RF(134) +RF(135) +RF(136) +RF(137) +RF(138) +RF(139) 
     *  +RF(140) +RF(141) +RF(142) +RF(143) +RF(144) +RF(174) 
     *  +RB( 59) +RB( 64) +RB( 75) +RB( 92) 
      A4_0 = ( +RB(  8) +RB( 89) +RF( 92) +RB(134) +RB(135) 
     *  +RB(136) +RB(137) +RB(139) +RB(140) +RB(143) +RB(174) )
     *  /MAX(DEN, SMALL)
      A4_2 = ( +RB( 48) )/MAX(DEN, SMALL)
      A4_3 = ( +RB(133) +RB(138) +RB(141) +RB(142) )/MAX(DEN, SMALL)
      A4_5 = ( +RB(  9) )/MAX(DEN, SMALL)
      A4_6 = ( +RF( 59) )/MAX(DEN, SMALL)
      A4_7 = ( +RF( 64) )/MAX(DEN, SMALL)
      A4_9 = ( +RB(144) )/MAX(DEN, SMALL)
      A4_10 = ( +RF( 75) )/MAX(DEN, SMALL)
C     HCO
      DEN = +RF( 13) +RF( 14) +RF( 51) +RF( 52) +RF( 95) 
     *  +RF(150) +RF(156) +RF(157) +RF(158) +RB(  7) +RB(  9) 
     *  +RB( 15) +RB( 23) +RB( 30) +RB( 55) +RB( 86) +RB( 96) 
     *  +RB(113) +RB(116) +RB(123) +RB(151) +RB(161) +RB(181) 
      A5_0 = ( +RB( 13) +RB( 14) +RF( 15) +RF( 23) +RF( 30) 
     *  +RB( 51) +RB( 52) +RF( 55) +RB( 95) +RF( 96) +RF(113) 
     *  +RB(150) +RF(151) +RB(156) +RB(157) +RB(158) )/MAX(DEN, SMALL)
      A5_2 = ( +RF( 86) +RF(116) +RF(123) )/MAX(DEN, SMALL)
      A5_3 = ( +RF(  7) )/MAX(DEN, SMALL)
      A5_4 = ( +RF(  9) )/MAX(DEN, SMALL)
      A5_8 = ( +RF(161) )/MAX(DEN, SMALL)
      A5_11 = ( +RF(180) +RF(180) +RF(181) +RF(184) )/MAX(DEN, SMALL)
C     CH2OH
      DEN = +RF( 16) +RF( 56) +RF( 57) +RF( 58) +RF( 59) 
     *  +RF( 97) +RF(159) +RB( 18) +RB( 53) +RB( 61) +RB( 65) 
     *  +RB( 99) +RB(152) 
      A6_0 = ( +RB( 16) +RF( 18) +RF( 53) +RB( 56) +RB( 57) 
     *  +RB( 58) +RF( 65) +RB( 97) +RF( 99) +RF(152) +RB(159) )
     *  /MAX(DEN, SMALL)
      A6_4 = ( +RB( 59) )/MAX(DEN, SMALL)
      A6_7 = ( +RF( 61) )/MAX(DEN, SMALL)
      A6_11 = ( +RF(184) )/MAX(DEN, SMALL)
C     CH3O
      DEN = +RF( 17) +RF( 60) +RF( 61) +RF( 62) +RF( 63) 
     *  +RF( 64) +RF( 98) +RF(160) +RB( 19) +RB( 54) +RB( 66) 
     *  +RB(100) +RB(111) +RB(145) +RB(153) 
      A7_0 = ( +RB( 17) +RF( 19) +RF( 54) +RB( 60) +RB( 62) 
     *  +RB( 63) +RF( 66) +RB( 98) +RF(100) +RF(111) +RF(145) 
     *  +RF(153) +RB(160) )/MAX(DEN, SMALL)
      A7_4 = ( +RB( 64) )/MAX(DEN, SMALL)
      A7_6 = ( +RB( 61) )/MAX(DEN, SMALL)
C     C2H3
      DEN = +RF( 22) +RF( 68) +RF( 69) +RF(103) +RF(161) 
     *  +RF(175) +RF(176) +RB( 67) +RB( 71) +RB(104) +RB(120) 
     *  +RB(132) +RB(154) 
      A8_0 = ( +RB( 22) +RF( 67) +RB( 68) +RB( 69) +RF( 71) 
     *  +RB(103) +RF(104) +RF(154) +RB(176) )/MAX(DEN, SMALL)
      A8_2 = ( +RF(120) )/MAX(DEN, SMALL)
      A8_5 = ( +RB(161) )/MAX(DEN, SMALL)
      A8_11 = ( +RB(175) )/MAX(DEN, SMALL)
C     C2H5
      DEN = +RF( 24) +RF( 72) +RF( 73) +RF(163) +RB( 25) 
     *  +RB( 70) +RB( 74) +RB(105) +RB(144) +RB(149) +RB(155) 
      A9_0 = ( +RB( 24) +RF( 25) +RF( 70) +RB( 72) +RB( 73) 
     *  +RF( 74) +RF(105) +RF(149) +RF(155) +RB(163) )/MAX(DEN, SMALL)
      A9_4 = ( +RF(144) )/MAX(DEN, SMALL)
C     HCCO
      DEN = +RF( 26) +RF( 75) +RF(164) +RB( 20) +RB( 27) 
     *  +RB( 76) +RB(106) +RB(122) 
      A10_0 = ( +RF( 20) +RB( 26) +RF( 27) +RF( 76) +RF(106) 
     *  +RB(125) +RB(164) +RB(165) +RB(165) )/MAX(DEN, SMALL)
      A10_2 = ( +RF(122) )/MAX(DEN, SMALL)
      A10_4 = ( +RB( 75) )/MAX(DEN, SMALL)
      A10_8 = ( +RB(132) )/MAX(DEN, SMALL)
C     CH2CHO
      DEN = +RF(178) +RF(179) +RF(180) +RF(181) +RF(182) 
     *  +RF(183) +RF(184) +RB(167) +RB(175) +RB(177) 
      A11_0 = ( +RF(167) +RF(177) +RB(179) +RB(182) +RB(183) )
     *        /MAX(DEN, SMALL)
      A11_5 = ( +RB(181) )/MAX(DEN, SMALL)
      A11_8 = ( +RF(175) )/MAX(DEN, SMALL)
C
      A2_0 = A2_0 + A2_1*A1_0
      DEN = 1 -A2_1*A1_2
      A2_0 = A2_0/MAX(DEN, SMALL)
      A2_4 = A2_4/MAX(DEN, SMALL)
      A2_5 = A2_5/MAX(DEN, SMALL)
      A2_3 = A2_3/MAX(DEN, SMALL)
      A2_10 = A2_10/MAX(DEN, SMALL)
      A2_8 = A2_8/MAX(DEN, SMALL)
      A4_0 = A4_0 + A4_9*A9_0
      DEN = 1 -A4_9*A9_4
      A4_0 = A4_0/MAX(DEN, SMALL)
      A4_2 = A4_2/MAX(DEN, SMALL)
      A4_5 = A4_5/MAX(DEN, SMALL)
      A4_3 = A4_3/MAX(DEN, SMALL)
      A4_6 = A4_6/MAX(DEN, SMALL)
      A4_10 = A4_10/MAX(DEN, SMALL)
      A4_7 = A4_7/MAX(DEN, SMALL)
      A5_0 = A5_0 + A5_11*A11_0
      A5_8 = A5_8 + A5_11*A11_8
      DEN = 1 -A5_11*A11_5
      A5_0 = A5_0/MAX(DEN, SMALL)
      A5_4 = A5_4/MAX(DEN, SMALL)
      A5_2 = A5_2/MAX(DEN, SMALL)
      A5_3 = A5_3/MAX(DEN, SMALL)
      A5_8 = A5_8/MAX(DEN, SMALL)
      A3_0 = A3_0 + A3_11*A11_0
      A3_5 = A3_5 + A3_11*A11_5
      A3_8 = A3_8 + A3_11*A11_8
      A6_0 = A6_0 + A6_11*A11_0
      A6_5 = A6_11*A11_5
      A6_8 = A6_11*A11_8
      A8_0 = A8_0 + A8_11*A11_0
      A8_5 = A8_5 + A8_11*A11_5
      DEN = 1 -A8_11*A11_8
      A8_0 = A8_0/MAX(DEN, SMALL)
      A8_2 = A8_2/MAX(DEN, SMALL)
      A8_5 = A8_5/MAX(DEN, SMALL)
      A2_0 = A2_0 + A2_8*A8_0
      A2_5 = A2_5 + A2_8*A8_5
      DEN = 1 -A2_8*A8_2
      A2_0 = A2_0/MAX(DEN, SMALL)
      A2_4 = A2_4/MAX(DEN, SMALL)
      A2_5 = A2_5/MAX(DEN, SMALL)
      A2_3 = A2_3/MAX(DEN, SMALL)
      A2_10 = A2_10/MAX(DEN, SMALL)
      A5_0 = A5_0 + A5_8*A8_0
      A5_2 = A5_2 + A5_8*A8_2
      DEN = 1 -A5_8*A8_5
      A5_0 = A5_0/MAX(DEN, SMALL)
      A5_4 = A5_4/MAX(DEN, SMALL)
      A5_2 = A5_2/MAX(DEN, SMALL)
      A5_3 = A5_3/MAX(DEN, SMALL)
      A3_0 = A3_0 + A3_8*A8_0
      A3_2 = A3_2 + A3_8*A8_2
      A3_5 = A3_5 + A3_8*A8_5
      A6_0 = A6_0 + A6_8*A8_0
      A6_2 = A6_8*A8_2
      A6_5 = A6_5 + A6_8*A8_5
      A10_0 = A10_0 + A10_8*A8_0
      A10_2 = A10_2 + A10_8*A8_2
      A10_5 = A10_8*A8_5
      A4_0 = A4_0 + A4_7*A7_0
      A4_6 = A4_6 + A4_7*A7_6
      DEN = 1 -A4_7*A7_4
      A4_0 = A4_0/MAX(DEN, SMALL)
      A4_2 = A4_2/MAX(DEN, SMALL)
      A4_5 = A4_5/MAX(DEN, SMALL)
      A4_3 = A4_3/MAX(DEN, SMALL)
      A4_6 = A4_6/MAX(DEN, SMALL)
      A4_10 = A4_10/MAX(DEN, SMALL)
      A6_0 = A6_0 + A6_7*A7_0
      A6_4 = A6_4 + A6_7*A7_4
      DEN = 1 -A6_7*A7_6
      A6_0 = A6_0/MAX(DEN, SMALL)
      A6_4 = A6_4/MAX(DEN, SMALL)
      A6_2 = A6_2/MAX(DEN, SMALL)
      A6_5 = A6_5/MAX(DEN, SMALL)
      A4_0 = A4_0 + A4_10*A10_0
      A4_2 = A4_2 + A4_10*A10_2
      A4_5 = A4_5 + A4_10*A10_5
      DEN = 1 -A4_10*A10_4
      A4_0 = A4_0/MAX(DEN, SMALL)
      A4_2 = A4_2/MAX(DEN, SMALL)
      A4_5 = A4_5/MAX(DEN, SMALL)
      A4_3 = A4_3/MAX(DEN, SMALL)
      A4_6 = A4_6/MAX(DEN, SMALL)
      A2_0 = A2_0 + A2_10*A10_0
      A2_4 = A2_4 + A2_10*A10_4
      A2_5 = A2_5 + A2_10*A10_5
      DEN = 1 -A2_10*A10_2
      A2_0 = A2_0/MAX(DEN, SMALL)
      A2_4 = A2_4/MAX(DEN, SMALL)
      A2_5 = A2_5/MAX(DEN, SMALL)
      A2_3 = A2_3/MAX(DEN, SMALL)
      A4_0 = A4_0 + A4_6*A6_0
      A4_2 = A4_2 + A4_6*A6_2
      A4_5 = A4_5 + A4_6*A6_5
      DEN = 1 -A4_6*A6_4
      A4_0 = A4_0/MAX(DEN, SMALL)
      A4_2 = A4_2/MAX(DEN, SMALL)
      A4_5 = A4_5/MAX(DEN, SMALL)
      A4_3 = A4_3/MAX(DEN, SMALL)
      A4_0 = A4_0 + A4_3*A3_0
      A4_2 = A4_2 + A4_3*A3_2
      A4_5 = A4_5 + A4_3*A3_5
      DEN = 1 -A4_3*A3_4
      A4_0 = A4_0/MAX(DEN, SMALL)
      A4_2 = A4_2/MAX(DEN, SMALL)
      A4_5 = A4_5/MAX(DEN, SMALL)
      A2_0 = A2_0 + A2_3*A3_0
      A2_4 = A2_4 + A2_3*A3_4
      A2_5 = A2_5 + A2_3*A3_5
      DEN = 1 -A2_3*A3_2
      A2_0 = A2_0/MAX(DEN, SMALL)
      A2_4 = A2_4/MAX(DEN, SMALL)
      A2_5 = A2_5/MAX(DEN, SMALL)
      A5_0 = A5_0 + A5_3*A3_0
      A5_4 = A5_4 + A5_3*A3_4
      A5_2 = A5_2 + A5_3*A3_2
      DEN = 1 -A5_3*A3_5
      A5_0 = A5_0/MAX(DEN, SMALL)
      A5_4 = A5_4/MAX(DEN, SMALL)
      A5_2 = A5_2/MAX(DEN, SMALL)
      A4_0 = A4_0 + A4_5*A5_0
      A4_2 = A4_2 + A4_5*A5_2
      DEN = 1 -A4_5*A5_4
      A4_0 = A4_0/MAX(DEN, SMALL)
      A4_2 = A4_2/MAX(DEN, SMALL)
      A2_0 = A2_0 + A2_5*A5_0
      A2_4 = A2_4 + A2_5*A5_4
      DEN = 1 -A2_5*A5_2
      A2_0 = A2_0/MAX(DEN, SMALL)
      A2_4 = A2_4/MAX(DEN, SMALL)
      A4_0 = A4_0 + A4_2*A2_0
      DEN = 1 -A4_2*A2_4
      A4_0 = A4_0/MAX(DEN, SMALL)
      XQ(4) = A4_0
      XQ(2) = A2_0 +A2_4*XQ(4)
      XQ(5) = A5_0 +A5_4*XQ(4) +A5_2*XQ(2)
      XQ(3) = A3_0 +A3_4*XQ(4) +A3_2*XQ(2) +A3_5*XQ(5)
      XQ(6) = A6_0 +A6_4*XQ(4) +A6_2*XQ(2) +A6_5*XQ(5)
      XQ(10) = A10_0 +A10_4*XQ(4) +A10_2*XQ(2) +A10_5*XQ(5)
      XQ(7) = A7_0 +A7_4*XQ(4) +A7_6*XQ(6)
      XQ(8) = A8_0 +A8_2*XQ(2) +A8_5*XQ(5)
      XQ(11) = A11_0 +A11_5*XQ(5) +A11_8*XQ(8)
      XQ(9) = A9_0 +A9_4*XQ(4)
      XQ(1) = A1_0 +A1_2*XQ(2)
C
      RF(  6) = RF(  6)*XQ( 2)
      RF(  7) = RF(  7)*XQ( 3)
      RB(  7) = RB(  7)*XQ( 5)
      RF(  8) = RF(  8)*XQ( 4)
      RF(  9) = RF(  9)*XQ( 4)
      RB(  9) = RB(  9)*XQ( 5)
      RF( 13) = RF( 13)*XQ( 5)
      RF( 14) = RF( 14)*XQ( 5)
      RB( 15) = RB( 15)*XQ( 5)
      RF( 16) = RF( 16)*XQ( 6)
      RF( 17) = RF( 17)*XQ( 7)
      RB( 18) = RB( 18)*XQ( 6)
      RB( 19) = RB( 19)*XQ( 7)
      RB( 20) = RB( 20)*XQ(10)
      RB( 21) = RB( 21)*XQ( 3)
      RF( 22) = RF( 22)*XQ( 8)
      RB( 23) = RB( 23)*XQ( 5)
      RF( 24) = RF( 24)*XQ( 9)
      RB( 25) = RB( 25)*XQ( 9)
      RF( 26) = RF( 26)*XQ(10)
      RB( 27) = RB( 27)*XQ(10)
      RB( 28) = RB( 28)*XQ( 3)
      RB( 30) = RB( 30)*XQ( 5)
      RF( 46) = RF( 46)*XQ( 2)
      RB( 46) = RB( 46)*XQ( 1)
      RF( 47) = RF( 47)*XQ( 3)
      RF( 48) = RF( 48)*XQ( 4)
      RB( 48) = RB( 48)*XQ( 2)
      RF( 51) = RF( 51)*XQ( 5)
      RF( 52) = RF( 52)*XQ( 5)
      RB( 53) = RB( 53)*XQ( 6)
      RB( 54) = RB( 54)*XQ( 7)
      RB( 55) = RB( 55)*XQ( 5)
      RF( 56) = RF( 56)*XQ( 6)
      RF( 57) = RF( 57)*XQ( 6)
      RF( 58) = RF( 58)*XQ( 6)
      RF( 59) = RF( 59)*XQ( 6)
      RB( 59) = RB( 59)*XQ( 4)
      RF( 60) = RF( 60)*XQ( 7)
      RF( 61) = RF( 61)*XQ( 7)
      RB( 61) = RB( 61)*XQ( 6)
      RF( 62) = RF( 62)*XQ( 7)
      RF( 63) = RF( 63)*XQ( 7)
      RF( 64) = RF( 64)*XQ( 7)
      RB( 64) = RB( 64)*XQ( 4)
      RB( 65) = RB( 65)*XQ( 6)
      RB( 66) = RB( 66)*XQ( 7)
      RB( 67) = RB( 67)*XQ( 8)
      RF( 68) = RF( 68)*XQ( 8)
      RF( 69) = RF( 69)*XQ( 8)
      RB( 70) = RB( 70)*XQ( 9)
      RB( 71) = RB( 71)*XQ( 8)
      RF( 72) = RF( 72)*XQ( 9)
      RF( 73) = RF( 73)*XQ( 9)
      RB( 74) = RB( 74)*XQ( 9)
      RF( 75) = RF( 75)*XQ(10)
      RB( 75) = RB( 75)*XQ( 4)
      RB( 76) = RB( 76)*XQ(10)
      RF( 85) = RF( 85)*XQ( 1)
      RF( 86) = RF( 86)*XQ( 2)
      RB( 86) = RB( 86)*XQ( 5)
      RF( 87) = RF( 87)*XQ( 3)
      RF( 88) = RF( 88)*XQ( 3)
      RB( 88) = RB( 88)*XQ( 2)
      RF( 89) = RF( 89)*XQ( 4)
      RB( 91) = RB( 91)*XQ( 3)
      RB( 92) = RB( 92)*XQ( 4)
      RF( 95) = RF( 95)*XQ( 5)
      RB( 96) = RB( 96)*XQ( 5)
      RF( 97) = RF( 97)*XQ( 6)
      RF( 98) = RF( 98)*XQ( 7)
      RB( 99) = RB( 99)*XQ( 6)
      RB(100) = RB(100)*XQ( 7)
      RF(103) = RF(103)*XQ( 8)
      RB(104) = RB(104)*XQ( 8)
      RB(105) = RB(105)*XQ( 9)
      RB(106) = RB(106)*XQ(10)
      RF(109) = RF(109)*XQ( 3)
      RB(111) = RB(111)*XQ( 7)
      RB(113) = RB(113)*XQ( 5)
      RF(114) = RF(114)*XQ( 1)
      RF(115) = RF(115)*XQ( 1)
      RF(116) = RF(116)*XQ( 2)
      RB(116) = RB(116)*XQ( 5)
      RF(117) = RF(117)*XQ( 2)
      RB(117) = RB(117)*XQ( 3)
      RF(118) = RF(118)*XQ( 2)
      RF(120) = RF(120)*XQ( 2)
      RB(120) = RB(120)*XQ( 8)
      RF(121) = RF(121)*XQ( 2)
      RF(122) = RF(122)*XQ( 2)
      RB(122) = RB(122)*XQ(10)
      RF(123) = RF(123)*XQ( 2)
      RB(123) = RB(123)*XQ( 5)
      RF(124) = RF(124)*XQ( 2)
      RF(126) = RF(126)*XQ( 3)
      RF(127) = RF(127)*XQ( 3)
      RF(129) = RF(129)*XQ( 3)
      RF(130) = RF(130)*XQ( 3)
      RF(131) = RF(131)*XQ( 3)
      RB(132) = RB(132)*XQ( 8)
      RF(133) = RF(133)*XQ( 4)
      RB(133) = RB(133)*XQ( 3)
      RF(134) = RF(134)*XQ( 4)
      RF(135) = RF(135)*XQ( 4)
      RF(136) = RF(136)*XQ( 4)
      RF(137) = RF(137)*XQ( 4)
      RF(138) = RF(138)*XQ( 4)
      RB(138) = RB(138)*XQ( 3)
      RF(139) = RF(139)*XQ( 4)
      RF(140) = RF(140)*XQ( 4)
      RF(141) = RF(141)*XQ( 4)
      RB(141) = RB(141)*XQ( 3)
      RF(142) = RF(142)*XQ( 4)
      RB(142) = RB(142)*XQ( 3)
      RF(143) = RF(143)*XQ( 4)
      RF(144) = RF(144)*XQ( 4)
      RB(144) = RB(144)*XQ( 9)
      RB(145) = RB(145)*XQ( 7)
      RB(149) = RB(149)*XQ( 9)
      RF(150) = RF(150)*XQ( 5)
      RB(151) = RB(151)*XQ( 5)
      RB(152) = RB(152)*XQ( 6)
      RB(153) = RB(153)*XQ( 7)
      RB(154) = RB(154)*XQ( 8)
      RB(155) = RB(155)*XQ( 9)
      RF(156) = RF(156)*XQ( 5)
      RF(157) = RF(157)*XQ( 5)
      RF(158) = RF(158)*XQ( 5)
      RF(159) = RF(159)*XQ( 6)
      RF(160) = RF(160)*XQ( 7)
      RF(161) = RF(161)*XQ( 8)
      RB(161) = RB(161)*XQ( 5)
      RF(163) = RF(163)*XQ( 9)
      RF(164) = RF(164)*XQ(10)
      RB(167) = RB(167)*XQ(11)
      RF(170) = RF(170)*XQ( 2)
      RF(171) = RF(171)*XQ( 3)
      RF(172) = RF(172)*XQ( 3)
      RF(174) = RF(174)*XQ( 4)
      RF(175) = RF(175)*XQ( 8)
      RB(175) = RB(175)*XQ(11)
      RF(176) = RF(176)*XQ( 8)
      RB(177) = RB(177)*XQ(11)
      RF(178) = RF(178)*XQ(11)
      RB(178) = RB(178)*XQ( 3)
      RF(179) = RF(179)*XQ(11)
      RF(180) = RF(180)*XQ(11)
      RF(181) = RF(181)*XQ(11)
      RB(181) = RB(181)*XQ( 5)
      RF(182) = RF(182)*XQ(11)
      RF(183) = RF(183)*XQ(11)
      RF(184) = RF(184)*XQ(11)
C
      END
C                                                                      C
C----------------------------------------------------------------------C
C                                                                      C
      SUBROUTINE RDOT(RF, RB, WDOT)
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      DIMENSION RF(*), RB(*), WDOT(*), ROP(184)
C
      DO I = 1, 184
         ROP(I) = RF(I) - RB(I)
      ENDDO
C
C     H2
      WDOT( 1) =   -ROP(  3)   +ROP(  8)   +ROP( 36)   +ROP( 37) 
     *             +ROP( 38)   +ROP( 39)   +ROP( 42)   +ROP( 44) 
     *             +ROP( 46)   +ROP( 48)   +ROP( 50)   +ROP( 52) 
     *             +ROP( 55)   +ROP( 57)   +ROP( 62)   +ROP( 65) 
     *             +ROP( 66)   +ROP( 69)   +ROP( 71)   +ROP( 73) 
     *             +ROP( 74)   +ROP( 76)   -ROP( 78)   -ROP( 79) 
     *             -ROP(117)   -ROP(127)   +ROP(128)   -ROP(136) 
     *             +ROP(162)   +ROP(166)   +ROP(169)   -ROP(170) 
     *             +ROP(174)   +ROP(182) 
C     H
      WDOT( 2) =   -ROP(  2)   +ROP(  3)   +ROP(  6)   -ROP(  7) 
     *           -2*ROP(  8)   -ROP(  9)   +ROP( 10)   +ROP( 14) 
     *             -ROP( 17)   +ROP( 19)   +ROP( 20) +2*ROP( 21) 
     *             -ROP( 24)   +ROP( 25)   +ROP( 26) +2*ROP( 28) 
     *             -ROP( 31)   -ROP( 32)   -ROP( 33)   -ROP( 34) 
     *             -ROP( 35) -2*ROP( 36) -2*ROP( 37) -2*ROP( 38) 
     *           -2*ROP( 39)   -ROP( 40)   -ROP( 41)   -ROP( 42) 
     *             -ROP( 43)   -ROP( 44)   -ROP( 45)   -ROP( 46) 
     *           -3*ROP( 47) -3*ROP( 48)   -ROP( 49)   -ROP( 50) 
     *             -ROP( 51)   -ROP( 52)   -ROP( 53)   -ROP( 55) 
     *             -ROP( 56)   -ROP( 57)   -ROP( 58)   +ROP( 59) 
     *           -2*ROP( 60)   -ROP( 61) -2*ROP( 62) -2*ROP( 63) 
     *             -ROP( 65) -2*ROP( 68) -2*ROP( 69) -2*ROP( 72) 
     *           -2*ROP( 73)   +ROP( 75)   -ROP( 76)   -ROP( 77) 
     *             +ROP( 79)   +ROP( 85)   +ROP( 86)   -ROP( 87) 
     *           -2*ROP( 88)   -ROP( 89) +2*ROP( 91) +2*ROP( 92) 
     *             +ROP( 94)   -ROP( 98)   +ROP(100)   +ROP(101) 
     *             -ROP(103)   +ROP(104)   +ROP(105) -2*ROP(109) 
     *             +ROP(111)   +ROP(115) +3*ROP(117)   +ROP(118) 
     *             -ROP(119) +2*ROP(120)   +ROP(121)   +ROP(124) 
     *             -ROP(126)   -ROP(127) -4*ROP(128)   -ROP(129) 
     *           -2*ROP(130) -2*ROP(131)   -ROP(132)   -ROP(134) 
     *           -2*ROP(135)   -ROP(136) -2*ROP(137)   -ROP(139) 
     *           -2*ROP(140) -2*ROP(143)   -ROP(144)   +ROP(145) 
     *           +2*ROP(149)   +ROP(153)   +ROP(154)   +ROP(155) 
     *             +ROP(156)   +ROP(157)   -ROP(160)   -ROP(161) 
     *             -ROP(163)   +ROP(166) +2*ROP(167) -2*ROP(172) 
     *           -2*ROP(173) -2*ROP(174)   -ROP(176) +2*ROP(178) 
     *             -ROP(179)   -ROP(180) -2*ROP(181) -2*ROP(182) 
     *             -ROP(183)   -ROP(184) 
C     O
      WDOT( 3) = -2*ROP(  1)   -ROP(  2)   -ROP(  3)   -ROP(  4) 
     *             -ROP(  5) -2*ROP(  6)   -ROP(  7)   -ROP(  8) 
     *             -ROP(  9)   -ROP( 10)   -ROP( 11)   -ROP( 12) 
     *             -ROP( 13)   -ROP( 14)   -ROP( 15)   -ROP( 16) 
     *             -ROP( 17)   -ROP( 18)   -ROP( 19)   -ROP( 20) 
     *             -ROP( 21) -2*ROP( 22)   -ROP( 23)   -ROP( 24) 
     *             -ROP( 25)   -ROP( 26)   -ROP( 27)   -ROP( 28) 
     *             +ROP( 29)   +ROP( 35)   +ROP( 41)   +ROP( 48) 
     *             +ROP( 67)   -ROP( 68)   -ROP( 69)   +ROP( 71) 
     *             +ROP( 81)   -ROP( 85)   -ROP( 86)   +ROP( 88) 
     *             -ROP(103)   +ROP(104)   -ROP(115)   -ROP(117) 
     *             -ROP(118)   -ROP(119)   -ROP(121)   -ROP(122) 
     *             -ROP(123)   -ROP(124)   -ROP(125)   +ROP(132) 
     *             +ROP(145)   +ROP(154)   -ROP(161)   -ROP(166) 
     *             -ROP(167)   -ROP(170)   +ROP(172)   -ROP(176) 
     *             -ROP(178) 
C     O2
      WDOT( 4) =   +ROP(  1)   +ROP(  4) +2*ROP(  6)   +ROP(  8) 
     *             +ROP( 13)   +ROP( 14)   -ROP( 15)   +ROP( 16) 
     *             -ROP( 18)   -ROP( 20)   -ROP( 21)   +ROP( 22) 
     *             -ROP( 23)   +ROP( 26)   -ROP( 27)   -ROP( 28) 
     *             -ROP( 29) -2*ROP( 30)   -ROP( 31)   -ROP( 32) 
     *             -ROP( 33)   -ROP( 34)   -ROP( 35)   +ROP( 42) 
     *             +ROP( 46)   +ROP( 47)   -ROP( 48)   +ROP( 51) 
     *             +ROP( 52)   -ROP( 53)   -ROP( 55)   +ROP( 56) 
     *             +ROP( 57)   +ROP( 58)   -ROP( 61)   -ROP( 64) 
     *             -ROP( 65)   -ROP( 67)   +ROP( 68)   +ROP( 69) 
     *             -ROP( 71)   -ROP( 76)   +ROP( 82)   +ROP( 85) 
     *             +ROP( 86)   +ROP( 87)   -ROP( 88)   +ROP( 89) 
     *             -ROP( 91)   -ROP( 92)   +ROP( 95)   -ROP( 96) 
     *             +ROP( 97)   -ROP( 99)   +ROP(103)   -ROP(104) 
     *             -ROP(106)   +ROP(107)   +ROP(108)   +ROP(109) 
     *             +ROP(110)   -ROP(113)   +ROP(115)   +ROP(117) 
     *           +2*ROP(118) +3*ROP(119)   +ROP(120) +2*ROP(121) 
     *             +ROP(122)   +ROP(123) +2*ROP(124) +3*ROP(125) 
     *             +ROP(127) +2*ROP(128)   +ROP(129)   +ROP(130) 
     *             +ROP(131)   +ROP(132)   +ROP(136)   +ROP(137) 
     *             +ROP(139)   +ROP(140)   +ROP(143)   +ROP(144) 
     *             -ROP(145)   -ROP(146)   +ROP(150)   -ROP(151) 
     *             -ROP(152)   -ROP(154)   +ROP(156)   +ROP(157) 
     *             -ROP(160)   -ROP(161)   -ROP(163) +2*ROP(165) 
     *             +ROP(168) +2*ROP(170) +2*ROP(173)   +ROP(174) 
     *             -ROP(178)   -ROP(179) -3*ROP(180)   -ROP(181) 
     *           -2*ROP(184) 
C     OH
      WDOT( 5) =   +ROP(  2)   +ROP(  3)   +ROP(  4)   +ROP(  5) 
     *             +ROP( 11)   +ROP( 13)   +ROP( 15)   +ROP( 16) 
     *             +ROP( 17)   +ROP( 18)   +ROP( 19)   +ROP( 20) 
     *             +ROP( 25)   -ROP( 26) +2*ROP( 27)   +ROP( 35) 
     *             -ROP( 40) +2*ROP( 43)   +ROP( 45)   +ROP( 58) 
     *             +ROP( 63)   -ROP( 75)   +ROP( 76)   -ROP( 79) 
     *           -2*ROP( 80) -2*ROP( 81)   -ROP( 82)   -ROP( 83) 
     *             -ROP( 84)   -ROP( 85)   -ROP( 86)   -ROP( 87) 
     *             -ROP( 88)   -ROP( 89)   -ROP( 90)   -ROP( 91) 
     *             -ROP( 92)   -ROP( 93)   -ROP( 94)   -ROP( 95) 
     *             -ROP( 96)   -ROP( 97)   -ROP( 98)   -ROP( 99) 
     *             -ROP(100)   -ROP(101)   -ROP(102)   -ROP(103) 
     *             -ROP(104)   -ROP(105)   +ROP(109)   +ROP(111) 
     *             +ROP(112)   +ROP(122)   -ROP(125)   +ROP(126) 
     *             -ROP(132)   +ROP(134)   +ROP(146) -2*ROP(165) 
     *             -ROP(168)   -ROP(169)   +ROP(179)   +ROP(180) 
     *             -ROP(183)   -ROP(184) 
C     H2O
      WDOT( 6) =   +ROP( 40)   +ROP( 41)   +ROP( 45)   +ROP( 59) 
     *             +ROP( 64)   +ROP( 79)   +ROP( 81)   +ROP( 82) 
     *             +ROP( 83)   +ROP( 84)   +ROP( 88)   +ROP( 91) 
     *             +ROP( 92)   +ROP( 93)   +ROP( 95)   +ROP( 96) 
     *             +ROP( 97)   +ROP( 98)   +ROP( 99)   +ROP(100) 
     *             +ROP(103)   +ROP(104)   +ROP(105)   +ROP(106) 
     *             -ROP(118)   +ROP(135)   -ROP(137)   +ROP(168) 
     *             -ROP(174)   +ROP(183) 
C     HO2
      WDOT( 7) =   -ROP(  4)   +ROP(  5)   -ROP(  6)   +ROP(  7) 
     *             +ROP(  9)   -ROP( 13)   -ROP( 14)   +ROP( 15) 
     *             -ROP( 16)   +ROP( 18)   +ROP( 23) +2*ROP( 30) 
     *             +ROP( 31)   +ROP( 32)   +ROP( 33)   +ROP( 34) 
     *             -ROP( 41)   -ROP( 42)   -ROP( 43)   +ROP( 44) 
     *             -ROP( 46)   +ROP( 48)   -ROP( 51)   -ROP( 52) 
     *             +ROP( 53)   +ROP( 55)   -ROP( 56)   -ROP( 57) 
     *             -ROP( 58)   -ROP( 59)   +ROP( 61)   +ROP( 65) 
     *             -ROP( 82)   +ROP( 83)   +ROP( 84)   +ROP( 88) 
     *             -ROP( 95)   +ROP( 96)   -ROP( 97)   +ROP( 99) 
     *           -2*ROP(107) -2*ROP(108)   -ROP(109)   -ROP(110) 
     *             -ROP(111)   -ROP(112)   -ROP(117)   -ROP(118) 
     *             -ROP(119)   -ROP(120)   -ROP(121)   -ROP(122) 
     *             -ROP(124)   -ROP(125)   +ROP(147)   -ROP(150) 
     *             +ROP(151)   +ROP(152)   -ROP(156)   -ROP(157) 
     *             +ROP(160)   +ROP(161)   +ROP(163)   -ROP(168) 
     *             -ROP(170)   +ROP(176) +2*ROP(180)   +ROP(181) 
     *           +2*ROP(184) 
C     H2O2
      WDOT( 8) =   -ROP(  5)   -ROP( 44)   -ROP( 45)   +ROP( 80) 
     *             -ROP( 83)   -ROP( 84)   +ROP(107)   +ROP(108) 
     *             +ROP(113)   -ROP(147) 
C     CH3
      WDOT( 9) =   -ROP( 10)   +ROP( 11)   +ROP( 23)   +ROP( 24) 
     *             +ROP( 47)   -ROP( 49)   +ROP( 50)   +ROP( 58) 
     *             +ROP( 63)   +ROP( 77)   -ROP( 90)   -ROP( 91) 
     *             -ROP( 92)   +ROP( 93)   +ROP(102)   -ROP(110) 
     *             -ROP(111)   -ROP(115)   -ROP(120)   +ROP(127) 
     *             -ROP(129) +2*ROP(130)   +ROP(136)   -ROP(139) 
     *           +2*ROP(140)   +ROP(144)   -ROP(145)   -ROP(146) 
     *             -ROP(147) -2*ROP(148) -2*ROP(149)   -ROP(150) 
     *             -ROP(151)   -ROP(152)   -ROP(153)   -ROP(154) 
     *             -ROP(155)   -ROP(166)   -ROP(169)   +ROP(170) 
     *             +ROP(181) 
C     CH4
      WDOT(10) =   -ROP( 11)   +ROP( 49)   -ROP( 50)   -ROP( 93) 
     *             +ROP(110)   -ROP(121)   -ROP(130)   -ROP(140) 
     *             +ROP(147)   +ROP(150)   +ROP(151)   +ROP(152) 
     *             +ROP(153)   +ROP(154)   +ROP(155) 
C     CO
      WDOT(11) =   +ROP(  7)   +ROP(  8)   +ROP(  9)   -ROP( 12) 
     *             -ROP( 14)   +ROP( 15) +2*ROP( 20)   +ROP( 21) 
     *             +ROP( 23) +2*ROP( 27)   -ROP( 29)   +ROP( 30) 
     *             +ROP( 48)   -ROP( 51)   +ROP( 55)   -ROP( 75) 
     *           +2*ROP( 76)   +ROP( 77)   -ROP( 78)   +ROP( 88) 
     *             -ROP( 94)   +ROP( 96)   +ROP(102) +2*ROP(106) 
     *             -ROP(112)   +ROP(113)   -ROP(115)   -ROP(117) 
     *             -ROP(118)   -ROP(119)   -ROP(120)   -ROP(121) 
     *             +ROP(123)   -ROP(124) -2*ROP(125)   +ROP(126) 
     *             -ROP(131)   -ROP(132)   +ROP(134)   +ROP(135) 
     *             +ROP(143)   +ROP(151)   +ROP(161) -2*ROP(165) 
     *             +ROP(166)   -ROP(170)   +ROP(179) +2*ROP(180) 
     *             +ROP(181)   +ROP(184) 
C     CO2
      WDOT(12) =   -ROP(  7)   -ROP(  8)   -ROP(  9)   +ROP( 12) 
     *             +ROP( 14)   +ROP( 21) +2*ROP( 28)   +ROP( 29) 
     *             -ROP( 47)   -ROP( 48)   +ROP( 59)   +ROP( 64) 
     *             +ROP( 75)   -ROP( 87)   -ROP( 88)   -ROP( 89) 
     *             +ROP( 91)   +ROP( 92)   +ROP( 94)   -ROP(109) 
     *             +ROP(112)   +ROP(117)   -ROP(119)   -ROP(123) 
     *             -ROP(126)   -ROP(127) -2*ROP(128)   -ROP(129) 
     *             -ROP(130)   -ROP(131)   -ROP(132)   -ROP(134) 
     *             -ROP(135)   -ROP(136)   -ROP(137)   -ROP(139) 
     *             -ROP(140) -2*ROP(143)   -ROP(144)   -ROP(172) 
     *           -2*ROP(173)   -ROP(174) +2*ROP(178) 
C     CH2O
      WDOT(13) =   +ROP( 10)   -ROP( 15)   +ROP( 18)   +ROP( 19) 
     *             +ROP( 24)   -ROP( 30)   +ROP( 51)   -ROP( 55) 
     *             -ROP( 56)   -ROP( 58)   -ROP( 59)   -ROP( 60) 
     *             -ROP( 63)   -ROP( 64)   +ROP( 65)   +ROP( 66) 
     *             +ROP( 78)   +ROP( 87)   +ROP( 89)   -ROP( 96) 
     *             +ROP( 99)   +ROP(100)   +ROP(109)   +ROP(111) 
     *             -ROP(113)   +ROP(118)   -ROP(124)   +ROP(143) 
     *             +ROP(145)   +ROP(146)   -ROP(151)   +ROP(152) 
     *             +ROP(153)   +ROP(161)   +ROP(169)   +ROP(172) 
     *             +ROP(174)   +ROP(179)   +ROP(184) 
C     CH3OH
      WDOT(14) =   -ROP( 18)   -ROP( 19)   +ROP( 56)   +ROP( 60) 
     *             -ROP( 65)   -ROP( 66)   +ROP( 90)   -ROP( 99) 
     *             -ROP(100)   +ROP(137)   -ROP(152)   -ROP(153) 
C     C2H2
      WDOT(15) =   -ROP( 20)   -ROP( 21)   -ROP( 67)   +ROP( 69) 
     *             -ROP(101)   -ROP(102)   +ROP(103)   +ROP(115) 
     *             +ROP(119)   +ROP(125)   +ROP(128)   +ROP(162) 
     *             +ROP(165)   +ROP(173)   +ROP(176) 
C     C2H4
      WDOT(16) =   -ROP( 23)   -ROP( 24)   +ROP( 25)   +ROP( 68) 
     *             -ROP( 71)   -ROP( 72)   +ROP( 74)   -ROP(104) 
     *             +ROP(105)   +ROP(121)   +ROP(129)   +ROP(139) 
     *             +ROP(144)   +ROP(149)   -ROP(154)   +ROP(155) 
     *             -ROP(162)   -ROP(167) 
C     C2H6
      WDOT(17) =   -ROP( 25)   +ROP( 72)   -ROP( 74)   -ROP(105) 
     *             -ROP(144)   +ROP(148)   -ROP(155) 
C     CH2CO
      WDOT(18) =   -ROP( 27)   -ROP( 28)   +ROP( 67)   -ROP( 68) 
     *             -ROP( 69)   +ROP( 71)   -ROP( 76)   -ROP( 77) 
     *             +ROP(101)   -ROP(103)   +ROP(104)   -ROP(106) 
     *             +ROP(120)   +ROP(124)   +ROP(131)   +ROP(132) 
     *             +ROP(154)   -ROP(161)   +ROP(167)   -ROP(176) 
     *             -ROP(178)   -ROP(179)   -ROP(180)   -ROP(181) 
     *             -ROP(184) 
C     N2
      WDOT(19) =   0.0
C
      END
