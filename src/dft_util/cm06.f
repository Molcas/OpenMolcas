************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2010, Yan Zhao                                         *
************************************************************************
      Subroutine CM06(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                 CoeffA,iSpin,F_xc,T_X,ijzy)
************************************************************************
*                                                                      *
*  CM06 evaluates the correlation part of the M06 suite of             *
*  functionals on a grid.                                              *
*  !!! Second derivatives are not available yet.                       *
*                                                                      *
*  Ref: (a) Zhao, Y.  and Truhlar, D. G. J. Chem. Phys. 125,           *
*    194101 (2006).                                                    *
*       (b) Y. Zhao and D. G. Truhlar, J. Phys. Chem. A (2006),        *
*    110(49),13126-13130.                                              *
*                                                                      *
*       ijzy - 1 M06-L                                                 *
*       ijzy - 2 M06-HF                                                *
*       ijzy - 3 M06                                                   *
*       ijzy - 4 M06-2X                                                *
*                                                                      *
*  YZ (10/07)                                                          *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)
      Integer mGrid

      REAL*8 F6, F43, Pi34, F13,
     &RS,RSP,Zeta,dZdA,dZdB,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ,
     &  P, EUEG
      REAL*8 PA,GAA,TauA,FA,FPA,FGA,FTA,EUA,ChiA,EUPA,ChiAP,ChiAG
      REAL*8 PB,GBB,TauB,FB,FPB,FGB,FTB,EUB,ChiB,EUPB,ChiBP,ChiBG
      REAL*8 sopp0, sopp1,sopp2, sopp3, sopp4
      REAL*8 U, W, dUdChiA,dUdChiB,dUdPA,dUdPB,dUdGA,dUdGB,
     &dWdU,dWdPA,dWdPB, dWdGA,dWdGB,EUEGPA,EUEGPB

      REAL*8   F1, F2, F3, F4, COpp
      Data COpp/0.0031d0/,F1/1.0d0/,F2/2.0d0/,F3/3.0d0/,F4/4.0d0/

      INTEGER ijzy

      if (ijzy.eq.1) then
C     Parameters for M06-L Correlation
         sopp0= 6.042374D-01
         sopp1= 1.776783D+02
         sopp2= -2.513252D+02
         sopp3= 7.635173D+01
         sopp4= -1.255699D+01
      elseif (ijzy.eq.2) then
C     Parameters for M06-HF Correlation
         sopp0= 1.674634D+00
         sopp1= 5.732017D+01
         sopp2= 5.955416D+01
         sopp3= -2.311007D+02
         sopp4= 1.255199D+02
      elseif (ijzy.eq.3) then
C     Parameters for M06 Correlation
         sopp0= 3.741539D+00
         sopp1= 2.187098D+02
         sopp2= -4.531252D+02
         sopp3= 2.936479D+02
         sopp4= -6.287470D+01
C     elseif (ijzy.eq.4) then
      else
C     Parameters for M06-2X Correlation
         sopp0= 8.833596D-01
         sopp1= 3.357972D+01
         sopp2= -7.043548D+01
         sopp3= 4.978271D+01
         sopp4= -1.852891D+01
      endif
*
      F6=6.0d0
      F43 = F4 / F3
      Pi34 = F3 / (F4*Pi)
      F13 = F1 / F3
*
      Ta=0.5D0*T_X
      If (iSpin.eq.1) Then
         Do iGrid = 1, mGrid
            PA=max(1.0D-24,Rho(ipR,iGrid))
            If (PA.lt.Ta) Go To 110
            grdrhoa_x=Rho(ipdRx,iGrid)
            grdrhoa_y=Rho(ipdRy,iGrid)
            grdrhoa_z=Rho(ipdRz,iGrid)
            GAA=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2
            TauA = Rho(ipTau,iGrid)
            Call m06css(T_X,PA,GAA,TauA,FA,FPA,FGA,FTA,EUA,
     &                  ChiA,EUPA,ChiAP,ChiAG,ijzy)

            F_xc(iGrid)=F_xc(iGrid)+(2.0D0*FA)
*           dF/dRho
            dF_dRho(ipR,iGrid)=dF_dRho(ipR,iGrid)+FPA
*           dF/dGamma, no Gamma_ab term
            dF_dRho(ipGxx,iGrid)=dF_dRho(ipGxx,iGrid)+ FGA
*           dF/dTau
            dF_dRho(ipT,iGrid)=dF_dRho(ipT,iGrid)+FTA
*
            P=2.D0*PA
            RS = (Pi34/P) ** F13
            RSP = -RS/(F3*P)
            Zeta = 0.0D0
            dZdA = F1/P
            Call lsdac(RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     &                 d2LdZZ)
            EUEG = P*PotLC - 2.0D0*EUA
            U = COpp*(2.0D0*ChiA)/(F1 + COpp*(2.0D0*ChiA))
            W = sopp0+U*(sopp1+U*(sopp2+U*(sopp3+U*sopp4)))
            F_xc(igrid) = F_xc(iGrid) + EUEG*W
*
            dUdChiA =COpp/(F1 + COpp*(2.0D0*ChiA))**2
            dUdPA= dUdChiA*ChiAP
            dUdGA= dUdChiA*ChiAG
            dWdU =sopp1+U*(F2*sopp2+U*(F3*sopp3+U*F4*sopp4))
            dWdPA= dWdU*dUdPA
            dWdGA= dWdU*dUdGA
            EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
*           dF/dRho
            dF_dRho(ipR,iGrid)=dF_dRho(ipR,iGrid)+EUEGPA*W
     &                                             +EUEG*dWdPA
*           dF/dGamma
            dF_dRho(ipGxx,iGrid)=dF_dRho(ipGxx,iGrid)+EUEG*dWdGA
110         continue
         End Do
      Else
         Do iGrid = 1, mGrid
            PA=max(1.0D-24,Rho(ipRa,iGrid))
            If (PA.lt.Ta) Go To 100
            grdrhoa_x=Rho(ipdRxa,iGrid)
            grdrhoa_y=Rho(ipdRya,iGrid)
            grdrhoa_z=Rho(ipdRza,iGrid)
            GAA=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2
            TauA = Rho(ipTaua,iGrid)
            Call m06css(T_X,PA,GAA,TauA,FA,FPA,FGA,FTA,EUA,
     &                  ChiA,EUPA,ChiAP,ChiAG,ijzy)
*
            F_xc(iGrid)=F_xc(iGrid)+FA
*           dF/dRhoa
            dF_dRho(ipRa,igrid)=dF_dRho(ipRa,igrid)+FPA
*           dF/dGammaaa
            dF_dRho(ipGaa,igrid)=dF_dRho(ipGaa,iGrid)+ FGA
*           dF/dTaua
            dF_dRho(ipTa,igrid)=dF_dRho(ipTa,iGrid)+FTA

100         continue
            PB=max(1.0D-24,Rho(ipRb,iGrid))
            if(PB.lt.Ta) goto 111
            grdrhob_x=Rho(ipdRxb,iGrid)
            grdrhob_y=Rho(ipdRyb,iGrid)
            grdrhob_z=Rho(ipdRzb,iGrid)
            GBB=grdrhob_x**2+grdrhob_y**2+grdrhob_z**2
            TauB = Rho(ipTaub,iGrid)
            call m06css(T_X,PB,GBB,TauB,FB,FPB,FGB,FTB,EUB,
     &                  ChiB,EUPB,ChiBP,ChiBG,ijzy)

            F_xc(iGrid)=F_xc(iGrid)+FB
*           dF/dRhob
            dF_dRho(ipRb,igrid)=dF_dRho(ipRb,igrid)+FPB
*           dF/dGammabb
            dF_dRho(ipGbb,igrid)=dF_dRho(ipGbb,iGrid)+ FGB
*           dF/dTaub
            dF_dRho(ipTb,igrid)=dF_dRho(ipTb,iGrid)+FTB
111         continue
*
            P=PA+PB
            IF (PB.lt.T_X.or.PA.lt.T_X) Go To 112
            RS = (Pi34/P) ** F13
            RSP = -RS/(F3*P)
            Zeta = (PA-PB)/P
            dZdA = (F1-Zeta)/P
            dZdB = (-F1-Zeta)/P
            Call lsdac(RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $                 d2LdZZ)
            EUEG = P*PotLC - EUA - EUB
            U = COpp*(ChiA+ChiB)/(F1 + COpp*(ChiA+ChiB))
            W = sopp0+U*(sopp1+U*(sopp2+U*(sopp3+U*sopp4)))
            F_xc(iGrid) = F_xc(iGrid) + EUEG*W
            dUdChiA =COpp/(F1 + COpp*(ChiA+ChiB))**2
            dUdChiB =COpp/(F1 + COpp*(ChiA+ChiB))**2
            dUdPA= dUdChiA*ChiAP
            dUdPB= dUdChiB*ChiBP
            dUdGA= dUdChiA*ChiAG
            dUdGB= dUdChiB*ChiBG
            dWdU =sopp1+U*(F2*sopp2+U*(F3*sopp3+U*F4*sopp4))
            dWdPA= dWdU*dUdPA
            dWdPB= dWdU*dUdPB
            dWdGA= dWdU*dUdGA
            dWdGB= dWdU*dUdGB
            EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
            EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
*           dF/dRho
            dF_dRho(ipRa,iGrid)=dF_dRho(ipRa,iGrid)+EUEGPA*W
     &                                             +EUEG*dWdPA
            dF_dRho(ipRb,iGrid)=dF_dRho(ipRb,iGrid)+EUEGPB*W
     &                                             +EUEG*dWdPB
*           dF/dGamma
            dF_dRho(ipGaa,iGrid)=dF_dRho(ipGaa,iGrid)+EUEG*dWdGA
            dF_dRho(ipGbb,iGrid)=dF_dRho(ipGbb,iGrid)+EUEG*dWdGB
 112        Continue
         End Do
      End If
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(CoeffA)
      End
      Subroutine m06css(DTol,PX,GX,TX,F,FP,FG,FT,EUEG,Chi,EUEGP,
     &                   ChiP,ChiG,ijzy)
      Implicit none
C
C     Compute the same-spin part of the m06 correlation functional for one grid
C     point and one spin-case.
C
C
      integer ijzy
      REAL*8 PX, GX, TX, F, FP, FG, FT, DTol
      REAL*8 EUEG, Chi, EUEGP, ChiP, ChiG
      REAL*8 Zero, Pt25, F1, F2, F3, F4, F5, F6, F8, F11
      REAL*8 ss, sss0,sss1, sss2, sss3, sss4, Css
      REAL*8 Pi, Pi34, F13, F23, F43, F53, F83, F113
      REAL*8 RS, FDUEG, D, Fscc, RSP, dFsccP, dFsccG
      REAL*8 E, W, U, dFsccT, dUdChi, dWdU, dWdP, dWdG
      REAL*8 d2LdSS,d2LdSZ,d2LdZZ,PotLC,dLdS,dLdZ


      Data Zero/0.0d0/, Pt25/0.25d0/, F1/1.0d0/, F2/2.0d0/, F3/3.0d0/,
     $  F4/4.0d0/, F5/5.0d0/, F6/6.0d0/, F8/8.0d0/, F11/11.0d0/,
     $  Css/0.06d0/
C
c     DTol=1.0D-7
      ss=1.0D0
      sss0=0.0D0
      sss1=0.0D0
      sss2=0.0D0
      sss3=0.0D0
      sss4=0.0D0

      if (ijzy.eq.1) then
C     Parameters for M06-L Correlation
         sss0=  5.349466D-01
         sss1=  5.396620D-01
         sss2=  -3.161217D+01
         sss3=  5.149592D+01
         sss4=  -2.919613D+01
      elseif (ijzy.eq.2) then
C     Parameters for M06-HF Correlation
         sss0=  1.023254D-01
         sss1=  -2.453783D+00
         sss2=  2.913180D+01
         sss3=  -3.494358D+01
         sss4=  2.315955D+01
      elseif (ijzy.eq.3) then
C     Parameters for M06 Correlation
         sss0=  5.094055D-01
         sss1=  -1.491085D+00
         sss2=  1.723922D+01
         sss3=  -3.859018D+01
         sss4=  2.845044D+01
      elseif (ijzy.eq.4) then
C     Parameters for M06-2X Correlation
         sss0=  3.097855D-01
         sss1=  -5.528642D+00
         sss2=  1.347420D+01
         sss3=  -3.213623D+01
         sss4=  2.846742D+01
      endif

      If ((PX.le.DTol))  then
        EUEG = Zero
        Chi = Zero
        EUEGP = Zero
        ChiP = Zero
        ChiG = Zero
C       PX = Zero
C       GX = Zero
C       TX = Zero
        F  = Zero
        FP = Zero
        FG = Zero
        FT = Zero
      else
        Pi = F4*ATan(F1)
        Pi34 = F3 / (F4*Pi)
        F13 = F1 / F3
        F23 = F2 / F3
        F43 = F2 * F23
        F53 = F5 / F3
        F83 = F8 / F3
        F113 = F11 / F3
        FDUEG = (F3/F5)*(F6*Pi*Pi)**F23
        RS = (Pi34/PX) ** F13
        Call lsdac(RS,F1,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
        EUEG = PX*PotLC
        D = TX - Pt25*GX/PX
C        DUEG = FDUEG*PX**F53
        Chi = GX/PX**F83
        U = Css*Chi/(F1 + Css*Chi)
        W = sss0+U*(sss1+U*(sss2+U*(sss3+U*sss4)))
        Fscc=D/TX
        E = Fscc*W*EUEG
        F = E*ss
        RSP = -RS/(F3*Px)
        ChiG = F1/PX**F83
        ChiP = -F83*Chi/PX
        dFsccP=Pt25*GX/(TX*PX**2)
        dFsccG=-Pt25/(TX*PX)
        dFsccT=Pt25*GX/(PX*TX**2)
        dUdChi=Css/((F1+Css*Chi)**2)
        dWdU=sss1+U*(F2*sss2+U*(F3*sss3+U*F4*sss4))
        dWdP=dWdU*dUdChi*ChiP
        dWdG=dWdU*dUdChi*ChiG
        EUEGP = PotLC + PX*dLdS*RSP
        FP = ss*(dFsccP*W*EUEG
     $                 + Fscc*dWdP*EUEG
     $                 + Fscc*W*EUEGP)
        FG = ss*(dFsccG*W*EUEG
     $                 + Fscc*dWdG*EUEG)

        FT = ss*(dFsccT*W*EUEG)
       Endif

       Return
       End

c
c Perdew 91 local correlation functional at one grid point
c

      Subroutine lsdac(rs,zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
     &
      Implicit Real*8 (A-H,O-Z)
C     Evaluate the Perdew 91 local correlation functional and its
C     derivatives at one point.

      REAL*8 rs
      REAL*8 PotLC,dLdS, dLdZ
      REAL*8 eps0c(6), eps1c(6), epsc(6)
      REAL*8 F1, F2, F3, F4, F6, F8, F9, F12
      REAL*8 GammaI,Zeta,FZeta,dfZdz,d2fZdz
      REAL*8 EU,dEUdRS,d2UdRS
      REAL*8 EP,dEPdRS,d2PdRS
      REAL*8 AlphaM,dAMdRS,d2AdRS
      REAL*8 GZ, HZ, dGZ, dHZ, d2GZ, d2HZ

      data eps0c/0.03109070D0,0.21370D0, 7.5957D0,3.5876D0,1.6382D0,
     &         0.49294D0/
      data eps1c/0.01554535D0,0.20548D0,14.1189D0,6.1977D0,3.3662D0,
     &         0.62517D0/
      data epsc/0.01688690D0,0.11125D0,10.3570D0,3.6231D0,0.88026D0,
     &         0.49671D0/
      data F1/1.0d0/, F2/2.0d0/, F3/3.0d0/, F4/4.0d0/, F6/6.0d0/,
     &  F8/8.0d0/, F9/9.0d0/, F12/12.0d0/

      Pi = F4*ATan(F1)
      Third = F1 / F3

      epsc(1) = F1 / (F6*Pi*Pi)
      FZZI = F9*(F2**Third-F1) / F4
      GammaI = F1 / (F2*F2**Third-F2)

      Call EvFZet(GammaI,Zeta,FZeta,dfZdz,d2fZdz)
      Call EvPWLC(eps0c(1),eps0c(2),eps0c(3),eps0c(4),eps0c(5),
     &            eps0c(6), RS,EU,dEUdRS,d2UdRS)
      Call EvPWLC(eps1c(1),eps1c(2),eps1c(3),eps1c(4),eps1c(5),
     &            eps1c(6),RS,EP,dEPdRS,d2PdRS)
      Call EvPWLC(epsc(1),epsc(2),epsc(3),epsc(4),epsc(5),
     &            epsc(6),RS,AlphaM,dAMdRS,d2AdRS)
      Z2 = Zeta*Zeta
      Z3 = Zeta*Z2
      Z4 = Zeta*Z3
      GZ = FZeta*Z4
      HZ = FZZI*(FZeta-GZ)
      PotLC = EU*(F1-GZ) + EP*GZ - AlphaM*HZ
      dLdS = dEUdRS*(F1-GZ) + dEPdRS*GZ - dAMdRS*HZ
      dGZ = dfZdz*Z4 + F4*FZeta*Z3
      dHZ = FZZI*(dFZdz-dGZ)
      dLdz = (EP-EU)*dGZ - AlphaM*dHZ
      d2GZ = d2fZdz*Z4 + F8*Z3*dfZdz + F12*FZeta*Z2
      d2HZ = FZZI*(d2FZdz-d2GZ)
      d2LdSS = d2UdRS*(F1-GZ) + d2PdRS*GZ - d2AdRS*HZ
      d2LdSZ = (dEPdRS-dEUdRS)*dGZ - dAMdRS*dHZ
      d2LdZZ = (EP-EU)*d2GZ - AlphaM*d2HZ
      Return
      End

c
c   f(zeta)
c
      Subroutine EvFZet(S,Zeta,FZeta,dfZdz,d2fZdz)
      Implicit Real*8(A-H,O-Z)
c
c     evaluate f(Zeta) and its derivatives for lsdac.
c
      REAL*8 Small
      REAL*8 S, Zeta, FZeta,dfZdz,d2fZdz
      REAL*8 Zero, One, Two, Three, Four, Nine, F8, F27
      REAL*8 OMZ, OPZ, OMZ2, OPZ2, OMZ3, OPZ3
      REAL*8 F13, F43, F49, F827
      data Zero/0.0d0/, One/1.0d0/, Two/2.0d0/, Three/3.0d0/,
     $  Four/4.0d0/, Nine/9.0d0/, F8/8.0D0/, F27/27.0D0/
C
      Small = 1.0d-14
      FZeta = -Two
      dfZdz = Zero
      d2fZdz = Zero
      OMZ = One - Zeta
      OPZ = One + Zeta
      OMZ2 = OMZ**2
      OPZ2 = OPZ**2
      F13 = One / Three
      F43 = Four / Three
      F49 = Four / Nine
      F827 = F8 / F27
      If(OMZ.gt.Small) then
        OMZ3 = OMZ ** F13
        fZeta = fZeta + OMZ*OMZ3
        dfZdz = dfZdz - OMZ3
        d2fZdz = d2fZdz + OMZ3/OMZ
        endIf
      If(OPZ.gt.Small) then
        OPZ3 = OPZ ** F13
        fZeta = fZeta + OPZ*OPZ3
        dfZdz = dfZdz + OPZ3
        d2fZdz = d2fZdz + OPZ3/OPZ
        endIf
      fZeta = fZeta * S
      dfZdz = dfZdz * F43 * S
      d2fZdz = d2fZdz * F49 * S
      Return
      End

c
c  pw91 local correlation
c
      Subroutine EvPWLC(A,A1,B1,B2,B3,B4,RS,V,dVdRS,d2VdRS)
      Implicit none
C
C     Evaluate the interpolation function for PW91 local correlation.
C
      REAL*8 A,A1,B1,B2,B3,B4,RS,V,dVdRS,d2VdRS
      REAL*8 F1,F2, F3, F4
      REAL*8 Q0, RS12, RS32,Q1,Q2
      REAL*8 dQ0dRS,dQ1dRS,dQ2dRS
      REAL*8 d2Q1dS, d2Q2dS
      data F1/1.0d0/, F2/2.0d0/, F3/3.0d0/, F4/4.0d0/
C
      Q0 = -F2*A*(F1+A1*RS)
      RS12 = Sqrt(RS)
      RS32 = RS*RS12
      Q1 = F2*A*(B1*RS12+B2*RS+B3*RS32+B4*RS*RS)
      Q2 = Log(F1+F1/Q1)
      V = Q0*Q2
C
      dQ0dRS = -F2*A*A1
      dQ1dRS = A*(B1/RS12+F2*B2+F3*B3*RS12+F4*B4*RS)
      dQ2dRS = -dQ1dRS/(Q1+Q1**2)
      dVdRS = dQ0dRS*Q2 + Q0*dQ2dRS
C
      d2Q1dS = A*(F3*B3/(RS12*F2)-B1/(RS32*F2)+F4*B4)
      d2Q2dS = (F2*Q1+F1)*(dQ1dRS/(Q1+Q1**2))**2 - d2Q1dS/(Q1+Q1**2)
      d2VdRS = F2*dQ0dRS*dQ2dRS + Q0*d2Q2dS
      Return
      End
