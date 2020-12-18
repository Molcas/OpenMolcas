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
      Subroutine CVS98(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                 CoeffA,iSpin,F_xc,T_X,ijzy)
************************************************************************
*                                                                      *
*  CVS98 evaluates the correlation part of the VS98 and M06 suite of   *
*  functionals on a grid.                                              *
*  !!! Second derivatives are not available yet.                       *
*                                                                      *
*  Ref:  T. V. Voorhis and G. E. Scuseria, J. Chem. Phys. 109,         *
*        400 (1998).                                                   *
*       Y. Zhao and D. G. Truhlar, J. Chem. Phys. 125,                 *
*        194101 (2006).                                                *
*                                                                      *
*                                                                      *
*       ijzy - 1 correlation functional in VS98                        *
*       ijzy - 2 correlation functional in M06-L                       *
*       ijzy - 3 correlation functional in M06-HF                      *
*       ijzy - 4 correlation functional in M06                         *
*       ijzy - 5 correlation functional in M06-2X                      *
*                                                                      *
*  YZ (10/07)                                                          *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)
      Integer mGrid

      REAL*8 Pi34, F13,
     &RS,RSP,Zeta,dZdA,dZdB,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ,
     &  P, EUEG, EUEGPA,EUEGPB
      REAL*8 PA,GAA,TauA,FA,FPA,FGA,FTA,EUA,ChiA,EUPA,ChiAP,
     &ChiAG,ZA,ZAP,ZAT
      REAL*8 PB,GBB,TauB,FB,FPB,FGB,FTB,EUB,ChiB,EUPB,ChiBP,
     &ChiBG,ZB,ZBP,ZBT
      REAL*8 ZAB, XAB, kab, xk, zk
      REAL*8 dgdx,dgdz,dgdPA,dgdGA,dgdTA,dgdPB,dgdGB,dgdTB
      REAL*8 gcab,gab
      REAL*8 r7, r8, r9, r10, r11, r12


      REAL*8  DTol,F1, F3, F4
      Data F1/1.0d0/,F3/3.0d0/,F4/4.0d0/
     $,gab/0.00304966d0/

C     Parameters for VS98
      if (ijzy.eq.1) then
              r7=   7.035010d-01
              r8=   7.694574d-03
              r9=   5.152765d-02
              r10=   3.394308d-05
              r11=  -1.269420d-03
              r12=   1.296118d-03
C     Parameters for M06-L
      elseif (ijzy.eq.2) then
              r7=      3.957626D-01
              r8=      -5.614546D-01
              r9=      1.403963D-02
              r10=     9.831442D-04
              r11=     -3.577176D-03
              r12=     0.000000D+00
C     Parameters for M06-HF
      elseif (ijzy.eq.3) then
              r7=    -6.746338D-01
              r8=    -1.534002D-01
              r9=    -9.021521D-02
              r10=   -1.292037D-03
              r11=   -2.352983D-04
              r12=   0.000000D+00

C     Parameters for M06
      elseif (ijzy.eq.4) then
               r7= -2.741539D+00
               r8= -6.720113D-01
               r9= -7.932688D-02
               r10=1.918681D-03
               r11=-2.032902D-03
               r12=0.000000D+00

C     Parameters for M06-2X
      elseif (ijzy.eq.5) then
              r7=  1.166404D-01
              r8=  -9.120847D-02
              r9=  -6.726189D-02
              r10= 6.720580D-05
              r11= 8.448011D-04
              r12= 0.000000D+00
      endif

      Pi34 = F3 / (F4*Pi)
      F13 = F1 / F3
      DTol = T_X
      Ta=0.5D0*T_X
*                                                                      *
************************************************************************
*                                                                      *
      If (iSpin.eq.1) Then
*
         Do iGrid = 1, mGrid
            PA=max(1.0D-24,Rho(ipR,iGrid),Ta)
            If (Rho(ipR,iGrid).lt.Ta) goto 110
            grdrhoa_x=Rho(ipdRx,iGrid)
            grdrhoa_y=Rho(ipdRy,iGrid)
            grdrhoa_z=Rho(ipdRz,iGrid)
            GAA=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2
            TauA = Rho(ipTau,iGrid)
            Call vs98ss(PA,GAA,TauA,FA,FPA,FGA,FTA,EUA,ZA,
     &                  ChiA,EUPA,ChiAP,ChiAG,ZAP,ZAT,ijzy)

            F_xc(iGrid)=F_xc(iGrid)+(2.0D0*FA)
*           dF/dRho
            dF_dRho(ipR,iGrid)=dF_dRho(ipR,iGrid)+FPA
*           dF/dGamma, no Gamma_ab term
            dF_dRho(ipGxx,iGrid)=dF_dRho(ipGxx,iGrid)+ FGA
*           dF/dTau
            dF_dRho(ipT,iGrid)=dF_dRho(ipT,iGrid)+FTA
110         Continue
*
            P=2.0D0*PA
            IF (PA.gt.DTol) THEN
               RS = (Pi34/P) ** F13
               RSP = -RS/(F3*P)
               Zeta = 0.0D0
               dZdA = F1/P
               Call lsdac(RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $                    d2LdZZ)
               EUEG = P*PotLC - 2.0D0*EUA
               ZAB = 2.0D0*ZA
               XAB = 2.0D0*ChiA
               kab = F1 + gab*(XAB+ZAB)
               xk = XAB/kab
               zk = ZAB/kab
               call gvt4(gcab,dgdx,dgdz,xk,zk,kab,gab,gab,r7,r8,r9,r10,
     &                   r11,r12)
               F_xc(iGrid)=F_xc(iGrid)+gcab*EUEG
               dgdPA = dgdx*ChiAP + dgdz*ZAP
               dgdGA = dgdx*ChiAG
               dgdTA = dgdz*ZAT
               EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
               EUEGPB = PotLC + P*dLdS*RSP - P*dLdZ*dZdA - EUPA
*              dF/dRho
               dF_dRho(ipR,iGrid)=dF_dRho(ipR,iGrid)
     &                           +(EUEGPA*gcab+EUEG*dgdPA)
*              dF/dGamma
               dF_dRho(ipGxx,iGrid)=dF_dRho(ipGxx,iGrid)+EUEG*dgdGA
*              dF/dTau
               dF_dRho(ipT,iGrid)=dF_dRho(ipT,iGrid)+EUEG*dgdTA
            ENDIF
         End Do
*                                                                      *
************************************************************************
*                                                                      *
      Else
         Do iGrid = 1, mGrid
            PA=max(1.0D-24,Rho(ipRa,iGrid),Ta)
            If(Rho(ipRa,iGrid).lt.Ta) goto 100
            grdrhoa_x=Rho(ipdRxa,iGrid)
            grdrhoa_y=Rho(ipdRya,iGrid)
            grdrhoa_z=Rho(ipdRza,iGrid)
            GAA=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2
            TauA = Rho(ipTaua,iGrid)
            Call vs98ss(PA,GAA,TauA,FA,FPA,FGA,FTA,EUA,ZA,
     &                  ChiA,EUPA,ChiAP,ChiAG,ZAP,ZAT,ijzy)

            F_xc(iGrid)=F_xc(iGrid)+FA
*           dF/dRhoa
            dF_dRho(ipRa,iGrid)=dF_dRho(ipRa,iGrid)+FPA
*           dF/dGammma_aa
            dF_dRho(ipGaa,iGrid)=dF_dRho(ipGaa,iGrid)+ FGA
*           dF/dTaua
            dF_dRho(ipTa,iGrid)=dF_dRho(ipTa,iGrid)+FTA
100         Continue
*
            PB=max(1.0D-24,Rho(ipRb,iGrid),Ta)
            if(Rho(ipRb,iGrid).lt.Ta) goto 111
            grdrhob_x=Rho(ipdRxb,iGrid)
            grdrhob_y=Rho(ipdRyb,iGrid)
            grdrhob_z=Rho(ipdRzb,iGrid)
            GBB=grdrhob_x**2+grdrhob_y**2+grdrhob_z**2
            TauB = Rho(ipTaub,iGrid)
            Call vs98ss(PB,GBB,TauB,FB,FPB,FGB,FTB,EUB,ZB,
     &                  ChiB,EUPB,ChiBP,ChiBG,ZBP,ZBT,ijzy)

            F_xc(iGrid)=F_xc(iGrid)+FB
*           dF/dRhob
            dF_dRho(ipRb,iGrid)=dF_dRho(ipRb,iGrid)+FPB
*           dF/dGammma_bb
            dF_dRho(ipGbb,iGrid)=dF_dRho(ipGbb,iGrid)+ FGB
*           dF/dTaub
            dF_dRho(ipTb,iGrid)=dF_dRho(ipTb,iGrid)+FTB
111         Continue

            P=PA+PB
            IF (PB.gt.DTol.and.PA.gt.DTol) THEN
               RS = (Pi34/P) ** F13
               RSP = -RS/(F3*P)
               Zeta = (PA-PB)/P
               dZdA = (F1-Zeta)/P
               dZdB = (-F1-Zeta)/P
               Call lsdac(RS,Zeta,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,
     $                    d2LdZZ)
               EUEG = P*PotLC - EUA - EUB
               ZAB = ZA + ZB
               XAB = ChiA+ChiB
               kab = F1 + gab*(XAB+ZAB)
               xk = XAB/kab
               zk = ZAB/kab
               call gvt4(gcab,dgdx,dgdz,xk,zk,kab,gab,gab,r7,r8,r9,r10,
     &                   r11,r12)
               F_xc(iGrid)=F_xc(iGrid)+gcab*EUEG
               dgdPA = dgdx*ChiAP + dgdz*ZAP
               dgdGA = dgdx*ChiAG
               dgdTA = dgdz*ZAT
               dgdPB = dgdx*ChiBP + dgdz*ZBP
               dgdGB = dgdx*ChiBG
               dgdTB = dgdz*ZBT
               EUEGPA = PotLC + P*dLdS*RSP + P*dLdZ*dZdA - EUPA
               EUEGPB = PotLC + P*dLdS*RSP + P*dLdZ*dZdB - EUPB
*              dF/dRho
               dF_dRho(ipRa,iGrid)=dF_dRho(ipRa,iGrid)
     &                            +(EUEGPA*gcab+EUEG*dgdPA)
               dF_dRho(ipRb,iGrid)=dF_dRho(ipRb,iGrid)
     &                            +(EUEGPB*gcab+EUEG*dgdPB)
*              dF/dGamma
               dF_dRho(ipGaa,iGrid)=dF_dRho(ipGaa,iGrid)+EUEG*dgdGA
               dF_dRho(ipGbb,iGrid)=dF_dRho(ipGbb,iGrid)+EUEG*dgdGB
*              dF/dTau
               dF_dRho(ipTa,iGrid)=dF_dRho(ipTa,iGrid)+EUEG*dgdTA
               dF_dRho(ipTb,iGrid)=dF_dRho(ipTb,iGrid)+EUEG*dgdTB
            END IF
         End Do
*
      End If
*                                                                      *
************************************************************************
*                                                                      *
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(CoeffA)
      End

      Subroutine vs98ss(PX,GX,TX,F,FP,FG,FT,EUEG,Z,Chi,EUEGP,
     &                   ChiP,ChiG,ZP,ZT,ijzy)
      Implicit none
C
C     Compute the same-spin part of the vs98 correlation functional for one grid
C     point and one spin-case.
C

      integer ijzy
      REAL*8 r13, r14, r15, r16, r17, r18
      REAL*8 PX, GX, TX, F, FP, FG, FT, DTol, Z, ZP, ZT
      REAL*8 EUEG, Chi, EUEGP, ChiP, ChiG, cf, gcc
      REAL*8 Zero, F1, F2, F3, F4, F5, F8, F11
      REAL*8 Pi, Pi34, F13, F23, F53, F83
      REAL*8 RS, D, RSP, PotLC, DX, DZ, dgdP, dgdG, dgdT
      REAL*8 E,DP, DG, DT, rhoo, rho43, rho53, rho83
      REAL*8 rrho, F4o3, kc, xk, zk, gc, dgdx, dgdz
      REAL*8 d2LdSS, d2LdSZ, d2LdZZ, dLdS, dLdZ

      Data Zero/0.0d0/, F1/1.0d0/, F2/2.0d0/, F3/3.0d0/,
     $  F4/4.0d0/, F5/5.0d0/, F8/8.0d0/, F11/11.0d0/,
     $  gcc/0.00515088d0/,cf/9.115599720d0/


      F4o3 = 4.0d0/3.0d0
C     Parameters for VS98
      if (ijzy.eq.1) then
              r13=   3.270912d-01
              r14=  -3.228915d-02
              r15=  -2.942406d-02
              r16=   2.134222d-03
              r17=  -5.451559d-03
              r18=   1.577575d-02
C     Parameters for M06-L
      elseif (ijzy.eq.2) then
              r13=   4.650534D-01
              r14=   1.617589D-01
              r15=   1.833657D-01
              r16=   4.692100D-04
              r17=  -4.990573D-03
              r18=   0.000000D+00
C     Parameters for M06-HF
      elseif (ijzy.eq.3) then
              r13=   8.976746D-01
              r14=  -2.345830D-01
              r15=   2.368173D-01
              r16=  -9.913890D-04
              r17=  -1.146165D-02
              r18=   0.000000D+00
C     Parameters for M06
      elseif (ijzy.eq.4) then
               r13=  4.905945D-01
               r14= -1.437348D-01
               r15=  2.357824D-01
               r16=  1.871015D-03
               r17= -3.788963D-03
               r18=  0.000000D+00
C     Parameters for M06-2X
      elseif (ijzy.eq.5) then
              r13=  6.902145D-01
              r14=  9.847204D-02
              r15=  2.214797D-01
              r16= -1.968264D-03
              r17= -6.775479D-03
              r18=  0.000000D+00
      endif

      DTol =1.0d-8
      If(PX.le.DTol) then
        EUEG = Zero
        Chi = Zero
        EUEGP = Zero
        ChiP = Zero
        ChiG = Zero
        PX = Zero
        GX = Zero
        TX = Zero
        F  = Zero
        FP = Zero
        FG = Zero
        FT = Zero
        Z  = Zero
        ZP = Zero
        ZT = Zero
      else
        Pi = F4*ATan(F1)
        Pi34 = F3 / (F4*Pi)
        F13 = F1 / F3
        F23 = F2 / F3
        F53 = F5 / F3
        F83 = F8 / F3
        rhoo = PX
        rrho = 1.0D0/rhoo
        rho43 = rhoo**F4o3
        rho53 = rhoo**F53
        rho83 = rho53*rhoo

        RS = (Pi34/PX) ** F13
        Call lsdac(RS,F1,PotLC,dLdS,dLdZ,d2LdSS,d2LdSZ,d2LdZZ)
        EUEG = PX*PotLC
        Chi = GX/rho83
        Z = (TX/rho53) - cf
        kc = F1 + gcc*(Chi + Z)
        xk = Chi/kc
        zk = Z/kc
*       D = F1 - Chi/(F4*(Z + cf))
        D = F1 - Chi*rho53/(F4*TX)
        call gvt4(gc,dgdx,dgdz,xk,zk,kc,gcc,gcc,r13,r14,r15,r16,r17,r18)
        E = D*EUEG*gc
c         write (*,*) "Chi, Z, gc", CHi, Z, gc
        F = E
c
        RSP = -RS/(F3*Px)
        ChiG = F1/PX**F83
        ChiP = -F83*Chi/PX
        ZP = -F53 * TX/rho83
        ZT =  F1/rho53
*       DZ = Chi/(F4*(Z + cf)*(Z + cf))
        DZ = Chi*rho53**2/(F4*TX**2)
*       DX = -F1/(F4*(Z + cf))
        DX = -F1*rho53/(F4*TX)
        DP = DZ*ZP + DX*ChiP
        DG = DX*ChiG
        DT = DZ*ZT
        dgdP = dgdx*ChiP + dgdz*ZP
        dgdG = dgdx*ChiG
        dgdT = dgdz*ZT
        EUEGP = PotLC + PX*dLdS*RSP
        FP = DP*EUEG*gc + D*EUEGP*gc + D*EUEG*dgdP
        FG = DG*EUEG*gc + D*EUEG*dgdG
        FT = DT*EUEG*gc + D*EUEG*dgdT
       Endif
       Return
       End
