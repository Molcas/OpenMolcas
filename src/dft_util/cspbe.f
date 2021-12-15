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
* Copyright (C) 2005, Per Ake Malmqvist                                *
*               2017, Giovanni Li Manni                                *
*               2017, Aron Cohen                                       *
************************************************************************
      Subroutine CSPBE(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                Coeff,iSpin,F_xc,T_X)
************************************************************************
*                                                                      *
* Object: To compute the functional called c_spbe in the Density       *
* Functional Repository (http://www.cse.clrc.ac.uk/qcg/dft)            *
* Original reference article: J.P. Perdew, K. Burke, and M. Ernzerhof, *
*  "Generalized gradient approximation made simple,"                   *
*      Phys. Rev. Lett. 77 (1996) 3865-3868.                           *
*                                                                      *
*      Author: G. Li Manni & A. Cohen, Max Planck Institute Stuttgart  *
*              Summer 2017, edited in Cambridge (UK) & Palermo (Sicily)*
************************************************************************
      use KSDFT_Info, only: tmpB
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
#include "ksdft.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)
* Local arrays:
      Real*8 func1(3),func2(3,3)
* Call arguments:
* Weights(mGrid) (input) integration weights.
* Rho(nRho,mGrid) (input) Density and density derivative values,
*   Rho(1,iGrid) is rho_alpha values, Rho(2,iGrid) is rho_beta values
*   Rho(i,iGrid) is grad_rho_alpha (i=3..5 for d/dx, d/dy, d/dz)
*   Rho(i,iGrid) is grad_rho_beta  (i=6..8 for d/dx, d/dy, d/dz)
* dF_dRho (inout) are (I believe) values of derivatives of the
*   DFT functional (*NOT* derivatives of Fock matrix contributions).
* F_xc is values of the DFT energy density functional (surprised?)

* IDORD=Order of derivatives to request from CsPBE:
      idord=1

* cspbe has three input variables apart from idord=1:
*  The density (rho_in), its gradient (abs value!) grdrho_in,
*  and spin polarization zeta (zeta_in).
* The result is returned as follows:
*   func0 = Correlation energy density
*   func1(1) = Its first derivative wrt rho
*   func1(2) = Its first derivative wrt gamma
*   func1(3) = Its first derivative wrt zeta
*   func2(1,1) = Second derivative wrt rho
*   func2(2,2) = Second derivative wrt gamma
*   func2(3,3) = Second derivative wrt zeta
*   func2(1,2) = Mixed derivative wrt rho and gamma
*   func2(1,3) = Mixed derivative wrt rho and zeta
*   func2(2,3) = Mixed derivative wrt gamma and zeta
* Here, gamma=(grad rho)**2; AKA sigma in some formulas.
* The derivatives w.r.t. spin components rhoa, rhob, and to the
* gradients grdrhoa_x, etcetc, ..., grdrhob_z are then:
*  Let F(rhoa,rhob,grdrhoa_x,...,grdrhob_z) = CPBE(rho,gamma,zeta),
* then dF_drhoa = func1(1)+2*rhob*func1(3)/rho**2
*      dF_drhob = func1(1)-2*rhoa*func1(3)/rho**2
*      dF_dgrdrhoa_x = dF_dgrdrhob_x = 2*func1(2)*grdrho_x and so on.
      if (ispin.eq.1) then
* ispin=1 means spin zero.
        do iGrid=1,mgrid
         rhoa=Rho(ipR,iGrid)
         rho_in=2.0D0*rhoa
         if(rho_in.lt.T_X) goto 110
         grdrho_x=2.D0*Rho(ipdRx,iGrid)
         grdrho_y=2.D0*Rho(ipdRy,iGrid)
         grdrho_z=2.D0*Rho(ipdRz,iGrid)

         gamma=(grdrho_x**2+grdrho_y**2+grdrho_z**2)
         grdrho_in=sqrt(gamma)
         zeta_in=0.0D0
         call cspbe_(idord,rho_in,grdrho_in,zeta_in,func0,func1,func2)
         F_xc(iGrid)=F_xc(iGrid)+Coeff*func0
* dF_drhoa:
         dF_dRho(ipR,iGrid)=dF_dRho(ipR,iGrid)+Coeff*func1(1)
* Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
         dF_dRho(ipGxx,iGrid)=dF_dRho(ipGxx,iGrid)+Coeff*func1(2)
         dF_dRho(ipGxy,iGrid)=dF_dRho(ipGxy,iGrid)+Coeff*2.0D0*func1(2)
 110     continue
        end do
      else
* ispin .ne. 1, use both alpha and beta components.
        do iGrid=1,mgrid
         rhoa=max(1.0D-24,Rho(ipRa,iGrid))
         rhob=max(1.0D-24,Rho(ipRb,iGrid))
         rho_in=rhoa+rhob
         if(rho_in.lt.T_X) goto 210
         grdrhoa_x=Rho(ipdRxa,iGrid)
         grdrhoa_y=Rho(ipdRya,iGrid)
         grdrhoa_z=Rho(ipdRza,iGrid)
         grdrhob_x=Rho(ipdRxb,iGrid)
         grdrhob_y=Rho(ipdRyb,iGrid)
         grdrhob_z=Rho(ipdRzb,iGrid)

         grdrho_x=grdrhoa_x+grdrhob_x
         grdrho_y=grdrhoa_y+grdrhob_y
         grdrho_z=grdrhoa_z+grdrhob_z
         gamma=(grdrho_x**2+grdrho_y**2+grdrho_z**2)
         grdrho_in=sqrt(gamma)
         zeta_in=(rhoa-rhob)/rho_in
         call cspbe_(idord,rho_in,grdrho_in,zeta_in,func0,func1,func2)
         F_xc(iGrid)=F_xc(iGrid)+Coeff*func0
         tmpB(iGrid)=F_xc(iGrid)-tmpB(iGrid)
* dF_drhoa:
         dF_dRho(ipRa,iGrid)=dF_dRho(ipRa,iGrid)+
     &            Coeff*(func1(1)+(2.0D0*func1(3))*(rhob/rho_in**2))
         dF_dRho(ipRb,iGrid)=dF_dRho(ipRb,iGrid)+
     &            Coeff*(func1(1)-(2.0D0*func1(3))*(rhoa/rho_in**2))
* Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
         dF_dRho(ipGaa,iGrid)=dF_dRho(ipGaa,iGrid)+Coeff*func1(2)
         dF_dRho(ipGab,iGrid)=dF_dRho(ipGab,iGrid)+Coeff*2.0D0*func1(2)
         dF_dRho(ipGbb,iGrid)=dF_dRho(ipGbb,iGrid)+Coeff*func1(2)
 210     continue
        end do
      end if

      Return
      End

      subroutine cspbe_(idord,rho_in,grdrho_in,zeta_in
     &,func0,func1,func2)
      implicit none

C
C     Parameter variables
C
      REAL*8        ONE
      PARAMETER     (ONE = 1.0D0)
      REAL*8        TWO
      PARAMETER     (TWO = 2.0D0)
      REAL*8        THREE
      PARAMETER     (THREE = 3.0D0)
      REAL*8        FOUR
      PARAMETER     (FOUR = 4.0D0)
      REAL*8        HALF
      PARAMETER     (HALF = 0.5D0)
      REAL*8        THIRD
      PARAMETER     (THIRD = ONE/THREE)
      REAL*8        TOTHRD
      PARAMETER     (TOTHRD = TWO/THREE)
      REAL*8        FOTHRD
      PARAMETER     (FOTHRD = FOUR/THREE)
      REAL*8        PI
      PARAMETER     (PI = 3.1415926535897932384626433832795D0)
      REAL*8        KFCNST
      PARAMETER     (KFCNST = 1.9191582926775130066248203262468D0)
      REAL*8        FZCNST
      PARAMETER     (FZCNST = 1.9236610509315363197594581232755D0)
      REAL*8        D2FZ0
      PARAMETER     (D2FZ0 = 8.0D0*FZCNST/9.0D0)
      REAL*8        RHO2RS
      PARAMETER     (RHO2RS = 0.62035049088842779322331643499916D0)
      REAL*8        BETA
      PARAMETER     (BETA = 0.06672455060314922D0)
C
C     Argument variables
C
      REAL*8        FUNC0,       FUNC1(3),    FUNC2(3,3)
      REAL*8        GRDRHO_IN,   RHO_IN,      ZETA_IN
      INTEGER       IDORD
C
C     Local variables
C
      REAL*8        AC,          AI,          CFF1
      REAL*8        CFF2,        CFF3,        CFF4,        CFF5
      REAL*8        CFF6,        D2ACDR2,     D2AIDR2
      REAL*8        D2AIDRDZ,    D2AIDZ2
      REAL*8        D2ECDR2
      REAL*8        D2ECDRDZ,    D2ECDZ2,     D2ECFDR2
      REAL*8        D2ECPDR2,    D2FZDZ2,     D2GDX2
      REAL*8        D2PDG2
      REAL*8        D2PDGDZ
      REAL*8        D2PDR2
      REAL*8        D2PDRDG
      REAL*8        D2PDRDZ,     D2PDX2,      D2PDZ2,      D2QDG2
      REAL*8        D2QDGDZ
      REAL*8        D2QDR2,      D2QDRDG,     D2QDRDZ,     D2QDX2
      REAL*8        D2QDZ2,      D2TDG2,      D2TDGDZ
      REAL*8        D2TDR2,      D2TDRDG,     D2TDRDZ,     D2TDZ2
      REAL*8        D2XDG2
      REAL*8        D2XDGDZ
      REAL*8        D2XDR2,      D2XDRDG,     D2XDRDZ,     D2XDZ2
      REAL*8        D2YDG2
      REAL*8        D2YDGDZ
      REAL*8        D2YDR2
      REAL*8        D2YDRDG
      REAL*8        D2YDRDZ,     D2YDZ2,      DACDR,       DAIDR
      REAL*8        DAIDZ
      REAL*8        DECDR,       DECDZ,       DECFDR,      DECPDR
      REAL*8        DFZDZ,       DGDX
      REAL*8        DPDG,        DPDR,        DPDX,        DPDZ
      REAL*8        DQDG,        DQDR,        DQDX,        DQDZ
      REAL*8        DTDG,        DTDR,        DTDZ
      REAL*8        DXDG,        DXDR,        DXDZ
      REAL*8        DYDG,        DYDR,        DYDZ,        EC
      REAL*8        ECF,         ECP
      REAL*8        FZ,          G,           GRDRHO,      PHI
      REAL*8        PHI3G,       KF
      REAL*8        KS,          P,           Q,           RHO
      REAL*8        RHOI,        RS,          T,           X
      REAL*8        Y,           Z2,          Z3
      REAL*8        Z4,          ZETA
      REAL*8 D2HDG2
      REAL*8 D2HDGDZ
      REAL*8 D2HDR2
      REAL*8 D2HDRDG
      REAL*8 D2HDRDZ
      REAL*8 D2HDZ2
      REAL*8 D2PHIDZ2
      REAL*8 DHDG
      REAL*8 DHDR
      REAL*8 DHDZ
      REAL*8 DPHIDZ
      REAL*8 GAMMAS
      REAL*8 H
      REAL*8 s, dHds, d2Hdrds, d2Hdsdz, d2Hds2

      rho=Max(1.0D-24,rho_in)
      grdrho=grdrho_in
      if(zeta_in.gt. 1.0D0) zeta= 1.0d0
      if(zeta_in.lt.-1.0D0) zeta=-1.0d0
      zeta=zeta_in*(1.0d0-2.3d-16)

      rhoi=1.0D0/rho
      rs=rho2rs*(rhoi**third)
      x=sqrt(rs)
      kF=kFcnst/rs
      ks=sqrt(four/pi)*sqrt(kF)
*--------------------------------------------------------------
* PAM 2009: Some compilers doing code checking erroneously conclude that
* a lot of variables may get used uninitialized (due to if-statements).
* So they must get initialized here:
      d2acdr2=9.9D99
      d2aidr2=9.9D99
      d2aidrdz=9.9D99
      d2aidz2=9.9D99
      d2ecdr2=9.9D99
      d2ecdrdz=9.9D99
      d2ecdz2=9.9D99
      d2ecfdr2=9.9D99
      d2ecpdr2=9.9D99
      d2fzdz2=9.9D99
      d2pdg2=9.9D99
      d2pdgdz=9.9D99
      d2pdr2=9.9D99
      d2pdrdg=9.9D99
      d2pdrdz=9.9D99
      d2pdz2=9.9D99
      d2qdg2=9.9D99
      d2qdgdz=9.9D99
      d2qdr2=9.9D99
      d2qdrdg=9.9D99
      d2qdrdz=9.9D99
      d2qdz2=9.9D99
      d2tdgdz=9.9D99
      d2tdr2=9.9D99
      d2tdrdg=9.9D99
      d2tdrdz=9.9D99
      d2tdz2=9.9D99
      d2xdg2=9.9D99
      d2xdgdz=9.9D99
      d2xdr2=9.9D99
      d2xdrdg=9.9D99
      d2xdrdz=9.9D99
      d2xdz2=9.9D99
      d2ydg2=9.9D99
      d2ydgdz=9.9D99
      d2ydr2=9.9D99
      d2ydrdg=9.9D99
      d2ydrdz=9.9D99
      d2ydz2=9.9D99
      dacdr=9.9D99
      daidr=9.9D99
      daidz=9.9D99
      decdr=9.9D99
      decdz=9.9D99
      decfdr=9.9D99
      decpdr=9.9D99
      dfzdz=9.9D99
      dgdx=9.9D99
      dpdg=9.9D99
      dpdr=9.9D99
      dpdx=9.9D99
      dpdz=9.9D99
      dqdg=9.9D99
      dqdr=9.9D99
      dqdx=9.9D99
      dqdz=9.9D99
      dtdg=9.9D99
      dtdr=9.9D99
      dtdz=9.9D99
      dxdg=9.9D99
      dxdr=9.9D99
      dxdz=9.9D99
      dydg=9.9D99
      dydr=9.9D99
      dydz=9.9D99
      d2hdg2=9.9D99
      d2hdgdz=9.9D99
      d2hdr2=9.9D99
      d2hdrdg=9.9D99
      d2hdrdz=9.9D99
      d2hdz2=9.9D99
      d2phidz2=9.9D99
      dhdg=9.9D99
      dhdr=9.9D99
      dhdz=9.9D99
      dphidz=9.9D99
      dhds=9.9D99
      d2hdrds=9.9D99
      d2hdsdz=9.9D99
      d2hds2=9.9D99
*--------------------------------------------------------------
* Ac, and its derivatives w.r.t. rho:
      cff1=2.D0*0.016887D0
      cff2=0.11125D0
      cff3=10.357d0
      cff4=3.6231d0
      cff5=0.88026D0
      cff6=0.49671d0
      Q=cff1*x*(cff3+x*(cff4+x*(cff5+x*cff6)))
      P=log(one+one/Q)
      G=-cff1*(one+cff2*x**2)*P
      Ac=G
      if(idord.ge.1) then
       dxdr=-x*rhoi/6.0d0
       dQdx=cff1*(cff3+x*(2.D0*cff4+x*(3.D0*cff5+x*4.D0*cff6)))
       dPdx=-one/(Q*(Q+one))*dQdx
       dGdx=-cff1*(2.D0*cff2*x*P+(one+cff2*x**2)*dPdx)
       dAcdr=dGdx*dxdr
      end if
      if(idord.ge.2) then
       d2xdr2=-(7.0d0*dxdr*rhoi)/6.0d0
       d2Qdx2=cff1*(2.D0*cff4+x*(6.D0*cff5+x*12.D0*cff6))
       d2Pdx2=(dQdx**2.D0*(two-one/(Q+one))/Q-d2Qdx2)/(Q*(Q+one))
       d2Gdx2=-cff1*(2.D0*cff2*P+4.D0*cff2*x*dPdx
     &       +(one+cff2*x**2)*d2Pdx2)
       d2Acdr2=d2Gdx2*(dxdr**2)+dGdx*d2xdr2
      end if
*--------------------------------------------------------------
* EcP, and its derivatives:
      cff1=2.D0*0.0310907D0
      cff2=0.21370D0
      cff3=7.5957d0
      cff4=3.5876d0
      cff5=1.6382D0
      cff6=0.49294d0
      Q=cff1*x*(cff3+x*(cff4+x*(cff5+x*cff6)))
      P=log(one+one/Q)
      G=-cff1*(one+cff2*x**2)*P
      EcP=G
      if(idord.ge.1) then
       dQdx=cff1*(cff3+x*(2.D0*cff4+x*(3.D0*cff5+x*4.D0*cff6)))
       dPdx=-one/(Q*(Q+one))*dQdx
       dGdx=-cff1*(2.D0*cff2*x*P+(one+cff2*x**2)*dPdx)
       dEcPdr=dGdx*dxdr
      end if
      if(idord.ge.2) then
       d2Qdx2=cff1*(2.D0*cff4+x*(6.D0*cff5+x*12.D0*cff6))
       d2Pdx2=(dQdx**2.D0*(two-one/(Q+one))/Q-d2Qdx2)/(Q*(Q+one))
       d2Gdx2=-cff1*(2.D0*cff2*P+4.D0*cff2*x*dPdx
     &       +(one+cff2*x**2)*d2Pdx2)
       d2EcPdr2=d2Gdx2*(dxdr**2)+dGdx*d2xdr2
      end if
*--------------------------------------------------------------
* EcF, and its derivatives:
      cff1=2.D0*0.015545D0
      cff2=0.20548D0
      cff3=14.1189D0
      cff4=6.1977d0
      cff5=3.3662D0
      cff6=0.62517d0
      Q=cff1*x*(cff3+x*(cff4+x*(cff5+x*cff6)))
      P=log(one+one/Q)
      G=-cff1*(one+cff2*x**2)*P
      EcF=G
      if(idord.ge.1) then
       dQdx=cff1*(cff3+x*(2.D0*cff4+x*(3.D0*cff5+x*4.D0*cff6)))
       dPdx=-one/(Q*(Q+one))*dQdx
       dGdx=-cff1*(2.D0*cff2*x*P+(one+cff2*x**2)*dPdx)
       dEcFdr=dGdx*dxdr
      end if
      if(idord.ge.2) then
       d2Qdx2=cff1*(2.D0*cff4+x*(6.D0*cff5+x*12.D0*cff6))
       d2Pdx2=(dQdx**2.D0*(two-one/(Q+one))/Q-d2Qdx2)/(Q*(Q+one))
       d2Gdx2=-cff1*(2.D0*cff2*P+4.D0*cff2*x*dPdx
     &       +(one+cff2*x**2)*d2Pdx2)
       d2EcFdr2=d2Gdx2*(dxdr**2)+dGdx*d2xdr2
      end if
*--------------------------------------------------------------
* fz, and its derivatives:
      P=(one+zeta)**(fothrd)
      Q=(one-zeta)**(fothrd)
      fz=fzcnst*(P+Q-two)
      if(idord.ge.1) then
       dPdz= fothrd*P/(one+zeta)
       dQdz=-fothrd*Q/(one-zeta)
       dfzdz=fzcnst*(dPdz+dQdz)
      end if
      if(idord.ge.2) then
       d2Pdz2=third*dPdz/(one+zeta)
       d2Qdz2=-third*dQdz/(one-zeta)
       d2fzdz2=fzcnst*(d2Pdz2+d2Qdz2)
      end if
*--------------------------------------------------------------
* Ec(rho,zeta), and its derivatives:
      z2=zeta**2
      z3=z2*zeta
      z4=z3*zeta
      P=(Ac/d2fz0)
      Q=(P+EcF-EcP)
      Ec=EcP+(Q*z4-P)*fz
      if(idord.ge.1) then
       dPdr=(dAcdr/d2fz0)
       dQdr=dPdr+dEcFdr-dEcPdr
       dEcdr=dEcPdr+(dQdr*z4-dPdr)*fz
       dEcdz=four*Q*z3*fz+(Q*z4-P)*dfzdz
      end if
      if(idord.ge.2) then
       d2Pdr2=(d2Acdr2/d2fz0)
       d2Qdr2=d2Pdr2+d2EcFdr2-d2EcPdr2
       d2Ecdr2=d2EcPdr2+(d2Qdr2*z4-d2Pdr2)*fz
       d2Ecdrdz=four*dQdr*z3*fz+(dQdr*z4-dPdr)*dfzdz
       d2Ecdz2=12.0D0*Q*z2*fz+8.0d0*Q*z3*dfzdz+(Q*z4-P)*d2fzdz2
      end if
*--------------------------------------------------------------
* phi(zeta), and derivatives:
      kF=kFcnst/rs
      ks=sqrt(four/pi)*sqrt(kF)
      P=half*(one+zeta)**(tothrd)
      Q=half*(one-zeta)**(tothrd)
      phi=P+Q
      if(idord.ge.1) then
       dPdz= tothrd*P/(one+zeta)
       dQdz=-tothrd*Q/(one-zeta)
       dphidz=dPdz+dQdz
      end if
      if(idord.ge.2) then
       d2Pdz2=-third*dPdz/(one+zeta)
       d2Qdz2= third*dQdz/(one-zeta)
       d2phidz2=d2Pdz2+d2Qdz2
      end if
*--------------------------------------------------------------
* A(rho,zeta)=(beta/gammas)/(exp(-Ec/(phi**3*gammas))-1.0D0)
      gammas=0.031090690869654895034940863712730D0
      phi3g=phi**3*gammas

* and then used below in the combination "A*t**2".
* Here, we use the inverse Ai instead.

* P=Ec/(phi**3*gammas)

      P=Ec/phi3g
      Q=(gammas/beta)*exp(-P)
      Ai=Q-(gammas/beta)
      y=dphidz/phi
      if(idord.ge.1) then
       dPdr=dEcdr/phi3g
       dPdz=(dEcdz-3.D0*Ec*y)/phi3g
       dAidr=-Q*dPdr
       dAidz=-Q*dPdz
      end if
      if(idord.ge.2) then
       d2Pdr2=d2Ecdr2/phi3g
       d2Pdrdz=(d2Ecdrdz-3.D0*dEcdr*y)/phi3g
       d2Pdz2=(d2Ecdz2-6.D0*dEcdz*y+(12.D0*y**2-3.D0*(d2phidz2/phi))*Ec)
     &       /phi3g
       d2Aidr2=Q*(dPdr**2-d2Pdr2)
       d2Aidrdz=Q*(dPdr*dPdz-d2Pdrdz)
       d2Aidz2=Q*(dPdz**2-d2Pdz2)
      end if
* Ai and its derivatives have been checked OK.

* T in the code is T**2 in the original formula
      s=grdrho**2
      d2Tdg2=half/(phi*ks*rho)**2
      T=half*d2Tdg2*s
      if(idord.ge.1) then
       dTdg=d2Tdg2*grdrho
       dTdr=(-7.0d0*third)*T*rhoi
       dTdz=-2.D0*T*y
      end if
      if(idord.ge.2) then
       d2Tdrdg=(-7.0d0*third)*dTdg*rhoi
       d2Tdr2=(-10.0D0*third)*dTdr*rhoi
       d2Tdgdz=-2.D0*dTdg*y
       d2Tdz2=-3.D0*y*dTdz-2.D0*d2phidz2*T/phi
       d2Tdrdz=-2.D0*dTdr*y
      end if
* T and its derivatives have been checked OK.
!      x=Ai*(Ai+T)
!      if(idord.ge.1) then
!       dxdr=Ai*(2.D0*dAidr+dTdr)+T*dAidr
!       dxdz=Ai*(2.D0*dAidz+dTdz)+T*dAidz
!       dxdg=Ai*dTdg
!      end if
!      if(idord.ge.2) then
!       d2xdr2=2.D0*dAidr*(dAidr+dTdr)+Ai*(2.D0*d2Aidr2+d2Tdr2)+T*d2Aidr2
!       d2xdrdg=Ai*d2Tdrdg+dTdg*dAidr
!       d2xdrdz=dAidz*(2.D0*dAidr+dTdr)+Ai*(2.D0*d2Aidrdz+d2Tdrdz)+
!     &         dTdz*dAidr+T*d2Aidrdz
!       d2xdg2=Ai*d2Tdg2
!       d2xdgdz=Ai*d2Tdgdz+dAidz*dTdg
!       d2xdz2=2.D0*dAidz*(dAidz+dTdz)+Ai*(2.D0*d2Aidz2+d2Tdz2)+T*d2Aidz2
!      end if
!*
!      P=T*x
!      if(idord.ge.1) then
!       dPdr=dTdr*x+T*dxdr
!       dPdz=dTdz*x+T*dxdz
!       dPdg=dTdg*x+T*dxdg
!      end if
!      if(idord.ge.2) then
!       d2Pdr2=d2Tdr2*x+2.D0*dTdr*dxdr+T*d2xdr2
!       d2Pdrdz=d2Tdrdz*x+dTdr*dxdz+dTdz*dxdr+T*d2xdrdz
!       d2Pdz2=d2Tdz2*x+2.D0*dTdz*dxdz+T*d2xdz2
!       d2Pdg2=d2Tdg2*x+2.D0*dTdg*dxdg+T*d2xdg2
!       d2Pdrdg=d2Tdrdg*x+dTdr*dxdg+dTdg*dxdr+T*d2xdrdg
!       d2Pdgdz=d2Tdgdz*x+dTdg*dxdz+dTdz*dxdg+T*d2xdgdz
!      end if
!      Q=x+T**2
!      if(idord.ge.1) then
!       dQdr=dxdr+2.D0*T*dTdr
!       dQdz=dxdz+2.D0*T*dTdz
!       dQdg=dxdg+2.D0*T*dTdg
!      end if
!      if(idord.ge.2) then
!       d2Qdr2=d2xdr2+2.D0*(dTdr**2+T*d2Tdr2)
!       d2Qdrdz=d2xdrdz+2.D0*(dTdz*dTdr+T*d2Tdrdz)
!       d2Qdz2=d2xdz2+2.D0*(dTdz**2+T*d2Tdz2)
!       d2Qdg2=d2xdg2+2.D0*(T*d2Tdg2+dTdg**2)
!       d2Qdrdg=d2xdrdg+2.D0*(dTdg*dTdr+T*d2Tdrdg)
!       d2Qdgdz=d2xdgdz+2.D0*(dTdg*dTdz+T*d2Tdgdz)
!      end if
       y = Ai*T/(Ai+T)
!      y=P/Q
      if(idord.ge.1) then
       dydr= (dAidr*T+Ai*dTdR)/(Ai+T)-Ai*T*(dAidr+dTdr)/(Ai+T)**2
       dydz= (dAidz*T+Ai*dTdz)/(Ai+T)-Ai*T*(dAidz+dTdz)/(Ai+T)**2
       dydg= (Ai*dTdg)/(Ai+T)-Ai*T*(dTdg)/(Ai+T)**2
      end if
!      if(idord.ge.2) then
!       d2ydr2=(d2Pdr2-2.D0*dydr*dQdr-y*d2Qdr2)/Q
!       d2ydrdz=(d2Pdrdz-dydr*dQdz-dydz*dQdr-y*d2Qdrdz)/Q
!       d2ydz2=(d2Pdz2-2.D0*dydz*dQdz-y*d2Qdz2)/Q
!       d2ydg2=(d2Pdg2-2.D0*dydg*dQdg-y*d2Qdg2)/Q
!       d2ydrdg=(d2Pdrdg-dydr*dQdg-dydg*dQdr-y*d2Qdrdg)/Q
!       d2ydgdz=(d2Pdgdz-dydg*dQdz-dydz*dQdg-y*d2Qdgdz)/Q
!      end if
      P=log(one+(beta/gammas)*y)
      Q=(beta/gammas)/(one+(beta/gammas)*y)
      if(idord.ge.1) then
       dPdr=Q*dydr
       dPdz=Q*dydz
       dPdg=Q*dydg
      end if
      if(idord.ge.2) then
       d2Pdr2=Q*(d2ydr2-Q*dydr**2)
       d2Pdrdz=Q*(d2ydrdz-Q*dydr*dydz)
       d2Pdz2=Q*(d2ydz2-Q*dydz**2)
       d2Pdg2=Q*(d2ydg2-Q*dydg**2)
       d2Pdrdg=Q*(d2ydrdg-Q*dydr*dydg)
       d2Pdgdz=Q*(d2ydgdz-Q*dydg*dydz)
      end if

      H=phi3g*P
      if(idord.ge.1) then
       dHdr=phi3g*dPdr
       dHdg=phi3g*dPdg
       y=dphidz/phi
       dHdz=phi3g*(dPdz+3.D0*y*P)
      end if
      if(idord.ge.2) then
       d2Hdr2=phi3g*d2Pdr2
       d2Hdrdz=phi3g*(d2Pdrdz+3.D0*y*dPdr)
       d2Hdz2=phi3g*(d2Pdz2+6.D0*y*dPdz+
     &          (6.D0*y**2+3.D0*(d2phidz2/phi))*P)
       d2Hdg2=phi3g*d2Pdg2
       d2Hdrdg=phi3g*d2Pdrdg
       d2Hdgdz=phi3g*(d2Pdgdz+3.D0*y*dPdg)
      end if

* Change derivatives w.r.t. grad_rho into derivatives w.r.t. sigma:
      if(idord.ge.1) then
       dHds=dHdg/(2.D0*grdrho)
      end if
      if(idord.ge.2) then
       d2Hdrds=d2Hdrdg/(2.D0*grdrho)
       d2Hdsdz=d2Hdgdz/(2.D0*grdrho)
       d2Hds2 =(d2Hdg2-2.D0*dHds) /(4.D0*s)
      end if
* Finally, the functional is defined as the integrand in the
* expression for the correlation energy, i.e. there is a factor
* rho:
      func0=rho_in*(Ec+H)
      if(idord.ge.1) then
       func1(1)=rho_in*(dEcdr+dHdr)+(Ec+H)
       func1(2)=rho_in*(      dHds)
       func1(3)=rho_in*(dEcdz+dHdz)
      end if
      if(idord.ge.2) then
       func2(1,1)=rho_in*(d2Ecdr2 +d2Hdr2)+2.0D0*(dEcdr+dHdr)
       func2(1,2)=rho_in*(         d2Hdrds)+dHds
       func2(1,3)=rho_in*(d2Ecdrdz+d2Hdrdz)+(dEcdz+dHdz)
       func2(2,2)=rho_in*(         d2Hds2)
       func2(2,3)=rho_in*(         d2Hdsdz)
       func2(3,3)=rho_in*(d2Ecdz2 +d2Hdz2)
       func2(2,1)=func2(1,2)
       func2(3,1)=func2(1,3)
       func2(3,2)=func2(2,3)
      end if

      return
      end
