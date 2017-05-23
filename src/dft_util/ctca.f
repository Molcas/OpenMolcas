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
*               2010, Grigory A. Shamov                                *
************************************************************************
      Subroutine CTCA(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                Coeff,iSpin,F_xc,T_X)
************************************************************************
*                                                                      *
* Object: To compute the Tognetti-Cortona-Adamo TCA functional         *
*         ref: J Chem Phys 128, 034101 (2008)                          *
*         template from CPBE functional, derived  with Maxima          *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
*             University of Lund, SWEDEN. December 2005 (CPBE)         *
*             Grigory Shamov, U of Manitoba (TCA) March 2010           *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
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

* IDORD=Order of derivatives to request from CPBE:
      idord=1

* cpbe has three input variables apart from idord=1:
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

         gamma_in=(grdrho_x**2+grdrho_y**2+grdrho_z**2)
C         grdrho_in=sqrt(gamma)
         zeta_in=0.0D0
         call ctca_(idord,rho_in,gamma_in,zeta_in,func0,func1,func2)
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
         gamma_in=(grdrho_x**2+grdrho_y**2+grdrho_z**2)
C         grdrho_in=sqrt(gamma)
         zeta_in=(rhoa-rhob)/rho_in
         call ctca_(idord,rho_in,gamma_in,zeta_in,func0,func1,func2)
         F_xc(iGrid)=F_xc(iGrid)+Coeff*func0
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

      subroutine ctca_(idord,rho_in,gamma_in,zeta_in,func0,func1,func2)
      implicit none

C
C     Parameter variables
C
      REAL*8       PI, alpha, sigma, one
      PARAMETER   (PI = 3.1415926535897932384626433832795D0, one=1.0d0)
      parameter (alpha = 2.30D0, sigma=1.43d0)
C
C     Argument variables
C
      REAL*8        FUNC0,       FUNC1(3),    FUNC2(3,3)
      REAL*8        GAMMA_IN,   RHO_IN,      ZETA_IN
      INTEGER       IDORD
C
C     Local variables
C
      real*8 rho, zeta, gamma

      real*8 T1
      real*8 T2
      real*8 T3
      real*8 T4
      real*8 T5
      real*8 T6
      real*8 T7
      real*8 T8
      real*8 T9
      real*8 T10
      real*8 T11
      real*8 T12
      real*8 T13
      real*8 T14
      real*8 T15
      real*8 T16
      real*8 T17
      real*8 T18
      real*8 T19
      real*8 T20
      real*8 T21
      real*8 T22
      real*8 T23
      real*8 T24


      rho=Max(1.0D-24,rho_in)
      gamma=gamma_in
      if(zeta_in.gt. 1.0D0) zeta= 1.0d0
      if(zeta_in.lt.-1.0D0) zeta=-1.0d0
      zeta=zeta_in*(1.0d0-2.3d-16)


* Finally, the functional is defined as the integrand in the
* expression for the correlation energy, i.e. there is a factor
* rho:

      T1=rho**(1.d+0/3.d+0)
      T2=1.970876462555557d+0/T1+4.88827d+0
      T3=8.97889d-1-6.55868d-1*atan(T2)
      T4=rho**(4.d+0/3.d+0)
      T5=one/T4
      T6=sqrt(gamma)
      T7=1.616204596739955d-1*T5*sigma*T6+1.d+0
      T8=one/T7**alpha
      T9=1.d+0-1.d+0*zeta
      T10=zeta+1.d+0
      T11=T10**(2.d+0/3.d+0)+T9**(2.d+0/3.d+0)
      T12=T11**3
      T13=T2**2+1.d+0
      T14=one/T13
      T15=one/rho
      T16=-alpha
      T17=T7**(T16-1)
      T18=one/T6
      T19=6.666666666666666d-1/T10**(1.d+0/3.d+0)-6.66666666666
     8   6666d-1/T9**(1.d+0/3.d+0)
      T20=T11**2
      T21=one/rho**(7.d+0/3.d+0)
      T22=-1.d+0*alpha-1.d+0
      T23=sigma**2
      T24=T7**(T16-2)

      func0 = 2.0149899425205864d-1*T3*T4*T8*T12


      if(idord.ge.1) then

      func1(1) = 4.342181343315399d-2*a
     ;  lpha*T3*T15*sigma*T17*T6*T12+2.686653256694115d-1*T3*T1*T8*T12+
     <  8.682153762983333d-2*T14*T8*T12
      func1(2) = -1.6283180037432754d-2*alpha*T3*sigma*T17*T18*T12
      func1(3) = 6.04496982756176d-1*T3*T4*T8*T19*T20

      end if

      if(idord.ge.2) then

      func2(1,1) = -9.357137929259761d-3*T22*alpha*T3*T23*T24
     ?  *gamma*T12/rho**(1.d+1/3.d+0)+1.4473937811051357d-2*alpha*T3*si
     @  gma*T17*T6*T12/rho**2+3.741903152356469d-2*alpha*T14*T21*sigma*
     1  T17*T6*T12+8.955510855647049d-2*T3*T8*T12/rho**(2.d+0/3.d+0)+1.
     2  1576205017311111d-1*T14*T15*T8*T12+1.140763499716801d-1*T2*T5*T
     3  8*T12/T13**2

      func2(2,2) = 8.141590018716377d-3*alpha*T3*sigma*T17*
     4  T12/T6**3-1.3158475213021542d-3*T22*alpha*T3*T5*T23*T24*T12/gam
     5  ma

      func2(3,3) = 6.04496982756176d-1*T3*T4*T8*(-2.2222222222222224d
     6  -1/T10**(4.d+0/3.d+0)-2.2222222222222224d-1/T9**(4.d+0/3.d+0))*
     7  T20+1.208993965512352d+0*T3*T4*T8*T19**2*T11

      func2(1,2) = 3.508926
     8  723472411d-3*T22*alpha*T3*T21*T23*T24*T12-7.01606841066838d-3*a
     9  lpha*T14*T5*sigma*T17*T18*T12

      func2(1,2) = 1.30265440299462d-1*alp
     =  ha*T3*T15*sigma*T17*T6*T19*T20+8.059959770082344d-1*T3*T1*T8*T1
     ;  9*T20+2.6046461288949996d-1*T14*T8*T19*T20

      func2(2,3) = -4.8849540
     <  112298256d-2*alpha*T3*sigma*T17*T18*T19*T20

      end if

      return
      end
