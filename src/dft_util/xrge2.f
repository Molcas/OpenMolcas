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
* Copyright (C) 2006, Per Ake Malmqvist                                *
*               2010, Grigory A. Shamov                                *
************************************************************************
      Subroutine XRGE2(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                   Coeff,iSpin,F_xc,T_X)
************************************************************************
*                                                                      *
* Object: To compute the exchange part of RGE2, by                     *
*  Ruzsinszky, Csonka, Scuseria, JCTC 2009, 5, 763-769                 *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
*             University of Lund, SWEDEN. June 2006                    *
*             Grigory A Shamov, U of Manitoba, spring 2010             *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),
     &       F_xc(mGrid)

* IDORD=Order of derivatives to request from XPBE:
      idord=1


      if (ispin.eq.1) then
* ispin=1 means spin zero.
* T_X: Screening threshold of total density.
        Ta=0.5D0*T_X
        do iGrid=1,mgrid
         rhoa=max(1.0D-24,Rho(ipR,iGrid))
         if(rhoa.lt.Ta) goto 110
         grdrhoa_x=Rho(ipdRx,iGrid)
         grdrhoa_y=Rho(ipdRy,iGrid)
         grdrhoa_z=Rho(ipdRz,iGrid)
         sigmaaa=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2

         call testRGE2_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2)
         F_xc(iGrid)=F_xc(iGrid)+Coeff*(2.0D0*Fa)
         dF_dRho(ipR,iGrid)=dF_dRho(ipR,iGrid)+Coeff*dFdrhoa
* Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
         dF_dRho(ipGxx,iGrid)=dF_dRho(ipGxx,iGrid)+Coeff*dFdgammaaa
* Note: For xpbe, dFdgammaab is zero.
 110     continue
        end do
      else
* ispin .ne. 1, use both alpha and beta components.
        do iGrid=1,mgrid
         rhoa=max(1.0D-24,Rho(ipRa,iGrid))
         rhob=max(1.0D-24,Rho(ipRb,iGrid))
         rho_tot=rhoa+rhob
         if(rho_tot.lt.T_X) goto 210
         grdrhoa_x=Rho(ipdRxa,iGrid)
         grdrhoa_y=Rho(ipdRya,iGrid)
         grdrhoa_z=Rho(ipdRza,iGrid)
         sigmaaa=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2
         call testRGE2_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         grdrhob_x=Rho(ipdRxb,iGrid)
         grdrhob_y=Rho(ipdRyb,iGrid)
         grdrhob_z=Rho(ipdRzb,iGrid)
         sigmabb=grdrhob_x**2+grdrhob_y**2+grdrhob_z**2
         call testRGE2_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,
     &          d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2)

         F_xc(iGrid)=F_xc(iGrid)+Coeff*(Fa+Fb)
         dF_dRho(ipRa,iGrid)=dF_dRho(ipRa,iGrid)+Coeff*dFdrhoa
         dF_dRho(ipRb,iGrid)=dF_dRho(ipRb,iGrid)+Coeff*dFdrhob
* Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
* Note: For xpbe, dFdgammaab is zero.
         dF_dRho(ipGaa,iGrid)=dF_dRho(ipGaa,iGrid)+Coeff*dFdgammaaa
         dF_dRho(ipGbb,iGrid)=dF_dRho(ipGbb,iGrid)+Coeff*dFdgammabb
 210     continue
        end do
      end if

      Return
      End

C Cmu to be modified to 10/81 to make it PBEsol
C      Data Ckp, Cmu / 0.804D0, 0.12345679012346D0/


      Subroutine testXPBE_(idord,rho_s,sigma_s,
     &                        f,dFdr,dFdg,d2Fdr2,d2Fdrdg,d2Fdg2)
************************************************************************
*                                                                      *
* Object: To compute the PBE exchange functional from MOLPRO manual    *
*         using Maxima derivatives and optimization/fortran generator  *
*  Ref: PBE,"Generalized gradient approximation made simple,"          *
*      Phys. Rev. Lett. 77 (1996) 3865-3868.                           *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author: Grigory Shamov, U of Manitoba, Feb-March 2010           *
************************************************************************
      Implicit None
      Integer idord
      real*8 rho_s,sigma_s, rho, gamma, k, mu, Tpi, Cx, one
      real*8 f,dFdr,dFdg,d2Fdr2,d2Fdrdg,d2Fdg2
      parameter(Tpi=3.141592653589793d0,
     &  Cx=0.9305257363491d0, one=1.0d0 )

      Data k, mu / 0.804D0, 0.2195149727645171D0/
C      Data CkF / 3.0936677262801359310D0/
C      Data CeX /-0.73855876638202240588D0/
      real*8 T1,T2,T3,T4,T5,T6 ,T7, T8
C
      rho=max(1.0D-24,rho_s)
      gamma=max(1.0D-24,sigma_s)


      T1=rho**(4.d+0/3.d+0)
      T2=one/k
      T3=1.6455307846020562d-2*T2*mu*gamma/rho**(8.d+0/3.d+0)+1.d+0
      T4=-1.d+0*k/T3+k+1.d+0
      T5=one/rho**(7.d+0/3.d+0)
      T6=one/T3**2
      T7=mu**2
      T8=one/T3**3

      F = -9.305257363491001d-1*T1*T4

      dFdr = 4.083223320071
     4   841d-2*mu*T5*gamma*T6-1.2407009817988002d+0*rho**(1.d+0/3.d+0)*
     5   T4

      dFdg = -1.5312087450269407d-2*mu*T6/T1

      d2Fdr2 = -4.135669939
     6   329333d-1*T4/rho**(2.d+0/3.d+0)-4.0832233200718443d-2*mu*gamma*
     7   T6/rho**(1.d+1/3.d+0)+3.583503825911055d-3*T2*T7*gamma**2*T8/rh
     8   o**6

      d2Fdrdg = 2.0416116600359208d-2*mu*T5*T6-1.343813934716646
     9   d-3*T2*T7*gamma*T8/rho**5

      d2Fdg2 = 5.03930225518742d-4*T2*T7*T8
     =   /rho**4

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(idord)
      End

      Subroutine testRGE2_(idord,rho_s,sigma_s,
     &                        f,dFdr,dFdg,d2Fdr2,d2Fdrdg,d2Fdg2)
************************************************************************
*                                                                      *
* Object: To compute the RGE exchange functional. Procedure: the PBE   *
*  exchange functional from MOLPRO manual with Fx changed to RGE, the  *
*  mu parameter changed to be as of PBEsol. Code was made              *
*         using Maxima derivatives and optimization/fortran generator  *
*                                                                      *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author: Grigory Shamov, U of Manitoba, Feb-March 2010           *
************************************************************************

      Implicit None
      Integer idord
      real*8 rho_s,sigma_s, rho, gamma, k, mu, Tpi, Cx, one
      real*8 f,dFdr,dFdg,d2Fdr2,d2Fdrdg,d2Fdg2
      parameter(Tpi=3.141592653589793d0,
     &  Cx=0.9305257363491d0, one=1.0d0 )

C Cmu to be modified to 10/81 to make it PBEsol
      Data k, mu / 0.804D0, 0.12345679012346D0/

C      Data k, mu / 0.804D0, 0.2195149727645171D0/
C      Data CkF / 3.0936677262801359310D0/
C      Data CeX /-0.73855876638202240588D0/
      real*8 T1, T2, T3, T4, T5, T6
      real*8 T7, T8, T9, T10, T11
      real*8 T12, T13, T14, T15, T16
C
      rho=max(1.0D-24,rho_s)
      gamma=max(1.0D-24,sigma_s)

      T1= rho**(4.d+0/3.d+0)
      T2=one/k
      T3=one/rho**(8.d+0/3.d+0)
      T4=one/k**2
      T5=  mu**2
      T6=one/rho**(1.6d+1/3.d+0)
      T7=gamma**2
      T8=2.707771563073058
     3   7d-4*T4*T5*T6*T7+1.6455307846020562d-2*T2*mu*T3*gamma+1.d+0
      T9=  -1.d+0*k/T8+k+1.d+0
      T10=one/rho**(1.1000000000000001d+1/3.d+0)
      T11=one/rho**(1.9d+1/3.d+0)
      T12=-1.4441448336389645d-3*T4*T5*T11*T7
     6   -4.388082092272149d-2*T2*mu*T10*gamma
      T13=one/T8**2
      T14=rho**(1.d+0/3.d+0)
      T15=5.415543126146117d-4*T4*T5*T6*gamma+1.64553078460
     8   20562d-2*T2*mu*T3
      T16=one/T8**3

      F = -9.305257363491001d-1*T1*T9

      dFdr = -1.2407009817988002d+0*T14*T9-9.305257363491001d-1*k*T1*
     =   T12*T13
      dFdg = -9.305257363491001d-1*k*T1*T15*T13

      d2Fdr2 = -4.1
     ;   35669939329333d-1*T9/rho**(2.d+0/3.d+0)-2.4814019635975995d+0*k
     <   *T14*T12*T13-9.305257363491001d-1*k*T1*(9.146250613046776d-3*T4
     =   *T5*T7/rho**(2.2000000000000003d+1/3.d+0)+1.6089634338331216d-1
     >   *T2*mu*gamma/rho**(1.3999999999999999d+1/3.d+0))*T13+1.86105147
     ?   26982003d+0*k*T1*T12**2*T16

      d2Fdrdg = -1.2407009817988002d+0*k*
     @   T14*T15*T13-9.305257363491001d-1*k*T1*(-2.888289667277929d-3*T4
     1   *T5*T11*gamma-4.388082092272149d-2*T2*mu*T10)*T13+1.86105147269
     2   82003d+0*k*T1*T15*T12*T16

      d2Fdg2 = 1.8610514726982003d+0*k*T1*T
     3   15**2*T16-5.03930225518742d-4*T2*T5*T13/rho**4

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(idord)
      End
