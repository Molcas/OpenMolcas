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
************************************************************************
      Subroutine XPBE_ofe(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                    Coeff,iSpin,F_xc,T_X)
************************************************************************
*                                                                      *
* Object: To compute the functional called x_pbe in the Density        *
* Functional Repository (http://www.cse.clrc.ac.uk/qcg/dft)            *
* Original reference article: J.P. Perdew, K. Burke, and M. Ernzerhof, *
*  "Generalized gradient approximation made simple,"                   *
*      Phys. Rev. Lett. 77 (1996) 3865-3868.                           *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
*             University of Lund, SWEDEN. June 2006                    *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),
     &       F_xc(mGrid)
* Call arguments:
* Weights(mGrid) (input) integration weights.
* Rho(nRho,mGrid) (input) Density and density derivative values,
*   Rho(1,iGrid) is rho_alpha values, Rho(2,iGrid) is rho_beta values
*   Rho(i,iGrid) is grad_rho_alpha (i=3..5 for d/dx, d/dy, d/dz)
*   Rho(i,iGrid) is grad_rho_beta  (i=6..8 for d/dx, d/dy, d/dz)
* dF_dRho (inout) are (I believe) values of derivatives of the
*   DFT functional (*NOT* derivatives of Fock matrix contributions).
* F_xc is values of the DFT energy density functional (surprised?)

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

         call xpbe_1(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2)
         F_xc(iGrid)=F_xc(iGrid)+2.0D0*Fa
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
         call xpbe_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         grdrhob_x=Rho(ipdRxb,iGrid)
         grdrhob_y=Rho(ipdRyb,iGrid)
         grdrhob_z=Rho(ipdRzb,iGrid)
         sigmabb=grdrhob_x**2+grdrhob_y**2+grdrhob_z**2
         call xpbe_1(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,
     &          d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2)

         F_xc(iGrid)=F_xc(iGrid)+Fa+Fb
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

      Subroutine XPBE_1(idord,rho_s,sigma_s,
     &                        f,dFdr,dFdg,d2Fdr2,d2Fdrdg,d2Fdg2)
************************************************************************
*                                                                      *
* Object: To compute the functional called x_pbe in the Density        *
* Functional Repository (http://www.cse.clrc.ac.uk/qcg/dft)            *
* Original reference article: J.P. Perdew, K. Burke, and M. Ernzerhof, *
*  "Generalized gradient approximation made simple,"                   *
*      Phys. Rev. Lett. 77 (1996) 3865-3868.                           *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
*             University of Lund, SWEDEN. December 2006                *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Data Ckp, Cmu / 0.804D0, 0.2195149727645171D0/
      Data CkF / 3.0936677262801359310D0/
      Data CeX /-0.73855876638202240588D0/

      rho=max(1.0D-24,rho_s)
      sigma=max(1.0D-24,sigma_s)

      rthrd=(2.D0*rho)**(1.0D0/3.0D0)
      XkF=CkF*rthrd

* FX, and its derivatives wrt S:
      s2=sigma/((2.D0*rho)*XkF)**2
      s=sqrt(s2)
      Cmus2=Cmu*s2
      t=1.0D0/(Ckp+Cmus2)
      fx=(Cmus2+Ckp*(1.0D0+Cmus2))*t
      a=2.0D0*Cmu*(Ckp*t)**2
      dfxds=a*s
      d2fxds2=-a*(3.0D0*Cmus2-Ckp)*t

* The derivatives of S wrt rho (r)  and sigma (g)
      a=1.0D0/(3.0D0*rho)
      b=1.0D0/(2.0D0*sigma)
      dsdr=-4.0D0*s*a
      dsdg=s*b
      d2sdr2=-7.0D0*dsdr*a
      d2sdrdg=dsdr*b
      d2sdg2=-dsdg*b

* Thus, derivatives of fx wrt rho and sigma
      dfxdr=dsdr*dfxds
      dfxdg=dsdg*dfxds
      d2fxdr2=d2sdr2*dfxds+dsdr**2*d2fxds2
      d2fxdrdg=d2sdrdg*dfxds+dsdr*dsdg*d2fxds2
      d2fxdg2=d2sdg2*dfxds+dsdg**2*d2fxds2

* rho*XeX, and its derivatives wrt rho
      rX=rho*CeX*rthrd
      drXdr=4.0d0*rX*a
      d2rXdr2=drXdr*a

* Put it together:
      F=rX*fx
      dFdr=drXdr*fx+rX*dfxdr
      dFdg=rX*dfxdg
      d2Fdr2=d2rXdr2*fx+2.0D0*drXdr*dfxdr+rX*d2fxdr2
      d2Fdrdg=drXdr*dfxdg +rX*d2fxdrdg
      d2Fdg2=rX*d2fxdg2

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(idord)
      End
