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
      Subroutine xB88(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                Coeff,iSpin,F_xc,T_X)
************************************************************************
*                                                                      *
* Object: To compute the functional called x_B88 in the Density        *
* Functional Repository (http://www.cse.clrc.ac.uk/qcg/dft)            *
* Following older code by Roland Lindh, this routine computes only     *
* the GGA addition to the LDA part.                                    *
* Original reference article:                                          *
*                                                                      *
*                                                                      *
*                                                                      *
* Called from:                                                         *
*                                                                      *
* Calling    :                                                         *
*                                                                      *
*      Author:Per Ake Malmqvist, Department of Theoretical Chemistry,  *
*             University of Lund, SWEDEN. June 2006                    *
************************************************************************
      use KSDFT_GLM, only: F_xca, F_xcb, tmpB
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
#include "ksdft.fh"
      Real*8 Rho(nRho,mGrid),
     &       dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)
cGLM     &       F_xca(mGrid),F_xcb(mGrid),tmpB(mGrid)

* IDORD=Order of derivatives to request from XPBE:
      idord=1
*
      Rho_Min=T_X*1.0D-2
*
      if (ispin.eq.1) then
* ispin=1 means spin zero.

* T_X: Screening threshold of total density.
        Ta=0.5D0*T_X
        do iGrid=1,mgrid
         Rhoa=Rho(ipR,iGrid)
         if(Rhoa.lt.Ta) goto 110
         grdrhoa_x=Rho(ipdRx,iGrid)
         grdrhoa_y=Rho(ipdRy,iGrid)
         grdrhoa_z=Rho(ipdRz,iGrid)
         sigmaaa=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2

         call xB88_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
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
         rhoa=Max(Rho_Min,rho(ipRa,iGrid))
         rhob=Max(Rho_Min,rho(ipRb,iGrid))
         rho_tot=rhoa+rhob
         if(rho_tot.lt.T_X) goto 210
         grdrhoa_x=rho(ipdRxa,iGrid)
         grdrhoa_y=rho(ipdRya,iGrid)
         grdrhoa_z=rho(ipdRza,iGrid)
         sigmaaa=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2
         call xB88_(idord,rhoa,sigmaaa,Fa,dFdrhoa,dFdgammaaa,
     &          d2Fdra2,d2Fdradgaa,d2Fdgaa2)

         grdrhob_x=rho(ipdRxb,iGrid)
         grdrhob_y=rho(ipdRyb,iGrid)
         grdrhob_z=rho(ipdRzb,iGrid)
         sigmabb=grdrhob_x**2+grdrhob_y**2+grdrhob_z**2
         call xB88_(idord,rhob,sigmabb,Fb,dFdrhob,dFdgammabb,
     &          d2Fdrb2,d2Fdrbdgbb,d2Fdgbb2)

         F_xc(iGrid) =F_xc(iGrid) +Coeff*(Fa+Fb)
         F_xca(iGrid)=F_xca(iGrid)+Coeff*(Fa)
         F_xcb(iGrid)=F_xcb(iGrid)+Coeff*(   Fb)
         dF_dRho(ipRa,iGrid)=dF_dRho(ipRa,iGrid)+Coeff*dFdrhoa
         dF_dRho(ipRb,iGrid)=dF_dRho(ipRb,iGrid)+Coeff*dFdrhob
* Maybe derivatives w.r.t. gamma_aa, gamma_ab, gamma_bb should be used instead.
* Note: For xpbe, dFdgammaab is zero.
         dF_dRho(ipGaa,iGrid)=dF_dRho(ipGaa,iGrid)+Coeff*dFdgammaaa
         dF_dRho(ipGbb,iGrid)=dF_dRho(ipGbb,iGrid)+Coeff*dFdgammabb
 210     continue
        end do
        tmpB(:)=F_xc(:)
      end if

      Return
      End

      subroutine xB88_(idord,rho_s,gamma_s,B88,dB88dr,dB88dg,d2B88dr2,
     &                d2B88drdg,d2B88dg2)
      implicit real*8 (a-h,o-z)
      parameter(third=1.0d0/3.0d0)
      parameter(four3=4.0d0/3.0d0)
      parameter(seven3=7.0d0/3.0d0)
      parameter(dcoef=0.0042d0)
C     parameter(xldacff=0.930525736349100025D0)

      rho=rho_s+1.0D-16
      gamma=gamma_s+1.0D-16
      r43 = rho**four3
      rhoinv=1.0d0/rho
* lda part:
*     xlda=-xldacff*r43
* Note: Use x=sqrt(gamma)/rho**four3
      x = sqrt(gamma_s)/r43
      hgi = 0.5D0/gamma
      p =sqrt(1.0D0+x**2)
      ash = log(x+p)
      d6 = 6.0D0*dcoef
      a = 1.0D0+d6*x*ash
      f = x**2/a

* Let b88(rho,gamma)=b(rho,x)
      dr43=-dcoef*r43
* The LDA part has been removed, just GGA part left.
*     b=dr43*f+xlda
      b=dr43*f
      b88 = b

      if(idord.lt.1) goto 99
      dxdr = -four3*x*rhoinv
      dxdg = hgi*x
      dadx = d6*(ash+x/p)
      dfdx = (2.0D0*x-f*dadx)/a
      dbdr = four3*b*rhoinv
      dbdx = dr43*dfdx
      db88dr = dbdr+dxdr*dbdx
      db88dg = dxdg*dbdx

      if(idord.lt.2) goto 99
      d2xdr2 = -seven3*dxdr*rhoinv
      d2xdrdg = hgi*dxdr
      d2xdg2 = -hgi*dxdg
      d2adx2 = d6*(1.0D0+p**2)/(p**3)
      d2fdx2 = (2.0D0-2.0D0*dadx*dfdx-d2adx2*f)/a
      d2bdr2 = third*dbdr*rhoinv
      d2bdx2 = dr43*d2fdx2
      d2bdrdx = four3*dbdx*rhoinv
      d2b88dr2 = d2bdr2+2.D0*dxdr*d2bdrdx+d2xdr2*dbdx+dxdr**2*d2bdx2
      d2b88dg2 = d2xdg2*dbdx+dxdg**2*d2bdx2
      d2b88drdg = dxdg*d2bdrdx+d2xdrdg*dbdx+dxdr*dxdg*d2bdx2

  99  continue

      return
      end
