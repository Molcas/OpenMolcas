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
      Subroutine XVS98(Rho,nRho,mGrid,dF_dRho,ndF_dRho,
     &                 CoeffA,iSpin,F_xc,T_X,ijzy)
************************************************************************
*                                                                      *
*  xvs98 evaluates the exchange part of the VS98 and M06 suite of      *
*  functionals on a grid.                                              *
*  !!! Second derivatives are not available yet.                       *
*                                                                      *
*  Ref:  T. V. Voorhis and G. E. Scuseria, J. Chem. Phys. 109,         *
*        400 (1998).                                                   *
*       Y. Zhao and D. G. Truhlar, J. Chem. Phys. 125,                 *
*        194101 (2006).                                                *
*                                                                      *
*       ijzy - 1 VS98                                                  *
*       ijzy - 2 M06-L                                                 *
*       ijzy - 3 M06-HF                                                *
*       ijzy - 4 M06                                                   *
*                                                                      *
*  YZ (10/07)                                                          *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "nq_index.fh"
      Real*8 Rho(nRho,mGrid),dF_dRho(ndF_dRho,mGrid),F_xc(mGrid)
      Integer mGrid

      integer ijzy
      REAL*8 rrho, rho43, rho13, rhoo, rho53, rho83
      REAL*8  Gamma
c
c     kinetic energy density or tau
c
      REAL*8 tauu,DTol
      REAL*8 f13, f43, f53, f83, f113, f4o3
      REAL*8 gx, gg, x, z, kx,xk,zk
      REAL*8 Nine, F10, F11
      REAL*8 cf, Axlsda, r1, r2, r3, r4, r5, r6

c      functional derivatives below FFFFFFFFFFFF

       REAL*8 dxdr, dxdg, dzdr, dzdt, dgdx, dgdz

c      functional derivatives above FFFFFFFFFFFF

       parameter (cf = 9.115599720d0, Axlsda = -0.9305257363491d0 )
       parameter (gg  = 0.00186726d0)
       parameter (f13=1.d0/3.d0,f43=4.0d0/3.0d0,f53=5.0d0/3.0d0)
       parameter (f83=8.d0/3.0d0, F113=11.0d0/3.d0,f4o3=4.0D0/3.D0)
       parameter (Nine=9.0d0,F10=10.d0, F11=11.d0)

      if (ijzy.eq.1) then
c
c     Parameters for VS98
c
        r1=  -9.800683d-01
        r2=  -3.556788d-03
        r3=   6.250326d-03
        r4=  -2.354518d-05
        r5=  -1.282732d-04
        r6=   3.574822d-04
      elseif (ijzy.eq.2) then
c
c     Parameters for M06-L
c
        r1 =   6.012244D-01*Axlsda
        r2 =   4.748822D-03*Axlsda
        r3 =  -8.635108D-03*Axlsda
        r4 =  -9.308062D-06*Axlsda
        r5 =   4.482811D-05*Axlsda
        r6 =   0.000000D+00
      elseif (ijzy.eq.3) then
c
c     Parameters for M06-HF
c
        r1 =   -1.179732D-01*Axlsda
        r2 =   -2.500000D-03*Axlsda
        r3 =   -1.180065D-02*Axlsda
        r4 =   0.000000D+00
        r5 =   0.000000D+00
        r6 =   0.000000D+00
      elseif (ijzy.eq.4) then
c
c     Parameters for M06
c
c
c     Parameters for M06
c
        r1 =   1.422057D-01*Axlsda
        r2 =   7.370319D-04*Axlsda
        r3 =   -1.601373D-02*Axlsda
        r4 =   0.000000D+00
        r5 =   0.000000D+00
        r6 =   0.000000D+00
      endif

      DTol = T_X

*                                                                      *
************************************************************************
*                                                                      *
*     iSpin=1
*
      If (iSpin.eq.1) Then
*                                                                      *
************************************************************************
*                                                                      *
        Ta=0.5D0*T_X
        Do iGrid = 1, mGrid
         rhoo=max(1.0D-24,rho(ipR,igrid))
         if(rhoo.lt.Ta) goto 110
         grdrhoa_x=rho(ipdRx,igrid)
         grdrhoa_y=rho(ipdRy,igrid)
         grdrhoa_z=rho(ipdRz,igrid)

         rho43 = rhoo**F4o3
         rrho = 1.0d0/rhoo       ! reciprocal of rho
         rho13 = rho43*rrho
         rho53 = rhoo**F53
         rho83 = rho53*rhoo
         tauu=Rho(ipTau,iGrid)
         Gamma=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2
         x = gamma/rho83
         dxdr = -f83*x*rrho
         dxdg = One/rho83
         z = tauu/rho53 - cf
         dzdr = -f53 * tauu/rho83
         dzdt = One/rho53
         kx = One + gg*x + gg*z
         xk = x/kx
         zk = z/kx
         call gvt4(gx,dgdx,dgdz,xk,zk,kx,gg,gg,r1,r2,r3,r4,r5,r6)
         F_xc(iGrid)=F_xc(iGrid)+(2.0D0*rho43*gx)
c
c     functional derivatives
c
*        dF/dRho
         dF_dRho(ipR,iGrid)=dF_dRho(ipR,iGrid)+ f43*rho13*gx +
     &                  rho43*(dgdx*dxdr + dgdz*dzdr)
*        dF/dGamma
         dF_dRho(ipGxx,iGrid)=dF_dRho(ipGxx,iGrid)+rho43*(dgdx*dxdg)
*        dF/dTau
         dF_dRho(ipT,iGrid)=dF_dRho(ipT,iGrid)+ rho43*(dgdz*dzdt)
110     continue
        Enddo
       else
* ispin .ne. 1, use both alpha and beta components.
        Ta=0.5D0*T_X
        do igrid=1,mgrid
c
* alpha component
c
         rhoo=max(1.0D-24,rho(ipRa,iGrid))
         if(rhoo.lt.Ta) goto 210
         grdrhoa_x=rho(ipdRxa,iGrid)
         grdrhoa_y=rho(ipdRya,iGrid)
         grdrhoa_z=rho(ipdRza,iGrid)

         rho43 = rhoo**F4o3
         rrho = 1.0d0/rhoo       ! reciprocal of rho
         rho13 = rho43*rrho
         rho53 = rhoo**F53
         rho83 = rho53*rhoo
         tauu=Rho(ipTaua,iGrid)
         Gamma=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2
         x = gamma/rho83
         dxdr = -f83*x*rrho
         dxdg = One/rho83
         z = tauu/rho53 - cf
         dzdr = -f53 * tauu/rho83
         dzdt = One/rho53
         kx = One + gg*x + gg*z
         xk = x/kx
         zk = z/kx
         call gvt4(gx,dgdx,dgdz,xk,zk,kx,gg,gg,r1,r2,r3,r4,r5,r6)
         F_xc(iGrid)=F_xc(iGrid)+rho43*gx
c
c     functional derivatives
c
*        dF/dRhoa
         dF_dRho(ipRa,iGrid)=dF_dRho(ipRa,iGrid)+ f43*rho13*gx +
     &                  rho43*(dgdx*dxdr + dgdz*dzdr)
*        dF/dGammaaa
         dF_dRho(ipGaa,iGrid)=dF_dRho(ipGaa,iGrid)+rho43*(dgdx*dxdg)
*        dF/dTaua
         dF_dRho(ipTa,iGrid)=dF_dRho(ipTa,iGrid)+ rho43*(dgdz*dzdt)
210      continue
c
c beta component
c
         rhoo=max(1.0D-24,rho(ipRb,iGrid))
         if(rhoo.lt.Ta) goto 310
         grdrhoa_x=rho(ipdRxb,iGrid)
         grdrhoa_y=rho(ipdRyb,iGrid)
         grdrhoa_z=rho(ipdRzb,iGrid)

         rho43 = rhoo**F4o3
         rrho = 1.0d0/rhoo       ! reciprocal of rho
         rho13 = rho43*rrho
         rho53 = rhoo**F53
         rho83 = rho53*rhoo
         tauu=Rho(ipTaub,iGrid)
         Gamma=grdrhoa_x**2+grdrhoa_y**2+grdrhoa_z**2
         x = gamma/rho83
         dxdr = -f83*x*rrho
         dxdg = One/rho83
         z = tauu/rho53 - cf
         dzdr = -f53 * tauu/rho83
         dzdt = One/rho53
         kx = One + gg*x + gg*z
         xk = x/kx
         zk = z/kx
         call gvt4(gx,dgdx,dgdz,xk,zk,kx,gg,gg,r1,r2,r3,r4,r5,r6)
         F_xc(iGrid)=F_xc(iGrid)+rho43*gx
c
c     functional derivatives
c
*        dF/dRhob
         dF_dRho(ipRb,iGrid)=dF_dRho(ipRb,iGrid)+ f43*rho13*gx +
     &                  rho43*(dgdx*dxdr + dgdz*dzdr)
*        dF/dGammabb
         dF_dRho(ipGbb,iGrid)=dF_dRho(ipGbb,iGrid)+rho43*(dgdx*dxdg)
*        dF/dTaub
         dF_dRho(ipTb,iGrid)=dF_dRho(ipTb,iGrid)+ rho43*(dgdz*dzdt)
310      continue
        enddo
       endif
*
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_real(CoeffA)
      End

      Subroutine gvt4(g,dgdx,dgdz,xk,zk,k,c,ct,r1,r2,r3,r4,r5,r6)
      Implicit none
c
c     Evaluate the GVT4 form in the VS98 functional
c
c
c    some working variables
      REAL*8 g,dgdx,dgdz,xk,zk,k,c,ct,r1,r2,r3,r4,r5,r6
      REAL*8 sxk,szk,sk,sc,sct,sr1,sr2,sr3,sr4,sr5,sr6,sk2
      REAL*8 One, Two, Three
      Data One/1.0d0/, Two/2.0d0/, Three/3.0d0/
C
      sxk = xk
      szk = zk
      sk = k
      sc = c
      sct = ct
      sr1 = r1
      sr2 = r2
      sr3 = r3
      sr4 = r4
      sr5 = r5
      sr6 = r6
      sk2 = sk*sk
      g =  (sr1 + sr2*sxk + sr3*szk
     $  +sr4*sxk*sxk + sr5*szk*sxk + sr6*szk*szk)/sk
      dgdx =   (-sr1*sc
     $  +sr2*(One-Two*sc*sxk)
     $  -Two*sr3*szk*sc
     $  +sr4*(Two*sxk-Three*sxk*sxk*sc)
     $  +sr5*(szk -Three*szk*sxk*sc)
     $  -Three*sr6*szk*szk*sc )/sk2
      dgdz =   (-sr1*sct
     $  -Two*sr2*sxk*sct
     $  +sr3*(One-Two*szk*sct)
     $  -Three*sr4*sxk*sxk*sct
     $  +sr5*(sxk-Three*sxk*szk*sct)
     $  +sr6*(Two*szk-Three*szk*szk*sct))/sk2

      return
      end
