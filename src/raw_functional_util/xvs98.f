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
      Subroutine XVS98(mGrid,
     &                 CoeffA,iSpin,F_xc,ijzy)
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
      use nq_Grid, only: Rho, Sigma, Tau
      use nq_Grid, only: vRho, vSigma, vTau
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
      Real*8 F_xc(mGrid)
      Real*8, Parameter:: T_X=1.0D-20
      Integer mGrid

      integer ijzy
      REAL*8 rrho, rho43, rho13, rhoo, rho53, rho83
      REAL*8  Gamma
c
c     kinetic energy density or tau
c
      REAL*8 tauu
      REAL*8 f43, f53, f83
      REAL*8 gx, gg, x, z, kx,xk,zk
      REAL*8 cf, Axlsda, r1, r2, r3, r4, r5, r6

c      functional derivatives below FFFFFFFFFFFF

       REAL*8 dxdr, dxdg, dzdr, dzdt, dgdx, dgdz

c      functional derivatives above FFFFFFFFFFFF

       parameter (cf = 9.115599720d0, Axlsda = -0.9305257363491d0 )
       parameter (gg  = 0.00186726d0)
       parameter (f43=4.0d0/3.0d0, f53=5.0d0/3.0d0, f83=8.d0/3.0d0)

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
         rhoo=max(1.0D-24,rho(1,igrid))
         if(rhoo.lt.Ta) goto 110

         rho43 = rhoo**F43
         rrho = 1.0d0/rhoo       ! reciprocal of rho
         rho13 = rho43*rrho
         rho53 = rhoo**F53
         rho83 = rho53*rhoo
         tauu=Tau(1,iGrid)
         Gamma=Sigma(1,iGrid)
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
         vRho(1,iGrid)=vRho(1,iGrid)+ f43*rho13*gx +
     &                  rho43*(dgdx*dxdr + dgdz*dzdr)
*        dF/dGamma
         vSigma(1,iGrid)=vSigma(1,iGrid)+rho43*(dgdx*dxdg)
*        dF/dTau
         vTau(1,iGrid)=vTau(1,iGrid)+ rho43*(dgdz*dzdt)
110     continue
        Enddo
       else
* ispin .ne. 1, use both alpha and beta components.
        Ta=0.5D0*T_X
        do igrid=1,mgrid
c
* alpha component
c
         rhoo=max(1.0D-24,rho(1,iGrid))
         if(rhoo.lt.Ta) goto 210

         rho43 = rhoo**F43
         rrho = 1.0d0/rhoo       ! reciprocal of rho
         rho13 = rho43*rrho
         rho53 = rhoo**F53
         rho83 = rho53*rhoo
         tauu=Tau(1,iGrid)
         Gamma=Sigma(1,iGrid)
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
         vRho(1,iGrid)=vRho(1,iGrid)+ f43*rho13*gx +
     &                  rho43*(dgdx*dxdr + dgdz*dzdr)
*        dF/dGammaaa
         vSigma(1,iGrid)=vSigma(1,iGrid)+rho43*(dgdx*dxdg)
*        dF/dTaua
         vTau(1,iGrid)=vTau(1,iGrid)+ rho43*(dgdz*dzdt)
210      continue
c
c beta component
c
         rhoo=max(1.0D-24,rho(2,iGrid))
         if(rhoo.lt.Ta) goto 310

         rho43 = rhoo**F43
         rrho = 1.0d0/rhoo       ! reciprocal of rho
         rho13 = rho43*rrho
         rho53 = rhoo**F53
         rho83 = rho53*rhoo
         tauu=Tau(2,iGrid)
         Gamma=Sigma(3,iGrid)
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
         vRho(2,iGrid)=vRho(2,iGrid)+ f43*rho13*gx +
     &                  rho43*(dgdx*dxdr + dgdz*dzdr)
*        dF/dGammabb
         vSigma(3,iGrid)=vSigma(3,iGrid)+rho43*(dgdx*dxdg)
*        dF/dTaub
         vTau(2,iGrid)=vTau(2,iGrid)+ rho43*(dgdz*dzdt)
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
