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
* Copyright (C) 1993,1999, Roland Lindh                                *
************************************************************************
      Subroutine TnsCtl(Wrk,nWrk,Coora,
     &                  nijkl,mabMax,mabMin,mcdMax,mcdMin,
     &                  HMtrxAB,HMtrxCD,la,lb,lc,ld,
     &                  iCmpa,jCmpb,kCmpc,lCmpd,
     &                  iShlla,jShllb,kShllc,lShlld,i_out)
************************************************************************
*                                                                      *
* Object: to transform the intermediate integral set directly to       *
*         the final integral set.                                      *
*         Note that the position in memory of the final set is not     *
*         fixed.                                                       *
*                                                                      *
*     Author: Roland Lindh, Dept. of Theoretical Chemistry,            *
*             University of Lund, SWEDEN                               *
*             Modified by R.L Februrary, 1999.                         *
************************************************************************
      use Real_Spherical
      use Basis_Info
      use Breit, only: nComp
      Implicit Real*8 (A-H,O-Z)
#include "itmax.fh"
#include "print.fh"
#include "real.fh"
      Parameter(lab=iTabMx*2+1,npMax=lab*(lab+1)*(lab+2)/6)
      Real*8 HMtrxAB(*),HMtrxCD(*)
      Real*8, Intent(inout) :: Wrk(nWrk)
      Real*8 Coora(3,4)
      Integer, Intent(out) :: i_out
      ![all others are intent(in)]

      Integer :: nDim
!
!     If nComp==1
!---- Integral are stored as e,f,IJKL in Wrk
!     If nComp/=1
!---- Integral are stored as ncomp,e,f,IJKL in Wrk
!
!---- Observe that Transf is false for s and p functions.

      ne=(mabMax-mabMin+1)
      nf=(mcdMax-mcdMin+1)
      nab=iCmpa*jCmpb
      ncd=kCmpc*lCmpd
      nDim=Max(ne*nf,nab*nf,nab*ncd)

      iW2=1
      iW3=1+nijkl*nDim

      If (nComp/=1) Then
         Wrk(iW3:iW3+ne*nf*nijkl-1)=Wrk(iW2:iW2+ne*nf*nijkl-1)
         Call DGetMO(Wrk(iW3),nComp,
     &               nComp,ne*nf*(nijkl/nComp),
     &               Wrk(iW2),ne*nf*(nijkl/nComp))
      End If

!---- If (ss|ss) integral exit.

      If (la+lb+lc+ld==0) Then
         i_out=1
         Return
      End If
*
*---- Transpose if no transformation is needed.
*
      If ((la*lb.eq.0).and.(lc*ld.eq.0).and.
     &    .Not.Shells(iShlla)%Transf .and.
     &    .Not.Shells(jShllb)%Transf .and.
     &    .Not.Shells(kShllc)%Transf .and.
     &    .Not.Shells(lShlld)%Transf ) Then
         Call DGeTMO(Wrk(iW2),ne*nf,ne*nf,nijkl,Wrk(iW3),nijkl)
         i_out=iW3
         Return
      End If
*
*---- Form matrix corresponding to the transfer equation, le,la,lb.
*     The matrix transforms directly from e0 cartesians to real
*     spherical harmonics.
*
      If (la+lb.eq.0) Then
         i_in=iW2
         i_out=iW3
      Else If ((la*lb.eq.0).and.
     &    .Not.Shells(iShlla)%Transf .and.
     &    .Not.Shells(jShllb)%Transf) Then
         Call DGeTMO(Wrk(iW2),ne,ne,nf*nijkl,Wrk(iW3),nf*nijkl)
         i_in=iW3
         i_out=iW2
      Else
*
*------- Now transform directly (e,[f,IJKL]) to ([f,IJKL],AB)
*        Int(lf*IJKL,lA*lB)=Int(le,lf*IJKL)*HMtrx(le,lA*lB)
*
         nfijkl = nf*nijkl
         Call Sp_Mlt(Wrk(iW2),ne,Wrk(iW3),nfijkl,HMtrxAB,
     &               iCmpa*jCmpb)
         i_in=iW3
         i_out=iW2
      End If
*
*
*---- Form matrix corresponding to the transfer equation, lf,lC,lD
*
      If (lc+ld.eq.0) Then
         i_out=i_in
      Else If ((lc*ld.eq.0).and.
     &    .Not.Shells(kShllc)%Transf .and.
     &    .Not.Shells(lShlld)%Transf) Then
         Call DGeTMO(Wrk(i_in),nf,nf,nijkl*iCmpa*jCmpb,Wrk(i_out),
     &               nijkl*iCmpa*jCmpb)
      Else
*
*------- Now transform directly (f,[IJKL,AB]) to ([IJKL,AB],CD)
*        Int(IJKL*lA*lB,lC*lD)=Int(lf,IJKL*lA*lB)*HMtrx(lf,lC*lD)
*
         nijklAB = nijkl*iCmpa*jCmpb
         Call Sp_Mlt(Wrk(i_in),nf,Wrk(i_out),nijklAB,HMtrxCD,
     &               kCmpc*lCmpd)
      End If
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(Coora)
      End If
      End
