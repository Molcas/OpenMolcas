!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
!!-----------------------------------------------------------------------!
!!
      Subroutine SetUpDipMat2(DipMat,max_term,C1,W1,det1,r01,           &
     &       C2,W2,det2,r02,                                            &
     &       max_mOrd,max_nOrd,max_nOrd2,                               &
     &       max_mInc,max_nInc,max_nInc2,                               &
     &       mMat,nMat,mInc,nInc,mDec,nDec,                             &
     &       det0,r0,r1,r2,base,TranDip,TranDipGrad,nnsiz,              &
     &  nOsc,nDimTot)
!!
!!  Purpose:
!!    Performs a least squares fit of the transition dipole at the two
!!    centra and at the intermediate oscillator. Calculates the matrix
!!    elements of the transition dipole at these centra.
!!
!!  Input:
!!    ipow       : Two dimensional integer array - terms of the
!!                 polynomial.
!!    var        : Real*8 two dimensional array - coordinates
!!                 to be used in the fit.
!!    dip        : Real*8 array - values of dipole at the
!!                 coordinates contained in var.
!!    trfName    : Character array - transformation associated with each
!!                 internal coordinate.
!!    use_weight : Logical
!!    max_term   : Integer - maximum order of the transition dipole terms.
!!    W1,W2      : Real*8 two dimensional arrays - eigenvectors
!!                 scaled by the square root of the eigenvalues.
!!    C1,C2      : Real*8 two dimensional arrays - inverses
!!                 of W1 and W2.
!!    det1,det2  : Real*8 variables - determinants of C1 and C2.
!!    r01,r02    : Real*8 arrays - coordinates of the two
!!                 oscillators.
!!    max_mOrd,
!!    max_nOrd,
!!    max_nOrd2
!!    max_mInc,
!!    max_nInc,
!!    max_nInc2  : Integer variables
!!    mMat,nMat,
!!    mInc,nInc,
!!    mDec,nDec  : Two dimensional integer arrays
!!
!!
!!  Output:
!!    DipMat     : Real*8 two dimensional array - contains the
!!                 matrix elements of the transition dipole.
!!
      Implicit Real*8 ( a-h,o-z )
#include "dims.fh"
      Real*8 DipMat( 0:nDimTot,0:nDimTot )
      Integer mMat(0:mdim1,mdim2),mInc(0:mdim1,mdim2),                  &
     &  mDec(0:mdim1,mdim2)
      Integer nMat(0:ndim1,ndim2),nInc(0:ndim1,ndim2),                  &
     &  nDec(0:ndim1,ndim2)
      Real*8 C1(nOsc,nOsc),C2(nOsc,nOsc),W1(nOsc,nOsc),W2(nOsc,nOsc)
      Real*8 r01(nOsc),r02(nOsc)
      Real*8 r0(nOsc),r1(nOsc),r2(nOsc)
      Real*8 Base (nOsc,nOsc)
      Real*8 TranDip(3), TranDipGrad (nOsc)
      Real*8 D0(3)
!       Real*8    max_err,stand_dev
#include "WrkSpc.fh"
!!
!!---- Initialize.
      call GetMem('Dij','Allo','Real',                                  &
     &  ipDij,(max_mOrd+1)*(max_mOrd+1))

      call GetMem('C','Allo','Real',ipC,nOsc*nOsc)
      call GetMem('W','Allo','Real',ipW,nOsc*nOsc)

      call GetMem('L','Allo','Real',ipL,(max_mOrd+1)*(max_mOrd+1))
      call GetMem('U','Allo','Real',ipU,(max_nOrd2+1)*(max_nOrd2+1))
      call GetMem('Sij','Allo','Real',ipSij,(max_mOrd+1)*(max_nOrd+1))
      call GetMem('r0vec','Allo','Real',ipr0vec,nOsc)
      call GetMem('alpha1','Allo','Real',ipalpha1,nOsc*nOsc)
      call GetMem('alpha2','Allo','Real',ipalpha2,nOsc*nOsc)
      call GetMem('beta','Allo','Real',ipbeta,nOsc*nOsc)
      call GetMem('D1','Allo','Real',ipD1,nOsc)
      call GetMem('D2','Allo','Real',ipD2,nOsc*nOsc)
      call GetMem('D3','Allo','Real',ipD3,nOsc*nOsc*nOsc)
      call GetMem('D4','Allo','Real',ipD4,nOsc*nOsc*nOsc*nOsc)
!!
!!---- Calculate terms of type
!!
!!                        dM  |
!!                        --  | < i | Q  | j >,
!!                        dQ  |        k
!!                          k  Q
!!                              0
!!     where M is the (transition) dipole moment, |i> and |j> are harmonic
!!     oscillator states, Q_0 is the equilibrium geometry and Q_k is the
!!     k:th normal coordinate.
!!
      l_C1=nOsc
      Call Calc_r00(C1,C2,W1,W2,Work(ipC),Work(ipW),                    &
     &  Work(ipalpha1),Work(ipalpha2),Work(ipr0vec),r01,                &
     &       r02,det0,det1,det2,FC00,l_C1)
      Call FCval(C1,W1,det1,r01,C2,W2,det2,r02,Work(ipSij),             &
     &  max_mOrd,                                                       &
     &       max_nOrd,max_nOrd2,                                        &
     &     max_mInc,max_nInc,max_nInc2,mMat,                            &
     &     nMat,mInc,nInc,mDec,nDec,                                    &
     &     Work(ipC),Work(ipW),det0,Work(ipr0vec),Work(ipL),            &
     & Work(ipU),FC00,Work(ipalpha1),Work(ipalpha2),Work(ipbeta),       &
     &  l_C1,nnsiz)
      D0(1) = TranDip(1)
      D0(2) = TranDip(2)
      D0(3) = TranDip(3)
      call dcopy_(nOsc,TranDipGrad,1,Work(ipD1),1)
      call dcopy_(nOsc*nOsc,[0.0d0],0,Work(ipD2),1)
      call dcopy_(nOsc*nOsc*nOsc,[0.0d0],0,Work(ipD3),1)
      call dcopy_(nOsc*nOsc*nOsc*nOsc,[0.0d0],0,Work(ipD4),1)

      Call DipMatEl(Work(ipDij),Work(ipW),Work(ipL),Work(ipU),          &
     &  FC00,nMat,nInc,nDec,D0(1),Work(ipD1),Work(ipD2),                &
     &       Work(ipD3),Work(ipD4),max_term,Base,ndim1,                 &
     & ndim2,max_mOrd,max_nOrd2)
!       DipMat(0:max_mOrd,0:max_mOrd) = Dij
      call dcopy_((max_mOrd+1)*(max_mOrd+1),Work(ipDij),1,DipMat,1)
!!
      call GetMem('Sij','Free','Real',ipSij,(max_mOrd+1)*(max_nOrd+1))
      call GetMem('L','Free','Real',ipL,(max_mOrd+1)*(max_mOrd+1))
      call GetMem('U','Free','Real',ipU,(max_nOrd2+1)*(max_nOrd2+1))
      call GetMem('C','Free','Real',ipC,nOsc*nOsc)
      call GetMem('W','Free','Real',ipW,nOsc*nOsc)
      call GetMem('Dij','Free','Real',ipDij,(max_mOrd+1)*(max_mOrd+1))
      call GetMem('D1','Free','Real',ipD1,nOsc)
      call GetMem('D2','Free','Real',ipD2,nOsc*nOsc)
      call GetMem('D3','Free','Real',ipD3,nOsc*nOsc*nOsc)
      call GetMem('D4','Free','Real',ipD4,nOsc*nOsc*nOsc*nOsc)
      call GetMem('r0vec','Free','Real',ipr0vec,nOsc)
      call GetMem('alpha1','Free','Real',ipalpha1,nOsc*nOsc)
      call GetMem('alpha2','Free','Real',ipalpha2,nOsc*nOsc)
      call GetMem('beta','Free','Real',ipbeta,nOsc*nOsc)
!!
! Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(r0)
         Call Unused_real_array(r1)
         Call Unused_real_array(r2)
      End If
      End
