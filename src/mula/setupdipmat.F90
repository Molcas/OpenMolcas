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

subroutine SetUpDipMat(DipMat,max_term,ipow,var,dip,trfName,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2,max_mInc, &
                       max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec,nDec,det0,r0,r1,r2,base,nnsiz,nOsc,nDimTot,nPolyTerm,ndata, &
                       nvar,MaxNumAt)
!  Purpose:
!    Performs a least squares fit of the transition dipole at the two
!    centra and at the inPolyTermediate oscillator. Calculates the matrix
!    elements of the transition dipole at these centra.
!
!  Input:
!    ipow       : Two dimensional integer array - terms of the
!                 polynomial.
!    var        : Real*8 two dimensional array - coordinates
!                 to be used in the fit.
!    dip        : Real*8 array - values of dipole at the
!                 coordinates contained in var.
!    trfName    : Character array - transformation associated with each
!                 internal coordinate.
!    max_term   : Integer - maximum order of the transition dipole terms.
!    W1,W2      : Real*8 two dimensional arrays - eigenvectors
!                 scaled by the square root of the eigenvalues.
!    C1,C2      : Real*8 two dimensional arrays - inverses
!                 of W1 and W2.
!    det1,det2  : Real*8 variables - determinants of C1 and C2.
!    r01,r02    : Real*8 arrays - coordinates of the two
!                 oscillators.
!    max_mOrd,
!    max_nOrd,
!    max_nOrd2
!    max_mInc,
!    max_nInc,
!    max_nInc2  : Integer variables
!    mMat,nMat,
!    mInc,nInc,
!    mDec,nDec  : Two dimensional integer arrays
!
!
!  Output:
!    DipMat     : Real*8 two dimensional array - contains the
!                 matrix elements of the transition dipole.

implicit real*8(a-h,o-z)
#include "dims.fh"
real*8 DipMat(0:nDimTot,0:nDimTot)
integer ipow(nPolyTerm,nvar)
real*8 var(ndata,nvar)
real*8 dip(ndata)
character*80 trfName(MaxNumAt)
integer mMat(0:mdim1,mdim2), mInc(0:mdim1,mdim2), mDec(0:mdim1,mdim2)
integer nMat(0:ndim1,ndim2), nInc(0:ndim1,ndim2), nDec(0:ndim1,ndim2)
real*8 C1(nOsc,nOsc), C2(nOsc,nOsc), W1(nOsc,nOsc), W2(nOsc,nOsc)
real*8 r01(nOsc), r02(nOsc)
real*8 r0(nOsc), r1(nOsc), r2(nOsc)
real*8 Base(nOsc,nOsc)
!real*8 max_err, stand_dev
#include "WrkSpc.fh"

! Initialize.
call GetMem('Dij','Allo','Real',ipDij,(max_mOrd+1)*(max_mOrd+1))
call GetMem('DijTrans','Allo','Real',ipDijTrans,(max_mOrd+1)*(max_mOrd+1))
call GetMem('C','Allo','Real',ipC,nOsc*nOsc)
call GetMem('W','Allo','Real',ipW,nOsc*nOsc)
call GetMem('L','Allo','Real',ipL,(max_mOrd+1)*(max_mOrd+1))
call GetMem('U','Allo','Real',ipU,(max_nOrd2+1)*(max_nOrd2+1))
call GetMem('FC','Allo','Real',ipFC,(max_mOrd+1)*(max_mOrd+1))
call GetMem('Sij','Allo','Real',ipSij,(max_mOrd+1)*(max_mOrd+1))
call GetMem('r0vec','Allo','Real',ipr0vec,nOsc)
call GetMem('alpha1','Allo','Real',ipalpha1,nOsc*nOsc)
call GetMem('alpha2','Allo','Real',ipalpha2,nOsc*nOsc)
call GetMem('beta','Allo','Real',ipbeta,nOsc*nOsc)
call GetMem('D1','Allo','Real',ipD1,nOsc)
call GetMem('D2','Allo','Real',ipD2,nOsc*nOsc)
call GetMem('D3','Allo','Real',ipD3,nOsc*nOsc*nOsc)
call GetMem('D4','Allo','Real',ipD4,nOsc*nOsc*nOsc*nOsc)
!----------------------------------------------------------------------!

call SetUpDipMat_a(DipMat,max_term,ipow,var,dip,trfName,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2,max_mInc, &
                   max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec,nDec,det0,r0,r1,r2,base,nnsiz,nOsc,nDimTot,nPolyTerm,ndata,nvar, &
                   MaxNumAt,Work(ipL),Work(ipU),Work(ipC),Work(ipSij),Work(ipr0vec),Work(ipalpha1),Work(ipalpha2),Work(ipbeta), &
                   Work(ipW),Work(ipDij),Work(ipDijTrans),Work(ipD1),Work(ipD2),Work(ipD3),Work(ipD4),Work(ipFC))

call GetMem('Dij','Free','Real',ipDij,(max_mOrd+1)*(max_mOrd+1))
call GetMem('DijTrans','Free','Real',ipDijTrans,(max_mOrd+1)*(max_mOrd+1))
call GetMem('C','Free','Real',ipC,nOsc*nOsc)
call GetMem('W','Free','Real',ipW,nOsc*nOsc)
call GetMem('L','Free','Real',ipL,(max_mOrd+1)*(max_mOrd+1))
call GetMem('U','Free','Real',ipU,(max_nOrd2+1)*(max_nOrd2+1))
call GetMem('FC','Free','Real',ipFC,(max_mOrd+1)*(max_mOrd+1))
call GetMem('Sij','Free','Real',ipSij,(max_mOrd+1)*(max_mOrd+1))
call GetMem('r0vec','Free','Real',ipr0vec,nOsc)
call GetMem('alpha1','Free','Real',ipalpha1,nOsc*nOsc)
call GetMem('alpha2','Free','Real',ipalpha2,nOsc*nOsc)
call GetMem('beta','Free','Real',ipbeta,nOsc*nOsc)
call GetMem('D1','Free','Real',ipD1,nOsc)
call GetMem('D2','Free','Real',ipD2,nOsc*nOsc)
call GetMem('D3','Free','Real',ipD3,nOsc*nOsc*nOsc)
call GetMem('D4','Free','Real',ipD4,nOsc*nOsc*nOsc*nOsc)

end subroutine SetUpDipMat
!####
subroutine SetUpDipMat_a(DipMat,max_term,ipow,var,dip,trfName,C1,W1,det1,r01,C2,W2,det2,r02,max_mOrd,max_nOrd,max_nOrd2,max_mInc, &
                         max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec,nDec,det0,r0,r1,r2,base,nnsiz,nOsc,nDimTot,nPolyTerm,ndata, &
                         nvar,MaxNumAt,L,U,C,Sij,r0vec,alpha1,alpha2,beta,W,Dij,DijTrans,D1,D2,D3,D4,FC)
!  Purpose:
!    Performs a least squares fit of the transition dipole at the two
!    centra and at the inPolyTermediate oscillator. Calculates the matrix
!    elements of the transition dipole at these centra.
!
!  Input:
!    ipow       : Two dimensional integer array - terms of the
!                 polynomial.
!    var        : Real*8 two dimensional array - coordinates
!                 to be used in the fit.
!    dip        : Real*8 array - values of dipole at the
!                 coordinates contained in var.
!    trfName    : Character array - transformation associated with each
!                 internal coordinate.
!    use_weight : Logical
!    max_term   : Integer - maximum order of the transition dipole terms.
!    W1,W2      : Real*8 two dimensional arrays - eigenvectors
!                 scaled by the square root of the eigenvalues.
!    C1,C2      : Real*8 two dimensional arrays - inverses
!                 of W1 and W2.
!    det1,det2  : Real*8 variables - determinants of C1 and C2.
!    r01,r02    : Real*8 arrays - coordinates of the two
!                 oscillators.
!    max_mOrd,
!    max_nOrd,
!    max_nOrd2
!    max_mInc,
!    max_nInc,
!    max_nInc2  : Integer variables
!    mMat,nMat,
!    mInc,nInc,
!    mDec,nDec  : Two dimensional integer arrays
!
!
!  Output:
!    DipMat     : Real*8 two dimensional array - contains the
!                 matrix elements of the transition dipole.

implicit real*8(a-h,o-z)
#include "dims.fh"
real*8 DipMat(0:nDimTot,0:nDimTot)
integer ipow(nPolyTerm,nvar)
real*8 var(ndata,nvar)
real*8 dip(ndata)
character*80 trfName(MaxNumAt)
logical use_weight
logical find_minimum
integer mMat(0:mdim1,mdim2), mInc(0:mdim1,mdim2), mDec(0:mdim1,mdim2)
integer nMat(0:ndim1,ndim2), nInc(0:ndim1,ndim2), nDec(0:ndim1,ndim2)
real*8 C1(nOsc,nOsc), C2(nOsc,nOsc), W1(nOsc,nOsc), W2(nOsc,nOsc)
real*8 r01(nOsc), r02(nOsc)
real*8 r0(nOsc), r1(nOsc), r2(nOsc)
real*8 Base(nOsc,nOsc)
real*8 max_err, stand_dev
real*8 Dij(0:max_mOrd,0:max_mOrd)
real*8 DijTrans(0:max_mOrd,0:max_mOrd)
real*8 C(nOsc,nOsc), W(nOsc,nOsc)
real*8 L(0:max_mOrd,0:max_mOrd)
real*8 U(0:max_nOrd2,0:max_nOrd2)
real*8 FC(0:max_mOrd,0:max_mOrd)
real*8 Sij(0:max_mOrd,0:max_mOrd)
real*8 r0vec(nOsc)
real*8 alpha1(nOsc,nOsc), alpha2(nOsc,nOsc), beta(nOsc,nOsc)
real*8 D1(nOsc)
real*8 D2(nOsc,nOsc)
real*8 D3(nOsc,nOsc,nOsc)
real*8 D4(nOsc,nOsc,nOsc,nOsc)
#include "WrkSpc.fh"

! Initialize.
find_minimum = .false.
use_weight = .false.

! Calculate terms of type
!
!                    dM  |
!                    --  | < i | Q  | j >,
!                    dQ  |        k
!                      k  Q
!                          0
! where M is the (transition) dipole moment, |i> and |j> are harmonic
! oscillator states, Q_0 is the equilibrium geometry and Q_k is the
! k:th normal coordinate.

l_C1 = nOsc
call Calc_r00(C1,C1,W1,W1,C,W,alpha1,alpha2,r0vec,r01,r01,det0,det1,det1,FC00,l_C1)
call FCval(C1,W1,det1,r01,C1,W1,det1,r01,Sij,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec, &
           nDec,C,W,det1,r0vec,L,U,FC00,alpha1,alpha2,beta,l_C1,nnsiz)
call GetMem('coef','Allo','Real',ipcoef,nPolyTerm)
l_r1 = nOsc
call PotFit(nPolyTerm,nvar,ndata,ipow,var,dip,Work(ipcoef),r1,l_r1,D0,D1,D2,D3,D4,trfName,stand_dev,max_err,find_minimum,max_term, &
            use_weight,nOsc,nOsc,nOsc)
call DipMatEl(Dij,W,L,U,FC00,nMat,ninc,ndec,D0,D1,D2,D3,D4,max_term,Base,ndim1,ndim2,max_mOrd,max_nOrd2)
!DipMat(0:max_mOrd,0:max_mOrd) = Dij
do iv=0,max_mOrd
  do jv=0,max_mOrd
    DipMat(iv,jv) = Dij(iv,jv)
  end do
end do

l_C2 = nOsc
call Calc_r00(C2,C2,W2,W2,C,W,alpha1,alpha2,r0vec,r02,r02,det0,det2,det2,FC00,l_C2)
call FCval(C2,W2,det2,r02,C2,W2,det2,r02,Sij,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec, &
           nDec,C,W,det2,r0vec,L,U,FC00,alpha1,alpha2,beta,l_C2,nnsiz)
l_r2 = nOsc
call PotFit(nPolyTerm,nvar,ndata,ipow,var,dip,Work(ipcoef),r2,l_r2,D0,D1,D2,D3,D4,trfName,stand_dev,max_err,find_minimum,max_term, &
            use_weight,nOsc,nOsc,nOsc)
call DipMatEl(Dij,W,L,U,FC00,nMat,ninc,ndec,D0,D1,D2,D3,D4,max_term,Base,ndim1,ndim2,max_mOrd,max_nOrd2)
!DipMat(max_mOrd+1:2*max_mOrd+1,max_mOrd+1:2*max_mOrd+1) = Dij
do iv=max_mOrd+1,2*max_mOrd+1
  do jv=max_mOrd+1,2*max_mOrd+1
    DipMat(iv,jv) = Dij(iv-max_mOrd-1,jv-max_mOrd-1)
  end do
end do

l_C1 = nOsc
call Calc_r00(C1,C2,W1,W2,C,W,alpha1,alpha2,r0vec,r01,r02,det0,det1,det2,FC00,l_C1)
call FCval(C1,W1,det1,r01,C2,W2,det2,r02,Sij,max_mOrd,max_nOrd,max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec, &
           nDec,C,W,det0,r0vec,L,U,FC00,alpha1,alpha2,beta,l_C1,nnsiz)
l_r0 = nOsc
call PotFit(nPolyTerm,nvar,ndata,ipow,var,dip,Work(ipcoef),r0,l_r0,D0,D1,D2,D3,D4,trfName,stand_dev,max_err,find_minimum,max_term, &
            use_weight,nOsc,nOsc,nOsc)
call GetMem('coef','Free','Real',ipcoef,nPolyTerm)
call DipMatEl(Dij,W,L,U,FC00,nMat,ninc,ndec,D0,D1,D2,D3,D4,max_term,Base,ndim1,ndim2,max_mOrd,max_nOrd2)
!DipMat(0:max_mOrd,max_mOrd+1:2*max_mOrd+1) = Dij
do iv=0,max_mOrd+1
  do jv=max_mOrd+1,2*max_mOrd+1
    DipMat(iv,jv) = Dij(iv,jv-max_mOrd-1)
  end do
end do

do iOrd=0,max_mOrd
  do jOrd=0,max_mOrd
    DijTrans(jOrd,iOrd) = Dij(iOrd,jOrd)
  end do
end do
do iv=0,max_mOrd
  do jv=0,max_mOrd
    DipMat(max_mOrd+1+iv,jv) = DijTrans(iv,jv)
  end do
end do

! Avoid unused argument warnings:
if (.false.) call Unused_real_array(FC)

end subroutine SetUpDipMat_a
