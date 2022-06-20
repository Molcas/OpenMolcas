************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
C!-----------------------------------------------------------------------!
C!
      Subroutine SetUpDipMat(DipMat,max_term,ipow,var,dip,trfName,
     &                                C1,W1,det1,r01,C2,W2,det2,r02,
     &                                max_mOrd,max_nOrd,max_nOrd2,
     &                                max_mInc,max_nInc,max_nInc2,
     &                                mMat,nMat,mInc,nInc,mDec,nDec,
     &                                det0,r0,r1,r2,base,nnsiz,nOsc,
     &  nDimTot,nPolyTerm,ndata,nvar,MaxNumAt)
C!
C!  Purpose:
C!    Performs a least squares fit of the transition dipole at the two
C!    centra and at the inPolyTermediate oscillator. Calculates the matrix
C!    elements of the transition dipole at these centra.
C!
C!  Input:
C!    ipow       : Two dimensional integer array - terms of the
C!                 polynomial.
C!    var        : Real*8 two dimensional array - coordinates
C!                 to be used in the fit.
C!    dip        : Real*8 array - values of dipole at the
C!                 coordinates contained in var.
C!    trfName    : Character array - transformation associated with each
C!                 internal coordinate.
C!    max_term   : Integer - maximum order of the transition dipole terms.
C!    W1,W2      : Real*8 two dimensional arrays - eigenvectors
C!                 scaled by the square root of the eigenvalues.
C!    C1,C2      : Real*8 two dimensional arrays - inverses
C!                 of W1 and W2.
C!    det1,det2  : Real*8 variables - determinants of C1 and C2.
C!    r01,r02    : Real*8 arrays - coordinates of the two
C!                 oscillators.
C!    max_mOrd,
C!    max_nOrd,
C!    max_nOrd2
C!    max_mInc,
C!    max_nInc,
C!    max_nInc2  : Integer variables
C!    mMat,nMat,
C!    mInc,nInc,
C!    mDec,nDec  : Two dimensional integer arrays
C!
C!
C!  Output:
C!    DipMat     : Real*8 two dimensional array - contains the
C!                 matrix elements of the transition dipole.
C!
      Implicit Real*8 ( a-h,o-z )
#include "dims.fh"
      Real*8 DipMat( 0:nDimTot,0:nDimTot )
      Integer ipow(nPolyTerm,nvar)
      Real*8 var(ndata,nvar)
      Real*8 dip (ndata)
      Character*80 trfName(MaxNumAt)
      Integer mMat(0:mdim1,mdim2),mInc(0:mdim1,mdim2),
     &  mDec(0:mdim1,mdim2)
      Integer nMat(0:ndim1,ndim2),nInc(0:ndim1,ndim2),
     &  nDec(0:ndim1,ndim2)
      Real*8 C1(nOsc,nOsc),C2(nOsc,nOsc),
     &  W1(nOsc,nOsc),W2(nOsc,nOsc)
      Real*8 r01(nOsc),r02(nOsc)
      Real*8 r0(nOsc),r1(nOsc),r2(nOsc)
      Real*8 Base (nOsc,nOsc)
c       Real*8   max_err,stand_dev
#include "WrkSpc.fh"

C!
C!---- Initialize.
      Call GetMem('Dij','Allo','Real',
     &  ipDij,(max_mOrd+1)*(max_mOrd+1))
      Call GetMem('DijTrans','Allo','Real',
     &  ipDijTrans,(max_mOrd+1)*(max_mOrd+1))
      Call GetMem('C','Allo','Real',ipC,nOsc*nOsc)
      Call GetMem('W','Allo','Real',ipW,nOsc*nOsc)
      Call GetMem('L','Allo','Real',ipL,(max_mOrd+1)*(max_mOrd+1))
      Call GetMem('U','Allo','Real',ipU,(max_nOrd2+1)*(max_nOrd2+1))
      Call GetMem('FC','Allo','Real',ipFC,(max_mOrd+1)*(max_mOrd+1))
      Call GetMem('Sij','Allo','Real',
     &  ipSij,(max_mOrd+1)*(max_mOrd+1))
      Call GetMem('r0vec','Allo','Real',ipr0vec,nOsc)
      Call GetMem('alpha1','Allo','Real',ipalpha1,nOsc*nOsc)
      Call GetMem('alpha2','Allo','Real',ipalpha2,nOsc*nOsc)
      Call GetMem('beta','Allo','Real',ipbeta,nOsc*nOsc)
      Call GetMem('D1','Allo','Real',ipD1,nOsc)
      Call GetMem('D2','Allo','Real',ipD2,nOsc*nOsc)
      Call GetMem('D3','Allo','Real',ipD3,nOsc*nOsc*nOsc)
      Call GetMem('D4','Allo','Real',ipD4,nOsc*nOsc*nOsc*nOsc)
C!-----------------------------------------------------------------------!
C!
      Call SetUpDipMat_a(DipMat,max_term,ipow,var,dip,trfName,
     &                                C1,W1,det1,r01,C2,W2,det2,r02,
     &                                max_mOrd,max_nOrd,max_nOrd2,
     &                                max_mInc,max_nInc,max_nInc2,
     &                                mMat,nMat,mInc,nInc,mDec,nDec,
     &                                det0,r0,r1,r2,base,nnsiz,nOsc,
     &  nDimTot,nPolyTerm,ndata,nvar,MaxNumAt,
     & Work(ipL),Work(ipU),Work(ipC),Work(ipSij),Work(ipr0vec),
     & Work(ipalpha1),Work(ipalpha2),Work(ipbeta),Work(ipW),
     & Work(ipDij),Work(ipDijTrans),Work(ipD1),Work(ipD2),
     & Work(ipD3),Work(ipD4),Work(ipFC))

      Call GetMem('Dij','Free','Real',
     &  ipDij,(max_mOrd+1)*(max_mOrd+1))
      Call GetMem('DijTrans','Free','Real',
     &  ipDijTrans,(max_mOrd+1)*(max_mOrd+1))
      Call GetMem('C','Free','Real',ipC,nOsc*nOsc)
      Call GetMem('W','Free','Real',ipW,nOsc*nOsc)
      Call GetMem('L','Free','Real',ipL,(max_mOrd+1)*(max_mOrd+1))
      Call GetMem('U','Free','Real',ipU,(max_nOrd2+1)*(max_nOrd2+1))
      Call GetMem('FC','Free','Real',ipFC,(max_mOrd+1)*(max_mOrd+1))
      Call GetMem('Sij','Free','Real',ipSij,(max_mOrd+1)*(max_mOrd+1))
      Call GetMem('r0vec','Free','Real',ipr0vec,nOsc)
      Call GetMem('alpha1','Free','Real',ipalpha1,nOsc*nOsc)
      Call GetMem('alpha2','Free','Real',ipalpha2,nOsc*nOsc)
      Call GetMem('beta','Free','Real',ipbeta,nOsc*nOsc)
      Call GetMem('D1','Free','Real',ipD1,nOsc)
      Call GetMem('D2','Free','Real',ipD2,nOsc*nOsc)
      Call GetMem('D3','Free','Real',ipD3,nOsc*nOsc*nOsc)
      Call GetMem('D4','Free','Real',ipD4,nOsc*nOsc*nOsc*nOsc)
      End

C!-----------------------------------------------------------------------!
C!
      Subroutine SetUpDipMat_a(DipMat,max_term,ipow,var,dip,trfName,
     &                                C1,W1,det1,r01,C2,W2,det2,r02,
     &                                max_mOrd,max_nOrd,max_nOrd2,
     &                                max_mInc,max_nInc,max_nInc2,
     &                                mMat,nMat,mInc,nInc,mDec,nDec,
     &                                det0,r0,r1,r2,base,nnsiz,nOsc,
     &  nDimTot,nPolyTerm,ndata,nvar,MaxNumAt,
     & L,U,C,Sij,r0vec,alpha1,alpha2,beta,W,
     & Dij,DijTrans,D1,D2,D3,D4,FC)
C!
C!  Purpose:
C!    Performs a least squares fit of the transition dipole at the two
C!    centra and at the inPolyTermediate oscillator. Calculates the matrix
C!    elements of the transition dipole at these centra.
C!
C!  Input:
C!    ipow       : Two dimensional integer array - terms of the
C!                 polynomial.
C!    var        : Real*8 two dimensional array - coordinates
C!                 to be used in the fit.
C!    dip        : Real*8 array - values of dipole at the
C!                 coordinates contained in var.
C!    trfName    : Character array - transformation associated with each
C!                 internal coordinate.
C!    use_weight : Logical
C!    max_term   : Integer - maximum order of the transition dipole terms.
C!    W1,W2      : Real*8 two dimensional arrays - eigenvectors
C!                 scaled by the square root of the eigenvalues.
C!    C1,C2      : Real*8 two dimensional arrays - inverses
C!                 of W1 and W2.
C!    det1,det2  : Real*8 variables - determinants of C1 and C2.
C!    r01,r02    : Real*8 arrays - coordinates of the two
C!                 oscillators.
C!    max_mOrd,
C!    max_nOrd,
C!    max_nOrd2
C!    max_mInc,
C!    max_nInc,
C!    max_nInc2  : Integer variables
C!    mMat,nMat,
C!    mInc,nInc,
C!    mDec,nDec  : Two dimensional integer arrays
C!
C!
C!  Output:
C!    DipMat     : Real*8 two dimensional array - contains the
C!                 matrix elements of the transition dipole.
C!
      Implicit Real*8 ( a-h,o-z )
#include "dims.fh"
      Real*8 DipMat( 0:nDimTot,0:nDimTot )
      Integer ipow(nPolyTerm,nvar)
      Real*8 var(ndata,nvar)
      Real*8 dip (ndata)
      Character*80 trfName(MaxNumAt)
      Logical      use_weight
      Logical      find_minimum
      Integer mMat(0:mdim1,mdim2),mInc(0:mdim1,mdim2),
     &  mDec(0:mdim1,mdim2)
      Integer nMat(0:ndim1,ndim2),nInc(0:ndim1,ndim2),
     &  nDec(0:ndim1,ndim2)
      Real*8 C1(nOsc,nOsc),C2(nOsc,nOsc),W1(nOsc,nOsc),W2(nOsc,nOsc)
      Real*8 r01(nOsc),r02(nOsc)
      Real*8 r0(nOsc),r1(nOsc),r2(nOsc)
      Real*8 Base (nOsc,nOsc)
      Real*8   max_err,stand_dev
      Real*8  Dij(0:max_mOrd,0:max_mOrd)
      Real*8  DijTrans(0:max_mOrd,0:max_mOrd)
      Real*8  C(nOsc,nOsc),W(nOsc,nOsc)
      Real*8  L(0:max_mOrd,0:max_mOrd)
      Real*8  U(0:max_nOrd2,0:max_nOrd2)
      Real*8  FC(0:max_mOrd,0:max_mOrd)
      Real*8  Sij(0:max_mOrd,0:max_mOrd)
      Real*8  r0vec(nOsc)
      Real*8  alpha1(nOsc,nOsc),alpha2(nOsc,nOsc),beta(nOsc,nOsc)
      Real*8  D1(nOsc)
      Real*8  D2(nOsc,nOsc)
      Real*8  D3(nOsc,nOsc,nOsc)
      Real*8  D4(nOsc,nOsc,nOsc,nOsc)
#include "WrkSpc.fh"

C!
C!---- Initialize.
      find_minimum = .false.
      use_weight = .false.
C!
C!---- Calculate terms of type
C!
C!                        dM  |
C!                        --  | < i | Q  | j >,
C!                        dQ  |        k
C!                          k  Q
C!                              0
C!     where M is the (transition) dipole moment, |i> and |j> are harmonic
C!     oscillator states, Q_0 is the equilibrium geometry and Q_k is the
C!     k:th normal coordinate.
C!
      l_C1=nOsc
      Call Calc_r00(C1,C1,W1,W1,C,W,alpha1,alpha2,r0vec,r01,r01,
     &       det0,det1,det1,FC00,l_C1)
      Call FCval(C1,W1,det1,r01,C1,W1,det1,r01,Sij,max_mOrd,
     &       max_nOrd,max_nOrd2,
     &                    max_mInc,max_nInc,max_nInc2,mMat,nMat,
     &     mInc,nInc,mDec,nDec,
     &                C,W,det1,r0vec,L,U,FC00,alpha1,alpha2,beta,
     &     l_C1,nnsiz)
      Call GetMem('coef','Allo','Real',ipcoef,nPolyTerm)
      l_r1=nOsc
      Call PotFit(nPolyTerm,nvar,ndata,ipow,var,dip,
     &  Work(ipcoef),r1,l_r1,D0,
     &       D1,D2,D3,D4,trfName,
     &                     stand_dev,max_err,find_minimum,max_term,
     &     use_weight,
     &    nOsc,nOsc,nOsc)
      Call DipMatEl(Dij,W,L,U,FC00,nMat,ninc,ndec,D0,D1,D2,D3,D4,
     &       max_term,Base,ndim1,ndim2,max_mOrd,max_nOrd2)
c       DipMat(0:max_mOrd,0:max_mOrd) = Dij
      do iv=0,max_mOrd
      do jv=0,max_mOrd
      DipMat(iv,jv)=Dij(iv,jv)
      enddo
      enddo
C!
      l_C2=nOsc
      Call Calc_r00(C2,C2,W2,W2,C,W,alpha1,alpha2,r0vec,r02,r02,
     &       det0,det2,det2,FC00,l_C2)
      Call FCval(C2,W2,det2,r02,C2,W2,det2,r02,Sij,max_mOrd,
     &     max_nOrd,max_nOrd2,
     &     max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,nInc,mDec,nDec,
     &     C,W,det2,r0vec,L,U,FC00,alpha1,alpha2,beta,
     &     l_C2,nnsiz)
      l_r2=nOsc
      Call PotFit(nPolyTerm,nvar,ndata,ipow,var,dip,
     &  Work(ipcoef),r2,l_r2,D0,
     &       D1,D2,D3,D4,trfName,
     &       stand_dev,max_err,find_minimum,max_term,use_weight,
     &    nOsc,nOsc,nOsc)
      Call DipMatEl(Dij,W,L,U,FC00,nMat,ninc,ndec,D0,D1,D2,
     &       D3,D4,max_term,Base,ndim1,ndim2,max_mOrd,max_nOrd2)
c       DipMat(max_mOrd+1:2*max_mOrd+1,max_mOrd+1:2*max_mOrd+1) = Dij
      do iv=max_mOrd+1,2*max_mOrd+1
      do jv=max_mOrd+1,2*max_mOrd+1
      DipMat(iv,jv)=Dij(iv-max_mOrd-1,jv-max_mOrd-1)
      enddo
      enddo
C!
      l_C1=nOsc
      Call Calc_r00(C1,C2,W1,W2,C,W,alpha1,alpha2,r0vec,r01,r02,
     &       det0,det1,det2,FC00,l_C1)
      Call FCval(C1,W1,det1,r01,C2,W2,det2,r02,Sij,max_mOrd,
     &       max_nOrd,max_nOrd2,
     &                    max_mInc,max_nInc,max_nInc2,mMat,nMat,
     &     mInc,nInc,mDec,nDec,
     &     C,W,det0,r0vec,L,U,FC00,alpha1,alpha2,beta,
     &  l_C1,nnsiz)
      l_r0=nOsc
      Call PotFit(nPolyTerm,nvar,ndata,ipow,var,dip,Work(ipcoef),
     &  r0,l_r0,D0,
     &       D1,D2,D3,D4,trfName,
     &       stand_dev,max_err,find_minimum,max_term,use_weight,
     &   nOsc,nOsc,nOsc)
      Call GetMem('coef','Free','Real',ipcoef,nPolyTerm)
      Call DipMatEl(Dij,W,L,U,FC00,nMat,ninc,ndec,D0,D1,D2,
     &       D3,D4,max_term,Base,ndim1,ndim2,max_mOrd,max_nOrd2)
c       DipMat(0:max_mOrd,max_mOrd+1:2*max_mOrd+1) = Dij
      do iv=0,max_mOrd+1
      do jv=max_mOrd+1,2*max_mOrd+1
      DipMat(iv,jv)=Dij(iv,jv-max_mOrd-1)
      enddo
      enddo

      Do iOrd = 0,max_mOrd
      Do jOrd = 0,max_mOrd
      DijTrans(jOrd,iOrd) = Dij(iOrd,jOrd)
      End Do
      End Do
      do iv=0,max_mOrd
      do jv=0,max_mOrd
      DipMat(max_mOrd+1+iv,jv) = DijTrans(iv,jv)
      enddo
      enddo
C!
c Avoid unused argument warnings:
      If (.False.) Call Unused_real_array(FC)
      End
