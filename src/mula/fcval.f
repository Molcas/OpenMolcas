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
* Copyright (C) 1995, Niclas Forsberg                                  *
************************************************************************
C!-----------------------------------------------------------------------!
C!
c       Module FCMod
C!
C!  Contains:
C!    FCval     (C1,W1,det1,r01,C2,W2,det2,r02,FC,max_mOrd,max_nOrd,
C!               max_nOrd2,max_mInc,max_nInc,max_nInc2,mMat,nMat,mInc,
C!               nInc,mDec,nDec,C,W,det0,r00,L,U,FC00,alpha1,alpha2,beta)
C!
C!  Uses:
C!    TabMod
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
C!-----------------------------------------------------------------------!
C!
c       Contains



C!-----------------------------------------------------------------------!
C!

C!-----------------------------------------------------------------------!
C!
      Subroutine FCval(
     &  C1,        W1,        det1,       r01,      C2,
     &  W2,        det2,      r02,        FC,       max_mOrd,
     &  max_nOrd,  max_nOrd2, max_mInc,   max_nInc, max_nInc2,
     &  mMat,      nMat,      mInc,       nInc,     mDec,
     &  nDec,      C,         W,          det0,     r00,
     &  L,         U,         FC00,       alpha1,   alpha2,
     &  beta,      nOsc,      nnsiz)
C!
C!  Purpose:
C!    Calculate multidimensional Franck Condon factors.
C!
C!  Input:
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
C!  Output:
C!    W          : Real*8 two dimensional array - eigenvectors
C!                 of intermediate oscillator scaled by the square root
C!                 of the eigenvalues.
C!    C          : Real*8 two dimensional array - inverse
C!                 of W.
C!    det0       : Real*8 variable - determinant of C.
C!    r00        : Real*8 array - coordinates of intermediate
C!                 oscillator.
C!    L,U        : Real*8 two dimensional arrays
C!    FC00       : Real*8 variable - zero-zero overlap.
C!    FC         : Real*8
C!    alpha1     : Real*8 two dimensional array -
C!                 0.5*C1(T)*C1.
C!    alpha2     : Real*8 two dimensional array -
C!                 0.5*C2(T)*C2.
C!    beta       : Real*8 two dimensional array -
C!                 0.5*alpha1*alpha^(-1)*alpha2
C!
C!  Calls:
C!    Cholesky    (LinAlg)
C!    Dcopy       (Essl)
C!    DGEMM_      (Essl)
C!    Dool        (LinAlg)
C!    Dscal       (Essl)
C!
C!  Uses:
C!    Linalg
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1995.
C!
      Implicit Real*8 ( a-h,o-z )
#include "dims.fh"
      Real*8 FC(0:max_mord,0:max_nord)
      Real*8 C1(nosc,nosc),C2(nosc,nosc),W1(nosc,nosc)
      Real*8 W2(nosc,nosc),C(nosc,nosc),W(nosc,nosc)
      Real*8 L (0:max_mord,0:max_ninc2)
      Real*8 U (0:max_nord,0:max_nord2)
      Real*8 alpha1(nosc,nosc),alpha2(nosc,nosc),beta(nosc,nosc)
      Real*8 r00(nosc),r01(nosc),r02(nosc)
      Integer mMat(0:mdim1,mdim2),mInc(0:mdim1,mdim2),
     &  mDec(0:mdim1,mdim2)
      Integer nMat(0:ndim1,ndim2),nInc(0:ndim1,ndim2),
     &  nDec(0:ndim1,ndim2)
#include "WrkSpc.fh"
C!
C!---- Initialize.
      nMaxMat = max(max_mord+1,max_nord+1)
      nTabDim = max(nMaxMat,8)
      nOscSqr = nOsc**2
      Call GetMem('temp','Allo','Real',iptemp,nOscSqr)
      Call GetMem('temp1','Allo','Real',iptemp1,nOscSqr)
      Call GetMem('temp2','Allo','Real',iptemp2,nOscSqr)
C!
C!---- Setup sqr table.
      n=nTabDim+1
      Call GetMem('sqr','Allo','Real',ipsqr,n+1)
      Do i = 0,nTabDim+1
      Work(ipsqr+i) = sqrt(dble(i))
      End Do
C!
C!---- Calculate alpha1, alpha2 and alpha.
      Call GetMem('alpha','Allo','Real',ipalpha,nOscSqr)
      Call DGEMM_('T','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,C1,nOsc,
     &            C1,nOsc,
     &            0.0d0,alpha1,nOsc)
      call dscal_(nOscSqr,0.5d0,alpha1,1)
      Call DGEMM_('T','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,C2,nOsc,
     &            C2,nOsc,
     &            0.0d0,alpha2,nOsc)
      call dscal_(nOscSqr,0.5d0,alpha2,1)
c       temp = alpha1+alpha2
      call dcopy_(nOscSqr,alpha1,1,Work(iptemp),1)
      Call Daxpy_(nOscSqr,1.0d0,alpha2,1,Work(iptemp),1)
c       alpha = 0.5d0*temp
      call dcopy_(nOscSqr,0.0d0,0,Work(ipalpha),1)
      Call Daxpy_(nOscSqr,0.5d0,Work(iptemp),1,Work(ipalpha),1)

C! Call xxDgemul(C,nOsc,'T',C,nOsc,'N',alpha,nOsc,nOsc,nOsc,nOsc)
C! call dscal_(nOscSqr,0.5d0,alpha,1)
C!
C!---- Calculate C using a Cholesky factorization of 2*alpha.
C! Call Cholesky(temp,C)
C!
C!---- Calculate W.
C! call dcopy_(nOscSqr,0.0d0,0,W,1) ; call dcopy_(nOsc,1.0d0,0,W,nOsc+1)
C! temp = C ; Call Dool(temp,W,det0)
C!
C!---- Calculate r00.
      Call GetMem('r_temp1','Allo','Real',ipr_temp1,nOsc)
      Call GetMem('r_temp2','Allo','Real',ipr_temp2,nOsc)
C!
C!---- Calculate beta.
      Do i = 1,nOsc
      Do j = 1,nOsc
      Work(iptemp1+j+nOsc*(i-1)-1) = C1(i,j)
      End Do
      End Do
c       temp1 = alpha1
      call dcopy_(nOscSqr,alpha1,1,Work(iptemp1),1)
c       temp  = 2.0d0*alpha
      call dcopy_(nOscSqr,0.0d0,0,Work(iptemp),1)
      Call Daxpy_(nOscSqr,2.0d0,Work(ipalpha),1,Work(iptemp),1)

      Call Dool(Work(iptemp),nOsc,nOsc,Work(iptemp1),nOsc,nOsc,det)
      Call DGEMM_('N','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,alpha2,nOsc,
     &            Work(iptemp1),nOsc,
     &            0.0d0,beta,nOsc)
C!
C!---- Calculate FC00.
c       r_temp1 = r02-r01
      call dcopy_(nOsc,r01,1,Work(ipr_temp1),1)
      Call Daxpy_(nOsc,-1.0d0,r02,1,Work(ipr_temp1),1)

      Call DGEMM_('N','N',
     &            nOsc,1,nOsc,
     &            1.0d0,beta,nOsc,
     &            Work(ipr_temp1),nOsc,
     &            0.0d0,Work(ipr_temp2),nOsc)
      FC00_exp = Ddot_(nOsc,Work(ipr_temp1),1,Work(ipr_temp2),1)
      FC00 = (sqrt(det1)*sqrt(det2)/det0)*exp(-FC00_exp)
C!
C!---- Calculate A, B and d matrices.
      Call GetMem('A1','Allo','Real',ipA1,nOscSqr)
      Call GetMem('B1','Allo','Real',ipB1,nOscSqr)
      Call GetMem('A2','Allo','Real',ipA2,nOscSqr)
      Call GetMem('B2','Allo','Real',ipB2,nOscSqr)
      Call GetMem('d1','Allo','Real',ipd1,nOsc)
      Call GetMem('d2','Allo','Real',ipd2,nOsc)
      Call DGEMM_('N','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,C1,nOsc,
     &            W,nOsc,
     &            0.0d0,Work(ipA1),nOsc)
      Call DGEMM_('T','T',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,W1,nOsc,
     &            C,nOsc,
     &            0.0d0,Work(iptemp),nOsc)
c       B1 = A1-temp
      call dcopy_(nOscSqr,Work(ipA1),1,Work(ipB1),1)
      Call Daxpy_(nOscSqr,-1.0d0,Work(iptemp),1,Work(ipB1),1)


      Call DGEMM_('T','N',
     &            nOsc,1,nOsc,
     &            1.0d0,W1,nOsc,
     &            Work(ipr_temp2),nOsc,
     &            0.0d0,Work(ipd1),nOsc)
      const = Work(ipsqr+8)
      call dscal_(nOsc,const,Work(ipd1),1)
      Call DGEMM_('N','N',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,C2,nOsc,
     &            W,nOsc,
     &            0.0d0,Work(ipA2),nOsc)
      Call DGEMM_('T','T',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,W2,nOsc,
     &            C,nOsc,
     &            0.0d0,Work(iptemp),nOsc)
c       B2 = A2-temp
      call dcopy_(nOscSqr,Work(ipA2),1,Work(ipB2),1)
      Call Daxpy_(nOscSqr,-1.0d0,Work(iptemp),1,Work(ipB2),1)

      Call DGEMM_('T','N',
     &            nOsc,1,nOsc,
     &            1.0d0,W2,nOsc,
     &            Work(ipr_temp2),nOsc,
     &            0.0d0,Work(ipd2),nOsc)
      const =-Work(ipsqr+8)
      call dscal_(nOsc,const,Work(ipd2),1)
C!
C!---- Calculate A1B1T and A2B2T.
      Call GetMem('A1B1T','Allo','Real',ipA1B1T,nOscSqr)
      Call GetMem('A2B2T','Allo','Real',ipA2B2T,nOscSqr)

      Call DGEMM_('N','T',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,Work(ipA1),nOsc,
     &            Work(ipB1),nOsc,
     &            0.0d0,Work(ipA1B1T),nOsc)
      Call DGEMM_('N','T',
     &            nOsc,nOsc,nOsc,
     &            1.0d0,Work(ipA2),nOsc,
     &            Work(ipB2),nOsc,
     &            0.0d0,Work(ipA2B2T),nOsc)
C!
      Call GetMem('temp','Free','Real',iptemp,nOscSqr)
      Call GetMem('temp1','Free','Real',iptemp1,nOscSqr)
      Call GetMem('temp2','Free','Real',iptemp2,nOscSqr)
      Call GetMem('alpha','Free','Real',ipalpha,nOscSqr)
      Call GetMem('r_temp1','Free','Real',ipr_temp1,nOsc)
      Call GetMem('r_temp2','Free','Real',ipr_temp2,nOsc)


C!
C!---- Initialize L matrix.
      call dcopy_((max_mord+1)*(max_ninc2+1),0.0d0,0,L,1)
c       L = 0.0d0
      L(0,0) = 1.0d0
C!

C!---- If max_mOrd.gt.0 then set up L(m,0).
      If ( max_mOrd.gt.0 ) Then
      Do kOsc = 1,nOsc
      L(mInc(0,kOsc),0) = Work(ipd1+kOsc-1)
      End Do
      If ( max_mInc.gt.0 ) Then
      Do iOrd = 1,max_mInc
      kOsc_start = nOsc
      Do While (( mMat(iOrd,kOsc_start).eq.0 ).and.
     &                       ( kOsc_start.gt.1 ))
      kOsc_start = kOsc_start-1
      End Do
      Do kOsc = kOsc_start,nOsc
      Do lOsc = 1,nOsc
      If ( mMat(iOrd,lOsc).gt.0 ) Then
      L(mInc(iOrd,kOsc),0) = L(mInc(iOrd,kOsc),0)+
     &                          Work(ipsqr+mMat(iOrd,lOsc))*
     &                          Work(ipA1B1T+kOsc+nOsc*(lOsc-1)-1)*
     &                          L(mDec(iOrd,lOsc),0)
      End If
      End Do
      L(mInc(iOrd,kOsc),0) =
     &                    (L(mInc(iOrd,kOsc),0)+Work(ipd1+kOsc-1)*
     &                    L(iOrd,0))/
     &                    Work(ipsqr+mMat(mInc(iOrd,kOsc),kOsc))
      End Do
      End Do
      End If
C!
C!---- Use recursion formula to obtain the rest of L.
C!    Do jOrd = 1,max_nInc2
      Write(6,*)'FCVAL Test prints:'
      Write(6,*)'In this loop, indices jOrd and iOrd will vary.'
      Write(6,*)'They are used to address L(iOrd,jOrd).'
      Write(6,'(1x,a,2i8)')'(1) jOrd goes from 1, up to max_mord=',
     &  max_mord
      Write(6,'(1x,a,2i8)')'(2) iOrd goes from 1, up to max_mord=',
     & max_mord
cvv      stop
      Do jOrd = 1,max_mOrd
      Write(6,'(1x,a,2i8)')' jOrd=',jOrd
      lOsc = nOsc
      Do While (( mMat(jOrd,lOsc).eq.0 ).and.( lOsc.gt.1 ))
      lOsc = lOsc-1
      End Do
      Do iOrd = 0,max_mOrd
      Do kOsc = 1,nOsc
      If ( mMat(iOrd,kOsc).gt.0 ) Then
      L(iOrd,jOrd) = L(iOrd,jOrd)+
     &              (Work(ipsqr+mMat(iOrd,kOsc))/
     &               Work(ipsqr+mMat(jOrd,lOsc)))*
     &               Work(ipA1+kOsc+nOsc*(lOsc-1)-1)*
     &               L(mDec(iOrd,kOsc),mDec(jOrd,lOsc))
      End If
      End Do
      End Do
      End Do
      End If
C!
C!---- Initialize U matrix.
c       U = 0.0d0
      call dcopy_((max_nord+1)*(max_nord2+1),0.0d0,0,U,1)
      U(0,0) = 1.0d0
C!
C!---- If max_nOrd.gt.0 then set up U(n,0).
      If ( max_nOrd.gt.0 ) Then
      Do kOsc = 1,nOsc
      U(nInc(0,kOsc),0) = Work(ipd2+kOsc-1)
      End Do
      If ( max_nInc.gt.0 ) Then
      Do iOrd = 1,max_nInc
      kOsc_start = nOsc
      Do While (( nMat(iOrd,kOsc_start).eq.0 ).and.
     &            ( kOsc_start.gt.1 ))
      kOsc_start = kOsc_start-1
      End Do
      Do kOsc = kOsc_start,nOsc
      Do lOsc = 1,nOsc
      If ( nMat(iOrd,lOsc).gt.0 ) Then

      U(nInc(iOrd,kOsc),0) = U(nInc(iOrd,kOsc),0)+
     &                          Work(ipsqr+nMat(iOrd,lOsc))*
     &                          Work(ipA2B2T+kOsc+nOsc*(lOsc-1)-1)*
     &                          U(nDec(iOrd,lOsc),0)
      End If
      End Do

      U(nInc(iOrd,kOsc),0) =
     &                    (U(nInc(iOrd,kOsc),0)+WOrk(ipd2+kOsc-1)*
     &                    U(iOrd,0))/
     &                    Work(ipsqr+nMat(nInc(iOrd,kOsc),kOsc))
      End Do
      End Do
      End If
C!
C!---- Use recursion formula to obtain the rest of U.
      Do jOrd = 1,max_nOrd2
      lOsc = nOsc
      Do While (( nMat(jOrd,lOsc).eq.0 ).and.( lOsc.gt.1 ))
      lOsc = lOsc-1
      End Do
      Do iOrd = 0,max_nOrd
      Do kOsc = 1,nOsc
      If ( nMat(iOrd,kOsc).gt.0 ) Then
      U(iOrd,jOrd) = U(iOrd,jOrd)+
     &         (Work(ipsqr+nMat(iOrd,kOsc))/
     &          Work(ipsqr+nMat(jOrd,lOsc)))*
     &          Work(ipA2+kOsc+nOsc*(lOsc-1)-1)*
     &         U(nDec(iOrd,kOsc),nDec(jOrd,lOsc))
      End If
      End Do
      End Do
      End Do
      End If
C!
      Call GetMem('A1','Free','Real',ipA1,nOscSqr)
      Call GetMem('B1','Free','Real',ipB1,nOscSqr)
      Call GetMem('A2','Free','Real',ipA2,nOscSqr)
      Call GetMem('B2','Free','Real',ipB2,nOscSqr)
      Call GetMem('A1B1T','Free','Real',ipA1B1T,nOscSqr)
      Call GetMem('A2B2T','Free','Real',ipA2B2T,nOscSqr)
      Call GetMem('d1','Free','Real',ipd1,nOsc)
      Call GetMem('d2','Free','Real',ipd2,nOsc)
C!
C!---- Calculate Franck-Condon factors.
c       FC = 0.0d0
      call dcopy_((max_mord+1)*(max_nord+1),0.0d0,0,FC,1)
      Do jOrd = 0,max_nOrd
      Do iOrd = 0,max_mOrd
      sum = 0.0d0
      Do kOrd = 0,min(max_mOrd,max_nOrd)

      sum = sum+L(iOrd,kOrd)*U(jOrd,kOrd)
      End Do
      FC(iOrd,jOrd) = FC00*sum
      End Do
      End Do
C!
      Call GetMem('sqr','Free','Real',ipsqr,n+1)
C!
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_real_array(r00)
         Call Unused_integer(nnsiz)
      End If
      End
