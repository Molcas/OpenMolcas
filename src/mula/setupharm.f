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
      Subroutine SetUpHarmDip(DipMat,max_term,m_max,n_max,
     &  mMat,mInc,mDec,
     &  nMat,nInc,nDec,C1,W1,det1,r01,C2,W2,det2,r02,
     &  C,W,det0,r00,TranDip,TranDipGrad,FC00,nnsiz,
     &   max_mOrd, max_nOrd, nOsc)
C!
C!  Purpose:
C!    Calculate the matrix elements of the transition dipole moment
C!    at the location of the intermediate oscillator.
C!
C!  Input:
C!    max_term   : Integer - maximum order of the transition dipole terms.
C!    W1,W2      : Real*8 two dimensional arrays - eigenvectors
C!                 scaled by the square root of the eigenvalues.
C!    C1,C2      : Real*8 two dimensional arrays - inverses
C!                 of W1 and W2.
C!    det1,det2  : Real*8 variables - determinants of C1 and C2.
C!    r01,r02    : Real*8 arrays - coordinates of the two
C!                 oscillators.
C!    Forcefield : Logical variable - whether or not to use transition
C!                 dipole from input.
C!
C!  Output:
C!    DipMat     : Real*8 two dimensional array - contains the
C!                 matrix elements of the transition dipole.
C!
C!  Uses:
C!    TabMod
C!    FCMod
C!
c       Use TabMod
cc      Use FCMod
      Implicit Real*8 ( a-h,o-z )
#include "dims.fh"
      Real*8 DipMat (0:max_mOrd,0:max_nOrd,0:3)
      Integer mMat(0:mdim1,mdim2),mInc(0:mdim1,mdim2),
     &  mDec(0:mdim1,mdim2)
      Integer nMat(0:ndim1,ndim2),nInc(0:ndim1,ndim2),
     &  nDec(0:ndim1,ndim2) !  (0:ndim1,ndim2)
      Real*8 C1(nosc,nosc),C2(nosc,nosc),
     &  W1(nosc,nosc),W2(nosc,nosc),C(nosc,nosc),W(nosc,nosc)
      Real*8 r01(nosc),r02(nosc),r00(nosc)
      Real*8  TranDipGrad (3,nosc)
      Real*8 TranDip (3)
      Integer  nvTabDim
#include "WrkSpc.fh"
C!
C!---- Initialize.
      Call TabDim2_drv(m_max,nosc,nvTabDim)
      m_max_ord  = nvTabDim-1
      Call TabDim2_drv(min(n_max,m_max+1),nosc,nvTabDim)

      mx_max_ord = nvTabDim-1
      Call TabDim2_drv(min(m_max,n_max+1),nosc,nvTabDim)
      nx_max_ord = nvTabDim-1
      Call TabDim2_drv(m_max-1,nosc,nvTabDim)
      max_mInc   = nvTabDim-1
      Call TabDim2_drv(n_max,nosc,nvTabDim)
      n_max_ord  = nvTabDim-1
      Call TabDim2_drv(n_max-1,nosc,nvTabDim)
      max_nInc   = nvTabDim-1
      Call GetMem('L','Allo','Real',
     &  ipL,(m_max_ord+1)*(nx_max_ord+1))
      Call GetMem('Sij','Allo','Real',
     &  ipSij,(m_max_ord+1)*(n_max_ord+1))

      Call GetMem('U','Allo','Real',ipU,(n_max_ord+1)*(mx_max_ord+1))
      Call GetMem('alpha1','Allo','Real',ipalpha1,nOsc*nOsc)
      Call GetMem('alpha2','Allo','Real',ipalpha2,nOsc*nOsc)
      Call GetMem('beta','Allo','Real',ipbeta,nOsc*nOsc)
C!
C!---- Calculate Franck-Condon factors.
      Call FCval(
     &  C1,        W1,         det1,     r01,      C2,
     &  W2,        det2,       r02,      Work(ipSij),      m_max_ord,
     &  n_max_ord, mx_max_ord, max_mInc, max_nInc, nx_max_ord,
     &  mMat,      nMat,       mInc,     nInc,     mDec,
     &  nDec,      C,          W,        det0,     r00,
     &  Work(ipL),         Work(ipU),          FC00,
     &  Work(ipalpha1),   Work(ipalpha2),
     &  Work(ipbeta),      nOsc,       nnsiz)
      Call GetMem('alpha1','Free','Real',ipalpha1,nOsc*nOsc)
      Call GetMem('alpha2','Free','Real',ipalpha2,nOsc*nOsc)
      Call GetMem('beta','Free','Real',ipbeta,nOsc*nOsc)
C!
C!---- Get the zeroth order contribution from transition dipole.
c       Dipmat = 0.0d0
      call dcopy_((max_mOrd+1)*(max_nOrd+1)*4, 0.0d0,0,Dipmat,1)
      Do iCar = 1,3
      do iv=0,max_mOrd
      do jv=0,max_nOrd
      DipMat(iv,jv,iCar) = TranDip(iCar)*
     &    Work(ipSij+iv+(m_max_ord+1)*jv)
      enddo
      enddo
      End Do
      Call GetMem('Sij','Free','Real',
     &  ipSij,(m_max_ord+1)*(n_max_ord+1))
C!
C! rt = sqrt(TranDip(1)**2+TranDip(2)**2+TranDip(3)**2)
C! If ( rt.gt.0.0d0 ) Then
C!    DipMat(0:,0:,0) = rt*Sij
C! End If
C!
C!---- Calculate LFU (just valid for tdm first derivatives)
      If ( max_term.eq.1 ) Then
      Call GetMem('temp1','Allo','Real',
     &    ipTemp1,(m_max_ord+1)*(mx_max_ord+1))
      Call GetMem('temp2','Allo','Real',
     &    ipTemp2,(m_max_ord+1)*(n_max_ord+1))
      If ( n_max.gt.m_max ) Then
      l_F=(m_max_ord+1)*(mx_max_ord+1)*3
      Call GetMem('F','Allo','Real',ipF,l_F)
      Call Fgenerator(
     &        nmat,Work(ipF),nInc,nDec,TranDipGrad,m_max_ord,
     &    mx_max_ord,nosc)
      Do iCar = 1,3
      Call DGEMM_('N','N',
     &            m_max_ord+1,mx_max_ord+1,m_max_ord+1,
     &            1.0d0,Work(ipL),m_max_ord+1,
     &            Work(ipF+(m_max_ord+1)*(mx_max_ord+1)*(iCar-1)),
     &            m_max_ord+1,
     &            0.0d0,Work(ipTemp1),m_max_ord+1)
      Call DGEMM_('N','T',
     &            m_max_ord+1,n_max_ord+1,mx_max_ord+1,
     &            1.0d0,Work(ipTemp1),m_max_ord+1,
     &            Work(ipU),n_max_ord+1,
     &            0.0d0,Work(ipTemp2),m_max_ord+1)

      do iv=0,max_mOrd
      do jv=0,max_nOrd
      DipMat(iv,jv,iCar) =
     &          DipMat(iv,jv,iCar)+
     &          Work(ipTemp2+iv+(m_max_ord+1)*jv)*FC00
      enddo
      enddo
      End Do
      Else
      l_F=(n_max_ord+1)*(nx_max_ord+1)*3
      Call GetMem('F','Allo','Real',ipF,l_F)
      Call Fgenerator(
     &          mmat,Work(ipF),mInc,mDec,trandipgrad,
     &        n_max_ord,nx_max_ord,nosc)
      Do iCar = 1,3
      Call DGEMM_('N','T',
     &            m_max_ord+1,n_max_ord+1,nx_max_ord+1,
     &            1.0d0,Work(ipL),m_max_ord+1,
     &            Work(ipF+(n_max_ord+1)*(nx_max_ord+1)*(iCar-1)),
     &            n_max_ord+1,
     &            0.0d0,Work(ipTemp1),m_max_ord+1)
      Call DGEMM_('N','T',
     &            m_max_ord+1,n_max_ord+1,n_max_ord+1,
     &            1.0d0,Work(ipTemp1),m_max_ord+1,
     &            Work(ipU),n_max_ord+1,
     &            0.0d0,Work(ipTemp2),m_max_ord+1)
c                DipMat(0:,0:,iCar) = DipMat(0:,0:,iCar)+Temp2*FC00
      do iv=0,max_mOrd
      do jv=0,max_nOrd
      DipMat(iv,jv,iCar) =
     &          DipMat(iv,jv,iCar)+
     &          Work(ipTemp2+iv+(m_max_ord+1)*jv)*FC00
      enddo
      enddo
      End Do
      End If
      Call GetMem('temp1','Free','Real',
     &    ipTemp1,(m_max_ord+1)*(mx_max_ord+1))
      Call GetMem('temp2','Free','Real',
     &    ipTemp2,(m_max_ord+1)*(n_max_ord+1))
      Call GetMem('F','Free','Real',ipF,l_F)

      End If
C!
      Call GetMem('U','Free','Real',ipU,(n_max_ord+1)*(mx_max_ord+1))
      Call GetMem('L','Free','Real',ipL,(m_max_ord+1)*(nx_max_ord+1))
C!
      End
