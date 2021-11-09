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

subroutine SetUpHarmDip(DipMat,max_term,m_max,n_max,mMat,mInc,mDec,nMat,nInc,nDec,C1,W1,det1,r01,C2,W2,det2,r02,C,W,det0,r00, &
                        TranDip,TranDipGrad,FC00,nnsiz,max_mOrd,max_nOrd,nOsc)
!  Purpose:
!    Calculate the matrix elements of the transition dipole moment
!    at the location of the intermediate oscillator.
!
!  Input:
!    max_term   : Integer - maximum order of the transition dipole terms.
!    W1,W2      : Real*8 two dimensional arrays - eigenvectors
!                 scaled by the square root of the eigenvalues.
!    C1,C2      : Real*8 two dimensional arrays - inverses
!                 of W1 and W2.
!    det1,det2  : Real*8 variables - determinants of C1 and C2.
!    r01,r02    : Real*8 arrays - coordinates of the two
!                 oscillators.
!    Forcefield : Logical variable - whether or not to use transition
!                 dipole from input.
!
!  Output:
!    DipMat     : Real*8 two dimensional array - contains the
!                 matrix elements of the transition dipole.
!
!  Uses:
!    TabMod
!    FCMod

use Constants, only: Zero, One

!use TabMod
!use FCMod
implicit real*8(a-h,o-z)
#include "dims.fh"
real*8 DipMat(0:max_mOrd,0:max_nOrd,0:3)
integer mMat(0:mdim1,mdim2), mInc(0:mdim1,mdim2), mDec(0:mdim1,mdim2)
integer nMat(0:ndim1,ndim2), nInc(0:ndim1,ndim2), nDec(0:ndim1,ndim2)
real*8 C1(nosc,nosc), C2(nosc,nosc), W1(nosc,nosc), W2(nosc,nosc), C(nosc,nosc), W(nosc,nosc)
real*8 r01(nosc), r02(nosc), r00(nosc)
real*8 TranDipGrad(3,nosc)
real*8 TranDip(3)
integer nvTabDim
#include "WrkSpc.fh"

! Initialize.
call TabDim2_drv(m_max,nosc,nvTabDim)
m_max_ord = nvTabDim-1
call TabDim2_drv(min(n_max,m_max+1),nosc,nvTabDim)

mx_max_ord = nvTabDim-1
call TabDim2_drv(min(m_max,n_max+1),nosc,nvTabDim)
nx_max_ord = nvTabDim-1
call TabDim2_drv(m_max-1,nosc,nvTabDim)
max_mInc = nvTabDim-1
call TabDim2_drv(n_max,nosc,nvTabDim)
n_max_ord = nvTabDim-1
call TabDim2_drv(n_max-1,nosc,nvTabDim)
max_nInc = nvTabDim-1
call GetMem('L','Allo','Real',ipL,(m_max_ord+1)*(nx_max_ord+1))
call GetMem('Sij','Allo','Real',ipSij,(m_max_ord+1)*(n_max_ord+1))

call GetMem('U','Allo','Real',ipU,(n_max_ord+1)*(mx_max_ord+1))
call GetMem('alpha1','Allo','Real',ipalpha1,nOsc*nOsc)
call GetMem('alpha2','Allo','Real',ipalpha2,nOsc*nOsc)
call GetMem('beta','Allo','Real',ipbeta,nOsc*nOsc)

! Calculate Franck-Condon factors.
call FCval(C1,W1,det1,r01,C2,W2,det2,r02,Work(ipSij),m_max_ord,n_max_ord,mx_max_ord,max_mInc,max_nInc,nx_max_ord,mMat,nMat,mInc, &
           nInc,mDec,nDec,C,W,det0,r00,Work(ipL),Work(ipU),FC00,Work(ipalpha1),Work(ipalpha2),Work(ipbeta),nOsc,nnsiz)
call GetMem('alpha1','Free','Real',ipalpha1,nOsc*nOsc)
call GetMem('alpha2','Free','Real',ipalpha2,nOsc*nOsc)
call GetMem('beta','Free','Real',ipbeta,nOsc*nOsc)

! Get the zeroth order contribution from transition dipole.
!Dipmat = Zero
call dcopy_((max_mOrd+1)*(max_nOrd+1)*4,[Zero],0,Dipmat,1)
do iCar=1,3
  do iv=0,max_mOrd
    do jv=0,max_nOrd
      DipMat(iv,jv,iCar) = TranDip(iCar)*Work(ipSij+iv+(m_max_ord+1)*jv)
    end do
  end do
end do
call GetMem('Sij','Free','Real',ipSij,(m_max_ord+1)*(n_max_ord+1))

!rt = sqrt(TranDip(1)**2+TranDip(2)**2+TranDip(3)**2)
!if (rt > Zero) DipMat(0:,0:,0) = rt*Sij
!end if

! Calculate LFU (just valid for tdm first derivatives)
if (max_term == 1) then
  call GetMem('temp1','Allo','Real',ipTemp1,(m_max_ord+1)*(mx_max_ord+1))
  call GetMem('temp2','Allo','Real',ipTemp2,(m_max_ord+1)*(n_max_ord+1))
  if (n_max > m_max) then
    l_F = (m_max_ord+1)*(mx_max_ord+1)*3
    call GetMem('F','Allo','Real',ipF,l_F)
    call Fgenerator(nmat,Work(ipF),nInc,nDec,TranDipGrad,m_max_ord,mx_max_ord,nosc)
    do iCar=1,3
      call DGEMM_('N','N',m_max_ord+1,mx_max_ord+1,m_max_ord+1,One,Work(ipL),m_max_ord+1, &
                  Work(ipF+(m_max_ord+1)*(mx_max_ord+1)*(iCar-1)),m_max_ord+1,Zero,Work(ipTemp1),m_max_ord+1)
      call DGEMM_('N','T',m_max_ord+1,n_max_ord+1,mx_max_ord+1,One,Work(ipTemp1),m_max_ord+1,Work(ipU),n_max_ord+1,Zero, &
                  Work(ipTemp2),m_max_ord+1)

      do iv=0,max_mOrd
        do jv=0,max_nOrd
          DipMat(iv,jv,iCar) = DipMat(iv,jv,iCar)+Work(ipTemp2+iv+(m_max_ord+1)*jv)*FC00
        end do
      end do
    end do
  else
    l_F = (n_max_ord+1)*(nx_max_ord+1)*3
    call GetMem('F','Allo','Real',ipF,l_F)
    call Fgenerator(mmat,Work(ipF),mInc,mDec,trandipgrad,n_max_ord,nx_max_ord,nosc)
    do iCar=1,3
      call DGEMM_('N','T',m_max_ord+1,n_max_ord+1,nx_max_ord+1,One,Work(ipL),m_max_ord+1, &
                  Work(ipF+(n_max_ord+1)*(nx_max_ord+1)*(iCar-1)),n_max_ord+1,Zero,Work(ipTemp1),m_max_ord+1)
      call DGEMM_('N','T',m_max_ord+1,n_max_ord+1,n_max_ord+1,One,Work(ipTemp1),m_max_ord+1,Work(ipU),n_max_ord+1,Zero, &
                  Work(ipTemp2),m_max_ord+1)
      !DipMat(0:,0:,iCar) = DipMat(0:,0:,iCar)+Temp2*FC00
      do iv=0,max_mOrd
        do jv=0,max_nOrd
          DipMat(iv,jv,iCar) = DipMat(iv,jv,iCar)+Work(ipTemp2+iv+(m_max_ord+1)*jv)*FC00
        end do
      end do
    end do
  end if
  call GetMem('temp1','Free','Real',ipTemp1,(m_max_ord+1)*(mx_max_ord+1))
  call GetMem('temp2','Free','Real',ipTemp2,(m_max_ord+1)*(n_max_ord+1))
  call GetMem('F','Free','Real',ipF,l_F)

end if

call GetMem('U','Free','Real',ipU,(n_max_ord+1)*(mx_max_ord+1))
call GetMem('L','Free','Real',ipL,(m_max_ord+1)*(nx_max_ord+1))

end subroutine SetUpHarmDip
