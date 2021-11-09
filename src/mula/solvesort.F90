!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************

subroutine SolveSort(A,C,S,D,W1,W2,W0,C1,C2,C0,r01,r02,r00,icre,iann,nMat,nd1,nd2,OccNumMat,nosc,nDimTot)
!  Purpose:
!    Solve the secular equation AC = SCD.
!
!  Input:
!    A         : Real*8 two dimensional array
!    S         : Real*8 two dimensional array
!    W0,W1,W2  : Real*8 two dimensional arrays - eigenvectors
!                scaled by square root of harmonic frequencies.
!    C0,C1,C2  : Real*8 two dimensional arrays - inverses
!                of W0,W1 and W2.
!    r00,
!    r01,
!    r02       : Real*8 arrays - geometry in internal
!                coordinates.
!    nMat      : Two dimensional integer array
!
!  Output:
!    C         : Real*8 two dimensional array
!    D         : Real*8 array
!    OccNumMat : Real*8 two dimensional array - occupation
!                numbers for the different modes.
!
!  Calls:
!    Jacob
!    Jacord
!    Linalg
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

use Constants, only: Zero, One, Two, Eight, Half, auTocm
use Definitions, only: wp, u6

!use Linalg
implicit real*8(a-h,o-z)
real*8 A(nDimTot,nDimTot)
real*8 C(nDimTot,nDimTot)
real*8 S(nDimTot,nDimTot)
real*8 D(nDimTot)
real*8 W1(nosc,nosc), W2(nosc,nosc), C1(nosc,nosc), C2(nosc,nosc), W0(nosc,nosc), C0(nosc,nosc)
real*8 r01(nosc), r02(nosc), r00(nosc)
integer nMat(0:nd1,nd2)
integer icre(0:nd1,nd2)
integer iann(0:nd1,nd2)
real*8 OccNumMat(0:nDimTot-1,nOsc)
#include "WrkSpc.fh"

! Initialize.
n = nDimTot
nOscSqr = nOsc**2
nSqr = n**2
nSqrTri = n*(n+1)/2
call GetMem('Scr','Allo','Real',ipScr,nSqrTri)
call GetMem('T','Allo','Real',ipT,n*n)
call GetMem('Tmp1','Allo','Real',ipTmp1,n*n)
call GetMem('Asymm','Allo','Real',ipAsymm,n*n)

! Cholesky decomposition of S.
do i=1,n
  S(i,i) = S(i,i)+1.0e-6_wp
end do
call Cholesky(S,Work(ipTmp1),n)
call dcopy_(nSqr,[Zero],0,Work(ipT),1)
call dcopy_(n,[One],0,Work(ipT),n+1)
call Dool_MULA(Work(ipTmp1),n,n,Work(ipT),n,n,det)

! Make A symmetric and transform it to lower packed storage in Scratch.
call DGEMM_('N','N',n,n,n,One,A,n,Work(ipT),n,Zero,Work(ipTmp1),n)
call DGEMM_('T','N',n,n,n,One,Work(ipT),n,Work(ipTmp1),n,Zero,Work(ipAsymm),n)
k = 1
do i=1,n
  do j=1,i
    Work(ipScr+k-1) = Work(ipAsymm+i+n*(j-1)-1)
    k = k+1
  end do
end do

! Diagonalize Scratch.
call dcopy_(nSqr,[Zero],0,Work(ipTmp1),1)
call dcopy_(n,[One],0,Work(ipTmp1),n+1)
call Jacob(Work(ipScr),Work(ipTmp1),n,n)
call Jacord(Work(ipScr),Work(ipTmp1),n,n)

! Store the eigenvalues in array D and calculate C.
do i=1,n
  ii = i*(i+1)/2
  D(i) = Work(ipScr+ii-1)
end do

call DGEMM_('N','N',n,n,n,One,Work(ipT),n,Work(ipTmp1),n,Zero,C,n)
call GetMem('Scr','Free','Real',ipScr,nSqrTri)
call GetMem('T','Free','Real',ipT,n*n)
call GetMem('Tmp1','Free','Real',ipTmp1,n*n)
call GetMem('Asymm','Free','Real',ipAsymm,n*n)

! Calculate alpha1, alpha2 and alpha.
call GetMem('temp1','Allo','Real',iptemp1,nOscSqr)
call GetMem('temp2','Allo','Real',iptemp2,nOscSqr)

call GetMem('alpha','Allo','Real',ipalpha,nOscSqr)
call GetMem('alpha1','Allo','Real',ipalpha1,nOscSqr)
call GetMem('alpha2','Allo','Real',ipalpha2,nOscSqr)

call DGEMM_('T','N',nOsc,nOsc,nOsc,One,C1,nOsc,C1,nOsc,Zero,Work(ipalpha1),nOsc)
call dscal_(nOscSqr,Half,Work(ipalpha1),1)
call DGEMM_('T','N',nOsc,nOsc,nOsc,One,C2,nOsc,C2,nOsc,Zero,Work(ipalpha2),nOsc)
call dscal_(nOscSqr,Half,Work(ipalpha2),1)
call DGEMM_('T','N',nOsc,nOsc,nOsc,One,C0,nOsc,C0,nOsc,Zero,Work(ipalpha),nOsc)
call dscal_(nOscSqr,Half,Work(ipalpha),1)

! Calculate beta.
call GetMem('beta','Allo','Real',ipbeta,nOscSqr)
call GetMem('r_temp1','Allo','Real',ipr_temp1,nOsc)
call GetMem('r_temp2','Allo','Real',ipr_temp2,nOsc)

!temp1 = alpha1
call dcopy_(nOscSqr,Work(ipalpha1),1,Work(iptemp1),1)
!temp2 = Two*alpha
call dcopy_(nOscSqr,Work(ipalpha),1,Work(iptemp2),1)
call dscal_(nOscSqr,Two,Work(iptemp2),1)

call Dool_MULA(Work(iptemp2),nOsc,nOsc,Work(iptemp1),nOsc,nOsc,det)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,Work(ipalpha2),nOsc,Work(iptemp1),nOsc,Zero,WOrk(ipbeta),nOsc)

! Calculate A, B and d matrices.
call GetMem('A1','Allo','Real',ipA1,nOscSqr)
call GetMem('B1','Allo','Real',ipB1,nOscSqr)
call GetMem('A2','Allo','Real',ipA2,nOscSqr)
call GetMem('B2','Allo','Real',ipB2,nOscSqr)
call GetMem('d1','Allo','Real',ipd1,nOsc)
call GetMem('d2','Allo','Real',ipd2,nOsc)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,C1,nOsc,W0,nOsc,Zero,Work(ipA1),nOsc)
call DGEMM_('T','T',nOsc,nOsc,nOsc,One,W1,nOsc,C0,nOsc,Zero,Work(iptemp1),nOsc)
call dcopy_(nOscSqr,Work(ipA1),1,Work(ipB1),1)
call daxpy_(nOscSqr,-One,Work(iptemp1),1,Work(ipB1),1)
!B1 = A1-temp1
!r_temp1 = r01-r00
call dcopy_(nOscSqr,r01,1,Work(ipr_temp1),1)
call daxpy_(nOscSqr,-One,r00,1,Work(ipr_temp1),1)
call DGEMM_('N','N',nOsc,1,nOsc,One,Work(ipbeta),nOsc,Work(ipr_temp1),nOsc,Zero,Work(ipr_temp2),nOsc)
call DGEMM_('T','N',nOsc,1,nOsc,One,W1,nOsc,Work(ipr_temp2),nOsc,Zero,Work(ipd1),nOsc)
call dscal_(nOsc,sqrt(Eight),Work(ipd1),1)
call DGEMM_('N','N',nOsc,nOsc,nOsc,One,C2,nOsc,W0,nOsc,Zero,Work(ipA2),nOsc)
call DGEMM_('T','T',nOsc,nOsc,nOsc,One,W2,nOsc,C0,nOsc,Zero,Work(iptemp2),nOsc)
call dcopy_(nOscSqr,Work(ipA2),1,Work(ipB2),1)
call daxpy_(nOscSqr,-One,Work(iptemp2),1,Work(ipB2),1)
!B2 = A2-temp2
!r_temp1 = r02-r00
call dcopy_(nOscSqr,r02,1,Work(ipr_temp1),1)
call daxpy_(nOscSqr,-One,r00,1,Work(ipr_temp1),1)
call DGEMM_('N','N',nOsc,1,nOsc,One,Work(ipbeta),nOsc,Work(ipr_temp1),nOsc,Zero,Work(ipr_temp2),nOsc)
call DGEMM_('T','N',nOsc,1,nOsc,One,W2,nOsc,Work(ipr_temp2),nOsc,Zero,Work(ipd2),nOsc)
call dscal_(nOsc,sqrt(Eight),Work(ipd2),1)

call GetMem('temp1','Free','Real',iptemp1,nOscSqr)
call GetMem('temp2','Free','Real',iptemp2,nOscSqr)

call GetMem('beta','Free','Real',ipbeta,nOscSqr)
call GetMem('r_temp1','Free','Real',ipr_temp1,nOsc)
call GetMem('r_temp2','Free','Real',ipr_temp2,nOsc)
call GetMem('alpha','Free','Real',ipalpha,nOscSqr)
call GetMem('alpha1','Free','Real',ipalpha1,nOscSqr)
call GetMem('alpha2','Free','Real',ipalpha2,nOscSqr)

m = (n/2)-1

call GetMem('c_col','Allo','Real',ipc_col,n)
call GetMem('SC_col','Allo','Real',ipSC_col,n)
call dcopy_(nDimTot*nOsc,[Zero],0,OccNumMat,1)
!OccNumMat = Zero
do k=1,n
  do jOsc=1,nosc
    !C_col = Zero
    call dcopy_(n,[Zero],0,Work(ipc_col),1)
    do l=0,m
      do iOsc=1,nOsc
        lc = iCre(l,iosc)
        la = iAnn(l,iosc)
        if (lc >= 0) Work(ipC_col+lc) = Work(ipC_col+lc)+sqrt(real(nmat(lc,iosc),kind=wp))*c(l+1,k)*Work(ipB1+iosc+nOsc*(josc-1)-1)
        if (la >= 0) Work(ipC_col+la) = Work(ipC_col+la)+sqrt(real(nmat(l,iosc),kind=wp))*c(l+1,k)*Work(ipA1+iosc+nOsc*(josc-1)-1)
      end do
      Work(ipC_col+l) = Work(ipC_col+l)-c(l+1,k)*Work(ipd1+josc-1)
    end do
    do l=0,m
      do iOsc=1,nosc
        lc = icre(l,iosc)
        la = iann(l,iosc)
        if (lc >= 0) Work(ipC_col+m+lc+1) = Work(ipC_col+m+lc+1)+sqrt(real(nmat(lc,iosc),kind=wp))*c(m+l+2,k)* &
                                            Work(ipB2+iosc+nOsc*(josc-1)-1)
        if (la >= 0) Work(ipC_col+m+la+1) = Work(ipC_col+m+la+1)+sqrt(real(nmat(l,iosc),kind=wp))*c(m+l+2,k)* &
                                            Work(ipA2+iosc+nOsc*(josc-1)-1)
      end do
      Work(ipC_col+m+l+1) = Work(ipC_col+m+l+1)-c(m+l+2,k)*Work(ipd2+josc-1)
    end do
    call DGEMM_('N','N',n,1,n,One,S,n,Work(ipC_col),n,Zero,Work(ipSC_col),n)
    OccNumMat(k-1,jOsc) = Ddot_(n,Work(ipC_col),1,Work(ipSC_col),1)
  end do
end do
write(u6,*)
write(u6,*) 'Frequencies'
write(u6,'(20i6)') (int((D(i+1)-D(1))*auTocm),i=1,m)
write(u6,*)

call GetMem('d1','Free','Real',ipd1,nOsc)
call GetMem('d2','Free','Real',ipd2,nOsc)
call GetMem('c_col','Free','Real',ipc_col,n)
call GetMem('SC_col','Free','Real',ipSC_col,n)
call GetMem('A1','Free','Real',ipA1,nOscSqr)
call GetMem('B1','Free','Real',ipB1,nOscSqr)
call GetMem('A2','Free','Real',ipA2,nOscSqr)
call GetMem('B2','Free','Real',ipB2,nOscSqr)

end subroutine SolveSort
