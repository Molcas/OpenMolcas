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

subroutine utmul(EXCH,N,Z,ML)
! the same as UTMU, except being that the input M is
! being transformed and only one projection is done at a time.

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer, intent(in) :: EXCH, N
! one projection is done
complex(kind=8), intent(inout) :: ML(EXCH,EXCH)
complex(kind=8), intent(in) :: Z(N,N)
! local variables:
integer :: I, J, i1, j1
logical :: DBG
real(kind=8) :: dznrm2_, R1, R2
external :: dznrm2_
complex(kind=8), allocatable :: TMP(:,:), MTMP(:,:)

DBG = .false.

if ((N <= 0) .or. (EXCH <= 0)) then
  write(6,'(A)') 'in UTMU2:   EXCH or N<=0 !!!'
  write(6,*) 'EXCH=',EXCH
  write(6,*) 'N   =',N
  call xquit(128)
end if

if (N > EXCH) then
  write(6,'(A)') 'in UTMU2:   EXCH < N !!!'
  write(6,*) 'EXCH=',EXCH
  write(6,*) 'N   =',N
  write(6,'(A)') 'Nothing is to be done >> Return'
  call xquit(128)
end if

R1 = dznrm2_(EXCH*EXCH,ML,1)
R2 = dznrm2_(N*N,Z,1)
if ((R1 < 1.0e-25_wp) .or. (R2 < 1.0e-25_wp)) then
  write(6,'(A)') 'in UTMU2:   M or Z are empty!!!'
  write(6,*) 'norm(M)=',R1
  write(6,*) 'norm(Z)=',R2
  call xquit(128)
end if

if (DBG) then
  write(6,'(A)') 'UTMU :: input moment'
  do i=1,EXCH
    do j=1,EXCH
      write(6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'|ML|',j,'>',ML(i,j)
    end do
  end do
end if

call mma_allocate(TMP,EXCH,EXCH,'TMP')
if (N < EXCH) then
  call mma_allocate(MTMP,(EXCH-N),(EXCH-N),'MTMP')
  ! save the part which is not altered
  do i=N+1,EXCH
    do j=N+1,EXCH
      i1 = i-N
      j1 = j-N
      MTMP(i1,j1) = ML(i,j)
    end do
  end do
end if

if (N == EXCH) then
  call zcopy_(EXCH*EXCH,[(0.0_wp,0.0_wp)],0,TMP,1)
  call zgemm_('C','N',EXCH,EXCH,EXCH,(1.0_wp,0.0_wp),Z,EXCH,ML,EXCH,(0.0_wp,0.0_wp),TMP,EXCH)
  call zcopy_(EXCH*EXCH,[(0.0_wp,0.0_wp)],0,ML(:,:),1)
  call zgemm_('N','N',EXCH,EXCH,EXCH,(1.0_wp,0.0_wp),TMP,EXCH,Z,EXCH,(0.0_wp,0.0_wp),ML,EXCH)

else

  call zcopy_(EXCH*EXCH,[(0.0_wp,0.0_wp)],0,TMP,1)
  call ZGEMM_('C','N',N,N,N,(1.0_wp,0.0_wp),Z(1:N,1:N),N,ML(1:N,1:N),N,(0.0_wp,0.0_wp),TMP,N)
  call zcopy_(N*N,[(0.0_wp,0.0_wp)],0,ML(1:N,1:N),1)

  call ZGEMM_('N','N',N,N,N,(1.0_wp,0.0_wp),TMP,N,Z(1:N,1:N),N,(0.0_wp,0.0_wp),ML(1:N,1:N),N)

  call zcopy_(EXCH*EXCH,[(0.0_wp,0.0_wp)],0,TMP,1)
  call ZGEMM_('C','N',N,EXCH,N,(1.0_wp,0.0_wp),Z(1:N,1:N),N,ML(1:N,1:EXCH),N,(0.0_wp,0.0_wp),TMP(1:N,1:EXCH),N)

  do I=1,N
    do J=N+1,EXCH
      ML(I,J) = TMP(I,J)
      ML(J,I) = conjg(TMP(I,J))
    end do
  end do
  do i=N+1,EXCH
    do j=N+1,EXCH
      i1 = i-N
      j1 = j-N
      ML(i,j) = MTMP(i1,j1)
    end do
  end do

end if !N == exch

if (DBG) then
  write(6,'(A)') 'UTMU :: unitary transformtion matrix'
  do i=1,N
    do j=1,N
      write(6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'| U |',j,'>',Z(i,j)
    end do
  end do
  write(6,'(A)') 'UTMU :: output moment'
  do i=1,EXCH
    do j=1,EXCH
      write(6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'|ML|',j,'>',ML(i,j)
    end do
  end do
end if

if (N < EXCH) call mma_deallocate(MTMP)
call mma_deallocate(TMP)

return

end subroutine utmul
