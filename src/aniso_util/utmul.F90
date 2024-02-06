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

use Constants, only: cZero, cOne
use Definitions, only: wp, u6

implicit none
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
  write(u6,'(A)') 'in UTMU2:   EXCH or N<=0 !!!'
  write(u6,*) 'EXCH=',EXCH
  write(u6,*) 'N   =',N
  call xquit(128)
end if

if (N > EXCH) then
  write(u6,'(A)') 'in UTMU2:   EXCH < N !!!'
  write(u6,*) 'EXCH=',EXCH
  write(u6,*) 'N   =',N
  write(u6,'(A)') 'Nothing is to be done >> Return'
  call xquit(128)
end if

R1 = dznrm2_(EXCH*EXCH,ML,1)
R2 = dznrm2_(N*N,Z,1)
if ((R1 < 1.0e-25_wp) .or. (R2 < 1.0e-25_wp)) then
  write(u6,'(A)') 'in UTMU2:   M or Z are empty!!!'
  write(u6,*) 'norm(M)=',R1
  write(u6,*) 'norm(Z)=',R2
  call xquit(128)
end if

if (DBG) then
  write(u6,'(A)') 'UTMU :: input moment'
  do i=1,EXCH
    do j=1,EXCH
      write(u6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'|ML|',j,'>',ML(i,j)
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
  call zgemm_('C','N',EXCH,EXCH,EXCH,cOne,Z,EXCH,ML,EXCH,cZero,TMP,EXCH)
  call zgemm_('N','N',EXCH,EXCH,EXCH,cOne,TMP,EXCH,Z,EXCH,cZero,ML,EXCH)

else

  call ZGEMM_('C','N',N,N,N,cOne,Z(1:N,1:N),N,ML(1:N,1:N),N,cZero,TMP,N)
  call zcopy_(N*N,[cZero],0,ML(1:N,1:N),1)
  call ZGEMM_('N','N',N,N,N,cOne,TMP,N,Z(1:N,1:N),N,cZero,ML(1:N,1:N),N)
  call ZGEMM_('C','N',N,EXCH,N,cOne,Z(1:N,1:N),N,ML(1:N,1:EXCH),N,cZero,TMP(1:N,1:EXCH),N)

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
  write(u6,'(A)') 'UTMU :: unitary transformtion matrix'
  do i=1,N
    do j=1,N
      write(u6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'| U |',j,'>',Z(i,j)
    end do
  end do
  write(u6,'(A)') 'UTMU :: output moment'
  do i=1,EXCH
    do j=1,EXCH
      write(u6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'|ML|',j,'>',ML(i,j)
    end do
  end do
end if

if (N < EXCH) call mma_deallocate(MTMP)
call mma_deallocate(TMP)

return

end subroutine utmul
