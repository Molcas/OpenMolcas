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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: EXCH, N
complex(kind=wp), intent(in) :: Z(N,N)
complex(kind=wp), intent(inout) :: ML(EXCH,EXCH)
integer(kind=iwp) :: I, i1, J, j1
real(kind=wp) :: R1, R2
complex(kind=wp), allocatable :: MTMP(:,:), TMP(:,:), TMP2(:,:)
real(kind=wp), external :: dznrm2_

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

#ifdef _DEBUGPRINT_
write(u6,'(A)') 'UTMU :: input moment'
do i=1,EXCH
  do j=1,EXCH
    write(u6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'|ML|',j,'>',ML(i,j)
  end do
end do
#endif

call mma_allocate(TMP,N,EXCH,'TMP')
if (N /= EXCH) then
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

  call mma_allocate(TMP2,N,EXCH,'TMP2')
  TMP2(:,:) = ML(1:N,:)
  call ZGEMM_('C','N',N,N,N,cOne,Z,N,TMP2(:,1:N),EXCH,cZero,TMP(:,1:N),N)
  call ZGEMM_('N','N',N,N,N,cOne,TMP(:,1:N),N,Z,N,cZero,TMP2(:,1:N),EXCH)
  ML(1:N,1:N) = TMP2(:,1:N)
  call ZGEMM_('C','N',N,EXCH,N,cOne,Z,N,TMP2,EXCH,cZero,TMP,N)
  call mma_deallocate(TMP2)

  do I=1,N
    ML(I,N+1:) = TMP(I,N+1:)
    ML(N+1:,I) = conjg(TMP(I,N+1:))
  end do
  ML(N+1:,N+1:) = MTMP(1:EXCH-N,1:EXCH-N)

end if !N == exch

#ifdef _DEBUGPRINT_
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
#endif

if (N < EXCH) call mma_deallocate(MTMP)
call mma_deallocate(TMP)

return

end subroutine utmul
