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

subroutine utmu(EXCH,N,Z,M1,M2)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: EXCH, N
complex(kind=wp), intent(in) :: Z(N,N), M1(3,EXCH,EXCH)
complex(kind=wp), intent(out) :: M2(3,EXCH,EXCH)
integer(kind=iwp) :: I, L
real(kind=wp) :: R1, R2
complex(kind=wp), allocatable :: M1_TMP(:,:), M2_TMP(:,:), TMP(:,:)
real(kind=wp), external :: dznrm2_
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: J
#endif

if ((N <= 0) .or. (EXCH <= 0)) then
  write(u6,'(A)') 'in UTMU:   EXCH or N<=0 !!!'
  write(u6,*) 'EXCH=',EXCH
  write(u6,*) 'N   =',N
  call xFlush(u6)
  call abend()
end if

if (N > EXCH) then
  write(u6,'(A)') 'in UTMU:   EXCH < N !!!'
  write(u6,*) 'EXCH=',EXCH
  write(u6,*) 'N   =',N
  write(u6,'(A)') 'Nothing is to be done >> Return'
  call xFlush(u6)
  call abend()
end if

R1 = dznrm2_(3*EXCH*EXCH,M1,1)
R2 = dznrm2_(N*N,Z,1)
if ((R1 < 1.0e-25_wp) .or. (R2 < 1.0e-25_wp)) then
  write(u6,'(A)') 'in UTMU:   M1 or Z are empty!!!'
  write(u6,*) 'norm(M1)=',R1
  write(u6,*) 'norm(Z )=',R2
  return
end if

#ifdef _DEBUGPRINT_
write(u6,'(A)') 'UTMU :: input moment'
do i=1,EXCH
  do j=1,EXCH
    write(u6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'|M1_l|',j,'>',(M1(l,i,j),l=1,3)
  end do
end do
#endif

call mma_allocate(M1_TMP,N,EXCH,'M1_TMP')
call mma_allocate(M2_TMP,N,N,'M2_TMP')
call mma_allocate(TMP,N,EXCH,'TMP')

if (N == EXCH) then

  do L=1,3
    M1_TMP(:,:) = M1(L,:,:)
    call ZGEMM_('C','N',EXCH,EXCH,EXCH,cOne,Z,EXCH,M1_TMP,EXCH,cZero,TMP,EXCH)
    call ZGEMM_('N','N',EXCH,EXCH,EXCH,cOne,TMP,EXCH,Z,EXCH,cZero,M2_TMP,EXCH)
    M2(L,:,:) = M2_TMP(:,:)
  end do !L

else

  M2(:,:,:) = cZero

  do L=1,3
    M1_TMP(:,:) = M1(L,1:N,:)
    call ZGEMM_('C','N',N,N,N,cOne,Z,N,M1_TMP(:,1:N),N,cZero,TMP(:,1:N),N)
    call ZGEMM_('N','N',N,N,N,cOne,TMP(:,1:N),N,Z,N,cZero,M2_TMP,N)
    call ZGEMM_('C','N',N,EXCH,N,cOne,Z,N,M1_TMP,N,cZero,TMP,N)
    M2(L,1:N,1:N) = M2_TMP(:,:)

    do I=1,N
      M2(L,I,N+1:) = TMP(I,N+1:)
      M2(L,N+1:,I) = conjg(TMP(I,N+1:))
    end do
    M2(L,N+1:,N+1:) = M1(L,N+1:,N+1:)
  end do !L

end if !N == exch

#ifdef _DEBUGPRINT_
write(u6,'(A)') 'UTMU :: input moment'
do i=1,EXCH
  do j=1,EXCH
    write(u6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'|M_l|',j,'>',(M1(l,i,j),l=1,3)
  end do
end do
write(u6,'(A)') 'UTMU :: unitary transformtion matrix'
do i=1,N
  do j=1,N
    write(u6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'| U |',j,'>',(Z(i,j),l=1,3)
  end do
end do
write(u6,'(A)') 'UTMU :: output moment'
do i=1,EXCH
  do j=1,EXCH
    write(u6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'|M_l|',j,'>',(M2(l,i,j),l=1,3)
  end do
end do
#endif
call mma_deallocate(M1_TMP)
call mma_deallocate(M2_TMP)
call mma_deallocate(TMP)

return

end subroutine utmu
