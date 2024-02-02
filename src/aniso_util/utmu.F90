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

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer, intent(in) :: EXCH, N
complex(kind=8), intent(in) :: M1(3,EXCH,EXCH)
complex(kind=8), intent(in) :: Z(N,N)
complex(kind=8), intent(out) :: M2(3,EXCH,EXCH)
! local variables:
integer :: L, I, J
logical :: DBG
real(kind=8) :: dznrm2_, R1, R2
external :: dznrm2_
complex(kind=8), allocatable :: TMP(:,:)

DBG = .false.

if ((N <= 0) .or. (EXCH <= 0)) then
  write(6,'(A)') 'in UTMU:   EXCH or N<=0 !!!'
  write(6,*) 'EXCH=',EXCH
  write(6,*) 'N   =',N
  call xFlush(6)
  call xquit(128)
end if

if (N > EXCH) then
  write(6,'(A)') 'in UTMU:   EXCH < N !!!'
  write(6,*) 'EXCH=',EXCH
  write(6,*) 'N   =',N
  write(6,'(A)') 'Nothing is to be done >> Return'
  call xFlush(6)
  call xquit(128)
end if

R1 = dznrm2_(3*EXCH*EXCH,M1,1)
R2 = dznrm2_(N*N,Z,1)
if ((R1 < 1.0e-25_wp) .or. (R2 < 1.0e-25_wp)) then
  write(6,'(A)') 'in UTMU:   M1 or Z are empty!!!'
  write(6,*) 'norm(M1)=',R1
  write(6,*) 'norm(Z )=',R2
  return
end if

if (DBG) then
  write(6,'(A)') 'UTMU :: input moment'
  do i=1,EXCH
    do j=1,EXCH
      write(6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'|M1_l|',j,'>',(M1(l,i,j),l=1,3)
    end do
  end do
end if

call mma_allocate(TMP,EXCH,EXCH,'TMP')
call zcopy_(3*EXCH*EXCH,[(0.0_wp,0.0_wp)],0,M2,1)

if (N == EXCH) then

  do L=1,3
    call zcopy_(EXCH*EXCH,[(0.0_wp,0.0_wp)],0,TMP,1)
    call ZGEMM_('C','N',EXCH,EXCH,EXCH,(1.0_wp,0.0_wp),Z,EXCH,M1(L,:,:),EXCH,(0.0_wp,0.0_wp),TMP,EXCH)
    call ZGEMM_('N','N',EXCH,EXCH,EXCH,(1.0_wp,0.0_wp),TMP,EXCH,Z,EXCH,(0.0_wp,0.0_wp),M2(L,:,:),EXCH)
  end do !L

else

  do L=1,3
    call zcopy_(EXCH*EXCH,[(0.0_wp,0.0_wp)],0,TMP,1)
    call ZGEMM_('C','N',N,N,N,(1.0_wp,0.0_wp),Z(1:N,1:N),N,M1(L,1:N,1:N),N,(0.0_wp,0.0_wp),TMP,N)
    call ZGEMM_('N','N',N,N,N,(1.0_wp,0.0_wp),TMP,N,Z(1:N,1:N),N,(0.0_wp,0.0_wp),M2(L,1:N,1:N),N)

    call zcopy_(EXCH*EXCH,[(0.0_wp,0.0_wp)],0,TMP,1)
    call ZGEMM_('C','N',N,EXCH,N,(1.0_wp,0.0_wp),Z(1:N,1:N),N,M1(L,1:N,1:EXCH),N,(0.0_wp,0.0_wp),TMP(1:N,1:EXCH),N)

    do I=1,N
      do J=N+1,EXCH
        M2(L,I,J) = TMP(I,J)
        M2(L,J,I) = conjg(TMP(I,J))
      end do
    end do
    do i=N+1,EXCH
      do j=N+1,EXCH
        M2(L,i,j) = M1(L,i,j)
      end do
    end do
  end do !L

end if !N == exch

if (DBG) then
  write(6,'(A)') 'UTMU :: input moment'
  do i=1,EXCH
    do j=1,EXCH
      write(6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'|M_l|',j,'>',(M1(l,i,j),l=1,3)
    end do
  end do
  write(6,'(A)') 'UTMU :: unitary transformtion matrix'
  do i=1,N
    do j=1,N
      write(6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'| U |',j,'>',(Z(i,j),l=1,3)
    end do
  end do
  write(6,'(A)') 'UTMU :: output moment'
  do i=1,EXCH
    do j=1,EXCH
      write(6,'(A,i3,A,i3,A,3(2ES16.8,2x))') '<',i,'|M_l|',j,'>',(M2(l,i,j),l=1,3)
    end do
  end do
end if
call mma_deallocate(TMP)

return

end subroutine utmu
