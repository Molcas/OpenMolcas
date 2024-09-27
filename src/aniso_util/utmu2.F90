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

subroutine utmu2(exch,n,z,m)
! the same as utmu, except being that the input m is being transformed.

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: exch, n
complex(kind=wp), intent(in) :: z(n,n)
complex(kind=wp), intent(inout) :: m(3,exch,exch)
integer(kind=iwp) :: i, i1, j, j1, l
real(kind=wp) :: r1, r2
complex(kind=wp), allocatable :: mtmp(:,:,:), tmp(:,:), tmp2(:,:)
real(kind=wp), external :: dznrm2_

if ((n <= 0) .or. (exch <= 0)) then
  write(u6,'(a)') 'in utmu2:   exch or n<=0 !!!'
  write(u6,*) 'exch=',exch
  write(u6,*) 'n   =',n
  call xFlush(u6)
  call abend()
end if

if (n > exch) then
  write(u6,'(a)') 'in utmu2:   exch < n !!!'
  write(u6,*) 'exch=',exch
  write(u6,*) 'n   =',n
  write(u6,'(a)') 'nothing is to be done >> return'
  call xFlush(u6)
  call abend()
end if

r1 = dznrm2_(3*exch*exch,m,1)
r2 = dznrm2_(n*n,z,1)
if ((r1 < 1.0e-25_wp) .or. (r2 < 1.0e-25_wp)) then
  write(u6,'(a)') 'in utmu2:   m or z are empty!!!'
  write(u6,*) 'norm(m)=',r1
  write(u6,*) 'norm(z)=',r2
  return
end if

#ifdef _DEBUGPRINT_
write(u6,'(a)') 'utmu2 :: input moment'
do i=1,exch
  do j=1,exch
    write(u6,'(a,i3,a,i3,a,3(2es16.8,2x))') '<',i,'|m_l|',j,'>',(m(l,i,j),l=1,3)
  end do
end do
#endif

call mma_allocate(tmp,n,exch,'tmp')
call mma_allocate(tmp2,n,exch,'tmp2')
if (n /= exch) then
  call mma_allocate(mtmp,3,(exch-n),(exch-n),'mtmp')
  ! save the part which is not altered
  do i=n+1,exch
    do j=n+1,exch
      i1 = i-n
      j1 = j-n
      mtmp(:,i1,j1) = m(:,i,j)
    end do
  end do
end if

if (n == exch) then

  do l=1,3
    tmp2(:,:) = m(l,:,:)
    call zgemm_('c','n',exch,exch,exch,cOne,z,exch,tmp2,exch,cZero,tmp,exch)
    call zgemm_('n','n',exch,exch,exch,cOne,tmp,exch,z,exch,cZero,tmp2,exch)
    m(l,:,:) = tmp2(:,:)
  end do !l

else

  do l=1,3
    tmp2(:,:) = m(l,1:n,:)
    call zgemm_('c','n',n,n,n,cOne,z,n,tmp2(:,1:n),n,cZero,tmp(:,1:n),n)
    call zgemm_('n','n',n,n,n,cOne,tmp(:,1:n),n,z,n,cZero,tmp2(:,1:n),n)
    m(l,1:n,1:n) = tmp2(:,1:n)
    call zgemm_('c','n',n,exch,n,cOne,z,n,tmp2,n,cZero,tmp,n)

    do i=1,n
      m(l,i,n+1:) = tmp(i,n+1:)
      m(l,n+1:,i) = conjg(tmp(i,n+1:))
    end do
    m(l,n+1:,n+1:) = mtmp(l,1:exch-n,1:exch-n)
  end do !l
end if !n == exch

#ifdef _DEBUGPRINT_
write(u6,'(a)') 'utmu2 :: unitary transformtion matrix'
do i=1,n
  do j=1,n
    write(u6,'(a,i3,a,i3,a,3(2es16.8,2x))') '<',i,'| u |',j,'>',z(i,j)
  end do
end do
write(u6,'(a)') 'utmu2 :: output moment'
do i=1,exch
  do j=1,exch
    write(u6,'(a,i3,a,i3,a,3(2es16.8,2x))') '<',i,'|m_l|',j,'>',(m(l,i,j),l=1,3)
  end do
end do
#endif

if (n < exch) call mma_deallocate(mtmp)
call mma_deallocate(tmp)
call mma_deallocate(tmp2)

#ifdef _DEBUGPRINT_
write(u6,*) 'at the end of utmu2'
call prmom('utmu2, moment',m,n)
#endif

return

end subroutine utmu2
