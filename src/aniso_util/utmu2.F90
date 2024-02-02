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

implicit none
integer, parameter :: wp = kind(0.d0)
#include "stdalloc.fh"
integer, intent(in) :: exch, n
complex(kind=8), intent(inout) :: m(3,exch,exch)
complex(kind=8), intent(in) :: z(n,n)
! local variables:
integer :: l, i, j, i1, j1
logical :: dbg
real(kind=8) :: dznrm2_, r1, r2
external :: dznrm2_
complex(kind=8), allocatable :: tmp(:,:), mtmp(:,:,:)

dbg = .false.

if ((n <= 0) .or. (exch <= 0)) then
  write(6,'(a)') 'in utmu2:   exch or n<=0 !!!'
  write(6,*) 'exch=',exch
  write(6,*) 'n   =',n
  call xflush(6)
  call xquit(128)
end if

if (n > exch) then
  write(6,'(a)') 'in utmu2:   exch < n !!!'
  write(6,*) 'exch=',exch
  write(6,*) 'n   =',n
  write(6,'(a)') 'nothing is to be done >> return'
  call xflush(6)
  call xquit(128)
end if

r1 = dznrm2_(3*exch*exch,m,1)
r2 = dznrm2_(n*n,z,1)
if ((r1 < 1.0e-25_wp) .or. (r2 < 1.0e-25_wp)) then
  write(6,'(a)') 'in utmu2:   m or z are empty!!!'
  write(6,*) 'norm(m)=',r1
  write(6,*) 'norm(z)=',r2
  return
end if

if (dbg) then
  write(6,'(a)') 'utmu2 :: input moment'
  do i=1,exch
    do j=1,exch
      write(6,'(a,i3,a,i3,a,3(2es16.8,2x))') '<',i,'|m_l|',j,'>',(m(l,i,j),l=1,3)
    end do
  end do
end if

call mma_allocate(tmp,exch,exch,'tmp')
if (n < exch) then
  call mma_allocate(mtmp,3,(exch-n),(exch-n),'mtmp')
  ! save the part which is not altered
  do i=n+1,exch
    do j=n+1,exch
      i1 = i-n
      j1 = j-n
      do l=1,3
        mtmp(l,i1,j1) = m(l,i,j)
      end do
    end do
  end do
end if

if (n == exch) then

  do l=1,3
    call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,tmp,1)
    call zgemm_('c','n',exch,exch,exch,(1.0_wp,0.0_wp),z,exch,m(l,:,:),exch,(0.0_wp,0.0_wp),tmp,exch)
    call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,m(l,:,:),1)
    call zgemm_('n','n',exch,exch,exch,(1.0_wp,0.0_wp),tmp,exch,z,exch,(0.0_wp,0.0_wp),m(l,:,:),exch)
  end do !l

else

  do l=1,3
    call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,tmp,1)
    call zgemm_('c','n',n,n,n,(1.0_wp,0.0_wp),z(1:n,1:n),n,m(l,1:n,1:n),n,(0.0_wp,0.0_wp),tmp,n)
    call zcopy_(n*n,[(0.0_wp,0.0_wp)],0,m(l,1:n,1:n),1)

    call zgemm_('n','n',n,n,n,(1.0_wp,0.0_wp),tmp,n,z(1:n,1:n),n,(0.0_wp,0.0_wp),m(l,1:n,1:n),n)

    call zcopy_(exch*exch,[(0.0_wp,0.0_wp)],0,tmp,1)
    call zgemm_('c','n',n,exch,n,(1.0_wp,0.0_wp),z(1:n,1:n),n,m(l,1:n,1:exch),n,(0.0_wp,0.0_wp),tmp(1:n,1:exch),n)

    do i=1,n
      do j=n+1,exch
        m(l,i,j) = tmp(i,j)
        m(l,j,i) = conjg(tmp(i,j))
      end do
    end do
    do i=n+1,exch
      do j=n+1,exch
        i1 = i-n
        j1 = j-n
        m(l,i,j) = mtmp(l,i1,j1)
      end do
    end do
  end do !l
end if !n == exch

if (dbg) then
  write(6,'(a)') 'utmu2 :: unitary transformtion matrix'
  do i=1,n
    do j=1,n
      write(6,'(a,i3,a,i3,a,3(2es16.8,2x))') '<',i,'| u |',j,'>',z(i,j)
    end do
  end do
  write(6,'(a)') 'utmu2 :: output moment'
  do i=1,exch
    do j=1,exch
      write(6,'(a,i3,a,i3,a,3(2es16.8,2x))') '<',i,'|m_l|',j,'>',(m(l,i,j),l=1,3)
    end do
  end do
end if

if (n < exch) call mma_deallocate(mtmp)
call mma_deallocate(tmp)

if (dbg) then
  write(6,*) 'at the end of utmu2'
  call prmom('utmu2, moment',m,n)
end if

return

end subroutine utmu2
