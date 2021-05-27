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

subroutine dkh_wspec(n,ord,ndk,ifodd,cdk,wr,rw,t1,t2,e,rer,or,ro,info,s1,s2,t3,t4,c)
! Calculate U(Word)O_{ord}U^{\dag}(Word) + U(Word)E_{0}U^{\dag}(Word)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, ord, ndk
logical(kind=iwp), intent(inout) :: ifodd
real(kind=wp), intent(in) :: cdk(ndk), wr(n,n), rw(n,n)
real(kind=wp), intent(inout) :: t1(n,n), t2(n,n), e(n,n,ndk), rer(n,n,ndk), or(n,n,ndk), ro(n,n,ndk)
integer(kind=iwp), intent(inout) :: info
real(kind=wp), intent(out) :: s1(n,n,*), s2(n,n,*), t3(n,n), t4(n,n), c(*)
integer(kind=iwp) :: m, i, j, k, L1, L2

m = ndk/ord
do L1=1,n
  do L2=1,n
    s1(L2,L1,1) = t1(L2,L1)
    s2(L2,L1,1) = t2(L2,L1)
  end do
end do
do i=1,m-1
  do L1=1,n
    do L2=1,n
      t1(L2,L1) = Zero
      t2(L2,L1) = Zero
    end do
  end do
  call dkh_cofu_spec(ndk,cdk,i+1,c)
  k = (i+1)*ord
  call dkh_woprig(n,ifodd,ord,i*ord,wr,rw,s1(1,1,i),s2(1,1,i),s1(1,1,i+1),s2(1,1,i+1),t3,t4)
  info = info+2
  do L1=1,n
    do L2=1,n
      t1(L2,L1) = t1(L2,L1)+s1(L2,L1,i+1)*c(i+1)
      t2(L2,L1) = t2(L2,L1)+s2(L2,L1,i+1)*c(i+1)
    end do
  end do
  do j=1,i
    call dkh_woplft(n,ifodd,ord,i*ord,wr,rw,s1(1,1,j),s2(1,1,j),s1(1,1,j),s2(1,1,j),t3,t4)
    info = info+2
    do L1=1,n
      do L2=1,n
        t1(L2,L1) = t1(L2,L1)+s1(L2,L1,j)*c(j)
        t2(L2,L1) = t2(L2,L1)+s2(L2,L1,j)*c(j)
      end do
    end do
  end do
  ifodd = .not. ifodd
  if (ifodd) then
    do L1=1,n
      do L2=1,n
        or(L2,L1,k) = or(L2,L1,k)+t1(L2,L1)
        ro(L2,L1,k) = ro(L2,L1,k)+t2(L2,L1)
      end do
    end do
  else
    do L1=1,n
      do L2=1,n
        e(L2,L1,k) = e(L2,L1,k)+t1(L2,L1)
        rer(L2,L1,k) = rer(L2,L1,k)+t2(L2,L1)
      end do
    end do
  end if
end do

end subroutine dkh_wspec
