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

#include "intent.fh"

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: n, ord, ndk
logical(kind=iwp), intent(inout) :: ifodd
real(kind=wp), intent(in) :: cdk(ndk), wr(n,n), rw(n,n)
real(kind=wp), intent(inout) :: t1(n,n), t2(n,n), e(n,n,ndk), rer(n,n,ndk), or(n,n,ndk), ro(n,n,ndk)
integer(kind=iwp), intent(inout) :: info
real(kind=wp), intent(_OUT_) :: s1(n,n,*), s2(n,n,*), c(*)
real(kind=wp), intent(out) :: t3(n,n), t4(n,n)
integer(kind=iwp) :: m, i, j, k

m = ndk/ord
s1(:,:,1) = t1(:,:)
s2(:,:,1) = t2(:,:)
do i=1,m-1
  t1(:,:) = Zero
  t2(:,:) = Zero
  call dkh_cofu_spec(ndk,cdk,i+1,c)
  k = (i+1)*ord
  call dkh_woprig(n,ifodd,wr,rw,s1(:,:,i),s2(:,:,i),s1(:,:,i+1),s2(:,:,i+1),t3,t4)
  info = info+2
  t1(:,:) = t1(:,:)+s1(:,:,i+1)*c(i+1)
  t2(:,:) = t2(:,:)+s2(:,:,i+1)*c(i+1)
  do j=1,i
    call dkh_woplft(n,ifodd,wr,rw,s1(:,:,j),s2(:,:,j),s1(:,:,j),s2(:,:,j),t3,t4)
    info = info+2
    t1(:,:) = t1(:,:)+s1(:,:,j)*c(j)
    t2(:,:) = t2(:,:)+s2(:,:,j)*c(j)
  end do
  ifodd = .not. ifodd
  if (ifodd) then
    or(:,:,k) = or(:,:,k)+t1(:,:)
    ro(:,:,k) = ro(:,:,k)+t2(:,:)
  else
    e(:,:,k) = e(:,:,k)+t1(:,:)
    rer(:,:,k) = rer(:,:,k)+t2(:,:)
  end if
end do

end subroutine dkh_wspec
