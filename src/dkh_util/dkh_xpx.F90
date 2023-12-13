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

subroutine dkh_xpx(n,dkord,xord,vord,EL,ES,OL,OS,Ep,E0,dkcof,cc,wr,rw,t1,t2,t3,t4,or,ro,e,rer,or_,ro_,e_,rer_,s1,s2,wsav)
! Evaluate DKH transformation in moment space for property operator

#include "intent.fh"

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
! Output : EL overwritten by the transformed property matrix
integer(kind=iwp), intent(in) :: n, dkord, xord, vord
real(kind=wp), intent(inout) :: EL(n,n)
real(kind=wp), intent(in) :: ES(n,n), OL(n,n), OS(n,n), Ep(n), E0(n), dkcof(*), cc(*), wsav(n,n,xord*2)
real(kind=wp), intent(out) :: wr(n,n), rw(n,n), t1(n,n), t2(n,n), t3(n,n), t4(n,n)
real(kind=wp), intent(_OUT_) :: or(n,n,*), ro(n,n,*), e(n,n,*), rer(n,n,*), or_(n,n,*), ro_(n,n,*), e_(n,n,*), rer_(n,n,*), &
                                s1(n,n,*), s2(n,n,*)
integer(kind=iwp) :: k, ord, cou, ks, ioe
logical(kind=iwp) :: ifodd

#include "macros.fh"
unused_var(dkord)
unused_var(Ep)
unused_var(E0)
unused_var(cc(1))

! Copy initial matrices

e(:,:,1) = EL(:,:)
rer(:,:,1) = ES(:,:)
or(:,:,1) = OL(:,:)
ro(:,:,1) = OS(:,:)

cou = 0
do ord=1,xord
  or_(:,:,1:vord) = Zero
  ro_(:,:,1:vord) = Zero
  e_(:,:,1:vord) = Zero
  rer_(:,:,1:vord) = Zero

  ! Set up W(ord) matrix, copy from saved set generated from dkh_ham

  wr(:,:) = wsav(:,:,ord*2-1)
  rw(:,:) = wsav(:,:,ord*2)

  ! Calculate [W(ord)->E(1-...)/O(1-...)]
  !   note that no odd term will be eliminated

  do ks=1,xord+1
    ! W1 only apply to E1/O1
    !DP if ((ord == 1) .and. (ks >= 2)) cycle
    if ((ord > 1) .or. (ks < 2)) then
      do ioe=1,2
        k = ks
        if (ioe == 1) then
          ifodd = .true.
        else
          ifodd = .false.
        end if
        ! copy initial matrix
        if (ifodd) then
          t1(:,:) = or(:,:,k)
          t2(:,:) = ro(:,:,k)
          or_(:,:,k) = or_(:,:,k)+t1(:,:)
          ro_(:,:,k) = ro_(:,:,k)+t2(:,:)
        else
          t1(:,:) = e(:,:,k)
          t2(:,:) = rer(:,:,k)
          e_(:,:,k) = e_(:,:,k)+t1(:,:)
          rer_(:,:,k) = rer_(:,:,k)+t2(:,:)
        end if
        call dkh_wgene(n,ord,k,xord+1,ifodd,dkcof,wr,rw,t1,t2,e_,rer_,or_,ro_,cou,s1,s2,t3,t4)
        ! cycle for odd/even operator
      end do
    end if
    ! cycle for ks
  end do
  ! copy
  or(:,:,1:xord+1) = or_(:,:,1:xord+1)
  ro(:,:,1:xord+1) = ro_(:,:,1:xord+1)
  e(:,:,1:xord+1) = e_(:,:,1:xord+1)
  rer(:,:,1:xord+1) = rer_(:,:,1:xord+1)
! cycle for ord
end do

! Sum over all k<=xord+1 terms
!   +1 appears because X itself is not counted as an order

do k=2,xord+1
  EL(:,:) = EL(:,:)+e(:,:,k)
end do
!write(u6,*) 'DKHX',xord,' Total matmul',cou

return

end subroutine dkh_xpx
