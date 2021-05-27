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

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
! Output : EL overwritten by the transformed property matrix
integer(kind=iwp), intent(in) :: n, dkord, xord, vord
real(kind=wp), intent(inout) :: EL(n,n)
real(kind=wp), intent(in) :: ES(n,n), OL(n,n), OS(n,n), Ep(n), E0(n), dkcof(*), cc(*), wsav(n,n,*)
real(kind=wp), intent(out) :: wr(n,n), rw(n,n), t1(n,n), t2(n,n), t3(n,n), t4(n,n), or(n,n,*), ro(n,n,*), e(n,n,*), rer(n,n,*), &
                              or_(n,n,*), ro_(n,n,*), e_(n,n,*), rer_(n,n,*), s1(n,n,*), s2(n,n,*)
integer(kind=iwp) :: i, j, k, ord, cou, ks, ioe
logical(kind=iwp) :: ifodd

! Copy initial matrices

do i=1,n
  do j=1,n
    e(j,i,1) = EL(j,i)
    rer(j,i,1) = ES(j,i)
    or(j,i,1) = OL(j,i)
    ro(j,i,1) = OS(j,i)
  end do
end do

cou = 0
do ord=1,xord
  do k=1,vord
    do i=1,n
      do j=1,n
        or_(j,i,k) = Zero
        ro_(j,i,k) = Zero
        e_(j,i,k) = Zero
        rer_(j,i,k) = Zero
      end do
    end do
  end do

  ! Set up W(ord) matrix, copy from saved set generated from dkh_ham

  do i=1,n
    do j=1,n
      wr(j,i) = wsav(j,i,ord*2-1)
      rw(j,i) = wsav(j,i,ord*2)
    end do
  end do

  ! Calculate [W(ord)->E(1-...)/O(1-...)]
  !   note that no odd term will be eliminated

  do ks=1,xord+1
    ! W1 only apply to E1/O1
    !DP if ((ord == 1) .and. (ks >= 2)) cycle
    if (ord > 1 .or. ks < 2) then
      do ioe=1,2
        k = ks
        if (ioe == 1) then
          ifodd = .true.
        else
          ifodd = .false.
        end if
        ! copy initial matrix
        if (ifodd) then
          do i=1,n
            do j=1,n
              t1(j,i) = or(j,i,k)
              t2(j,i) = ro(j,i,k)
              or_(j,i,k) = or_(j,i,k)+t1(j,i)
              ro_(j,i,k) = ro_(j,i,k)+t2(j,i)
            end do
          end do
        else
          do i=1,n
            do j=1,n
              t1(j,i) = e(j,i,k)
              t2(j,i) = rer(j,i,k)
              e_(j,i,k) = e_(j,i,k)+t1(j,i)
              rer_(j,i,k) = rer_(j,i,k)+t2(j,i)
            end do
          end do
        end if
        call dkh_wgene(n,ord,k,xord+1,ifodd,dkcof,wr,rw,t1,t2,e_,rer_,or_,ro_,cou,s1,s2,t3,t4)
        ! cycle for odd/even operator
      end do
    end if
    ! cycle for ks
  end do
  ! copy
  do k=1,xord+1
    do i=1,n
      do j=1,n
        or(j,i,k) = or_(j,i,k)
        ro(j,i,k) = ro_(j,i,k)
        e(j,i,k) = e_(j,i,k)
        rer(j,i,k) = rer_(j,i,k)
      end do
    end do
  end do
! cycle for ord
end do

! Sum over all k<=xord+1 terms
!   +1 appears because X itself is not counted as an order

do k=2,xord+1
  do i=1,n
    do j=1,n
      EL(j,i) = EL(j,i)+e(j,i,k)
    end do
  end do
end do
!write(u6,*) 'DKHX',xord,' Total matmul',cou

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_integer(dkord)
  call Unused_real_array(Ep)
  call Unused_real_array(E0)
  call Unused_real_array(cc)
end if

end subroutine dkh_xpx
