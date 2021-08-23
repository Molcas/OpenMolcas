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

subroutine dkh_ham(n,dkord,xord,vord,EL,ES,OL,OS,Ep,E0,dkcof,cc,wr,rw,t1,t2,t3,t4,or,ro,e,rer,or_,ro_,e_,rer_,s1,s2,wsav)
! Evaluate DKH Hamiltonian in moment space

#include "intent.fh"

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
! Input :
!   n       dimension of matrix
!   dkord   order of DKH Hamiltonian
!   xord    order for property integrals
!   vord    actual calculated order satisfy both dkord and xord
!   ( EL OL )
!   ( OS ES ) potential matrix in fpFW space
!   Ep,E0   diagonal kinetic matrix, E0=Ep-c^{2}
!   dkcof   expansion coefficient of unitary transformation in terms of anti-Hermitian matrix W
! Output :
!   EL      overwritten by the transformed Hamiltonian
!   wsav    store W matrices
integer(kind=iwp), intent(in) :: n, dkord, xord, vord
real(kind=wp), intent(inout) :: EL(n,n)
real(kind=wp), intent(in) :: ES(n,n), OL(n,n), OS(n,n), Ep(n), E0(n), dkcof(*)
real(kind=wp), intent(_OUT_) :: cc(*), or(n,n,*), ro(n,n,*), e(n,n,*), rer(n,n,*), or_(n,n,*), ro_(n,n,*), e_(n,n,*), rer_(n,n,*), &
                                s1(n,n,*), s2(n,n,*), wsav(n,n,*)
real(kind=wp), intent(out) :: wr(n,n), rw(n,n), t1(n,n), t2(n,n), t3(n,n), t4(n,n)
integer(kind=iwp) :: i, j, k, ord, cou, ks, ioe
logical(kind=iwp) :: ifodd

! Copy initial matrices

e(:,:,1) = EL(:,:)
rer(:,:,1) = ES(:,:)
or(:,:,1) = OL(:,:)
ro(:,:,1) = OS(:,:)
! counter of total number of matrix multiplications
cou = 0
do ord=1,vord/2
  or_(:,:,1:vord) = Zero
  ro_(:,:,1:vord) = Zero
  e_(:,:,1:vord) = Zero
  rer_(:,:,1:vord) = Zero

  ! Set up W(ord) matrix

  do i=1,n
    do j=1,n
      wr(j,i) = or(j,i,ord)/(Ep(j)+Ep(i))
      rw(j,i) = -ro(j,i,ord)/(Ep(j)+Ep(i))
    end do
  end do
  if (ord <= xord) then
    ! Save W matrix
    wsav(:,:,ord*2-1) = wr(:,:)
    wsav(:,:,ord*2) = rw(:,:)
  end if

  ! Calculate [W(ord)->O(ord)], the terms of W(ord) apply on O(ord)
  !   also include the contributions from [W(ord)->E(0)] via recalculated coefficients

  k = ord
  ifodd = .true.
  t1(:,:) = or(:,:,k)
  t2(:,:) = ro(:,:,k)
  call dkh_wspec(n,ord,vord,ifodd,dkcof,wr,rw,t1,t2,e_,rer_,or_,ro_,cou,s1,s2,t3,t4,cc)

  ! Calculate [W(ord)->E(1-...)/O(ord+1-...)]

  do ks=1,vord
    ! W1 only apply to E1
    !DP if ((ord == 1) .and. (ks >= 2)) cycle
    if ((ord > 1) .or. (ks < 2)) then
      do ioe=1,2
        ! only even operator for k<=ord survived, odd terms were eliminated
        !DP if ((ioe == 1) .and. (ks <= ord)) cycle
        if ((ioe == 2) .or. (ks > ord)) then
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
          call dkh_wgene(n,ord,k,vord,ifodd,dkcof,wr,rw,t1,t2,e_,rer_,or_,ro_,cou,s1,s2,t3,t4)
        end if
        ! cycle for odd/even operator
      end do
    end if
    ! cycle for ks
  end do
  ! copy
  or(:,:,1:vord) = or_(:,:,1:vord)
  ro(:,:,1:vord) = ro_(:,:,1:vord)
  e(:,:,1:vord) = e_(:,:,1:vord)
  rer(:,:,1:vord) = rer_(:,:,1:vord)
! cycle for ord
end do

! Sum over all k<=dkord terms

EL(:,:) = zERO
do i=1,n
  EL(i,i) = E0(i)
end do
do k=1,dkord
  EL(:,:) = EL(:,:)+e(:,:,k)
end do
!write(u6,*) 'DKH',vord,' Total matmul',cou

return

end subroutine dkh_ham
