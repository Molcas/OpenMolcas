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

integer n, dkord, xord, vord
real*8 EL(n,n), ES(n,n), OL(n,n), OS(n,n), Ep(n), E0(n)
real*8 dkcof(*), cc(*)
! Output :
!   EL      overwritten by the transformed Hamiltonian
!   wsav    store W matrices
real*8 wsav(n,n,*)
! Temp
real*8 wr(n,n), rw(n,n), t1(n,n), t2(n,n), t3(n,n), t4(n,n), or(n,n,*), ro(n,n,*), e(n,n,*), rer(n,n,*), or_(n,n,*), ro_(n,n,*), &
       e_(n,n,*), rer_(n,n,*), s1(n,n,*), s2(n,n,*)
integer i, j, k, ord, cou, ks, ioe
logical ifodd

! Copy initial matrices

do i=1,n
  do j=1,n
    e(j,i,1) = EL(j,i)
    rer(j,i,1) = ES(j,i)
    or(j,i,1) = OL(j,i)
    ro(j,i,1) = OS(j,i)
  end do
end do
! counter of total number of matrix multiplications
cou = 0
do ord=1,vord/2
  do k=1,vord
    do i=1,n
      do j=1,n
        or_(j,i,k) = 0.d0
        ro_(j,i,k) = 0.d0
        e_(j,i,k) = 0.d0
        rer_(j,i,k) = 0.d0
      end do
    end do
  end do

  ! Set up W(ord) matrix

  do i=1,n
    do j=1,n
      wr(j,i) = or(j,i,ord)/(Ep(j)+Ep(i))
      rw(j,i) = -ro(j,i,ord)/(Ep(j)+Ep(i))
      if (ord <= xord) then
        ! Save W matrix
        wsav(j,i,ord*2-1) = wr(j,i)
        wsav(j,i,ord*2) = rw(j,i)
      end if
    end do
  end do

  ! Calculate [W(ord)->O(ord)], the terms of W(ord) apply on O(ord)
  !   also include the contributions from [W(ord)->E(0)] via recalculated coefficients

  k = ord
  ifodd = .true.
  do i=1,n
    do j=1,n
      t1(j,i) = or(j,i,k)
      t2(j,i) = ro(j,i,k)
    end do
  end do
  call dkh_wspec(n,ord,vord,ifodd,dkcof,wr,rw,t1,t2,e_,rer_,or_,ro_,cou,s1,s2,t3,t4,cc)

  ! Calculate [W(ord)->E(1-...)/O(ord+1-...)]

  do ks=1,vord
    ! W1 only apply to E1
    !DP if ((ord == 1) .and. (ks >= 2)) cycle
    if (ord > 1 .or. ks < 2) then
      do ioe=1,2
        ! only even operator for k<=ord survived, odd terms were eliminated
        !DP if ((ioe == 1) .and. (ks <= ord)) cycle
        if (ioe == 2 .or. ks > ord) then
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
          call dkh_wgene(n,ord,k,vord,ifodd,dkcof,wr,rw,t1,t2,e_,rer_,or_,ro_,cou,s1,s2,t3,t4)
        end if
        ! cycle for odd/even operator
      end do
    end if
    ! cycle for ks
  end do
  ! copy
  do k=1,vord
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

! Sum over all k<=dkord terms

do i=1,n
  do j=1,n
    if (j == i) then
      EL(j,i) = E0(i)
    else
      EL(j,i) = 0.d0
    end if
  end do
end do
do k=1,dkord
  do i=1,n
    do j=1,n
      EL(j,i) = EL(j,i)+e(j,i,k)
    end do
  end do
end do
!write(6,*) 'DKH',vord,' Total matmul',cou

return

end subroutine dkh_ham
