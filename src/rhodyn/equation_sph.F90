!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2022-2023, Vladislav Kochetov                          *
!               2023, Thies Romig                                      *
!***********************************************************************

subroutine equation_sph(time,rhot,res)
!***********************************************************************
! Purpose : RHS of Liouville equation is obtained here in ITOs basis
!
!***********************************************************************
!
!  time   : current time
!  rhot   : density matrix at current time
!  res    : obtained RHS of Liouville equation d(rhot)/d(time)

use rhodyn_data, only: d, hamiltonian, hamiltoniant, ipglob, flag_pulse, k_max, len_sph, k_ranks, q_proj, &
                       Y1, Y2, irs1, irs2, ics1, ics2, q_max, mirr, lroots, threshold, n
use rhodyn_utils, only: get_kq_order, compare_matrices
use linalg_mod, only: mult
use Constants, only: Zero, One, cOne, cZero, Onei
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: time
complex(kind=wp), intent(in) :: rhot(len_sph,d,d)
complex(kind=wp), intent(out) :: res(len_sph,d,d)
!complex(kind=wp) :: z(d,d)
integer(kind=iwp) :: c, i, k, K_prime, l, l_prime, m, q, n1, n2
!integer(kind=iwp), save :: icount = 0

if (ipglob > 3) write(u6,*) 'solve equation, time: ',time

! number of states with ground state spin to distinguish contributions to different spin manifolds
! now works only with 2 spin manifolds
!ngs = lroots(1)

! if pulse is enabled, modify Hamiltonian at time t:
if (flag_pulse) call pulse(hamiltonian,hamiltoniant,time,-1)

i = 1
do l=1,len_sph
  k = k_ranks(l)
  q = q_proj(l)
  if ((q>=-q_max) .and. (q<=0)) then
  ! get right part of Liouville equation -i*(hamiltoniant*rhot - rhot*hamiltoniant)
  call zgemm_('N','N',d,d,d,-Onei,hamiltoniant,d,rhot(l,:,:),d,cZero,res(l,:,:),d)
  call zgemm_('N','N',d,d,d,Onei,rhot(l,:,:),d,hamiltoniant,d,cOne,res(l,:,:),d)
  ! Vsoc contribution
  do m=1,3
    do K_prime=k-1,k+1,1
      if ((K_prime < 0) .or. (K_prime > k_max)) cycle
      if ((q-m+2 > K_prime) .or. (q-m+2 < -K_prime)) cycle
      !if (icount == 0) write(u6,*) i
      l_prime = get_kq_order(K_prime,q-m+2)
      ! masked contributions, working
      !call zgemm_('N','N',d,d,d,-Onei,Y1(:,:,i),d,rhot(l_prime,:,:)*irs1,d,cOne,res(l,:,:),d)
      !call zgemm_('N','N',d,d,d,-Onei,Y1(:,:,i+1),d,rhot(l_prime,:,:)*irs2,d,cOne,res(l,:,:),d)
      !call zgemm_('N','N',d,d,d,-Onei,rhot(l_prime,:,:)*ics1,d,Y2(:,:,i),d,cOne,res(l,:,:),d)
      !call zgemm_('N','N',d,d,d,-Onei,rhot(l_prime,:,:)*ics2,d,Y2(:,:,i+1),d,cOne,res(l,:,:),d)
      ! sliced contributions, working
      !call zgemm_('N','N',d,lroots(1),d,-Onei,Y1(:,:,i),d,  rhot(l_prime,:,1:lroots(1)),d,cOne, res(l,:,1:lroots(1)),d)
      !call zgemm_('N','N',d,lroots(2),d,-Onei,Y1(:,:,i+1),d,rhot(l_prime,:,lroots(1)+1:), d,cOne, res(l,:,lroots(1)+1:), d)
      !call zgemm_('N','N',lroots(1),d,d,-Onei,rhot(l_prime,1:lroots(1),:),&
      !            lroots(1),Y2(:,:,i),  d,cOne,res(l,1:lroots(1),:),lroots(1))
      !call zgemm_('N','N',lroots(2),d,d,-Onei,rhot(l_prime,lroots(1)+1:,:),&
      !            lroots(2),Y2(:,:,i+1),d,cOne,res(l,lroots(1)+1:,:), lroots(2))
      !i = i+2
      do c=1,n
        ! set boundaries for slice multiplication
        if (c == 1) then
          n1 = 1
          n2 = lroots(1)
        else
          n1 = n1 + lroots(c-1)
          n2 = n2 + lroots(c)
        end if
        call zgemm_('N','N',d,lroots(c),d,-Onei,Y1(:,:,i),d,rhot(l_prime,:,n1:n2),d,cOne,res(l,:,n1:n2),d)
        call zgemm_('N','N',lroots(c),d,d,-Onei,rhot(l_prime,n1:n2,:),lroots(c),Y2(:,:,i),d,cOne,res(l,n1:n2,:),lroots(c))
        i = i+1
      end do
    end do
  end do
  else if (q == 1) then
      !write(u6,*) 'rank = ', k, q
      !write(u6,*) 'value = ', sum(res(l,:,:))
      !write(u6,*) 'error = ', sum(res(l,:,:)-transpose(conjg(res(l-2,:,:)))*mirr)
      !call compare_matrices(res(l,:,:),transpose(conjg(res(l-2,:,:)))*mirr,d,'Check mirror',threshold)
      res(l,:,:) = transpose(conjg(res(l-2,:,:)))*mirr
  end if
end do
!averaging
!do l=1, k_max
! call zgemm_('N','N',d,d,d,-Onei,hamiltoniant,d,rhot(l,:,:),d,cZero,res(l,:,:),d)
!  call zgemm_('N','N',d,d,d,Onei,rhot(l,:,:),d,hamiltoniant,d,cOne,res(l,:,:),d)
  ! Vsoc contribution
!  do m=1,3
!    do K_prime=k-1,k+1,1
!      if ((K_prime < 0) .or. (K_prime > k_max)) cycle
!      !if (icount == 0) write(u6,*) i
!      l_prime = get_kq_order(K_prime,q-m+2)
!      call zgemm_('N','N',d,lroots(1),d,-Onei,Y1(:,:,i),d,  rhot(l_prime,:,1:lroots(1)),d,cOne, res(l,:,1:lroots(1)),d)
!      call zgemm_('N','N',d,lroots(2),d,-Onei,Y1(:,:,i+1),d,rhot(l_prime,:,lroots(1)+1:), d,cOne, res(l,:,lroots(1)+1:), d)
!      call zgemm_('N','N',lroots(1),d,d,-Onei,rhot(l_prime,1:lroots(1),:),&
!                  lroots(1),Y2(:,:,i),  d,cOne,res(l,1:lroots(1),:),lroots(1))
!      call zgemm_('N','N',lroots(2),d,d,-Onei,rhot(l_prime,lroots(1)+1:,:),&
!                  lroots(2),Y2(:,:,i+1),d,cOne,res(l,lroots(1)+1:,:), lroots(2))
!      i = i+2
!    end do
!  end do

! ELEMENT-WISE PROPAGATION
! use rhodyn_data, only: d, E_SF, hamiltonian, ipglob, flag_pulse, flag_so, k_max, list_sf_spin, len_sph, k_ranks, q_proj, &
!                        threshold, V_SO_red
! use rhodyn_utils, only: W3J, W6J, get_kq_order
! use Constants, only: Zero, One, cZero, Onei
! use Definitions, only: wp, iwp, u6

! implicit none
! real(kind=wp), intent(in) :: time
! complex(kind=wp), intent(in) :: rhot(len_sph,d,d)
! complex(kind=wp), intent(out) :: res(len_sph,d,d)
! complex(kind=wp) :: z
! integer(kind=iwp) :: a, b, c, k, K_prime, l, l_prime, m, q, sa, sb, sc

! if (ipglob > 3) write(u6,*) 'solve equation, time: ',time

! res(:,:,:) = Zero

! ! a loop over rows
! do a=1,d
!   sa = nint(2*list_sf_spin(a))
!   ! b loop over columns
!   do b=1,d
!     sb = nint(2*list_sf_spin(b))
!     ! l loop over indices k,q of ITOs basis
!     do l=1,len_sph
!       k = k_ranks(l)
!       if (k > (sa+sb)/2) cycle ! delta {s_a k s_b}
!       !if (k > nint(sa+sb)) cycle ! delta {s_a k s_b}
!       q = q_proj(l)
!       z = cZero ! accumulates contribution from V_SO term
!       ! c loop over rows and columns for calculating V_SO contribution
!       do c=1,d
!         sc = nint(2*list_sf_spin(c))
!         ! dipole contribution
!         if (flag_pulse) then
!           if (sa == sc) z = z-rhot(l,c,b)*hamiltonian(a,c)
!           if (sb == sc) z = z+rhot(l,a,c)*hamiltonian(c,b)
!         end if
!         ! Vsoc contribution
!         if (flag_so) then
!           do m=1,3
!             if ((abs(V_SO_red(a,c,m)) <= threshold) .and. (abs(V_SO_red(c,b,m)) <= threshold)) cycle
!             do K_prime=k-1,k+1,1
!               if ((K_prime < 0) .or. (K_prime > k_max)) cycle
!               !do Q_prime=-K_prime,K_prime,1
!               !  l_prime = get_kq_order(K_prime,Q_prime)
!               l_prime = get_kq_order(K_prime,q-m+2)
!               z = z+(-1)**((sa+sb)/2-q+m)* &
!                   sqrt(real(3*(2*k+1)*(2*K_prime+1),kind=wp))* &
!                   W3J(real(K_prime,kind=wp),One,real(k,kind=wp),real(q-m+2,kind=wp),real(m-2,kind=wp),real(-q,kind=wp))* &
!                   (W6J(2*K_prime,2,2*k,sa,sb,sc)*rhot(l_prime,c,b)*V_SO_red(a,c,m)- &
!                    (-1)**(K_prime+k+1)*W6J(2,2*K_prime,2*k,sa,sb,sc)*rhot(l_prime,a,c)*V_SO_red(c,b,m))
!               !end do
!             end do
!           end do
!         end if
!       end do
!       res(l,a,b) = -Onei*((E_SF(a)-E_SF(b))*rhot(l,a,b)+z)
!     end do
!   end do
! end do

end subroutine equation_sph
