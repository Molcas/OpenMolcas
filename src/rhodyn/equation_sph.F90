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

use rhodyn_data, only: d, E_SF, threshold, hamiltonian, flag_pulse, flag_so, ipglob, V_SO_red, list_sf_spin, len_sph, k_max, &
                       k_ranks, q_proj
use rhodyn_utils, only : W3J, W6J, get_kq_order
use Constants, only: Zero, One, cZero, cOne, Onei
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: time
complex(kind=wp), intent(in) :: rhot(len_sph,d,d)
complex(kind=wp), intent(out) :: res(len_sph,d,d)
complex(kind=wp) :: z
!real(kind=wp) :: sa, sb, sc, test_sixj
integer(kind=iwp) :: sa, sb, sc
integer(kind=iwp) :: a, b, c, l, k, q, l_prime, K_prime, m

if (ipglob > 2) then
  write(u6,*) 'solve equation, time: ', time
end if

res(:,:,:) = Zero

! a loop over rows
do a=1,d
  sa = nint(2*list_sf_spin(a))
  !sa = list_sf_spin(a)
  ! b loop over columns
  do b=1,d
    sb = nint(2*list_sf_spin(b))
    !sb = list_sf_spin(b)
    ! l loop over indices k,q of ITOs basis
    do l=1,len_sph
      k = k_ranks(l)
      if (k>(sa+sb)/2) cycle ! delta function {s_a k s_b}
      !if (k>nint(sa+sb)) cycle ! delta function {s_a k s_b}
      q = q_proj(l)
      z = cZero ! accumulates contribution from V_SO term
      !c loop over rows and columns for calculating V_SO contribution
      do c=1,d
        sc = nint(2*list_sf_spin(c))
        !sc = list_sf_spin(c)
        !dipole contribution
        if (flag_pulse) then
          if (sa == sc) z = z - rhot(l,c,b)*hamiltonian(a,c)
          if (sb == sc) z = z + rhot(l,a,c)*hamiltonian(c,b)
        endif
        !Vsoc contribution
        if (flag_so) then
          do m=1,3
            if ((abs(V_SO_red(a,c,m)) <= threshold) .and. (abs(V_SO_red(c,b,m)) <= threshold)) cycle
            do K_prime=k-1,k+1,1
              if ((K_prime < 0) .or. (K_prime > k_max)) cycle
              ! try: remove Q_prime loop, non-zero only for Q_prime=q-m+2
              !do Q_prime=-K_prime,K_prime,1
              !l_prime = get_kq_order(K_prime,Q_prime)
              l_prime = get_kq_order(K_prime,q-m+2)
              ! sum with real arguments for 6j symbols
              ! z = z + (-1)**(nint(sa+sb)-Q_prime) * &
              !     sqrt(dble(3*(2*k+1)*(2*K_prime+1))) * &
              !     W3J(dble(K_prime),One,dble(k),dble(Q_prime),dble(m-2),dble(-q)) * &
              !     (W6J(dble(K_prime),One,dble(k),sa,sb,sc)*rhot(l_prime,c,b)*V_SO_red(a,c,m) - &
              !     (-1)**(K_prime+k+1) * W6J(One,dble(K_prime),dble(k),sa,sb,sc)*rhot(l_prime,a,c)*V_SO_red(c,b,m))
              z = z + (-1)**((sa+sb)/2-q+m) * &
                sqrt(dble(3*(2*k+1)*(2*K_prime+1))) * &
                W3J(dble(K_prime),One,dble(k),dble(q-m+2),dble(m-2),dble(-q)) * &
                (W6J(2*K_prime,2,2*k,sa,sb,sc)*rhot(l_prime,c,b)*V_SO_red(a,c,m) - &
                (-1)**(K_prime+k+1) * W6J(2,2*K_prime,2*k,sa,sb,sc)*rhot(l_prime,a,c)*V_SO_red(c,b,m))
              !enddo
            enddo
          enddo
        endif
      enddo
      res(l,a,b) = - Onei*((E_SF(a)-E_SF(b))*rhot(l,a,b) + z)
    enddo
  enddo
enddo

end subroutine equation_sph
