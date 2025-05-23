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

use rhodyn_data, only: d, flag_pulse, hamiltonian, hamiltoniant, ipglob, k_max, k_ranks, len_sph, lroots, mirr, n, q_max, q_proj, &
                       Y1, Y2
use rhodyn_utils, only: get_kq_order
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cOne, cZero, Onei
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: time
complex(kind=wp), intent(in) :: rhot(len_sph,d,d)
complex(kind=wp), intent(out) :: res(len_sph,d,d)
integer(kind=iwp) :: c, i, k, K_prime, l, l_prime, m, n1, n2, q
complex(kind=wp), allocatable :: rhot_tmp(:), res_tmp(:)

if (ipglob > 3) write(u6,*) 'solve equation, time: ',time

! number of states with ground state spin to distinguish contributions to different spin manifolds
! now works only with 2 spin manifolds
!ngs = lroots(1)

! if pulse is enabled, modify Hamiltonian at time t:
if (flag_pulse) call pulse(hamiltonian,hamiltoniant,time,-1)

call mma_allocate(rhot_tmp,d**2)
call mma_allocate(res_tmp,d**2)

i = 1
do l=1,len_sph
  k = k_ranks(l)
  q = q_proj(l)
  if ((q >= -q_max) .and. (q <= 0)) then
    ! get right part of Liouville equation -i*(hamiltoniant*rhot - rhot*hamiltoniant)
    rhot_tmp(:) = pack(rhot(l,:,:),.true.)
    call zgemm_('N','N',d,d,d,-Onei,hamiltoniant,d,rhot_tmp,d,cZero,res_tmp,d)
    call zgemm_('N','N',d,d,d,Onei,rhot_tmp,d,hamiltoniant,d,cOne,res_tmp,d)
    res(l,:,:) = reshape(res_tmp,[d,d])
    ! Vsoc contribution
    do m=1,3
      do K_prime=k-1,k+1,1
        if ((K_prime < 0) .or. (K_prime > k_max)) cycle
        if ((q-m+2 > K_prime) .or. (q-m+2 < -K_prime)) cycle
        l_prime = get_kq_order(K_prime,q-m+2)
        do c=1,n
          ! set boundaries for slice multiplication
          if (c == 1) then
            n1 = 1
            n2 = lroots(1)
          else
            n1 = n1+lroots(c-1)
            n2 = n2+lroots(c)
          end if
          rhot_tmp(1:d*lroots(c)) = pack(rhot(l_prime,:,n1:n2),.true.)
          call zgemm_('N','N',d,lroots(c),d,-Onei,Y1(:,:,i),d,rhot_tmp,d,cZero,res_tmp,d)
          res(l,:,n1:n2) = res(l,:,n1:n2)+reshape(res_tmp(1:d*lroots(c)),[d,lroots(c)])
          rhot_tmp(1:d*lroots(c)) = pack(rhot(l_prime,n1:n2,:),.true.)
          call zgemm_('N','N',lroots(c),d,d,-Onei,rhot_tmp,lroots(c),Y2(:,:,i),d,cZero,res_tmp,lroots(c))
          res(l,n1:n2,:) = res(l,n1:n2,:)+reshape(res_tmp(1:d*lroots(c)),[lroots(c),d])
          i = i+1
        end do
      end do
    end do
  ! mirror from q=-1 to q=1
  else if (q == 1) then
    res(l,:,:) = transpose(conjg(res(l-2,:,:)))*mirr
  end if
end do

call mma_deallocate(rhot_tmp)
call mma_deallocate(res_tmp)

end subroutine equation_sph
