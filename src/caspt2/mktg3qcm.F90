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
! Copyright (C) 2025, Stefano Battaglia                                *
!***********************************************************************
#ifdef _DMRG_

subroutine mktg3qcm(lsym1, lsym2, state1, state2, ovl, tg1, tg2, ntg3, tg3)

  use stdalloc, only: mma_allocate, mma_deallocate
  use qcmaquis_interface
  use definitions, only: wp, iwp
  use gugx, only: SGS

  implicit none

#include "caspt2.fh"
! #include "pt2_guga.fh"

  Integer(kind=iwp), intent(in)  :: lsym1, lsym2, state1, state2, ntg3
  Real(kind=wp), intent(out) :: tg1(nasht, nasht), tg2(nasht, nasht, nasht, nasht)
  Real(kind=wp), intent(out) :: tg3(ntg3), ovl

  Real(kind=wp), allocatable :: tg3_tmp(:, :, :, :, :, :)
  Real(kind=wp) :: val
  Integer(kind=iwp) :: t, u, v, w, x, y, z
  Integer(kind=iwp) :: ituvxyz, sym_sig1, sym_sig2, sym_tau

  ! This might be memory hungry
  call mma_allocate(tg3_tmp, nasht, nasht, nasht, nasht, nasht, nasht)
  tg3_tmp(:, :, :, :, :, :) = 0.0_wp

  ! TODO: maybe we should implement an interface to get the tdms that mirrors that of rdms
  ! call qcmaquis_interface_get_1tdm_full(tg1, state1, state2)
  ! call qcmaquis_interface_get_2tdm_full(tg2, state1, state2)
  ! call qcmaquis_interface_get_3tdm_compact(tg3, ntg3, state1, state2)

  ! TODO: compute the overlap, check in the original code how

  do z = 1, nasht
    do y = 1, nasht
      ! symmetry of sigma2 = E_yz|Psi2>
      sym_sig2 = mul(mul(SGS%ism(y), SGS%ism(z)), lsym2)
      do x = 1, nasht
        do v = 1, nasht
          ! symmetry of tau = E_vx|sigma2>
          sym_tau = mul(mul(SGS%ism(x), SGS%ism(v)), sym_sig2)
          do t = 1, nasht
            do u = 1, nasht
              ! symmetry of sigma1 = <Psi1|E_tu
              sym_sig1 = mul(mul(SGS%ism(t), SGS%ism(u)), lsym1)
              ! only for for matching symmetries we have an element different from 0
              if (sym_sig1 == sym_tau) then
                ! generate the flat index
                call get_tg3_index(t, u, v, x, y, z, nasht, ituvxyz)
                tg3(ituvxyz) = tg3_tmp(t, u, v, x, y, z)
              end if
            end do
          end do
        end do
      end do
    end do
  end do

  call mma_deallocate(tg3_tmp)
end subroutine mktg3qcm

#elif defined (NAGFOR)
! Some compilers do not like empty files
subroutine empty_mkfg3qcm()
end subroutine empty_mkfg3qcm
#endif


