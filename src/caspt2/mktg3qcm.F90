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

  ! state 1: bra
  ! state 2: ket
  Integer(kind=iwp), intent(in)  :: lsym1, lsym2, state1, state2, ntg3
  Real(kind=wp), intent(out) :: tg1(nasht, nasht), tg2(nasht, nasht, nasht, nasht)
  Real(kind=wp), intent(out) :: tg3(ntg3), ovl

  Real(kind=wp), allocatable :: tg3_tmp(:, :, :, :, :, :)
  Real(kind=wp) :: val
  Integer(kind=iwp) :: t, u, v, w, x, y, z
  Integer(kind=iwp) :: ituvxyz, sym_sig1, sym_sig2, sym_tau

  ! Detecting phase
  Real(kind=wp), allocatable :: tg1_tmp(:, :)


  write(*,*) "=== QCM: Building TRANSITION-RDM === "
  write(*,*) "between bra", state1, " and ket", state2

  ! This might be memory hungry
  call mma_allocate(tg3_tmp, nasht, nasht, nasht, nasht, nasht, nasht)
  tg3_tmp(:, :, :, :, :, :) = 0.0_wp
  tg2(:,:,:,:) = 0.0
  tg1(:,:) = 0.0

  ! Needed to detect phase
  ! We compute the transition RDM in two ways in order to detect if a phase 
  ! Flip occured whilst optimizing
  call mma_allocate(tg1_tmp, nasht, nasht)


  ! Remeasure TDM using rotated MPS
  call qcmaquis_interface_get_trans_1rdm_full(tg1, state2-1, state1-1)
  call qcmaquis_interface_get_trans_2rdm_full(tg2, state2-1, state1-1)
  call qcmaquis_interface_get_trans_3rdm_full(tg3_tmp, state2-1, state1-1)

  ! Fetch rotated TDM from unrotated MPS
  call qcmaquis_interface_read_rdm_full(int(state2-1, c_int), &
    int(state1-1, c_int), tg1_tmp, int(1, c_int))
  call qcmaquis_interface_read_rdm_full(int(state2-1, c_int), &
    int(state1-1, c_int), tg2, int(2, c_int))
  call qcmaquis_interface_read_rdm_full(int(state2-1, c_int), &
    int(state1-1, c_int), tg3_tmp, int(3, c_int))


  ! Detecting phase flip
  ! TODO: Implement this check, should abort
  ! currently fails as the two TDMs are not always the same
  if (dabs(tg1(1,1)) - dabs(tg1_tmp(1,1)) > 1e-9) then
    write(*,*) "Internal error: 1-TDM do not correspond to the same state"
     write(*,*) "TG1:"
     do u = 1, nasht
       do t = 1, nasht
        write(*, '(2I3,2F18.12)') t, u, tg1(t, u), tg1_tmp(t,u)
       end do
     end do
     call exit(1)
   end if

   ! Rotating MPS or TDM yields different phase different phase
   if (tg1(1, 1) * tg1_tmp(1,1) < 0.0) then
    write(*,*) "Phase flip detected: flipping phase"
    tg1(:,:) = -1.0 * tg1(:,:)
    tg2(:,:,:,:) = -1.0 * tg2(:,:,:,:) 
    tg3_tmp(:,:,:,:,:,:) = -1.0 * tg3_tmp(:,:,:,:,:,:) 
  end if
  call mma_deallocate(tg1_tmp)

  ! TODO: compute the overlap, check in the original code how
  ! ovl = qcmaquis_interface_get_overlap_with_ket_bra(int(state1-1, c_int), int(state2-1, c_int))

  do z = 1, nasht
    do y = 1, nasht
      ! symmetry of sigma2 = E_yz|Psi2>
      sym_sig2 = mul(mul(SGS%ism(y), SGS%ism(z)), lsym2)
      do x = 1, nasht
        do v = 1, nasht
          ! if (y + (z - 1) * nasht  < v + (x - 1) * nasht) then 
          !   cycle
          ! end if
          ! symmetry of tau = E_vx|sigma2>
          sym_tau = mul(mul(SGS%ism(x), SGS%ism(v)), sym_sig2)
          do u = 1, nasht
            do t = 1, nasht
            ! if (v + (x - 1) * nasht  < t + (u - 1) * nasht) then 
            !   cycle
            ! end if
              ! symmetry of sigma1 = <Psi1|E_tu
              sym_sig1 = mul(mul(SGS%ism(t), SGS%ism(u)), lsym1)
              ! only for for matching symmetries we have an element different from 0
              if (sym_sig1 == sym_tau) then
                ! generate the flat index
                call get_tg3_index(t, u, v, x, y, z, nasht, ituvxyz)
                tg3(ituvxyz) = tg3_tmp(t, v, y, u, x, z)
                ! if (tg3(ituvxyz) > 1e-3) then
                !   write(*,*) ituvxyz, tg3(ituvxyz)
                ! end if
              end if
            end do
          end do
        end do
      end do
    end do
  end do

   ! Debug output
   write(*,*) "TG1:"
   do u = 1, nasht
     do t = 1, nasht
      write(*, '(2I3,F18.12)') t, u, tg1(t, u)
     end do
   end do
   write(*,*) "TG2:"
   do x = 1, nasht
     do v = 1, nasht
       do u = 1, nasht
         do t = 1, nasht
            write(*, '(4I3,F18.12)') t, u, v, x, tg2(t, u, v, x)
         end do
       end do
     end do
   end do
   write(*,*) "TG3 (linearized):"
   do t = 1, ntg3
   end do

   ! write(*,*) "# TG1 as numpy array (for Python):"
   ! write(*,'(A)', advance='no') "dmrg1tdm = np.array(["
   ! do u = 1, nasht
   !   do t = 1, nasht
   !     write(*, '(F25.16, A)', advance='no') tg1(t,u), ", "
   !   end do
   ! end do
   ! write(*,'(A)') "]).reshape((4, 4), order='F')"
   ! 
   ! write(*,*) "# TG2 as numpy array (for Python):"
   ! write(*,'(A)', advance='no') "dmrg2tdm = np.array(["
   ! do x = 1, nasht
   !   do v = 1, nasht
   !     do u = 1, nasht
   !       do t = 1, nasht
   !     write(*, '(F25.16,A)', advance='no') tg2(t,u,v,x), ", "
   !       end do
   !     end do
   !   end do
   ! end do
   ! write(*,'(A)') "]).reshape((4, 4, 4, 4), order='F')"
   ! call exit(1)

  call mma_deallocate(tg3_tmp)
end subroutine mktg3qcm

#elif defined (NAGFOR)
! Some compilers do not like empty files
subroutine empty_mkfg3qcm()
end subroutine empty_mkfg3qcm
#endif


