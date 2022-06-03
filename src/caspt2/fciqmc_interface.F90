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
! Copyright (C) 2022, Arta Safari                                      *
!***********************************************************************
module fciqmc_interface

  use definitions, only: wp, u6
  ! use caspt2_data, only: epsa
  use linalg_mod, only: verify_, abort_
  use fortran_strings, only: str
  use stdalloc, only: mma_allocate, mma_deallocate
  use mh5, only: mh5_open_file_r, mh5_close_file, mh5_fetch_dset, &
                 mh5_open_group, mh5_close_group

  implicit none

  public :: DoFCIQMC
  logical :: DoFCIQMC = .false.

  contains

    !>  @brief
    !>    Wrapper to collect all required density and Fock matrices and feed
    !>    them into poly3.
    !>
    !>  @param[in]     iff       integer specifying whether contractions of 1-4RDMs with
    !>                           the Fock matrix should be computed. 1 = yes, else no.
    !>  @param[in]     idxG3     Table containing the active space indices
    !>  @param[out]    g1        dense, redundant 1RDM
    !>  @param[out]    g2        dense, redundant 2RDM
    !>  @param[out]    g3        sparse 3RDM
    !>  @param[out]    f1        dense contraction of Fockian with 2RDM
    !>  @param[out]    f2        dense contraction of Fockian with 3RDM
    !>  @param[out]    f3        sparse contraction of Fockian with 4RDM
    ! subroutine mkfg3neci(iff, idxG3, g1, g2, g3, f1, f2, f3)
    !      integer, intent(in) :: iff
    !      integer(kind=1), intent(in) :: idxG3(6, *)
    !      ! TODO: what is the difference between nlev and nact?
    !      real(wp), intent(out) :: g1(nlev, nlev), g2(nlev, nlev, nlev, nlev), g3(*), &
    !                               f1(nlev, nlev), f2(nlev, nlev, nlev, nlev), f3(*)

    !      if (nAc < 4) then
    !        call abort_('FCIQMC-CASPT2 requires at least 4 electrons.')
    !      else if (nAc >= 4) then
    !        call read_twordm(nlev, g2, mstate(jstate))
    !        call read_f3_or_g3("3RDM", nAc, idxG3, nG3, g3, epsa, f2, &
    !                           mstate(jstate))
    !        call read_f3_or_g3("cF4RDM", nAc, idxG3, nG3, f3, epsa, f2, &
    !                           mstate(jstate))
    !        g1 = calc_onerdm(nlev, nAc, g2)
    !        f1 = calc_f1_matrix(iff, nlev, f1, g2, epsa)
    !        f2 = calc_f2_matrix()
    !      end if
    ! end subroutine mkfg3neci


    !>  @brief
    !>    Read stochastically sampled 2RDMs stored in HDF5 format.
    ! subroutine read_twordm(nlev, g2, iroot)
    !   integer, intent(in) :: nlev
    !   real(wp), intent(_OUT_) :: g2(nlev, nlev, nlev, nlev)

    !   continue
    ! end subroutine


    !>  @brief
    !>    Read stochastically sampled 3RDMs and contracted Fock tensors
    !>    stored in HDF5 format. When 3RDMs are read, the contraction with
    !>    the Fockian (f2) is performed on the fly.
    !>
    !>  @param[in]     name       String, either "3RDM" or "contractedFock"
    !>  @param[in]     nAc        Number of active electrons
    !>  @param[in]     idxG3      Table containing the active space indices
    !>  @param[in]     nG3
    !>  @param[_OUT_]  storage    Array into which either a 3RDM or contracted
    !>                            Fock matrix is read
    !>  @param[in]     epsa       Orbital energies in the active space
    !>  @param[in]     f2         Contraction of Fock matrix with 3RDM
    !>  @param[in]     iroot      MCSCF root number.
    ! subroutine read_f3_or_g3(name, nAc, idxG3, nG3, storage, epsa, f2, iroot)
    !   character(len=*), intent(in) :: name
    !   integer, intent(in) :: nAc, idxG3(6, nG3), nG3
    !   real(wp), intent(_OUT_) :: storage(*)
    !   real(wp), intent(in) :: epsa(nAc)
    !   real(wp), intent(out) :: f2(nAc, nAc, nAc, nAc)
    !   integer, intent(in) :: iroot
    !   integer :: file_id, iG3, ip1, iq1, ip2, iq2, ip3, iq3, idx
    !   logical :: tExist
    !   real(wp), allocatable :: buffer(:)

    !   ! allocate the tensor with 6 active indices
    !   call mma_allocate(buffer, nAc**6)
    !   if (name == '3RDM') then
    !     call f_Inquire('spinfree-ThreeRDM.' // str(iroot) // '.h5', tExist)
    !     call verify_(tExist, 'spinfree-ThreeRDM.' // str(iroot) // '.h5 does not exist.')
    !     hdf5_file = mh5_open_file_r('spinfree-ThreeRDM.' // str(iroot) // '.h5')
    !     hdf5_group = mh5_open_group(hdf5_file, '3RDM')
    !   else if (name == 'contractedFock') then
    !     call f_Inquire('cF4RDM.' // str(iroot) // '.h5', tExist)
    !     call verify_(tExist, 'cF4RDM.' // str(iroot) // '.h5 does not exist.')
    !     hdf5_file = mh5_open_file_r('cF4RDM.' // str(iroot) // '.h5')
    !     hdf5_group = mh5_open_group(hdf5_file, 'cF4RDM')
    !   end if
    !   call mh5_fetch_dset(hdf5_group, 'elements', buffer)
    !   call mh5_close_group(hdf5_group)
    !   call mh5_close_file(hdf5_file)

    !   ! pardon these magic numbers, they are taken from the ChemPS2 interface...
    !   do iG3=1,nG3
    !     ip1 = idxG3(1, iG3) - 1
    !     iq1 = idxG3(2, iG3) - 1
    !     ip2 = idxG3(3, iG3) - 1
    !     iq2 = idxG3(4, iG3) - 1
    !     ip3 = idxG3(5, iG3) - 1
    !     iq3 = idxG3(6, iG3) - 1
    !     idx = ip1 + nAc * (ip2 + nAc * (ip3 + nAc * (iq1 + nAc * (iq2 + nAc * iq3))))
    !     storage(iG3) = buffer(1 + idx)
    !   end do

    !   if (name == "3RDM") then
    !     do iq2 = 1, nAc
    !       do ip2 = 1, nAc
    !         do iq1 = 1, nAc
    !           do ip1 = 1, nAc
    !             f2(ip1, iq1, ip2, iq2) = 0.0_wp
    !             do ip3 = 1, nAc
    !               idx = ip1 + nAc * (ip2 - 1 + nAc *(ip3 - 1 +nAc &
    !                            * (iq1 - 1 + nAC * (iq2 - 1 + nAc * (ip3 - 1)))))
    !               f2(ip1, iq1, ip2, iq2) = f2(ip1, iq1, ip2, iq2) &
    !                            + epsa(ip3) * buffer(idx)
    !             end do
    !           end do
    !         end do
    !       end do
    !     end do
    !   end if

    !   mma_deallocate(buffer)
    ! end subroutine read_f3_or_g3


    !>  @brief
    !>    Contracts one pair of indices from 2RDM to retrieve 1RDM.
    pure function calc_onerdm(g2, nAct) result(g1)
      real(wp), intent(in) :: g2(nAct, nAct, nAct, nAct)
      integer, intent(in) :: nAct
      real(wp) :: g1(nAct, nAct)
      integer :: p, q, r

      g1(:,:) = 0.0_wp
      do r = 1, nAct
        do q = 1, nAct
          do p = 1, nAct
            g1(p, r) = g1(p, r) + g2(p, q, r, q)
          end do
        end do
      end do
      g1 = g1 / (nAct - 1)
    end function calc_onerdm


    !>  @brief
    !>    Contracts one pair of indices from 2RDM with Fock matrix.
    ! pure function calc_f1_matrix(g2, epsa) result(f1)
    !   real(wp), intent(in) :: g2(nAct, nAct, nAct, nAct), epsa(nAct)
    !   integer, intent(in) :: nAct
    !   real(wp) :: f1(nAct, nAct)

    !   f1(:,:) = 0.0_wp
    !   do iz=1,nlev
    !     iysym=ism(iz)
    !     do iy=1,nlev
    !       ixysym=mul(ism(iy),iysym)
    !       if(iff.ne.0.and.ixysym.eq.1) then
    !         f1(iy,iz) = 0.0
    !         do iw=1,nlev
    !           f1(iy,iz)=f1(iy,iz)+g2(iw,iw,iy,iz)*epsa(iw)
    !         end do
    !       end if
    !     end do
    !   end do
    ! end function calc_f1_matrix


    !>  @brief
    !>    Contracts one pair of indices from 3RDM with Fock matrix.
    ! pure function calc_f2_matrix(g3, epsa) result(f2)
    !   real(wp), intent(in) :: g3(*), epsa(nAct)
    !   integer, intent(in) :: nAct
    !   real(wp) :: f2(nAct, nAct, nAct, nAct)

    !   f2(:,:,:,:) = 0.0_wp
    !   do iz=1,nlev
    !     iysym=ism(iz)
    !     do iy=1,nlev
    !       ixysym=mul(ism(iy),iysym)
    !       if(iff.ne.0.and.ixysym.eq.1) then
    !         f1(iy,iz) = 0.0
    !         do iw=1,nlev
    !           f1(iy,iz)=f1(iy,iz)+g2(iw,iw,iy,iz)*epsa(iw)
    !         end do
    !       end if
    !     end do
    !   end do
    ! end function calc_f2_matrix

end module fciqmc_interface
