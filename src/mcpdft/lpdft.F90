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
! Copyright (C) 2023, Matthew R. Hennefarth                            *
!***********************************************************************

module lpdft
  use stdalloc, only: mma_allocate, mma_deallocate
  use definitions, only: wp, iwp
  use mcpdft_output, only: lf
  implicit none
  private
  logical :: do_lpdft=.false.

  real(kind=wp) :: hconst = 0.0D0
  real(kind=wp), allocatable, dimension(:) :: inactD1
  real(kind=wp), allocatable, dimension(:) :: casd1_0, casd1s_0, casd2_0
  real(kind=wp), allocatable, dimension(:) :: veff1, veff2

  public :: do_lpdft
  public :: lpdft_kernel
  contains
    subroutine lpdft_kernel(CMO)
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
      real(kind=wp), dimension(ntot2), intent(in) :: CMO
      call init_vars(CMO)

      call set_iadr15()
      ! Want to first get the zero-order densities
      call get_zero_order_densities(casd1_0, casd1s_0, casd2_0)
      call compute_Eot_and_potentials(CMO)
      call compute_hconst()

      call release_vars()
    end subroutine lpdft_kernel

    subroutine get_zero_order_densities(casd1_0, casd1s_0, casd2_0)
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
      real(kind=wp), dimension(nacpar), intent(out) :: casd1_0, casd1s_0
      real(kind=wp), dimension(nacpr2), intent(out) :: casd2_0

      real(kind=wp), dimension(nacpar) :: casd1, casd1s
      real(kind=wp), dimension(nacpr2) :: casd2
      integer(kind=iwp) :: disk, root

      disk = iAdr15(3)
      do root=1, nroots
        call dDaFile(JOBOLD,2, casd1, nacpar, disk)
        call dDaFile(JOBOLD,2, casd1s, nacpar, disk)
        call dDaFile(JOBOLD,2, casd2, nacpr2, disk)
        ! Need a dummy read for some reason
        call dDaFile(JOBOLD,0, casd2, nacpr2, disk)

        call daxpy_(nacpar, weight(root), casd1, 1, casd1_0, 1)
        call daxpy_(nacpar, weight(root), casd1s, 1, casd1s_0, 1)
        call daxpy_(nacpr2, weight(root), casd2, 1, casd2_0, 1)
      end do
    end subroutine get_zero_order_densities

    subroutine compute_Eot_and_potentials(CMO)
      use KSDFT_Info, only: do_pdftpot
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
!For lRf
#include "rctfld.fh"

      real(kind=wp), dimension(ntot2), intent(in) :: CMO

      ! Some dummy vars that don't matter
      real(kind=wp), dimension(ntot1) :: htmp, gtmp
      ! The density matrices
      real(kind=wp), dimension(ntot2) :: casd1AO_0
      real(kind=wp), dimension(ntot1) :: d1ao_0, casD1AO_folded_0
      integer(kind=iwp) :: charge

      logical, parameter :: First=.True., Dff=.False., Do_DFT=.True.

      ! Convert CASD1 (MO) -> AO Basis
      call cas_mo_to_ao(CMO, casd1_0, casd1AO_0)

      call Fold(nSym, nBas, inactD1, d1ao_0)
      call Fold(nSym, nBas, casD1AO_0, casD1AO_folded_0)
      ! Create full 1RDM (inactive/frozen + active) and store in d1ao_0
      call daxpy_(nTot1, 1.0D0, casD1AO_folded_0, 1, d1ao_0, 1)

      ! Load in the nuclear potential
      Call Get_dScalar('PotNuc',potNuc)

      call get_charge(charge)

      ! We need the 1 and 2-electron potentials
      do_pdftpot = .True.

      veff1 = 0.0D0
      veff2 = 0.0D0
      call put_dArray('ONTOPO', veff1, ntot1)
      call put_dArray('ONTOPT', veff2, nfint)

      htmp = 0.0D0
      gtmp = 0.0D0
      call DrvXV(htmp, gtmp, d1ao_0, potnuc, nTot1, First, Dff, NonEq, &
                 lRF, KSDFT_TEMP, ExFac, charge, ispin, inactD1, &
                 casD1AO_0, nTot1, DFTFOCK, Do_DFT)

      ! Load in the on-top potentials
      call get_dArray("ONTOPO", veff1, ntot1)
      call get_dArray("ONTOPT", veff2, nfint)

    end subroutine compute_Eot_and_potentials

    subroutine compute_hconst()
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
      real(kind=wp) :: hnuc
      integer(kind=iwp) :: i

      Call Get_dScalar('PotNuc',hnuc)

      write(lf, *) "1E On Top Pot"
      do i=1, ntot1
        write(lf,*) veff1(i)
      end do
      write(lf, *) "2E On Top Pot"
      do i=1, nfint
        write(lf, *) veff2(i)
      end do

      hconst = hnuc
    end subroutine compute_hconst

    subroutine set_iadr15()
      ! Set the iAdr15 variable in rasscf.fh so I don't have to do any
      ! weird jobiph/jobold looking.
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
      logical :: found
      integer :: iad15

      iadr15 = 0
      iad15 = 0
      call f_inquire('JOBOLD', Found)
      if(.not. found) then
        call f_inquire('JOBIPH', Found)
        if(found) then
          JOBOLD=JOBIPH
        end if
      end if
      if(found) then
        if(JOBOLD <= 0) then
          JOBOLD = 20
          call DaName(JOBOLD, 'JOBOLD')
        end if
      end if
      call iDaFile(JOBOLD, 2, iadr15, 15, iad15)
    end subroutine set_iadr15

    subroutine init_vars(CMO)
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
      real(kind=wp), dimension(ntot2) :: CMO

      call mma_allocate(casd1_0, nacpar, "CASD1_0")
      call mma_allocate(casd1s_0, nacpar, "CASD1s_0")
      call mma_allocate(casd2_0, nacpr2, "CASD2_0")
      call mma_allocate(inactD1, ntot2, "inactD1")
      call mma_allocate(veff1, ntot1, "1E-OnTopPot")
      call mma_allocate(veff2, nfint, "2E-OnTopPot")

      call get_d1i_mcpdft(CMO, inactD1)

    end subroutine init_vars

    subroutine release_vars()
      call mma_deallocate(casd1_0)
      call mma_deallocate(casd1s_0)
      call mma_deallocate(casd2_0)
      call mma_deallocate(inactD1)
      call mma_deallocate(veff1)
      call mma_deallocate(veff2)
    end subroutine release_vars
end module lpdft