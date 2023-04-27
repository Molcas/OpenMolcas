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

module mspdft
  implicit none
  private

#include "mspdft.fh"

  character(len=8) :: mspdftmethod
  logical :: do_rotate = .False.
  integer :: iF1MS,iF2MS,iFxyMS,iFocMS,iIntS,iDIDA,IP2MOt
  integer :: D1AOMS,D1SAOMS

  ! CMS-NACS stuff
  logical :: DoNacMSPD, DoMeciMSPD, isCMSNAC
  integer :: cmsNACStates(2)

  public :: dogradmspd, mspdftmethod, do_rotate, iF1MS, iF2MS
  public :: iFxyMS, iFocMS, iIntS, iDIDA, IP2MOt, D1AOMS, D1SAOMS

  public :: DoNacMSPD, DoMeciMSPD, isCMSNAC, cmsNACStates

  public :: mspdft_finalize

  contains
  subroutine mspdft_finalize(heff, nroots, irlxroot, iadr19)
    ! Performs the final parts of MS-PDFT, namely:
    !   1. diagonalize heff
    !   2. print out the final results
    !   3. save to jobiph if requested
    !
    ! Args:
    !   heff: ndarray of length nroots*nroots
    !     MS-PDFT effective hamiltonian matrix. Diagonal elements should
    !     already be replaced with the correct energies.
    !
    !   nroots: integer
    !     number of roots, or dimension of heff
    !
    !   irlxroot: integer
    !     root to relax for some reason it is needed
    !
    !   iadr19: integer list of length 15
    !     Holds some information when writing out

    use definitions, only: wp
    use stdalloc, only: mma_allocate, mma_deallocate
    use mcpdft_output, only: lf
    use mspdft_util, only: print_effective_ham, print_final_energies, &
                           print_mspdft_vectors
    use write_pdft_job, only: iwjob, writejob
    use hybridpdft, only: do_hybrid

    integer, intent(in) :: nroots, irlxroot
    real(kind=wp), dimension(nroots**2), intent(in) :: heff
    integer, dimension(15), intent(in) :: IADR19

    real(kind=wp), dimension(nroots) :: e_mspdft
    real(kind=wp), dimension(nroots**2) :: si_pdft
    integer :: i, info, dim_scratch
    real(kind=wp), dimension(1) :: wgronk
    real(kind=wp), dimension(:), allocatable :: scratch

    write(lf, '(6X,80a)') ('*',i=1,80)
    write(lf,*)
    write(lf,'(34X,2A)')MSPDFTMethod,' FINAL RESULTS'
    write(lf,*)
    write(lf,'(6X,80a)') ('*',i=1,80)
    write(lf,*)

    if(do_hybrid) then
      write(lf,'(6X,3A)') 'Hybrid ', MSPDFTMethod, ' Effective Hamiltonian'
    else
      write(lf, '(6X,2A)') mspdftmethod, ' Effective Hamiltonian'
    end if
    call print_effective_ham(heff, nroots, 10)

    ! Now we diagonalize the final MS-PDFT H-matrix
    ! Eigenvectors will be stored in heff.
    ! Eigenvalues will be stored in e_mspdft

    ! Since the dsyev_ call will override the heff with the orthonormal
    ! eigenvectors, I am going to set the new variable now.
    si_pdft = heff

    ! This first call is to get how big of a scratch we need
    call dsyev_('V','U', nroots, si_pdft, nroots, e_mspdft, wgronk, -1, info)

    dim_scratch = int(wgronk(1))
    call mma_allocate(scratch, dim_scratch, label="XScratch")
    ! Now we actually do the diagonalization
    call dsyev_('V', 'U', nroots, si_pdft, nroots, e_mspdft, scratch, dim_scratch, info)
    call mma_deallocate(scratch)

    call print_final_energies(e_mspdft, nroots, mspdftmethod)

    ! Update information on the runfile for possible gradient calculations.
    call put_dArray('Last energies', e_mspdft, nroots)
    call Put_dScalar('Last energy', e_mspdft(iRlxRoot))

    if(do_hybrid) then
      write(lf, '(6X,3A)') 'Hybrid ',MSPDFTMethod,' Eigenvectors:'
    else
        write(lf,'(6X,2A)') MSPDFTMethod,' Eigenvectors:'
    end if
    call print_mspdft_vectors(si_pdft, nroots)

    ! Added by Chen to write energies and states of MS-PDFT into JOBIPH
    if (iwjob==1) call writejob(iadr19, e_mspdft=e_mspdft, si_pdft=si_pdft)

    write(lf,'(6X,80a)') ('*',i=1,80)

    if (dogradmspd) then
      if (doNACMSPD) then
        call MSPDFTNAC_Misc(si_pdft)
      else
        call mspdftgrad_misc(si_pdft)
      end if
    end if

  end subroutine mspdft_finalize
end module mspdft