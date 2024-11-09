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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************

!***********************************************************************
! History:                                                             *
! Jie J. Bao on May 09, 2022, created this file                        *
! Matthew R. Hennefarth Apr 05, 2023, upgraded to module               *
!***********************************************************************

module mspdft_util
  implicit none
  private

  public :: print_final_energies,print_mspdft_vectors,print_effective_ham
  public :: replace_diag

contains
  subroutine print_final_energies(e_mspdft,nroots)
    ! Prints the Final MS-PDFT Energies
    !
    ! Args:
    !   e_mspdft: ndarray of len nroots
    !     Final MS-PDFT energies
    !
    !   nroots: integer
    !     Number of roots in the calculation
    use definitions,only:iwp,wp,u6
    implicit none

    integer(kind=iwp),intent(in) :: nroots
    real(kind=wp),dimension(nroots),intent(in) :: e_mspdft

    integer(kind=iwp) :: state

    do state = 1,nroots
      call printresult(u6,'(6X,A,I3,A,F16.8)','MSPDFT root number',state,' Total energy:',e_mspdft(state),1)
    enddo
    write(u6,*)

  endsubroutine print_final_energies

  subroutine print_mspdft_vectors(si_pdft,nroots)
    ! Prints the final mspdft eigenvectors in the intermediate state
    ! basis and reference state basis
    !
    ! Args:
    !   si_pdft: ndarray of length nroots*nroots
    !     Contains the orthonormal eigenvectors in the intermediate
    !     state basis.
    !
    !   nroots: integer
    !     number of roots (lroots) or dimension of u and eigenvectors
    !
    use constants,only:zero
    use definitions,only:wp,iwp,u6
    implicit none

    integer(kind=iwp),intent(in) :: nroots
    real(kind=wp),dimension(nroots,nroots),intent(in) :: si_pdft

    logical :: refbas
    character(len=9),dimension(nroots) :: VecStat
    character(len=9) :: StatVec
    character(len=30) :: mspdftfmt
    character(Len=18) :: MatInfo

    integer(kind=iwp) :: root
    real(kind=wp),dimension(nroots,nroots) :: reference_vectors,eig_vecs_in_ref

    refbas = .false.

    do root = 1,nroots
      write(statvec,'(A5,I4)') 'Root ',root
      vecstat(root) = statvec
    enddo

    write(u6,*)

    write(u6,'(6X,A)') 'Intermediate-state Basis'
    write(mspdftfmt,'(A4,I5,A9)') '(6X,',nroots,'(A10,5X))'
    write(u6,mspdftfmt)((VecStat(root)),root=1,nroots)
    Call RecPrt(' ','(7X,10(F9.6,6X))',si_pdft,size(si_pdft,1),size(si_pdft,2))
    write(u6,*)

    call f_inquire('ROT_VEC',refbas)
    if(refbas) then
      ! Generate reference state basis
      eig_vecs_in_ref = zero
      call readmat2('ROT_VEC',MatInfo,reference_vectors,nroots,nroots,7,18,'T')
      call dgemm_('n','n',nroots,nroots,nroots,1.0d0,reference_vectors, &
                  nroots,si_pdft,nroots,0.0d0,eig_vecs_in_ref,nroots)
      write(u6,'(6X,A)') 'Reference-state Basis'
      write(u6,mspdftfmt)((VecStat(root)),root=1,nroots)
      call RecPrt(' ','(7X,10(F9.6,6X))',eig_vecs_in_ref,nroots,nroots)
      call PrintMat2('FIN_VEC',MatInfo,eig_vecs_in_ref,nroots,nroots,7,18,'T')
    endif

  endsubroutine print_mspdft_vectors

  subroutine print_effective_ham(heff,nroots,digit)
    ! Prints the effective Hamiltonian
    !
    ! Args:
    !   heff: ndarray of len nroots*nroots
    !     Effective MS-PDFT Hamiltonian
    !
    !   nroots: integer
    !     Dimension of heff, or number of roots
    !
    !   digit: integer
    !     Threshold value to shift diagonal elements by when printed
    use constants,only:zero
    use definitions,only:iwp,wp,u6
    use printlevel,only:silent
    use mcpdft_output,only:iPrLoc
    implicit none

    integer(kind=iwp),intent(in) :: nroots,digit
    real(kind=wp),dimension(nroots,nroots),intent(in) :: heff

    real(kind=wp),dimension(nroots,nroots) :: shifted_heff
    real(kind=wp) :: shift
    integer(kind=iwp) :: root

    if(iPrLoc(1) == silent) then
      return
    endif

    shifted_heff = heff

    call should_shift_diag(heff,nroots,digit,shift)

    if(shift /= zero) then
      write(u6,'(6X,A,F9.2,A)') '(diagonal values increased by',-shift,' hartree)'
      do root = 1,nroots
        shifted_heff(root,root) = shifted_heff(root,root)-shift
      enddo
    endif

    call recprt(' ','(7X,10(F9.6,1X))',shifted_heff,nroots,nroots)
    write(u6,*)
  endsubroutine print_effective_ham

  subroutine should_shift_diag(heff,nroots,digit,shift)
    ! Determines if diagonal elements should be shifted, and if so the number
    ! of decimal places to shift.
    !
    ! Args:
    !   heff: ndarray of len nroots*nroots
    !     Effective MS-PDFT Hamiltonian
    !
    !   nroots: integer
    !     Dimension of heff, or number of roots
    !
    !   digit: integer
    !     Threshold value to shift diagonal elements by when printed
    !
    ! Returns:
    !   shift: real
    !     Amount to shift diagonal elements by
    use constants,only:zero
    use definitions,only:iwp,wp
    implicit none

    integer(kind=iwp),intent(in) :: nroots,digit
    real(kind=wp),dimension(nroots,nroots),intent(in) :: heff

    real(kind=wp),intent(out) :: shift

    real(kind=wp),dimension(nroots) :: rdiag
    real(kind=wp) :: maxelem
    integer(kind=iwp) :: i,ishift

    shift = zero

    do i = 1,nroots
      rdiag(i) = heff(i,i)
    enddo

    maxelem = maxval(rdiag)

    if(abs(maxelem) < real(digit,8)) then
      return
    endif

    ishift = int(maxelem,8)/digit*digit
    shift = real(ishift,8)

  endsubroutine should_shift_diag

  subroutine replace_diag(mat,diag,ndim)
    use definitions,only:iwp,wp
    implicit none

    integer(kind=iwp),intent(in) :: ndim
    real(kind=wp),dimension(ndim),intent(in) :: diag
    real(kind=wp),dimension(ndim,ndim),intent(inout) :: mat

    integer(kind=iwp) :: i

    do i = 1,ndim
      mat(i,i) = diag(i)
    enddo

  endsubroutine replace_diag

endmodule mspdft_util
