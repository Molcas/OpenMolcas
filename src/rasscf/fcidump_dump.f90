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
! Copyright (C) 2014, Giovanni Li Manni                                *
!               2019, Oskar Weser                                      *
!***********************************************************************
module fcidump_dump
  use fcidump_tables
  implicit none
  private
  public :: dump_ascii, dump_hdf5
contains

!>  @brief
!>    Create FCIDUMP file
!>
!>  @author Oskar Weser
!>
!>  @details
!>  Create an ASCII formatted FCIDUMP with core energy,
!>  orbital energies, Fock matrix elements and two electron integrals.
!>  Contains information about \p nAsh, \p nSym, \p nActEl,
!>  \p iSpin, and \p lSym.
!>
!>  @param[in] EMY Core energy
!>  @param[in] orbital_table Orbital energies with index
!>  @param[in] fock_table
!>  @param[in] two_el_table
  subroutine dump_ascii(path, EMY, orbital_table, fock_table, &
                        two_el_table, orbsym)
    use general_data, only : nSym, nActEl, iSpin, lSym, nAsh
    implicit none
    character(len=*), intent(in) :: path
    real*8, intent(in) :: EMY
    type(OrbitalTable), intent(in) :: orbital_table
    type(FockTable), intent(in) :: fock_table
    type(TwoElIntTable), intent(in) :: two_el_table
    integer, intent(in) :: orbsym(:)
    integer :: i, j, ireturn, isFreeUnit, LuFCI

    call qEnter('dump_ascii')

    LuFCI = isFreeUnit(38)
    call molcas_open(LuFCI, path)

    write(LuFCI,'(1X,A11,I3,A7,I3,A5,I3,A)') ' &FCI NORB=',sum(nAsh), &
       ',NELEC=',nActEl,',MS2=',int((ISPIN-1.0d0)),','
    write(LuFCI,'(A,500(I2,","))')'  ORBSYM=',(orbsym(i),i=1,size(orbsym))
    write(LuFCI,'(2X,A5,I1)') 'ISYM=', LSYM -1
    write(LuFCI,'(A)') ' &END'

    do j = 1, length(two_el_table)
      write(LuFCI, '(1X,G20.11,4I5)') &
        two_el_table%values(j), (two_el_table%index(i, j), i = 1, 4)
    end do

    do j = 1, length(fock_table)
      write(LuFCI, '(1X,G20.11,4I5)') &
        fock_table%values(j), (fock_table%index(i, j), i = 1, 2), 0, 0
    end do

    do j = 1, length(orbital_table)
      write(LuFCI, '(1X,G20.11,4I5)') &
        orbital_table%values(j), orbital_table%index(j), 0, 0, 0
    end do

    write(LuFCI,'(1X,G20.11,4I5)') EMY, 0, 0, 0, 0

    close(LuFCI)

! ========== For testing purposes FROM HERE =============
    if (length(orbital_table) /= 0 .and. length(fock_table) /=0 &
        .and. length(two_el_table) /= 0) then
      call Add_Info('core energy', [EMY], 1, 8)
      call Add_Info('Orbital Energy', orbital_table%values(1), 1, 8)
      call Add_Info('Fock element', fock_table%values(1), 1, 8)
      call Add_Info('TwoEl Integral element', two_el_table%values(1), 1, 8)
    end if
! ========== For testing purposes TO HERE ===============

    call FastIO('STATUS')
    ireturn = 0
    call qExit('dump_ascii')

    return
  end subroutine dump_ascii

!>  @brief
!>    Create H5FCIDUMP file
!>
!>  @author Oskar Weser
!>
!>  @details
!>  Create an FCIDUMP in HDF5-file format with core energy,
!>  orbital energies, Fock matrix elements and two electron integrals.
!>  Contains information about \p nAsh, \p nSym, \p nActEl,
!>  \p iSpin, and \p lSym.
!>
!>  @param[in] EMY Core energy
!>  @param[in] orbital_table Orbital energies with index
!>  @param[in] fock_table
!>  @param[in] two_el_table
  subroutine dump_hdf5(path, EMY, orbital_table, fock_table, two_el_table, orbsym)
    use general_data, only : nSym, nActEl, multiplicity => iSpin, lSym, nAsh
    use gas_data, only : iDoGAS
    use gugx_data, only : IfCAS
#ifdef _HDF5_
    use mh5
#endif
    implicit none
    character(len=*), intent(in) :: path
    real*8, intent(in) :: EMY
    type(OrbitalTable), intent(in) :: orbital_table
    type(FockTable), intent(in) :: fock_table
    type(TwoElIntTable), intent(in) :: two_el_table
    integer, intent(in) :: orbsym(:)
#ifdef _HDF5_
    integer :: file_id, dset_id
    integer :: ireturn
    character :: lIrrep(24)

    call qEnter('dump_hdf5')

    file_id = mh5_create_file(path)

!   symmetry information
    call Get_cArray('Irreps', lIrrep, 24)
    call mh5_init_attr(file_id, 'IRREP_LABELS', 1, [nSym], lIrrep, 3)
    call mh5_init_attr(file_id, 'ORBSYM', 1, [size(orbsym)], orbsym)
    call mh5_init_attr(file_id, 'NORB', sum(nAsh(:nSym)))

    ! Set wavefunction type
    if (iDoGAS) then
      call mh5_init_attr(file_id, 'CI_TYPE', 'GAS')
    else if (IfCAS .eq. 0) then
      call mh5_init_attr(file_id, 'CI_TYPE', 'CAS')
    else
      call mh5_init_attr(file_id, 'CI_TYPE', 'RAS')
    end if

    call mh5_init_attr(file_id, 'MOLCAS_MODULE', 'RASSCF')
    call mh5_init_attr(file_id, 'CORE_ENERGY', EMY)
    call mh5_init_attr(file_id, 'NELEC', nActEl)
    call mh5_init_attr(file_id, 'MULTIPLICITY', multiplicity)
    call mh5_init_attr(file_id, 'ISYM', lSym - 1)

    dset_id = mh5_create_dset_int(file_id, 'ORBITAL_INDEX', &
      1, [length(orbital_table)])
    call mh5_init_attr(dset_id, 'DESCRIPTION', &
      'Index for the orbitals in active space.')
    call mh5_put_dset_array_int(dset_id, orbital_table%index)
    call mh5_close_dset(dset_id)

    dset_id = mh5_create_dset_real(file_id, 'ORBITAL_ENERGIES', &
      1, [length(orbital_table)])
    call mh5_init_attr(dset_id, 'DESCRIPTION', &
      'Energies of orbitals in active space.')
    call mh5_put_dset_array_real(dset_id, orbital_table%values)
    call mh5_close_dset(dset_id)

    dset_id = mh5_create_dset_int(file_id, 'FOCK_INDEX', &
      2, [2, length(fock_table)])
    call mh5_init_attr(dset_id, 'DESCRIPTION', &
      'The index i, j for the Fock matrix elements <i| F |j>.')
    call mh5_put_dset_array_int(dset_id, fock_table%index)
    call mh5_close_dset(dset_id)

    dset_id = mh5_create_dset_real(file_id, 'FOCK_VALUES', &
      1, [length(fock_table)])
    call mh5_init_attr(dset_id, 'DESCRIPTION', &
      'The Fock matrix elements <i| F |j>.')
    call mh5_init_attr(dset_id, 'CUTOFF', fock_table%cutoff)
    call mh5_put_dset_array_real(dset_id, fock_table%values)
    call mh5_close_dset(dset_id)

    dset_id = mh5_create_dset_int(file_id, 'TWO_EL_INT_INDEX', &
      2, [4, length(two_el_table)])
    call mh5_init_attr(dset_id, 'DESCRIPTION', &
      'The index i, j, k, l for the two electron integrals <i j | 1/r_{12} | k l >.')
    call mh5_put_dset_array_int(dset_id, two_el_table%index)
    call mh5_close_dset(dset_id)

    dset_id = mh5_create_dset_real(file_id, 'TWO_EL_INT_VALUES', &
      1, [length(two_el_table)])
    call mh5_init_attr(dset_id, 'DESCRIPTION', &
      'The two electron integrals <i j | 1/r_{12} | k l >.')
    call mh5_init_attr(dset_id, 'CUTOFF', two_el_table%cutoff)
    call mh5_put_dset_array_real(dset_id, two_el_table%values)
    call mh5_close_dset(dset_id)

    call mh5_close_file(file_id)

    call FastIO('STATUS')
    ireturn = 0
    call qExit('dump_hdf5')
#else
! Avoid unused argument warnings
    if (.false.) then
      call unused_real(EMY)
      call unused_character(path)
      call unused(orbital_table)
      call unused(fock_table)
      call unused(two_el_table)
      call unused_integer_array(orbsym)
    end if
#endif
  end subroutine dump_hdf5
end module fcidump_dump
