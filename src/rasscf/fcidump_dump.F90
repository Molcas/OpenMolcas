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
<<<<<<< HEAD
!               2026, Nike Dattani                                     *
!               2026, Jaafar Mehrez                                    *
=======
>>>>>>> upstream-openmolcas/master
!***********************************************************************

#include "macros.fh"

module fcidump_dump

use fcidump_tables, only: FockTable, length, OrbitalTable, TwoElIntTable
use Definitions, only: wp, iwp

implicit none
private

<<<<<<< HEAD
public :: dump_ascii, dump_fort55, dump_hdf5
=======
public :: dump_ascii, dump_hdf5
>>>>>>> upstream-openmolcas/master

contains

!>  @brief
!>    Create FCIDUMP file
!>
!>  @author Oskar Weser
!>
!>  @details
!>  Create an ASCII formatted FCIDUMP with core energy,
!>  orbital energies, Fock matrix elements and two electron integrals.
!>  Contains information about \p nAsh, \p nActEl,
!>  \p iSpin, and \p stSym.
!>
!>  @param[in] path
!>  @param[in] EMY Core energy
!>  @param[in] orbital_table Orbital energies with index
!>  @param[in] fock_table
!>  @param[in] two_el_table
!>  @param[in] orbsym
subroutine dump_ascii(path,EMY,orbital_table,fock_table,two_el_table,orbsym)

  use general_data, only: nActEl, iSpin, stSym, nAsh

  character(len=*), intent(in) :: path
  real(kind=wp), intent(in) :: EMY
  type(OrbitalTable), intent(in) :: orbital_table
  type(FockTable), intent(in) :: fock_table
  type(TwoElIntTable), intent(in) :: two_el_table
  integer(kind=iwp), intent(in) :: orbsym(:)
  integer(kind=iwp) :: i, j, LuFCI
  integer(kind=iwp), external :: isFreeUnit

  LuFCI = isFreeUnit(38)
  call molcas_open(LuFCI,path)

  write(LuFCI,'(1X,A11,I3,A7,I3,A5,I3,A)') ' &FCI NORB=',sum(nAsh),',NELEC=',nActEl,',MS2=',ISPIN-1,','
  write(LuFCI,'(A,500(I2,","))') '  ORBSYM=',(orbsym(i),i=1,size(orbsym))
  write(LuFCI,'(2X,A5,I1)') 'ISYM=',STSYM-1
  write(LuFCI,'(A)') ' &END'

  do j=1,length(two_el_table)
<<<<<<< HEAD
    write(LuFCI,'(1X,E27.20,4I5)') two_el_table%values(j),(two_el_table%idx(i,j),i=1,4)
  end do

  do j=1,length(fock_table)
    write(LuFCI,'(1X,E27.20,4I5)') fock_table%values(j),(fock_table%idx(i,j),i=1,2),0,0
  end do

  do j=1,length(orbital_table)
    write(LuFCI,'(1X,E27.20,4I5)') orbital_table%values(j),orbital_table%idx(j),0,0,0
  end do

  write(LuFCI,'(1X,E27.20,4I5)') EMY,0,0,0,0
=======
    write(LuFCI,'(1X,G20.11,4I5)') two_el_table%values(j),(two_el_table%idx(i,j),i=1,4)
  end do

  do j=1,length(fock_table)
    write(LuFCI,'(1X,G20.11,4I5)') fock_table%values(j),(fock_table%idx(i,j),i=1,2),0,0
  end do

  do j=1,length(orbital_table)
    write(LuFCI,'(1X,G20.11,4I5)') orbital_table%values(j),orbital_table%idx(j),0,0,0
  end do

  write(LuFCI,'(1X,G20.11,4I5)') EMY,0,0,0,0
>>>>>>> upstream-openmolcas/master

  close(LuFCI)

  ! ========== For testing purposes FROM HERE =============
  if ((length(orbital_table) /= 0) .and. (length(fock_table) /= 0) .and. (length(two_el_table) /= 0)) then
    call Add_Info('core energy',[EMY],1,8)
    call Add_Info('Orbital Energy',orbital_table%values(1),1,8)
    call Add_Info('Fock element',fock_table%values(1),1,8)
    call Add_Info('TwoEl Integral element',two_el_table%values(1),1,8)
  end if
  ! ========== For testing purposes TO HERE ===============

  call FastIO('STATUS')

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
!>  \p iSpin, and \p stSym.
!>
!>  @param[in] path
!>  @param[in] EMY Core energy
!>  @param[in] orbital_table Orbital energies with index
!>  @param[in] fock_table
!>  @param[in] two_el_table
!>  @param[in] orbsym
subroutine dump_hdf5(path,EMY,orbital_table,fock_table,two_el_table,orbsym)

# ifdef _HDF5_
  use general_data, only: nSym, nActEl, multiplicity => iSpin, stSym, nAsh
<<<<<<< HEAD
  use gas_data, only: iDoGAS
  use sguga, only: SGS
=======
  use sguga_states, only: SGS
  use gas_data, only: iDoGAS
>>>>>>> upstream-openmolcas/master
  use mh5, only: mh5_close_dset, mh5_close_file, mh5_create_file, mh5_create_dset_int, mh5_create_dset_real, mh5_init_attr, &
                 mh5_put_dset
# endif

  character(len=*), intent(in) :: path
  real(kind=wp), intent(in) :: EMY
  type(OrbitalTable), intent(in) :: orbital_table
  type(FockTable), intent(in) :: fock_table
  type(TwoElIntTable), intent(in) :: two_el_table
  integer(kind=iwp), intent(in) :: orbsym(:)
# ifdef _HDF5_
  integer(kind=iwp) :: dset_id, file_id
  character :: lIrrep(24)
<<<<<<< HEAD
=======
  integer(kind=iwp), parameter:: istate=1
>>>>>>> upstream-openmolcas/master

  file_id = mh5_create_file(path)

  ! symmetry information
  call Get_cArray('Irreps',lIrrep,24)
  call mh5_init_attr(file_id,'IRREP_LABELS',1,[nSym],lIrrep,3)
  call mh5_init_attr(file_id,'ORBSYM',1,[size(orbsym)],orbsym)
  call mh5_init_attr(file_id,'NORB',sum(nAsh(:nSym)))

  ! Set wavefunction type
  if (iDoGAS) then
    call mh5_init_attr(file_id,'CI_TYPE','GAS')
<<<<<<< HEAD
  else if (SGS%IFRAS == 0) then
=======
  else if (SGS(istate)%IFRAS == 0) then
>>>>>>> upstream-openmolcas/master
    call mh5_init_attr(file_id,'CI_TYPE','CAS')
  else
    call mh5_init_attr(file_id,'CI_TYPE','RAS')
  end if

  call mh5_init_attr(file_id,'MOLCAS_MODULE','RASSCF')
  call mh5_init_attr(file_id,'CORE_ENERGY',EMY)
  call mh5_init_attr(file_id,'NELEC',nActEl)
  call mh5_init_attr(file_id,'MULTIPLICITY',multiplicity)
  call mh5_init_attr(file_id,'ISYM',stSym-1)

  dset_id = mh5_create_dset_int(file_id,'ORBITAL_INDEX',1,[length(orbital_table)])
  call mh5_init_attr(dset_id,'DESCRIPTION','Index for the orbitals in active space.')
  call mh5_put_dset(dset_id,orbital_table%idx)
  call mh5_close_dset(dset_id)

  dset_id = mh5_create_dset_real(file_id,'ORBITAL_ENERGIES',1,[length(orbital_table)])
  call mh5_init_attr(dset_id,'DESCRIPTION','Energies of orbitals in active space.')
  call mh5_put_dset(dset_id,orbital_table%values)
  call mh5_close_dset(dset_id)

  dset_id = mh5_create_dset_int(file_id,'FOCK_INDEX',2,[2,length(fock_table)])
  call mh5_init_attr(dset_id,'DESCRIPTION','The index i, j for the Fock matrix elements <i| F |j>.')
  call mh5_put_dset(dset_id,fock_table%idx)
  call mh5_close_dset(dset_id)

  dset_id = mh5_create_dset_real(file_id,'FOCK_VALUES',1,[length(fock_table)])
  call mh5_init_attr(dset_id,'DESCRIPTION','The Fock matrix elements <i| F |j>.')
  call mh5_init_attr(dset_id,'CUTOFF',fock_table%cutoff)
  call mh5_put_dset(dset_id,fock_table%values)
  call mh5_close_dset(dset_id)

  dset_id = mh5_create_dset_int(file_id,'TWO_EL_INT_INDEX',2,[4,length(two_el_table)])
  call mh5_init_attr(dset_id,'DESCRIPTION','The index i, j, k, l for the two electron integrals <i j | 1/r_{12} | k l >.')
  call mh5_put_dset(dset_id,two_el_table%idx)
  call mh5_close_dset(dset_id)

  dset_id = mh5_create_dset_real(file_id,'TWO_EL_INT_VALUES',1,[length(two_el_table)])
  call mh5_init_attr(dset_id,'DESCRIPTION','The two electron integrals <i j | 1/r_{12} | k l >.')
  call mh5_init_attr(dset_id,'CUTOFF',two_el_table%cutoff)
  call mh5_put_dset(dset_id,two_el_table%values)
  call mh5_close_dset(dset_id)

  call mh5_close_file(file_id)

  call FastIO('STATUS')
# else
  unused_var(EMY)
  unused_var(path)
  unused_var(orbital_table)
  unused_var(fock_table)
  unused_var(two_el_table)
  unused_var(orbsym)
# endif

end subroutine dump_hdf5

<<<<<<< HEAD
!>  @brief
!>    Create fort.55 file in MRCC format
!>
!>  @author Jaafar Mehrez
!>
!>  @details
!>  Writes an FCIDUMP compatible with MRCC (fort.55).
!>  Unlike the standard ASCII FCIDUMP, this format uses:
!>    - Line 1: NMO NELEC
!>    - Line 2: space-separated ORBSYM
!>    - Line 3:  150000
!>  One-electron integrals are written from the Fock table only;
!>  orbital energies are omitted.
subroutine dump_fort55(path,EMY,orbital_table,fock_table,two_el_table,orbsym)

  use general_data, only: nActEl, nAsh, nSym
  use Symmetry_Info, only: SymLab, lIrrep, Symmetry_Info_Get
  use Constants, only: Zero
  use Definitions, only: u6
  use Index_Functions, only: iTri
  use fcidump_tables, only: cutoff_default, length

  character(len=*), intent(in) :: path
  real(kind=wp), intent(in) :: EMY
  type(OrbitalTable), intent(in) :: orbital_table
  type(FockTable), intent(in) :: fock_table
  type(TwoElIntTable), intent(in) :: two_el_table
  integer(kind=iwp), intent(in) :: orbsym(:)
  integer(kind=iwp) :: i, j, k, l, m, ij, kl, LuFCI, mrcc_orbsym(size(orbsym)), nmo, npair
  real(kind=wp) :: val
  real(kind=wp), allocatable :: eri_2d(:,:)
  integer(kind=iwp), external :: isFreeUnit
  character(len=3) :: label

  call Symmetry_Info_Get()

  do i=1,size(orbsym)
    label = adjustl(trim(lIrrep(orbsym(i)-1)))
    call to_lower(label)
    mrcc_orbsym(i) = molcas_irrep_to_mrcc(orbsym(i),trim(adjustl(SymLab)),label)
  end do

  nmo = sum(nAsh(:nSym))
  npair = nmo*(nmo+1)/2
  LuFCI = isFreeUnit(38)
  call molcas_open(LuFCI,path)

  write(LuFCI,'(I0,1X,I0)') nmo,nActEl

  if (size(mrcc_orbsym) > 0) then
    do i=1,size(mrcc_orbsym)-1
      write(LuFCI,'(I0,1X)',advance='NO') mrcc_orbsym(i)
    end do
    write(LuFCI,'(I0)') mrcc_orbsym(size(mrcc_orbsym))
  end if

  write(LuFCI,'(A)') ' 150000'

  ! 4-fold symmetry
  allocate(eri_2d(npair,npair))
  eri_2d = Zero

  do j=1,length(two_el_table)
    i = two_el_table%idx(1,j)
    k = two_el_table%idx(2,j)
    l = two_el_table%idx(3,j)
    m = two_el_table%idx(4,j)
    ij = iTri(i,k)
    kl = iTri(l,m)
    eri_2d(ij,kl) = two_el_table%values(j)
    eri_2d(kl,ij) = two_el_table%values(j)
  end do

  do i=1,nmo
    do j=1,i
      ij = iTri(i,j)
      do k=1,nmo
        do l=1,k
          kl = iTri(k,l)
          val = eri_2d(ij,kl)
          if (abs(val) > cutoff_default) then
            write(LuFCI,'(1X,E27.20,4I5)') val,i,j,k,l
          end if
        end do
      end do
    end do
  end do

  deallocate(eri_2d)

  do j=1,length(fock_table)
    write(LuFCI,'(1X,E27.20,4I5)') fock_table%values(j),(fock_table%idx(i,j),i=1,2),0,0
  end do

  write(LuFCI,'(1X,E27.20,4I5)') EMY,0,0,0,0

  close(LuFCI)

  call FastIO('STATUS')

end subroutine dump_fort55

!>  @brief
!>    Map irrep label to the MRCC numbering used in fort.55
!>
!>  @author Jaafar Mehrez
!>

function molcas_irrep_to_mrcc(molcas_irrep,group,label) result(mrcc_irrep)

  integer(kind=iwp), intent(in) :: molcas_irrep
  character(len=*), intent(in) :: group, label
  integer(kind=iwp) :: mrcc_irrep

  mrcc_irrep = molcas_irrep

  select case (group)
  case ('D2h')
    select case (label)
    case ('ag')
      mrcc_irrep = 1
    case ('b1g')
      mrcc_irrep = 2
    case ('b2g')
      mrcc_irrep = 3
    case ('b3g')
      mrcc_irrep = 4
    case ('au')
      mrcc_irrep = 5
    case ('b1u')
      mrcc_irrep = 6
    case ('b2u')
      mrcc_irrep = 7
    case ('b3u')
      mrcc_irrep = 8
    end select
  case ('C2v')
    select case (label)
    case ('a1')
      mrcc_irrep = 1
    case ('a2')
      mrcc_irrep = 2
    case ('b1')
      mrcc_irrep = 3
    case ('b2')
      mrcc_irrep = 4
    end select
  case ('C2h')
    select case (label)
    case ('ag')
      mrcc_irrep = 1
    case ('bg')
      mrcc_irrep = 2
    case ('au')
      mrcc_irrep = 3
    case ('bu')
      mrcc_irrep = 4
    end select
  case ('D2')
    select case (label)
    case ('a')
      mrcc_irrep = 1
    case ('b1')
      mrcc_irrep = 2
    case ('b2')
      mrcc_irrep = 3
    case ('b3')
      mrcc_irrep = 4
    end select
  case ('Cs')
    select case (label)
    case ("a'")
      mrcc_irrep = 1
    case ('a"')
      mrcc_irrep = 2
    end select
  case ('C2')
    select case (label)
    case ('a')
      mrcc_irrep = 1
    case ('b')
      mrcc_irrep = 2
    end select
  case ('Ci')
    select case (label)
    case ('ag')
      mrcc_irrep = 1
    case ('au')
      mrcc_irrep = 2
    end select
  case ('C1')
    mrcc_irrep = 1
  end select

end function molcas_irrep_to_mrcc

!>  @brief
!>    Convert a character string to lowercase in-place.
!>
!>  @author Jaafar Mehrez
!>
subroutine to_lower(s)

  character(len=*), intent(inout) :: s
  integer(kind=iwp) :: i, ic

  do i=1,len(s)
    ic = ichar(s(i:i))
    if ((ic >= ichar('A')) .and. (ic <= ichar('Z'))) then
      s(i:i) = char(ic+ichar('a')-ichar('A'))
    end if
  end do

end subroutine to_lower

=======
>>>>>>> upstream-openmolcas/master
end module fcidump_dump
