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
! Copyright (C) 2026, Dong Q. Le                                       *
!***********************************************************************
! Data derived from EasySpin (https://easyspin.org/)                   *
! Copyright (c) 2026                                                   *
! Licensed under the MIT License.                                      *
! Please see the file readme_easyspin_license.md in the                *
! Tools/create_easy_spin_data directory for licence details.           *
!***********************************************************************

module spin_data

use Constants, only: Zero, Half
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: iwp, wp, u6

implicit none
private

type :: isotope
  integer(kind=iwp) :: AtNumb, MassNumb
  real(kind=wp) :: Abundance, NucSpin, GNuc
  logical(kind=iwp) :: Stable
end type isotope

type :: element
  integer(kind=iwp) :: first_isotope, last_isotope
end type element

integer(kind=iwp) :: first_isotope, last_isotope
logical(kind=iwp) :: not_found_gnuc, not_found_spin

type(element), allocatable :: elements(:)
type(isotope), allocatable :: isotopes(:)

public :: free_spin_data, get_first_nonzero_GNUC, GNUC_by_nucspin, GNUC_NUCSPIN_by_nucmass, init_spin_data, NUCSPIN_by_gnuc

! Private extensions to mma interfaces

interface mma_allocate
  module procedure :: element_mma_allo_1D, element_mma_allo_1D_lim
  module procedure :: isotope_mma_allo_1D, isotope_mma_allo_1D_lim
end interface
interface mma_deallocate
  module procedure :: element_mma_free_1D
  module procedure :: isotope_mma_free_1D
end interface

contains

subroutine init_spin_data()

  integer(kind=iwp) :: AtNumb, Err, Err1, Err2, Err3, Err4, Err5, FirstIso, iRec, LastIso, NumbElem, NumbIso, SpinData
  logical(kind=iwp) :: Found
  character(len=180) :: Line
  character(len=*), parameter :: SPINDATA_NAME = 'SPINDATA'
  integer(kind=iwp), external :: IsFreeUnit

  call f_Inquire(SPINDATA_NAME,Found)

  if (.not. Found) then
    write(u6,*) 'File data/spin_data.txt does not exist'
    call abend()
  end if

  SpinData = IsFreeUnit(79)
  call molcas_open(SpinData,SPINDATA_NAME)

  do
    read(SpinData,'(A)',iostat=Err) Line
    if (Err < 0) then
      exit
    else if (Err > 0) then
      write(u6,*) 'spin_data.F90: Error reading spin_data.txt'
      call abend()
    end if
    !----------------------------------------------
    if ((Line(1:1) == '!') .or. (Line == '')) cycle
    !----------------------------------------------
    if (Line(1:1) == '@') then
      read(Line(22:24),*,iostat=Err1) NumbElem
      read(Line(57:59),*,iostat=Err2) NumbIso

      if (Err1 > 0) write(u6,*) 'spin_data.F90: Error reading NumbElem'
      if (Err2 > 0) write(u6,*) 'spin_data.F90: Error reading NumbIso'
      if (Err1+Err2 > 0) call AbEnd()

      call mma_allocate(elements,NumbElem,Label='elements')
      call mma_allocate(isotopes,NumbIso,Label='isotopes')
#     include "macros.fh"
      unused_proc(mma_allocate(elements,[0,0]))
      unused_proc(mma_allocate(isotopes,[0,0]))
    end if
    !----------------------------------------------
    if (Line(1:1) == '>') then
      read(Line(22:24),*,iostat=Err1) AtNumb
      read(Line(41:43),*,iostat=Err2) FirstIso
      read(Line(47:49),*,iostat=Err3) LastIso

      if (Err1 > 0) write(u6,*) "spin_data.F90: Error reading AtNumb"
      if (Err2 > 0) write(u6,*) "spin_data.F90: Error reading FirstIso"
      if (Err3 > 0) write(u6,*) "spin_data.F90: Error reading LastIso"
      if ((Err1+Err2+Err3) > 0) call AbEnd()

      elements(AtNumb)%first_isotope = FirstIso
      elements(AtNumb)%last_isotope = LastIso

      ! Skip header file
      read(SpinData,'(A)',iostat=Err) Line

      ! Reading all isotopes of this element
      do iRec=FirstIso,LastIso
        read(SpinData,'(A)',iostat=Err) Line
        read(Line(17:24),*,iostat=Err1) isotopes(iRec)%MassNumb
        read(Line(28:39),*,iostat=Err2) isotopes(iRec)%Abundance
        read(Line(43:49),*,iostat=Err3) isotopes(iRec)%NucSpin
        read(Line(54:64),*,iostat=Err4) isotopes(iRec)%GNuc
        read(Line(71:73),*,iostat=Err5) isotopes(iRec)%Stable

        if (Err > 0) write(u6,*) "spin_data.F90: Error reading Line in spin_data.txt"
        if (Err1 > 0) write(u6,*) "spin_data.F90: Error reading MassNum"
        if (Err2 > 0) write(u6,*) "spin_data.F90: Error reading Abundance"
        if (Err3 > 0) write(u6,*) "spin_data.F90: Error reading NucSpin"
        if (Err4 > 0) write(u6,*) "spin_data.F90: Error reading GNuc"
        if (Err5 > 0) write(u6,*) "spin_data.F90: Error reading Stable"
        if ((Err+Err1+Err2+Err3+Err4+Err5) > 0) call AbEnd()
      end do
    end if
    !----------------------------------------------
  end do

  close(SpinData)

end subroutine

subroutine GNUC_by_nucspin(AtNumb,MassNumb,GNuc,NucSpin,Stability)

  integer(kind=iwp), intent(in) :: AtNumb
  integer(kind=iwp), intent(out) :: MassNumb
  real(kind=wp), intent(out) :: GNuc
  real(kind=wp), intent(in) :: NucSpin
  character, intent(out), optional :: Stability
  integer(kind=iwp) :: iRec

  first_isotope = elements(AtNumb)%first_isotope
  last_isotope = elements(AtNumb)%last_isotope

  GNuc = Zero
  not_found_gnuc = .true.
  do iRec=first_isotope,last_isotope
    if (abs(isotopes(iRec)%NucSpin-NucSpin) < 1.0e-5_wp) then
      not_found_gnuc = .false.
      GNuc = isotopes(iRec)%GNuc
      MassNumb = isotopes(iRec)%MassNumb

      if (present(Stability)) then
        if (isotopes(iRec)%Stable) then
          Stability = ' '
        else
          Stability = '*'
        end if
      end if

      exit
    end if
  end do

  if (not_found_gnuc) call warning_not_found(1)
  call isot_info(AtNumb,GNuc,NucSpin)

end subroutine GNUC_by_nucspin

subroutine NUCSPIN_by_gnuc(AtNumb,MassNumb,GNuc,NucSpin,Stability)

  integer(kind=iwp), intent(in) :: AtNumb
  integer(kind=iwp), intent(out) :: MassNumb
  real(kind=wp), intent(in) :: GNuc
  real(kind=wp), intent(out) :: NucSpin
  character, intent(out), optional :: Stability
  integer(kind=iwp) :: iRec

  first_isotope = elements(AtNumb)%first_isotope
  last_isotope = elements(AtNumb)%last_isotope

  NucSpin = Half
  not_found_spin = .true.
  do iRec=first_isotope,last_isotope
    if (abs(isotopes(iRec)%GNuc-GNuc) < 0.001_wp) then
      not_found_spin = .false.
      NucSpin = isotopes(iRec)%NucSpin
      MassNumb = isotopes(iRec)%MassNumb

      if (present(Stability)) then
        if (isotopes(iRec)%Stable) then
          Stability = ' '
        else
          Stability = '*'
        end if
      end if

      exit
    end if
  end do

  if (not_found_spin) call warning_not_found(0)
  call isot_info(AtNumb,GNuc,NucSpin)

end subroutine NUCSPIN_by_gnuc

subroutine GNUC_NUCSPIN_by_nucmass(AtNumb,MassNumb,GNuc,NucSpin,Stability)

  integer(kind=iwp), intent(in) :: AtNumb, MassNumb
  real(kind=wp), intent(out) :: GNUc, NucSpin
  character, intent(out), optional :: Stability
  integer(kind=iwp) :: iRec

  first_isotope = elements(AtNumb)%first_isotope
  last_isotope = elements(AtNumb)%last_isotope

  NucSpin = Half
  GNuc = Zero
  not_found_gnuc = .true.
  not_found_spin = .true.
  do iRec=first_isotope,last_isotope
    if (isotopes(iRec)%MassNumb == MassNumb) then
      not_found_gnuc = .false.
      not_found_spin = .false.
      NucSpin = isotopes(iRec)%NucSpin
      GNuc = isotopes(iRec)%GNuc

      if (present(Stability)) then
        if (isotopes(iRec)%Stable) then
          Stability = ' '
        else
          Stability = '*'
        end if
      end if

      exit
    end if
  end do

  if (not_found_spin) call warning_not_found(0)
  if (not_found_gnuc) call warning_not_found(1)
  call isot_info(AtNumb,GNuc,NucSpin,MassNumb)

end subroutine GNUC_NUCSPIN_by_nucmass

subroutine get_first_nonzero_GNUC(AtNumb,MassNumb,GNuc,NucSpin,Stability)

  integer(kind=iwp), intent(in) :: AtNumb
  integer(kind=iwp), intent(out) :: MassNumb
  real(kind=wp), intent(out) :: GNuc, NucSpin
  character, intent(out), optional :: Stability
  integer(kind=iwp) :: iRec

  first_isotope = elements(AtNumb)%first_isotope
  last_isotope = elements(AtNumb)%last_isotope

  NucSpin = Half
  GNuc = Zero
  not_found_gnuc = .true.
  not_found_spin = .true.
  do iRec=first_isotope,last_isotope
    if (isotopes(iRec)%GNuc /= Zero) then
      not_found_gnuc = .false.
      not_found_spin = .false.
      GNuc = isotopes(iRec)%GNuc
      NucSpin = isotopes(iRec)%NucSpin
      MassNumb = isotopes(iRec)%MassNumb

      if (present(Stability)) then
        if (isotopes(iRec)%Stable) then
          Stability = ' '
        else
          Stability = '*'
        end if
      end if

      exit
    end if
  end do

  if (not_found_spin) call warning_not_found(0)
  if (not_found_gnuc) call warning_not_found(1)
  call isot_info(AtNumb,GNuc,NucSpin)

end subroutine get_first_nonzero_GNUC

subroutine isot_info(AtNumb,GNuc,NucSpin,MassNumb)

  integer(kind=iwp), intent(in) :: AtNumb
  real(kind=wp), intent(in) :: GNuc, NucSpin
  integer(kind=iwp), intent(in), optional :: MassNumb

  write(u6,'(11X,A18)') repeat('_',18)
  write(u6,'(11X,A18,I0)') 'INFO          Z = ',AtNumb
  if (present(MassNumb)) then
    write(u6,'(11X,A18,I0)') '              A = ',MassNumb
  end if
  write(u6,'(11X,A18,F12.8)') '       g-factor = ',GNuc
  write(u6,'(11X,A18,F5.2)') '        NucSpin = ',NucSpin
  write(u6,*) ''

end subroutine isot_info

subroutine warning_not_found(option)

  integer(kind=iwp), intent(in) :: option

  if (option == 0) then
    write(u6,'(11X,A37)') 'WARNING: Not found a matching record!'
    write(u6,'(20X,A61)') 'Defaulting to NucSpin = 0.5. This value might be fictitious. '
  else if (option == 1) then
    write(u6,'(11X,A37)') 'WARNING: Not found a matching record!'
    write(u6,'(20X,A61)') 'Defaulting to g-factor = 0.0. This value might be fictitious.'
  end if

end subroutine warning_not_found

subroutine free_spin_data()

  call mma_deallocate(isotopes)
  call mma_deallocate(elements)

end subroutine free_spin_data

! Private extensions to mma_interfaces, using preprocessor templates
! (see mma_util/stdalloc.F90)

! Define element_mma_allo_1D, element_mma_allo_1D_lim, element_mma_free_1D
#define _TYPE_ type(element)
#  define _SUBR_NAME_ element_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'elm_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

! Define isotope_mma_allo_1D, isotope_mma_allo_1D_lim, isotope_mma_free_1D
#define _TYPE_ type(isotope)
#  define _SUBR_NAME_ isotope_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'iso_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module spin_data
