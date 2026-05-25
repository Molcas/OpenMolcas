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
  use Definitions,   only: iwp, wp, u6
  implicit none

  private

  type :: isotope
    integer(kind=iwp)             :: AtNumb
    integer(kind=iwp)             :: MassNum
    real(kind=wp)                 :: Abundance
    real(kind=wp)                 :: NucSpin
    real(kind=wp)                 :: NucGfac
    logical(iwp)                  :: Stable
  end type isotope

  type :: element
    integer(kind=iwp)             :: first_isotope
    integer(kind=iwp)             :: last_isotope
  end type element

  type(element),allocatable   :: elements(:)
  type(isotope),allocatable   :: isotopes(:)

  public :: init_spin_data, free_spin_data, GNUC_NUCSPIN_by_nucmass, &
            GNUC_by_nucspin, NUCSPIN_by_gnuc, get_first_nonzero_GNUC


  contains

  subroutine init_spin_data()
    logical(kind=iwp)   :: Found
    character(len=180)  :: Line
    character(len=*), parameter :: SPINDATA_NAME = 'SPINDATA'
    integer(kind=iwp) :: NumbElem, NumbIso, SpinData, AtNumb, iRec, FirstIso, LastIso
    integer(kind=iwp) :: Err, Err1,Err2,Err3,Err4,Err5
    integer(kind=iwp), external :: IsFreeUnit


    call f_Inquire(SPINDATA_NAME,Found)

    if (.not.Found) then
      write(u6,*) 'File data/spin_data.txt does not exist'
      call abend()
    end if

  SpinData = IsFreeUnit(79)
  call molcas_open(SpinData,SPINDATA_NAME)

    do
      read(SpinData,'(A)',iostat=Err) Line
      if (Err < 0) then
       exit
      else if(Err > 0) then
        write(u6,*) 'spin_data.F90: Error reading spin_data.txt'
        call abend()
      endif
      !----------------------------------------------
      if ((Line(1:1) == '!') .or. (Line == '')) cycle
      !----------------------------------------------
      if (Line(1:1) == '@') then
        read(Line(22:24), *, iostat=Err1) NumbElem
        read(Line(57:59), *, iostat=Err2) NumbIso

        if (Err1 > 0) write(u6,*) 'spin_data.F90: Error reading NumbElem'
        if (Err2 > 0) write(u6,*) 'spin_data.F90: Error reading NumbIso'
        if (Err1 + Err2 >0) call AbEnd()

        allocate(elements(NumbElem),stat=Err1)
        allocate(isotopes(NumbIso),stat=Err2)
        if (Err1 > 0) write(u6,*) 'spin_data.F90: Error allocating elements(:)'
        if (Err2 > 0) write(u6,*) 'spin_data.F90: Error allocating isotopes(:)'
        if (Err1 + Err2 >0) call AbEnd()
      endif
      !----------------------------------------------
      if (Line(1:1) == '>') then
        read(Line(22:24), *, iostat=Err1) AtNumb
        read(Line(41:43), *, iostat=Err2) FirstIso
        read(Line(47:49), *, iostat=Err3) LastIso

        if (Err1 > 0) write(u6,*) "spin_data.F90: Error reading AtNumb"
        if (Err2 > 0) write(u6,*) "spin_data.F90: Error reading FirstIso"
        if (Err3 > 0) write(u6,*) "spin_data.F90: Error reading LastIso"
        if ((Err1+Err2+Err3) > 0) call AbEnd()

        elements(AtNumb)%first_isotope = FirstIso
        elements(AtNumb)%last_isotope  = LastIso

        ! Skip header file
        read(SpinData,'(A)',iostat=Err) Line

        ! Reading all isotopes of this element
        do iRec = FirstIso, LastIso
          read(SpinData,'(A)',iostat=Err) Line
          read(Line(17:24), *, iostat=Err1) isotopes(iRec)%MassNum
          read(Line(28:39), *, iostat=Err2) isotopes(iRec)%Abundance
          read(Line(43:49), *, iostat=Err3) isotopes(iRec)%NucSpin
          read(Line(54:64), *, iostat=Err4) isotopes(iRec)%NucGfac
          read(Line(71:73), *, iostat=Err5) isotopes(iRec)%Stable

          if (Err  > 0) write(u6,*) "spin_data.F90: Error reading Line in spin_data.txt"
          if (Err1 > 0) write(u6,*) "spin_data.F90: Error reading MassNum"
          if (Err2 > 0) write(u6,*) "spin_data.F90: Error reading Abundance"
          if (Err3 > 0) write(u6,*) "spin_data.F90: Error reading NucSpin"
          if (Err4 > 0) write(u6,*) "spin_data.F90: Error reading NucGfac"
          if (Err5 > 0) write(u6,*) "spin_data.F90: Error reading Stable"
          if ((Err+Err1+Err2+Err3+Err4+Err5) > 0) call AbEnd()
        end do
      endif
      !----------------------------------------------
    enddo

    close(SpinData)
  end subroutine


!======================================================================
  function GNUC_by_nucspin(AtNumb,NucSpin) result(NucGfac)
    integer(kind=iwp),intent(in) :: AtNumb
    real(kind=wp), intent(in)    :: NucSpin
    real(kind=wp)                :: NucGfac
    integer(kind=iwp)            :: iRec, first_isotope, last_isotope
    logical(kind=iwp)            :: not_found_gfac


    first_isotope = elements(AtNumb)%first_isotope
    last_isotope  = elements(AtNumb)%last_isotope

    NucGfac = 0.0_wp
    not_found_gfac = .true.
    do iRec=first_isotope,last_isotope
      if (abs(isotopes(iRec)%NucSpin - NucSpin) < 1e-5_wp) then
        NucGfac = isotopes(iRec)%NucGfac
        exit
      endif
    enddo

    if (not_found_gfac) call warning_not_found(1)
    call isot_info(AtNumb, NucGfac, NucSpin)

  end function GNUC_by_nucspin
!======================================================================


!======================================================================
  function NUCSPIN_by_gnuc(AtNumb,NucGfac) result(NucSpin)
    integer(kind=iwp),intent(in) :: AtNumb
    real(kind=wp), intent(in)    :: NucGfac
    real(kind=wp)                :: NucSpin
    integer(kind=iwp)            :: iRec, first_isotope, last_isotope
    logical(kind=iwp)            :: not_found_spin

    first_isotope = elements(AtNumb)%first_isotope
    last_isotope  = elements(AtNumb)%last_isotope

    NucSpin = 0.5_wp
    not_found_spin = .true.
    do iRec=first_isotope,last_isotope
      if (abs(isotopes(iRec)%NucGfac - NucGfac) < 0.001_wp) then
        not_found_spin = .false.
        NucSpin = isotopes(iRec)%NucSpin
        exit
      endif
    enddo

    if (not_found_spin) call warning_not_found(0)
    call isot_info(AtNumb, NucGfac, NucSpin)

  end function NUCSPIN_by_gnuc
!======================================================================


!======================================================================
  subroutine GNUC_NUCSPIN_by_nucmass(AtNumb,MassNum,NucGfac,NucSpin)
    integer(kind=iwp),intent(in) :: AtNumb, MassNum
    real(kind=wp), intent(out)   :: NucSpin, NucGfac
    integer(kind=iwp)            :: iRec, first_isotope, last_isotope
    logical(kind=iwp)            :: not_found_gfac , not_found_spin

    first_isotope = elements(AtNumb)%first_isotope
    last_isotope  = elements(AtNumb)%last_isotope

    NucSpin = 0.5_wp
    NucGfac = 0.0_wp
    not_found_gfac = .true.
    not_found_spin = .true.
    do iRec=first_isotope,last_isotope
      if (isotopes(iRec)%MassNum == MassNum) then
        not_found_gfac = .false.
        not_found_spin = .false.
        NucSpin = isotopes(iRec)%NucSpin
        NucGfac = isotopes(iRec)%NucGfac
        exit
      endif
    enddo


    if (not_found_spin) call warning_not_found(0)
    if (not_found_gfac) call warning_not_found(1)
    call isot_info(AtNumb, NucGfac, NucSpin, MassNum)

  end subroutine GNUC_NUCSPIN_by_nucmass
!======================================================================


!======================================================================
  subroutine get_first_nonzero_GNUC(AtNumb, NucGfac,NucSpin)
    integer(kind=iwp),intent(in)           :: AtNumb
    real(kind=wp), intent(out)             :: NucGfac, NucSpin
    integer(kind=iwp)                      :: iRec, first_isotope, last_isotope
    logical(kind=iwp)                      :: not_found_gfac, not_found_spin

    first_isotope = elements(AtNumb)%first_isotope
    last_isotope  = elements(AtNumb)%last_isotope

    NucSpin = 0.5_wp
    NucGfac = 0.0_wp
    not_found_gfac = .true.
    not_found_spin = .true.
    do iRec=first_isotope,last_isotope
      if (isotopes(iRec)%NucGfac /= 0.0_wp) then
        not_found_gfac = .false.
        not_found_spin = .false.
        NucGfac = isotopes(iRec)%NucGfac
        NucSpin = isotopes(iRec)%NucSpin
        exit
      endif
    enddo

    if (not_found_spin) call warning_not_found(0)
    if (not_found_gfac) call warning_not_found(1)
    call isot_info(AtNumb, NucGfac, NucSpin)

  end subroutine get_first_nonzero_GNUC
!======================================================================


!======================================================================
  subroutine isot_info(AtNumb, NucGfac, NucSpin,MassNum)
    integer(kind=iwp), intent(in) :: AtNumb
    integer(kind=iwp), intent(in), optional :: MassNum
    real(kind=wp), intent(in)     :: NucGfac, NucSpin

    write(6,'(11X,A18)')    repeat('_', 18)
    write(6,'(11X,A18,I0)')    'INFO          Z = ', AtNumb
    if (present(MassNum)) then
      write(6,'(11X,A18,I0)')    '              A = ', MassNum
    endif
    write(6,'(11X,A18,F12.8)') '       g-factor = ', NucGfac
    write(6,'(11X,A18,F5.2)')  '        NucSpin = ', NucSpin
    write(6,*) ''
  end subroutine isot_info
!======================================================================


!======================================================================
  subroutine warning_not_found(option)
    integer(kind=iwp), intent(in) :: option

    if(option == 0) then
      write(6,'(11X,A37)')  'WARNING: Not found a matching record!'
      write(6,'(20X,A61)')  'Defaulting to NucSpin = 0.5. This value might be fictitious. '
    else if(option == 1) then
      write(6,'(11X,A37)')  'WARNING: Not found a matching record!'
      write(6,'(20X,A61)')  'Defaulting to g-factor = 0.0. This value might be fictitious.'
    endif
  end subroutine warning_not_found
!======================================================================


!======================================================================
  subroutine free_spin_data()
    integer :: istat

    deallocate(isotopes, stat=istat)
    if (istat > 0) then
      write(6,*) ' Error in deallocating EasySpin isotopes.'
    end if

    deallocate(elements, stat=istat)
    if (istat > 0) then
      write(6,*) ' Error in deallocating EasySpin elements.'
    end if

  end subroutine free_spin_data

end module spin_data
