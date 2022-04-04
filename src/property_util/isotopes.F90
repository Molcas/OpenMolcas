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
! Copyright (C) 2017, Ignacio Fdez. Galvan                             *
!***********************************************************************
!
! Isotope numbers and masses
!
! Each item in the the array ElementList contains data for an element:
!  - %Symbol: symbol
!  - %Natural: number of natural ocurring isotopes
!  - %Isotopes: array with all the isotopes for the element, sorted
!               by order of abundance (most to least, artificial
!               sorted by increasing mass number); elements with no
!               natural isotopes have the most stable first. Each item
!               in this array contains:
!    - %A: mass number (protons + neutrons)
!    - %m: isotopic mass in Da
!
! The "default" isotope for each element is simply the first item in
! the %Isotopes member.

module Isotopes

use stdalloc, only: mma_Allocate, mma_Deallocate
use Definitions, only: wp, iwp, u6

implicit none
private
type Iso_t
  integer(kind=iwp) :: A
  real(kind=wp) :: m
end type Iso_t
type Element_t
  character(len=2) :: Symbol
  integer(kind=iwp) :: Natural
  type(Iso_t), allocatable :: Isotopes(:)
end type Element_t
integer(kind=iwp), parameter :: MaxAtomNum = 118
type(Element_t), allocatable :: ElementList(:)
character(len=2), parameter :: PTab(0:MaxAtomNum) = [' X', &
                                                     ' H','He','Li','Be',' B',' C',' N',' O',' F','Ne', &
                                                     'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca', &
                                                     'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
                                                     'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr', &
                                                     'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn', &
                                                     'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd', &
                                                     'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb', &
                                                     'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg', &
                                                     'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th', &
                                                     'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm', &
                                                     'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
                                                     'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og' &
                                                    ]
#include "constants2.fh"

interface Isotope
  module procedure Isotope_sym, Isotope_num
end interface Isotope

protected :: ElementList
public :: MaxAtomNum, Isotope, ElementList, Initialize_Isotopes, Free_Isotopes, NuclideMass, PTab

! Private extensions to mma interfaces

interface cptr2loff
  module procedure elm_cptr2loff
  module procedure iso_cptr2loff
end interface
interface mma_Allocate
  module procedure element_mma_allo_1D, element_mma_allo_1D_lim
  module procedure isotope_mma_allo_1D, isotope_mma_allo_1D_lim
end interface
interface mma_Deallocate
  module procedure element_mma_free_1D
  module procedure isotope_mma_free_1D
end interface

contains

! This subroutine allocates and fills the data in ElementList
! Since each array has a different size, it has to be done dynamically

subroutine Initialize_Isotopes()
  use Constants, only: Zero, One
  use Definitions, only: u6
  integer(kind=iwp) :: Err, i, Lu_iso, Most, NumElem, NumIso, NumNat, n1, n2
  logical(kind=iwp) :: Found
  character(len=180) :: Line
  character(len=20) :: Word
  real(kind=wp), allocatable :: Swap(:), Tab(:,:), Tmp(:,:)
  character(len=*), parameter :: ISODATA_NAME = 'ISODATA'
  integer(kind=iwp), external :: IsFreeUnit
# ifdef _GARBLE_
  interface
    subroutine c_null_alloc(A)
      import :: Iso_t
      type(Iso_t), allocatable :: A(:)
    end subroutine c_null_alloc
  end interface
# endif

  if (allocated(ElementList)) return
  call mma_Allocate(ElementList,MaxAtomNum,'ElmList')
# ifdef _GARBLE_
  ! Garbling corrupts the allocation status of allocatable components, use a hack to reset it
  do NumElem=1,size(ElementList,1)
    call c_null_alloc(ElementList(NumElem)%Isotopes)
  end do
# endif

# include "macros.fh"
  unused_proc(mma_allocate(ElementList,[0,0]))
  unused_proc(mma_allocate(ElementList(0)%Isotopes,[0,0]))

  call f_Inquire(ISODATA_NAME,Found)
  if (.not.Found) then
    write(u6,*) 'Isotope: The isotope data file does not exist'
    call abend()
  end if
  Lu_iso = IsFreeUnit(10)
  call molcas_open(Lu_iso,ISODATA_NAME)

  ! The data is read from the ISODATA file (see that file for source)

  ! guess for max isotopes per element: 10
  ! but don't worry, it will be increased when needed
  call mma_allocate(Tab,3,10,label='Tab')
  call mma_allocate(Swap,3,label='Swap')

  NumElem = 0
  NumIso = 0
  NumNat = 0
  Most = 0
  do
    read(Lu_iso,'(A)',iostat=Err) Line
    if (Err < 0) then
      exit
    else if (Err > 0) then
      write(u6,*) 'Isotope: Error reading the isotope data file'
      call abend()
    end if
    if ((Line(1:1) == '#') .or. (Line == '')) cycle
    ! A new element starts with the atomic number in columns 1:3
    read(Line(1:3),*,iostat=Err) i
    if (Err == 0) then
      if ((i < lbound(ElementList,1)) .or. (i > ubound(ElementList,1))) then
        write(u6,*) 'Isotope: Error reading the isotope data file'
        write(u6,*) 'Num: ',i
        call abend()
      end if
      NumElem = i
      if (allocated(ElementList(NumElem)%Isotopes)) then
        write(u6,*) 'Isotope: Error reading the isotope data file'
        write(u6,*) 'Duplicate: ',NumElem
        call abend()
      end if
      ! Columns 47:51 may contain the most stable isotope between brackets
      Word = Line(47:51)
      n1 = index(Word,'[')
      n2 = index(Word,']')
      if ((n1 > 0) .and. (n2 > 0)) then
        read(Word(n1+1:n2-1),*,iostat=Err) Most
        if (Err /= 0) then
          write(u6,*) 'Isotope: Error reading the isotope data file'
          write(u6,*) 'Word: ',Word
          call abend()
        end if
      else
        Most = 0
      end if
      NumIso = 0
      NumNat = 0
      Tab(:,:) = Zero
    else if (NumElem == 0) then
      cycle
    end if
    ! When all the isotopes for an element have been read, store them
    if (Line(1:1) == '_') then
      ElementList(NumElem)%Symbol = adjustl(PTab(NumElem))
      ElementList(NumElem)%Natural = NumNat
      call mma_Allocate(ElementList(NumElem)%Isotopes,NumIso)
      do i=1,NumIso
        ElementList(NumElem)%Isotopes(i) = Iso_t(nint(Tab(1,i)),Tab(2,i))
      end do
      cycle
    end if
    ! Now reading the data for each isotope
    NumIso = NumIso+1
    ! Increase the Tab size if needed
    if (NumIso > size(Tab,2)) then
      call mma_allocate(Tmp,size(Tab,1),2*size(Tab,2),label='Tmp')
      Tmp(:,:size(Tab,2)) = Tab
      Tmp(:,size(Tab,2)+1) = Zero
      call mma_deallocate(Tab)
      call move_alloc(Tmp,Tab)
    end if
    ! Columns 9:11 contain the mass number
    read(Line(9:11),*,iostat=Err) i
    if (Err /= 0) then
      write(u6,*) 'Isotope: Error reading the isotope data file'
      write(u6,*) 'Num: ',i
      call abend()
    end if
    Tab(1,NumIso) = real(i,kind=wp)
    ! Columns 14:31 contain the mass (with uncertainty)
    ! NB: nint(%m) should equal to %A, but you never know...
    Word = Line(14:31)
    n1 = index(Word,'(')
    if (n1 > 0) Word = Word(1:n1-1)
    read(Word,*,iostat=Err) Tab(2,NumIso)
    if (Err /= 0) then
      write(u6,*) 'Isotope: Error reading the isotope data file'
      write(u6,*) 'Word: ',Word
      call abend()
    end if
    ! Columns 33:46 may contain the abundance (with uncertainty)
    Word = Line(33:46)
    if (Word == '') then
      if (nint(Tab(1,NumIso)) == Most) then
        ! assign a small non-zero abundance to the most stable isotope
        Tab(3,NumIso) = tiny(Tab)
      else
        Tab(3,NumIso) = -One
      end if
    else
      n1 = index(Word,'(')
      if (n1 > 0) Word = Word(1:n1-1)
      read(Word,*,iostat=Err) Tab(3,NumIso)
      if (Err /= 0) then
        write(u6,*) 'Isotope: Error reading the isotope data file'
        write(u6,*) 'Word: ',Word
        call abend()
      end if
      if (Tab(3,NumIso) > Zero) NumNat = NumNat+1
    end if
    ! sort this isotope by order of abundance
    do i=1,NumIso-1
      if (Tab(3,NumIso) > Tab(3,i)) then
        Swap(:) = Tab(:,NumIso)
        Tab(:,i+1:NumIso) = Tab(:,i:NumIso-1)
        Tab(:,i) = Swap
        exit
      end if
    end do
  end do

  call mma_deallocate(Tab)
  call mma_deallocate(Swap)

  close(Lu_iso)

end subroutine Initialize_Isotopes

! This subroutine frees up the memory

subroutine Free_Isotopes()
  integer(kind=iwp) :: i
  if (.not. allocated(ElementList)) return
  do i=1,size(ElementList,1)
    call mma_Deallocate(ElementList(i)%Isotopes)
  end do
  call mma_Deallocate(ElementList)
end subroutine Free_Isotopes

! Subroutine(s) to get the Mass of the isotope IsNr belonging to the
! element Atom. If IsNr=0, the most abundant isotope (or the most
! stable if all are radioactive) is selected. The mass is returned
! in atomic units (m_e).
! Atom can be an atomic symbol or an atomic number.

subroutine Isotope_sym(IsNr,Atom,Mass)
  integer(kind=iwp), intent(inout) :: IsNr
  character(len=2), intent(in) :: Atom
  real(kind=wp), intent(out) :: Mass
  integer(kind=iwp) :: i, This
  character(len=2) :: Sym, Sym2

  call Initialize_Isotopes()

  Sym2 = adjustl(Atom)
  call UpCase(Sym2)
  if ((Sym2 == 'D') .or. (Sym2 == 'T')) Sym2 = 'H'
  This = 0
  do i=1,MaxAtomNum
    Sym = adjustl(ElementList(i)%Symbol)
    call UpCase(Sym)
    if (Sym == Sym2) then
      This = i
      exit
    end if
  end do

  if (This == 0) then
    write(u6,*) 'Isotope: Did not find atom!'
    write(u6,*) 'Atom=',Atom
    call Abend()
  end if

  if (IsNr == 0) IsNr = ElementList(This)%Isotopes(1)%A
  if (Sym2 == 'D') IsNr = 2
  if (Sym2 == 'T') IsNr = 3
  do i=1,size(ElementList(This)%Isotopes,1)
    if (ElementList(This)%Isotopes(i)%A == IsNr) then
      Mass = uToau*ElementList(This)%Isotopes(i)%m
      return
    end if
  end do

  write(u6,*) 'Isotope: Did not find isotope!'
  write(u6,*) 'IsNr=',IsNr
  write(u6,*) 'Atom=',Atom
  call Abend()

end subroutine Isotope_sym

subroutine Isotope_num(IsNr,Atom,Mass)
  integer(kind=iwp), intent(inout) :: IsNr
  integer(kind=iwp), intent(in) :: Atom
  real(kind=wp), intent(out) :: Mass
  integer(kind=iwp) :: i

  call Initialize_Isotopes()

  if ((Atom < 0) .or. (Atom > MaxAtomNum)) then
    write(u6,*) 'Isotope: Did not find atom!'
    write(u6,*) 'Atom=',Atom
    call Abend()
  end if

  if (IsNr == 0) IsNr = ElementList(Atom)%Isotopes(1)%A
  do i=1,size(ElementList(Atom)%Isotopes,1)
    if (ElementList(Atom)%Isotopes(i)%A == IsNr) then
      Mass = uToau*ElementList(Atom)%Isotopes(i)%m
      return
    end if
  end do

  write(u6,*) 'Isotope: Did not find isotope!'
  write(u6,*) 'IsNr=',IsNr
  write(u6,*) 'Atom=',Atom
  call Abend()

end subroutine Isotope_num

! Function that returns the mass in atomic units (m_e) of a particular
! nuclide with Z protons and A-Z neutrons. Returns -1.0 if the nuclide
! is unknown.

function NuclideMass(Z,A)
  use Constants, only: One
  real(kind=wp) :: NuclideMass
  integer(kind=iwp), intent(in) :: Z, A
  integer(kind=iwp) :: i

  call Initialize_Isotopes()

  NuclideMass = -One
  if ((Z < 1) .or. (Z > MaxAtomNum)) return
  do i=1,size(ElementList(Z)%Isotopes,1)
    if (ElementList(Z)%Isotopes(i)%A == A) then
      NuclideMass = uToau*ElementList(Z)%Isotopes(i)%m
      exit
    end if
  end do

end function NuclideMass

! Private extensions to mma_interfaces, using preprocessor templates
! (see src/mma_util/stdalloc.f)

! Define elm_cptr2loff, element_mma_allo_1D, element_mma_allo_1D_lim, element_mma_free_1D
#define _TYPE_ type(element_t)
#  define _FUNC_NAME_ elm_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ element_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'elm_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

! Define iso_cptr2loff, isotope_mma_allo_1D, isotope_mma_allo_1D_lim, isotope_mma_free_1D
#define _TYPE_ type(iso_t)
#  define _FUNC_NAME_ iso_cptr2loff
#  include "cptr2loff_template.fh"
#  undef _FUNC_NAME_
#  define _SUBR_NAME_ isotope_mma
#  define _DIMENSIONS_ 1
#  define _DEF_LABEL_ 'iso_mma'
#  include "mma_allo_template.fh"
#  undef _SUBR_NAME_
#  undef _DIMENSIONS_
#  undef _DEF_LABEL_
#undef _TYPE_

end module Isotopes
