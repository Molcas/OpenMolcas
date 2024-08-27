!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine DrvEFP(First)

#ifdef _EFP_
use EFP_Module, only: EFP_Instance, nEFP_FRAGMENTS, Coor_Type, FRAG_Type, EFP_Coors
use EFP, only: EFP_Add_Fragment, EFP_Add_Potential, EFP_Create, EFP_Get_Frag_Atom_Count, EFP_Prepare, &
               EFP_Set_Electron_Density_Field_FN, EFP_Set_Frag_Coordinates
use iso_c_binding, only: c_int, c_char, c_ptr, c_size_t, c_loc, c_funloc

implicit none
logical First
external Molcas_ELECTRON_DENSITY_FIELD_FN
character(len=180) :: CurrDir, MolDir
type(c_ptr) :: cptr1
integer(c_int) :: irc
character(kind=c_char) :: Name*180, PATH*180
integer(c_int) :: Molcas_ELECTRON_DENSITY_FIELD_FN
integer(c_size_t) :: frag_idx
integer(c_size_t), target :: n_atoms
integer :: iFrag, i, j, iLast
#ifdef _DEBUGPRINT_
type(efp_energy), target :: Energy
integer(c_int) :: do_gradient
#endif

! Initate the EFP object,

if (First) then

  ! Create a new EFP object

  EFP_Instance = EFP_Create()

  ! Now add the potentials

# ifdef _DEBUGPRINT_
  write(6,*) 'Initiation of EFP'
  write(6,*) 'nEFP_fragments=',nEFP_fragments
# endif
  iFrag = 0
  do i=1,nEFP_fragments

    ! Find a unique fragment potential.
    ! Procastinate.

    do j=i+1,nEFP_fragments

      ! Branch out if there is an entry later with the same potential
      if (FRAG_TYPE(i) == FRAG_TYPE(j)) Go To 999
    end do

    ! At this point we have a unique potential label in FRAG_TYPE(i)

    if (index(FRAG_TYPE(i),'_l') /= 0) then

      ! .efg file found in the library directory

      iLast = index(FRAG_TYPE(i),'_l')
      Name = FRAG_TYPE(i)(1:iLast-1)
      call GetEnvf('MOLCAS',MolDir)
      iLast = index(MolDir,' ')
      if (MolDir(iLast-1:iLast-1) /= '/') then
        MolDir(iLast:iLast) = '/'
        iLast = iLast+1
      end if
      Path = MolDir(1:iLast-1)//'External/efp/fraglib/'
    else

      ! .efg file found in the $CurrDir directory

      iLast = index(FRAG_TYPE(i),' ')
      Name = FRAG_TYPE(i)(1:iLast-1)
      call GetEnvf('CurrDir',CurrDir)
      iLast = index(CurrDir,' ')
      Path = CurrDir(1:iLast-1)//'/'
    end if

    iLast = index(Path,' ')
    Path = Path(1:iLast-1)//Name
    iLast = index(Path,' ')
    Path = Path(1:iLast-1)//'.efp'//char(0)

    irc = EFP_ADD_POTENTIAL(EFP_Instance,Path)
    if (irc /= 0) then
      write(6,*) 'EFP potential file error.'
      write(6,*) Path
      write(6,*) 'Return code:',irc
      call Abend()
    end if

    ! Loop over all fragment again an add those that are of the
    ! current type

    iLast = index(FRAG_TYPE(i),' ')
    Name = FRAG_TYPE(i)(1:iLast-1)
    iLast = index(NAME,' ')
    NAME = NAME(1:iLast-1)//char(0)
    do j=1,i
      if (FRAG_TYPE(j) /= FRAG_TYPE(i)) cycle
      irc = efp_add_fragment(EFP_Instance,Name)
      if (irc /= 0) then
        write(6,*) 'EFP_ADD_FRAGMET error.'
        write(6,*) 'Return code:',irc
        call Abend()
      end if
      cptr1 = c_loc(EFP_COORS(1,j))
      irc = efp_set_frag_Coordinates(EFP_Instance,iFrag,Coor_type,cptr1)
      if (irc /= 0) then
        write(6,*) 'EFP_SET_FRAG_COORDINATES error.'
        write(6,*) 'Return code:',irc
        call Abend()
      end if
      iFrag = iFrag+1
    end do

999 continue
  end do

  irc = EFP_PREPARE(EFP_Instance)
  if (irc /= 0) then
    write(6,*) 'EFP_PREPARE error.'
    write(6,*) 'Return code:',irc
    call Abend()
  end if

# ifdef _DEBUGPRINT_
  do_gradient = 0
  irc = EFP_COMPUTE(EFP_Instance,do_gradient)
  if (irc /= 0) then
    write(6,*) 'EFP_COMPUTE error.'
    write(6,*) 'Return code:',irc
    call Abend()
  end if

  irc = EFP_GET_ENERGY(EFP_Instance,c_loc(Energy))
  write(6,*) Energy%Total
# endif

  irc = EFP_SET_ELECTRON_DENSITY_FIELD_FN(EFP_Instance,c_funloc(Molcas_ELECTRON_DENSITY_FIELD_FN))

end if

! Add EFP charges to the nuclear repulsion term

do frag_idx=1,nEFP_fragments

  ! Pick up the number of atoms in the fragment

  irc = EFP_GET_FRAG_ATOM_COUNT(EFP_Instance,frag_idx,c_loc(n_atoms))

  ! Pick up the fragment coordinates and charges

  !irc = EFP_GET_FRAG_ATOMS(EFP_Instance,...

end do
#else
! Dummy routine
logical First
if (First .or. (.not. First)) return
#endif

return

end subroutine DrvEFP
