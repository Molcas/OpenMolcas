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

#include "compiler_features.h"
#ifdef _EFP_

!#define _DEBUGPRINT_
subroutine DrvEFP(First)

use, intrinsic :: iso_c_binding, only: c_char, c_funloc, c_int, c_loc, c_ptr, c_size_t
use EFP_Module, only: Coor_Type, EFP_Coors, EFP_Instance, FRAG_Type, nEFP_FRAGMENTS
use EFP, only: EFP_Add_Fragment, EFP_Add_Potential, EFP_Create, EFP_Get_Frag_Atom_Count, EFP_Prepare, &
               EFP_Set_Electron_Density_Field_FN, EFP_Set_Frag_Coordinates
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: First
integer(kind=iwp) :: i, iFrag, iLast, j
integer(c_int) :: irc
integer(c_size_t) :: frag_idx
integer(c_size_t), target :: n_atoms
character(len=180) :: CurrDir, MolDir
character(kind=c_char,len=180) :: FragName, PATH
type(c_ptr) :: cptr1
integer(c_int), external :: Molcas_ELECTRON_DENSITY_FIELD_FN
#ifdef _DEBUGPRINT_
integer(c_int) :: do_gradient
type(efp_energy), target :: Energy
#endif

! Initate the EFP object,

if (First) then

  ! Create a new EFP object

  EFP_Instance = EFP_Create()

  ! Now add the potentials

# ifdef _DEBUGPRINT_
  write(u6,*) 'Initiation of EFP'
  write(u6,*) 'nEFP_fragments=',nEFP_fragments
# endif
  iFrag = 0
  outer: do i=1,nEFP_fragments

    ! Find a unique fragment potential.
    ! Procastinate.

    do j=i+1,nEFP_fragments

      ! Branch out if there is an entry later with the same potential
      if (FRAG_TYPE(i) == FRAG_TYPE(j)) cycle outer
    end do

    ! At this point we have a unique potential label in FRAG_TYPE(i)

    if (index(FRAG_TYPE(i),'_l') /= 0) then

      ! .efg file found in the library directory

      iLast = index(FRAG_TYPE(i),'_l')
      FragName = FRAG_TYPE(i)(1:iLast-1)
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
      FragName = FRAG_TYPE(i)(1:iLast-1)
      call GetEnvf('CurrDir',CurrDir)
      iLast = index(CurrDir,' ')
      Path = CurrDir(1:iLast-1)//'/'
    end if

    iLast = index(Path,' ')
    Path = Path(1:iLast-1)//FragName
    iLast = index(Path,' ')
    Path = Path(1:iLast-1)//'.efp'//char(0)

    irc = EFP_ADD_POTENTIAL(EFP_Instance,Path)
    if (irc /= 0) then
      write(u6,*) 'EFP potential file error.'
      write(u6,*) Path
      write(u6,*) 'Return code:',irc
      call Abend()
    end if

    ! Loop over all fragment again an add those that are of the
    ! current type

    iLast = index(FRAG_TYPE(i),' ')
    FragName = FRAG_TYPE(i)(1:iLast-1)
    iLast = index(FragNAME,' ')
    FragName = FragName(1:iLast-1)//char(0)
    do j=1,i
      if (FRAG_TYPE(j) /= FRAG_TYPE(i)) cycle
      irc = efp_add_fragment(EFP_Instance,FragName)
      if (irc /= 0) then
        write(u6,*) 'EFP_ADD_FRAGMET error.'
        write(u6,*) 'Return code:',irc
        call Abend()
      end if
      cptr1 = c_loc(EFP_COORS(1,j))
      irc = efp_set_frag_Coordinates(EFP_Instance,iFrag,Coor_type,cptr1)
      if (irc /= 0) then
        write(u6,*) 'EFP_SET_FRAG_COORDINATES error.'
        write(u6,*) 'Return code:',irc
        call Abend()
      end if
      iFrag = iFrag+1
    end do

  end do outer

  irc = EFP_PREPARE(EFP_Instance)
  if (irc /= 0) then
    write(u6,*) 'EFP_PREPARE error.'
    write(u6,*) 'Return code:',irc
    call Abend()
  end if

# ifdef _DEBUGPRINT_
  do_gradient = 0
  irc = EFP_COMPUTE(EFP_Instance,do_gradient)
  if (irc /= 0) then
    write(u6,*) 'EFP_COMPUTE error.'
    write(u6,*) 'Return code:',irc
    call Abend()
  end if

  irc = EFP_GET_ENERGY(EFP_Instance,c_loc(Energy))
  write(u6,*) Energy%Total
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

return

end subroutine DrvEFP

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(DrvEFP)

#endif
