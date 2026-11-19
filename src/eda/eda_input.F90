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
! Copyright (C) 2026, Okan Koeksal                                     *
!***********************************************************************

subroutine EDA_Read_Input(iReturn)

use Definitions, only: iwp, u5, u6
use EDA_Data, only: EDA_nFrag, EDA_nAtom, EDA_Charge, EDA_Mult, &
                    EDA_RefMethod, EDA_FragBasis, EDA_Debug, EDA_Reset_Input, &
                    EDA_Allocate_Fragments
use Para_Info, only: nProcs, MyRank, Is_Real_Par

implicit none
#include "warnings.h"

integer(kind=iwp), intent(out) :: iReturn
integer(kind=iwp) :: i, AllocError
character(len=180) :: Line, Value
character(len=180), external :: Get_Ln
logical :: IsRoot

IsRoot = (nProcs <= 1) .or. (.not. Is_Real_Par()) .or. (MyRank == 0)

call EDA_Reset_Input()
iReturn = _RC_INPUT_ERROR_

call RdNLst(u5,'EDA')

read_loop: do
  Line = Get_Ln(u5)
  Value = adjustl(Line)
  call UpCase(Value)

  if ((Value(1:1) == ' ') .or. (Value(1:1) == '*')) cycle

  select case (Value(1:4))
    case ('NFRA')
      Line = Get_Ln(u5)
      call Get_I1(1,EDA_nFrag)
      if (EDA_nFrag < 1) then
        if (IsRoot) write(u6,'(6X,A,I0,A)') 'EDA: NFRAG must be positive; received ',EDA_nFrag,'.'
        return
      end if
      call EDA_Allocate_Fragments(EDA_nFrag,AllocError)
      if (AllocError /= 0) then
        if (IsRoot) write(u6,'(6X,A,I0,A)') 'EDA: unable to allocate fragment data for NFRAG=',EDA_nFrag,'.'
        return
      end if

    case ('NATO')
      if ((EDA_nFrag < 1) .or. (.not. allocated(EDA_nAtom))) then
        if (IsRoot) write(u6,'(6X,A)') 'EDA: NFRAG must appear before NATOM.'
        return
      end if
      Line = Get_Ln(u5)
      call Get_I(1,EDA_nAtom,EDA_nFrag)

    case ('CHAR')
      if ((EDA_nFrag < 1) .or. (.not. allocated(EDA_Charge))) then
        if (IsRoot) write(u6,'(6X,A)') 'EDA: NFRAG must appear before CHARGE.'
        return
      end if
      Line = Get_Ln(u5)
      call Get_I(1,EDA_Charge,EDA_nFrag)

    case ('MULT')
      if ((EDA_nFrag < 1) .or. (.not. allocated(EDA_Mult))) then
        if (IsRoot) write(u6,'(6X,A)') 'EDA: NFRAG must appear before MULT.'
        return
      end if
      Line = Get_Ln(u5)
      call Get_I(1,EDA_Mult,EDA_nFrag)

    case ('REFE')
      Line = Get_Ln(u5)
      EDA_RefMethod = '        '
      EDA_RefMethod(1:min(8,len_trim(Line))) = adjustl(Line(1:min(8,len_trim(Line))))
      call UpCase(EDA_RefMethod)

    case ('FBAS')
      ! Fragment-basis convention:
      !   OWN   : each fragment RunFile was generated in its own AO basis.
      !           Fragment AO blocks must be contiguous in the supermolecule.
      !   SUPER : each fragment RunFile was generated in the full supermolecule
      !           basis, e.g. with GATEWAY/BSSE ghost fragments.
      Line = Get_Ln(u5)
      EDA_FragBasis = '        '
      EDA_FragBasis(1:min(8,len_trim(Line))) = adjustl(Line(1:min(8,len_trim(Line))))
      call UpCase(EDA_FragBasis)

    case ('DEBU')
      EDA_Debug = .true.

    case ('END ')
      exit read_loop

    case default
      if (IsRoot) write(u6,'(6X,2A)') 'EDA: unrecognized keyword: ',trim(Value)
      return
  end select
end do read_loop

if ((EDA_nFrag < 1) .or. (.not. allocated(EDA_nAtom)) .or. &
    (.not. allocated(EDA_Charge)) .or. (.not. allocated(EDA_Mult))) then
  if (IsRoot) write(u6,'(6X,A)') 'EDA: NFRAG is required.'
  return
end if
if (any(EDA_nAtom <= 0)) then
  if (IsRoot) write(u6,'(6X,A)') 'EDA: NATOM must contain a positive atom count for every fragment.'
  return
end if
if (any(EDA_Mult /= 1)) then
  if (IsRoot) write(u6,'(6X,A)') 'EDA currently supports only closed-shell singlet fragments (MULT=1).'
  return
end if

if ((trim(EDA_RefMethod) /= 'AUTO') .and. (trim(EDA_RefMethod) /= 'RHF') .and. &
    (trim(EDA_RefMethod) /= 'MBPT2')) then
  if (IsRoot) write(u6,'(6X,2A)') 'EDA: REFERENCE must be AUTO, RHF, or MBPT2; received ',trim(EDA_RefMethod)
  return
end if

if ((trim(EDA_FragBasis) /= 'OWN') .and. (trim(EDA_FragBasis) /= 'SUPER')) then
  if (IsRoot) write(u6,'(6X,2A)') 'EDA: FBASIS must be OWN or SUPER; received ',trim(EDA_FragBasis)
  return
end if

if (EDA_Debug .and. IsRoot) then
  write(u6,'(/,6X,A)') 'EDA input summary:'
  write(u6,'(8X,A,I0)') 'Number of fragments: ',EDA_nFrag
  do i=1,EDA_nFrag
    write(u6,'(8X,A,I0,A,I0,A,I0,A,I0)') 'Fragment ',i,': NATOM=',EDA_nAtom(i), &
                                          ' CHARGE=',EDA_Charge(i),' MULT=',EDA_Mult(i)
  end do
  write(u6,'(8X,2A)') 'Reference: ',trim(EDA_RefMethod)
  write(u6,'(8X,2A)') 'Fragment basis: ',trim(EDA_FragBasis)
end if

iReturn = _RC_ALL_IS_WELL_

end subroutine EDA_Read_Input
