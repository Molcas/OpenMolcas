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
! Copyright (C) Valera Veryazov                                        *
!***********************************************************************
!  get_BasisType
!
!> @brief
!>   Logical function to check basis set
!> @author V. Veryazov
!>
!> @details
!> Function returns ``.true.`` if the basis set in the current
!> calculation has specific type.
!> Only 3 first characters are used, so
!> ``get_BasisType('segmented')`` is the same as ``get_BasisType('SEG')``.
!>
!> The list of available basis set types is available at
!> src/input_util/basistype.fh and basis_library/basistype.tbl.
!>
!> @param[in] BasisType Basis set type
!>
!> @return ``.true.`` if the basis set in the current calculation has specific type
!***********************************************************************

function get_BasisType(BasisType)

use BasisType_Mod, only: BasTypeAll, BasTypeCon, BasTypeRel
use Definitions, only: iwp

implicit none
logical(kind=iwp) :: get_BasisType
character(len=*), intent(in) :: BasisType
integer(kind=iwp) :: BasisTypes(4), i, nData
logical(kind=iwp) :: Found
character(len=3) :: temp, TypeCon, TypeAll, TypeRel

nData = 0
Found = .false.
call Qpg_iArray('BasType',Found,nData)
if (.not. Found) then
  get_BasisType = .false.
  return
end if
call get_bastype(BasisTypes,nData)
i = BasisTypes(1)
if (i <= 0) then
  TypeCon = 'UNK'
else
  TypeCon = BasTypeCon((i-1)*4+1:(i-1)*4+3)
end if
i = BasisTypes(2)
if (i <= 0) then
  TypeAll = 'UNK'
else
  TypeAll = BasTypeAll((i-1)*4+1:(i-1)*4+3)
end if
i = BasisTypes(3)
if (i <= 0) then
  TypeRel = 'UNK'
else
  TypeRel = BasTypeRel((i-1)*4+1:(i-1)*4+3)
end if

i = len(BasisType)
if (i > 3) i = 3
temp = '___'
temp(1:i) = BasisType(1:i)
call UpCase(temp)
get_BasisType = .false.
if ((temp == TypeCon) .or. (temp == TypeAll) .or. (temp == TypeRel)) get_BasisType = .true.

return

end function get_BasisType
