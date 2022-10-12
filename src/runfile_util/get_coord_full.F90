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
! Copyright (C) 2015, Ignacio Fdez. Galvan                             *
!***********************************************************************
!  Get_Coord_Full
!
!> @brief
!>   Get coordinates from RUNFILE
!> @author I. Fdez. Galv&aacute;n
!>
!> @details
!> Place Cartesian coordinates (in a.u.) into array \p Coord_Full(3,*).
!> Includes MM atoms otherwise invisible to gateway/slapaf.
!>
!> @param[out] Coord_Full  Array of coordinates
!> @param[in]  nAtoms_Full Number of atoms
!***********************************************************************

subroutine Get_Coord_Full(Coord_Full,nAtoms_Full)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms_Full
real(kind=wp), intent(out) :: Coord_Full(3,nAtoms_Full)
integer(kind=iwp) :: nAtoms_All, nAtoms_Fullx, nCoordMM
logical(kind=iwp) :: Found

call Get_nAtoms_Full(nAtoms_Fullx)
if (nAtoms_Full /= nAtoms_Fullx) then
  write(u6,*) 'Get_Coord_Full: nAtoms_Full /= nAtoms_Fullx'
  write(u6,*) 'nAtoms_Full=',nAtoms_Full
  write(u6,*) 'nAtoms_Fullx=',nAtoms_Fullx
  call Abend()
end if
call Get_nAtoms_All(nAtoms_All)
if (nAtoms_Full < nAtoms_All) then
  write(u6,*) 'Get_Coord_Full: nAtoms_Full < nAtoms_All'
  write(u6,*) 'nAtoms_Full=',nAtoms_Full
  write(u6,*) 'nAtoms_Fullx=',nAtoms_All
  call Abend()
end if
call Get_Coord_All(Coord_Full,nAtoms_All)
call Qpg_dArray('MMO Coords',Found,nCoordMM)
if (Found) call Get_dArray('MMO Coords',Coord_Full(:,nAtoms_All+1:),nCoordMM)

return

end subroutine Get_Coord_Full
