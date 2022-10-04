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
!  Get_Grad_Full
!
!> @brief
!>   Get gradient from RUNFILE
!> @author I. Fdez. Galv&aacute;n
!>
!> @details
!> Place Cartesian gradient (in a.u.) into array \p Grad_Full(3,*).
!> Includes MM atoms otherwise invisible to gateway/slapaf.
!>
!> @param[out] Grad_Full   Array of gradient
!> @param[in]  nAtoms_Full Number of atoms
!***********************************************************************

subroutine Get_Grad_Full(Grad_Full,nAtoms_Full)

implicit none
integer nAtoms_Full, nAtoms_Fullx, nAtoms_All, nGrad, nGradMM
real*8 Grad_Full(3,nAtoms_Full)
logical Found

call Get_nAtoms_Full(nAtoms_Fullx)
if (nAtoms_Full /= nAtoms_Fullx) then
  write(6,*) 'Get_Grad_Full: nAtoms_Full /= nAtoms_Fullx'
  write(6,*) 'nAtoms_Full=',nAtoms_Full
  write(6,*) 'nAtoms_Fullx=',nAtoms_Fullx
  call Abend()
end if
call Get_nAtoms_All(nAtoms_All)
if (nAtoms_Full < nAtoms_All) then
  write(6,*) 'Get_Coord_Full: nAtoms_Full < nAtoms_All'
  write(6,*) 'nAtoms_Full=',nAtoms_Full
  write(6,*) 'nAtoms_Fullx=',nAtoms_All
  call Abend()
end if
call Qpg_dArray('GRAD',Found,nGrad)
if ((.not. Found) .or. (nGrad == 0)) then
  write(6,*) 'Get_Grad_Full: Did not find GRAD'
  call Abend()
end if
call Get_dArray('GRAD',Grad_Full,nGrad)
call Qpg_dArray('MMO Grad',Found,nGradMM)
if (Found) call Get_dArray('MMO Grad',Grad_Full(1,nAtoms_All+1),nGradMM)

return

end subroutine Get_Grad_Full
