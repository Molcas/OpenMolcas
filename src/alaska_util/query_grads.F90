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
!  Query_Grads
!
!> @brief Query sizes from a gradients file
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Reads the number of roots and length of the gradients from the
!> gradients file (GRADS).
!>
!> @param[out] Exists Whether or not the gradients file exists
!> @param[out] nRoots Number of roots allowed in the gradients file
!> @param[out] nGrad  Length of the vectors in the gradients file
!***********************************************************************

subroutine Query_Grads(Exists,nRoots,nGrad)

use Definitions, only: iwp

implicit none
logical(kind=iwp), intent(out) :: Exists
integer(kind=iwp), intent(out) :: nRoots, nGrad
integer(kind=iwp) :: iAd, iDum(1), LuGrad, TOC(5)
logical(kind=iwp) :: Found
character(len=5), parameter :: Filename = 'GRADS'

call f_Inquire(Filename,Found)
if (.not. Found) then
  Exists = .false.
  nRoots = 0
  nGrad = 0
  return
end if
Exists = .true.

LuGrad = 20
call DaName(LuGrad,Filename)
iAd = 0
call iDaFile(LuGrad,2,TOC,size(TOC),iAd)
call iDaFile(LuGrad,2,iDum,1,iAd)
nRoots = iDum(1)
call iDaFile(LuGrad,2,iDum,1,iAd)
nGrad = iDum(1)
call DaClos(LuGrad)

end subroutine Query_Grads
