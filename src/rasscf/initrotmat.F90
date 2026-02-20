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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
!*****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 07, 2022, created this file.               *
!*****************************************************************

subroutine InitRotMat(RotMat,lRoots,CMSSFile,LenCMSS)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lRoots, LenCMSS
real(kind=wp), intent(out) :: RotMat(lRoots,lRoots)
character(len=LenCMSS), intent(in) :: CMSSFile
character(len=16) :: ScrChar

if (CMSSFile == 'XMS') then
  call ReadMat('ROT_VEC',ScrChar,RotMat,lroots,lroots,7,16,'T')
else
  call ReadMat(CMSSFile,ScrChar,RotMat,lroots,lroots,LenCMSS,16,'T')
end if

return

end subroutine InitRotMat
