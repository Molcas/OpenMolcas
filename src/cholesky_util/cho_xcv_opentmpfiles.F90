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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_XCV_OpenTmpFiles()

use Cholesky, only: LuTmp, nSym
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: iSym
character(len=6) :: Filename

do iSym=1,nSym
  LuTmp(iSym) = 7
  write(Filename,'(A4,I2.2)') 'VTMP',iSym
  call DAName_MF_WA(LuTmp(iSym),Filename)
end do

end subroutine Cho_XCV_OpenTmpFiles
