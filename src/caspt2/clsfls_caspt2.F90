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
! Copyright (C) 1993, Markus P. Fuelscher                              *
!***********************************************************************

subroutine ClsFls_CASPT2()
!***********************************************************************
!     Close files.                                                     *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!***********************************************************************

use definitions, only: iwp
use caspt2_global, only: iPrGlb
use PrintLevel, only: SILENT
use caspt2_global, only: LUCIEX, LUONEM, LUHLF1, LUHLF2, LUHLF3, LUINTM, LUDMAT, LUDRA, LUDRATOT, LURHS, LUH0T, LUSOLV, LUSBT
use caspt2_module, only: IfChol

implicit none
integer(kind=iwp) IMAT, iRc, IVEC

!----------------------------------------------------------------------*
!     Start                                                            *
!-------------------------------------- -------------------------------*

call DaClos(LUCIEX)
call DaClos(LUONEM)
call DaClos(LUINTM)
call DaClos(LUDRA)
call DaClos(LUDRATOT)
call DaClos(LUHLF1)
call DaClos(LUHLF2)
call DaClos(LUHLF3)
call DaClos(LUDMAT)
call DaClos(LUSOLV)
call DaClos(LUSBT)
do IVEC=1,8
  call DaClos(LURHS(IVEC))
end do
do IMAT=1,4
  call DaClos(LUH0T(IMAT))
end do
!---  close the ORDINT file -------------------------------------------*
if (.not. IfChol) then
  iRc = -1
  call ClsOrd(iRc)
  if ((IRC /= 0) .and. (IPRGLB > SILENT)) call WarningMessage(1,'Failed to close ORDINT file.')
end if
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*

end subroutine ClsFls_CASPT2
