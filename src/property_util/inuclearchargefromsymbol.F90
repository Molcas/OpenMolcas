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
! Copyright (C) 2000, Per-Olof Widmark                                 *
!***********************************************************************

function iNuclearChargeFromSymbol(Symbol)
!***********************************************************************
!                                                                      *
! Routine: iNuclearChargeFromSymbol                                    *
! Purpose: Wrapper for ixNuclearChargeFromSymbol, to give the nuclear  *
!          charge for atom with symbol 'Symbol'.                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Author:  Per-Olof Widmark                                            *
!          Lund University, Sweden.                                    *
! Written: March 2000                                                  *
! History: none.                                                       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
! Parameters:                                                          *
!                                                                      *
!***********************************************************************

implicit none
integer iNuclearChargeFromSymbol
!----------------------------------------------------------------------*
! Dummy parameters.                                                    *
!----------------------------------------------------------------------*
character*(*) Symbol
!----------------------------------------------------------------------*
! Local variables.                                                     *
!----------------------------------------------------------------------*
integer Z
integer Opt
integer Rc
!----------------------------------------------------------------------*
! External references                                                  *
!----------------------------------------------------------------------*
integer ixNuclearChargeFromSymbol
external ixNuclearChargeFromSymbol

Opt = 0
Rc = 0
Z = ixNuclearChargeFromSymbol(Symbol,Rc,Opt)
if (Rc /= 0) then
  call SysAbendMsg('inuclearchargefromsymbol','Fail to get nuclear charge',' ')
end if
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
iNuclearChargeFromSymbol = Z

return

end function iNuclearChargeFromSymbol
