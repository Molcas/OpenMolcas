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
!               2024, Matthew R. Hennefarth                            *
!***********************************************************************
subroutine close_files_mcpdft()
  use definitions,only:iwp,u6
  use Fock_util_global,only:docholesky
  use general_data,only:jobiph,jobold,luintm

  implicit none

#include "warnings.h"

  integer(kind=iwp) :: return_code,iOpt

  !---  close the JOBOLD file -------------------------------------------*
  If(JOBOLD > 0 .and. JOBOLD /= JOBIPH) Then
    Call DaClos(JOBOLD)
    JOBOLD = -1
  Else If(JOBOLD > 0) Then
    JOBOLD = -1
  EndIf
  !---  close the JOBIPH file -------------------------------------------*
  If(JOBIPH > 0) Then
    Call DaClos(JOBIPH)
    JOBIPH = -1
  EndIf
  !---  close the ORDINT file -------------------------------------------*
  If(.not. DoCholesky) then
    return_code = -1
    Call ClsOrd(return_code)
    If(return_code /= _RC_ALL_IS_WELL_) Then
      Call WarningMessage(1,'Failed to close the ORDINT file.')
    EndIf
  EndIf
  !---  close the file carrying the transformed two-electron integrals --*
  Call DaClos(LUINTM)

  !--- close the one-electorn integral file
  return_code = -1
  iOpt = 0
  call clsone(return_code,iOpt)
  if(return_code /= _RC_ALL_IS_WELL_) then
    write(u6,*) "Error when trying to close the one-electron"
    write(u6,*) "integral file."
    call abend()
  endif
  Return
End
