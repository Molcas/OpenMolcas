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
subroutine close_files_mcpdft()
!***********************************************************************
!                                                                      *
!     Close files.                                                     *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1993                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************
    use Fock_util_global, only: docholesky
    use mcpdft_output, only: lf

    implicit none

    integer :: return_code, iOpt

#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "warnings.h"

    !---  close the JOBOLD file -------------------------------------------*
    If(JOBOLD.gt.0.and.JOBOLD.ne.JOBIPH) Then
        Call DaClos(JOBOLD)
        JOBOLD=-1
    Else If (JOBOLD.gt.0) Then
        JOBOLD=-1
    End If
    !---  close the JOBIPH file -------------------------------------------*
    If(JOBIPH.gt.0) Then
        Call DaClos(JOBIPH)
        JOBIPH=-1
    End If
    !---  close the ORDINT file -------------------------------------------*
    If (.not.DoCholesky) then
        return_code = -1
        Call ClsOrd(return_code)
        If ( return_code.ne._RC_ALL_IS_WELL_) Then
            Call WarningMessage(1,'Failed to close the ORDINT file.')
        End If
    End If
    !---  close the file carrying the transformed two-electron integrals --*
    Call DaClos(LUINTM)

    !--- close the one-electorn integral file
    return_code = -1
    iOpt = 0
    call clsone(return_code, iOpt)
    if (return_code .ne. _RC_ALL_IS_WELL_) then
        write(lf, *) "Error when trying to close the one-electron"
        write(lf, *) "integral file."
        call abend()
    end if
    Return
End
