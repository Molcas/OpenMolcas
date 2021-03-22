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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************
      Subroutine RPA(rc)
!
!     Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
!     Driver routine for the calculation of correlation energies in the
!     random-phase approximation (RPA) using Cholesky/DF integrals.
!
!     NOTE: conventional integrals are not implemented!
!
      Implicit None
      Integer rc
#include "warnings.fh"

      Character*3 SecNam
      Parameter (SecNam='RPA')
      Character*80 string

      Integer irc

!=======================================================================
!     Enter
!=======================================================================

      rc=_RC_ALL_IS_WELL_ ! init return code
      irc=0 ! init internal return code

!=======================================================================
!     Setup: initialize data, read reference orbitals and process input
!=======================================================================

      Call StatusLine('RPA: ','Setup')
      Call RPA_Setup()

!=======================================================================
!     Exit after cleanup
!=======================================================================

!100  Continue ! failures jump to this point
      Call StatusLine('RPA: ','Cleanup')
      Call RPA_Cleanup(irc)
      If (irc.ne.0) Then
         Write(string,'(A,A,I4)')                                       &
     &   SecNam,': Cleanup failed! rc=',irc
         Call WarningMessage(2,string)
         If (rc.eq._RC_ALL_IS_WELL_) rc=_RC_INTERNAL_ERROR_
      End If

      End
