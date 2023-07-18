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
! Copyright (C) Thomas Bondo Pedersen                                  *
!***********************************************************************
!  Cho_X_Final
!
!> @brief
!>   Finalize Cholesky utilities
!> @author Thomas Bondo Pedersen
!>
!> @details
!> Deallocates memory, closes files, etc., as initialized
!> by ::Cho_X_Init. On exit, \p irc = ``0`` signals successful finalization.
!>
!> @param[out] irc Return code
!***********************************************************************
      Subroutine Cho_X_Final(irc)
      use ChoArr, only: MySP
      use ChoBkm, only: BkmVec, BkmThr, nRow_BkmVec, nCol_BkmVec,       &
     &                   nRow_BkmThr, nCol_BkmThr
      use ChoIni
      use stdalloc, only: mma_deallocate
      Implicit None
      Integer irc

      Character*11 SecNam
      Parameter (SecNam = 'Cho_X_Final')

      Integer ChoIsIni

!     Register entry.
!     ---------------


!     Set  error code.
!     ----------------

      irc = 0

!     Read initialization integer flag from runfile.
!     ----------------------------------------------

      Call Get_iScalar('ChoIni',ChoIsIni)

!     Finalize if needed.
!     -------------------

      If (ChoIsIni .eq. ChoIniCheck) Then

!        Close files.
!        ------------

         Call Cho_OpenVR(2,2)

!        Deallocate vector buffer.
!        -------------------------

         Call Cho_VecBuf_Final()

!        Deallocate memory.
!        ------------------

         Call Cho_X_Dealloc(irc)
         If (irc .ne. 0) Go To 1

         If (Allocated(MySP)) Call mma_deallocate(MySP)
         If (Allocated(BkmVec)) Then
            Call mma_deallocate(BkmVec)
            nRow_BkmVec=0
            nCol_BkmVec=0
         End If
         If (Allocated(BkmThr)) Then
            Call mma_deallocate(BkmThr)
            nRow_BkmThr=0
            nCol_BkmThr=0
         End If

!        Reset initialization integer on runfile to "not set".
!        -----------------------------------------------------

    1    ChoIsIni = ChoIniCheck + 1
         Call Put_iScalar('ChoIni',ChoIsIni)

      End If

!     Register exit and return.
!     -------------------------

      End
