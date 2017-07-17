************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
*  Cho_X_Final
*
*> @brief
*>   Finalize Cholesky utilities
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Deallocates memory, closes files, etc., as initialized
*> by ::Cho_X_Init. On exit, \p irc = ``0`` signals successful finalization.
*>
*> @param[out] irc Return code
************************************************************************
      Subroutine Cho_X_Final(irc)
      Implicit None
      Integer irc
#include "choini.fh"
#include "choptr2.fh"
#include "chobkm.fh"

      Character*11 SecNam
      Parameter (SecNam = 'Cho_X_Final')

      Integer ChoIsIni

C     Register entry.
C     ---------------

      Call Qenter(SecNam)

C     Set  error code.
C     ----------------

      irc = 0

C     Read initialization integer flag from runfile.
C     ----------------------------------------------

      Call Get_iScalar('ChoIni',ChoIsIni)

C     Finalize if needed.
C     -------------------

      If (ChoIsIni .eq. ChoIniCheck) Then

C        Close files.
C        ------------

         Call Cho_OpenVR(2,2)

C        Deallocate vector buffer.
C        -------------------------

         Call Cho_VecBuf_Final()

C        Deallocate memory.
C        ------------------

         Call Cho_X_Dealloc(irc)
         If (irc .ne. 0) Go To 1

         If (l_mySP .gt. 0) Then
            Call GetMem('mySP','Free','Inte',ip_mySP,l_mySP)
            l_mySP = 0
         End If
         If (l_BkmVec.gt.0) Then
            Call GetMem('BkmVec','Free','Inte',ip_BkmVec,l_BkmVec)
            ip_BkmVec=0
            l_BkmVec=0
            nRow_BkmVec=0
            nCol_BkmVec=0
         End If
         If (l_BkmThr.gt.0) Then
            Call GetMem('BkmThr','Free','Real',ip_BkmThr,l_BkmThr)
            ip_BkmThr=0
            l_BkmThr=0
            nRow_BkmThr=0
            nCol_BkmThr=0
         End If

C        Reset initialization integer on runfile to "not set".
C        -----------------------------------------------------

    1    ChoIsIni = ChoIniCheck + 1
         Call Put_iScalar('ChoIni',ChoIsIni)

      End If

C     Register exit and return.
C     -------------------------

      Call Qexit(SecNam)
      End
