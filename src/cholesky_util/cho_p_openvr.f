************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SubRoutine Cho_P_OpenVR(iOpt)
C
C     Purpose: open (iOpt=1) or close (iOpt=2) local and global storage
C              files.
C
      Use Para_Info, Only: nProcs, Is_Real_Par
      Implicit None
      Integer iOpt
#include "cholesky.fh"
      Integer ID
#include "cho_para_info.fh"
#include "choglob.fh"
      Integer iSym
      Character*5 FNRed
      Character*6 FNRst, FNVec(8)
      Character*12 SecNam
      Parameter (SecNam = 'Cho_P_OpenVR')

C     Local files.
C     ------------

      If (Cho_Real_Par) Then
         ID = 1
      Else
         ID = 2
      End If
      Call Cho_OpenVR(iOpt,ID)

C     Global files for restart info and reduced set indices.
C     ------------------------------------------------------

      If (Cho_Real_Par) Then
         If (iOpt .eq. 1) Then
            LuRed_G = 7
            FNRed = 'CHRED'
            Call DAName_MF_WA(LuRed_G,FNRed)
            LuRst_G = 7
            FNRst = 'CHORST'
            Call DAName_MF_WA(LuRst_G,FNRst)
            Do iSym = 1,nSym
               LuCho_G(iSym) = 7
               Write(FNVec(iSym),'(A5,I1)') 'CHVEC',iSym
               Call DaName_MF_WA(LuCho_G(iSym),FNVec(iSym))
            End Do
         Else If (iOpt .eq. 2) Then
            If (LuRed_G .gt. 0) Then
               Call DAClos(LuRed_G)
               LuRed_G = 0
            End If
            If (LuRst_G .gt. 0) Then
               Call DAClos(LuRst_G)
               LuRst_G = 0
            End If
            Do iSym = 1,nSym
               If (LuCho_G(iSym) .gt. 0) Then
                  Call DaClos(LuCho_G(iSym))
                  LuCho_G(iSym) = 0
               End If
            End Do
         Else
            Write(Lupri,*) SecNam,': iOpt out of bounds: ',iOpt
            Call Cho_Quit('Error in '//SecNam,104)
         End If
      Else
         If (CHO_FAKE_PAR .and. nProcs.gt.1 .and. Is_Real_Par()) Then
            If (iOpt .eq. 1) Then
               If (CHO_ADRVEC .eq. 1) Then
                  Do iSym = 1,nSym
                     LuCho_G(iSym) = 7
                     Write(FNVec(iSym),'(A5,I1)') 'CHVCL',iSym
                     Call DaName_MF_WA(LuCho_G(iSym),FNVec(iSym))
                  End Do
               Else If (CHO_ADRVEC .eq. 2) Then
                  Do iSym = 1,nSym
                     LuCho_G(iSym) = 7
                     Write(FNVec(iSym),'(A5,I1)') 'CHVCL',iSym
                     Call DaName_MF(LuCho_G(iSym),FNVec(iSym))
                  End Do
               Else
                  Call Cho_Quit('CHO_ADRVEC out of bounds in '//SecNam,
     &                          102)
                  Call iZero(LuCho_G,nSym)
               End If
C              Swap units so that
C                 LuCho_G points to 'CHVEC'
C                 LuCho   points to 'CHVCL'
               Call iSwap(nSym,LuCho,1,LuCho_G,1)
            Else If (iOpt .eq. 2) Then
               Do iSym = 1,nSym
                  If (LuCho_G(iSym) .gt. 0) Then
                     Call DaClos(LuCho_G(iSym))
                     LuCho_G(iSym) = 0
                  End If
               End Do
            Else
               Write(Lupri,*) SecNam,': iOpt out of bounds: ',iOpt
               Call Cho_Quit('Error in '//SecNam,104)
            End If
         End If
      End If

      End
