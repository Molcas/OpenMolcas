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
      SubRoutine Cho_1VecRd_SP(Vec,lVec,jVec,iSym,LstSP,nSP,iRedC,
     &                         iLoc)
C
C     Purpose: read vector jVec, sym. iSym, from disk. Read only
C              components from shell pairs in list LstSP. Use scratch
C              location iLoc to set index arrays. On input, iRedC
C              specifies the reduced set available at location iLoc:
C              specify -1 if not set (or unknown). On exit, iRedC
C              identifies the reduced set for which indices are
C              available at location iLoc. NOTE: only WA files!!
C
      use ChoSwp, only: nnBstRSh, iiBstRSh
      Implicit Real*8 (a-h,o-z)
      Real*8  Vec(lVec)
      Integer LstSP(nSP)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Character*13 SecNam
      Parameter (SecNam = 'Cho_1VecRd_SP')

      Parameter (N2 = InfVec_N2)

      Integer  Cho_P_LocalSP
      External Cho_P_LocalSP

      InfVec(i,j,k)=iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)

C     Return if no vectors are available on disk.
C     -------------------------------------------

      If (NumCho(iSym) .lt. 1) Return

C     Check that vector storage mode is word-addressable (WA).
C     --------------------------------------------------------

      If (Cho_AdrVec .ne. 1) Then
         Write(Lupri,*) SecNam,': WA address mode is required!'
         Write(Lupri,*) 'Cho_AdrVec is: ',Cho_AdrVec,
     &                  ' (should be 1)'
         Call Cho_Quit('WA address mode is required in '//SecNam,104)
      End If

C     Get reduced set of this vector.
C     -------------------------------

      If (jVec.gt.0 .and. jVec.le.NumCho(iSym)) Then
         iRed = InfVec(jVec,2,iSym)
      Else
         Call Cho_Quit('Red. set error in '//SecNam,104)
         iRed = -999999
      End If

C     Set reduced set (if needed).
C     ----------------------------

      If (iRedC .ne. iRed) Then
         Call Cho_X_SetRed(irc,iLoc,iRed)
         If (irc .ne. 0) Then
            Write(Lupri,*) SecNam,': Cho_X_SetRed returned ',irc
            Call Cho_Quit('Error in '//SecNam,104)
         End If
         iRedC = iRed
      End If

C     Read vector elements.
C     ---------------------

      iAdr0 = InfVec(jVec,3,iSym)
      kV = 1
      Do iSP = 1,nSP
         iShlAB = Cho_P_LocalSP(LstSP(iSP))
         iOpt = 2
         lTot = nnBstRSh(iSym,iShlAB,iLoc)
         iAdr = iAdr0 + iiBstRSh(iSym,iShlAB,iLoc)
         Call dDAFile(LuCho(iSym),iOpt,Vec(kV),lTot,iAdr)
         kV = kV + lTot
      End Do

      End
