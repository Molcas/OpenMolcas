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
      Subroutine InitB(nMult,natom,nAtQM,nGrdPt,ipCord,ipGrid,ipT,ipTT,
     &           ipTTT,ipExt,ipB,ipIsMM)
      Implicit Real*8 (A-H,O-Z)
*
*     Compute the electrostatic tensor matrix between QM atoms and
*     grid points
*     Then form the B = ExtPot[TtT^-1]Tt vector
*
#include "espf.fh"
*
      Call QEnter('initb')
      iPL = iPL_espf()
      nOrd = nMult/nAtQM
*
*     T
*
      Do iPnt = 1, nGrdPt
         iQM = 0
         Do jAt = 1, natom
            If (iWork(ipIsMM+jAt-1).eq.1) GoTo 1
            iQM = iQM + 1
            X = Work(ipGrid+(IPnt-1)*3  ) - Work(ipCord+(jAt-1)*3  )
            Y = Work(ipGrid+(IPnt-1)*3+1) - Work(ipCord+(jAt-1)*3+1)
            Z = Work(ipGrid+(IPnt-1)*3+2) - Work(ipCord+(jAt-1)*3+2)
            R = sqrt(X*X+Y*Y+Z*Z)
            iCur = ipT + (iPnt-1)*nMult+nOrd*(iQM-1)
            Work(iCur) = One / R
            If (nOrd.gt.1) Then
               R3 = R * R * R
               Work(iCur+1) = X / R3
               Work(iCur+2) = Y / R3
               Work(iCur+3) = Z / R3
            End If
1           Continue
         End Do
      End Do
      If (iQM.ne.nAtQM) Then
         Write(6,'(A,I4,A4,I4)') ' Error in espf/initb: iQM != nAtQM ',
     &                           iQM,' != ',nAtQM
         Call Abend()
      End If
*
*
*
*      Call RecPrt('T',' ',Work(ipT),nGrdPt,nMult)
*      nMax = Max(nGrdPt,nMult)
*      Call Allocate_Work(ipU,nMax*nMult)
*      Call Allocate_Work(ipW,nMult)
*      Call Allocate_Work(ipV,nMax*nMult)
*      Call Allocate_Work(ipScr,nMult)
*      Call SVD(nMax,nGrdPt,nMult,Work(ipT),Work(ipW),.True.,
*     &         Work(ipU),.True.,Work(ipV),iErr,Work(ipScr))
*      print*,'iErr=',iErr
*      Call RecPrt('U',' ',Work(ipU),nMax,nMult)
*      Call RecPrt('w',' ',Work(ipW),nMult,1)
*      Call RecPrt('V',' ',Work(ipV),nMax,nMult)
*
*      Call Free_Work(ipW)
*      Call Free_Work(ipU)
*      Call Free_Work(ipV)
*      Call Free_Work(ipScr)
*
*     TtT
*
      Do iMlt = 1, nMult
         Do jMlt = 1, nMult
            iCur = ipTT + (iMlt-1)*nMult + (jMlt-1)
            Work(iCur) = Zero
            Do kPnt = 1, nGrdPt
               Work(iCur) = Work(iCur)
     &                    + Work(ipT+(kPnt-1)*nMult+(iMlt-1))
     &                    * Work(ipT+(kPnt-1)*nMult+(jMlt-1))
            End Do
         End Do
      End Do
*
*     TtT^-1
*
      Call Allocate_Work(ipScr,nMult*nMult)
      Call minv(Work(ipTT),Work(ipScr),Ising,Det,nMult)
      Call dCopy_(nMult*nMult,Work(ipScr),1,Work(ipTT),1)
      Call Free_Work(ipScr)
*
*     [TtT^-1]Tt
*
      Do iMlt = 1, nMult
         Do jPnt = 1, nGrdPt
            iCur = ipTTT + (iMlt-1)*nGrdPt + (jPnt-1)
            Work(iCur) = Zero
            Do kMlt = 1, nMult
               Work(iCur) = Work(iCur)
     &                    + Work(ipTT+(iMlt-1)*nMult+(kMlt-1))
     &                    * Work(ipT+(jPnt-1)*nMult+(kMlt-1))
            End Do
         End Do
      End Do
      If (iPL.ge.4) Call RecPrt('(TtT)^(-1)Tt matrix in InitB',' ',
     &             Work(ipTTT),nMult,nGrdPt)
*
*     B = ExtPot[TtT^-1]Tt
*
      Do iPnt = 1, nGrdPt
         iQM = 0
         iCur = ipB + (iPnt-1)
         Work(iCur) = Zero
         Do jAt = 1, natom
            If (iWork(ipIsMM+jAt-1).eq.1) Goto 2
            iQM = iQM + 1
            Work(iCur) = Work(iCur) + Work(ipExt+(jAt-1)*10)
     &                 * Work(ipTTT+nOrd*(iQM-1)*nGrdPt+(iPnt-1))
            If (nOrd.gt.1) Then
               Work(iCur) = Work(iCur)
     &                    + Work(ipExt+(jAt-1)*10+1)
     &                    * Work(ipTTT+(nOrd*(iQM-1)+1)*nGrdPt+(iPnt-1))
     &                    + Work(ipExt+(jAt-1)*10+2)
     &                    * Work(ipTTT+(nOrd*(iQM-1)+2)*nGrdPt+(iPnt-1))
     &                    + Work(ipExt+(jAt-1)*10+3)
     &                    * Work(ipTTT+(nOrd*(iQM-1)+3)*nGrdPt+(iPnt-1))
            End If
2           Continue
         End Do
      End Do
      If (iPL.ge.4) Then
         Write(6,'(A)') ' In InitB (grid coordinates, B value)'
         Do iPnt = 1, nGrdPt
           Write(6,1234) iPnt,(Work(ipGrid+(iPnt-1)*3+J),J=0,2),
     &                        Work(ipB+iPnt-1)
1234       Format(I4,4F12.6)
         End do
      End If
*
      Call QExit('initb')
      Return
      End
