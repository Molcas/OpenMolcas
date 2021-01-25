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
      Subroutine InitDB(nMult,natom,nAtQM,nGrdPt,ipCord,ipGrid,ipT,ipTT,
     &           ipTTT,ipExt,ipDB,ipIsMM,iRMax,DeltaR,iGrdTyp,ipDGrd)
      Implicit Real*8 (A-H,O-Z)
*
*     Compute DB = d(ExtPot[TtT^-1]Tt)/dq
*
#include "espf.fh"
*
*     Various derivatives
*
      Call GetMem('DTTT','Allo','Real',ipDTTT,nMult*nGrdPt*nAtQM*3)
      Call GetMem('DT','Allo','Real',ipDT,nGrdPt*nMult*3*nAtQM)
      Call GetMem('DTT','Allo','Real',ipDTT,nMult*nMult*nAtQM*3)
      Call GetMem('DTTTT','Allo','Real',ipDTTTT,nMult*nMult*nAtQM*3)

      Call CalcDT(nMult,nGrdPt,natom,nAtQM,ipIsMM,iGrdTyp,Work(ipCord),
     &            Work(ipGrid),Work(ipDGrd),Work(ipT),Work(ipTT),
     &            Work(ipDT),Work(ipDTT),Work(ipDTTTT),Work(ipDTTT))

      Call GetMem('DTTTT','Free','Real',ipDTTTT,nMult*nMult*nAtQM*3)
      Call GetMem('DTT','Free','Real',ipDTT,nMult*nMult*nAtQM*3)
      Call GetMem('DT','Free','Real',ipDT,nGrdPt*nMult*3*nAtQM)
*
*     Finally dB = dTTT * V_ext + TTT * dV_ext
*
      Call CalcDB(nMult,nGrdPt,natom,nAtQM,ipIsMM,Work(ipTTT),
     &            Work(ipDTTT),Work(ipExt),Work(ipDB))
      Call GetMem('DTTT','Free','Real',ipDTTT,nMult*nGrdPt*nAtQM*3)
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer(iRMax)
        Call Unused_real(DeltaR)
      End If
      End
*
      Subroutine CalcDB(nMult,nGrdPt,natom,nAtQM,ipIsMM,TTT,DTTT,ExtPot,
     &                  DB)
      Implicit Real*8 (A-H,O-Z)
*
*     dB = dTTT * V_ext + TTT * dV_ext
*
#include "espf.fh"
*
      Dimension TTT(nGrdPt,nMult),DTTT(nMult,nGrdPt,3,nAtQM),
     &          ExtPot(10,natom),DB(nGrdPt,3,nAtQM)
*
      iPL = iPL_espf()
*
      If (iPL.ge.4) Call RecPrt('TTT in calcdb',' ',TTT,nMult,nGrdPt)
      nOrd = nMult/nAtQM
      Do iPnt = 1, nGrdPt
         iQM = 0
         Do iAt = 1, natom
            If (iWork(ipIsMM+iAt-1) .ne. 0) Goto 1
            iQM = iQM + 1
            DB(iPnt,1,iQM)= TTT(iPnt,nOrd*(iQM-1)+1)*ExtPot(2,iAt)
            DB(iPnt,2,iQM)= TTT(iPnt,nOrd*(iQM-1)+1)*ExtPot(3,iAt)
            DB(iPnt,3,iQM)= TTT(iPnt,nOrd*(iQM-1)+1)*ExtPot(4,iAt)
            If (nOrd.eq.4) Then
               DB(iPnt,1,iQM)= DB(iPnt,1,iQM)
     &                       + TTT(iPnt,nOrd*(iQM-1)+2)*ExtPot( 5,iAt)
     &                       + TTT(iPnt,nOrd*(iQM-1)+3)*ExtPot( 8,iAt)
     &                       + TTT(iPnt,nOrd*(iQM-1)+4)*ExtPot( 9,iAt)
               DB(iPnt,2,iQM)= DB(iPnt,2,iQM)
     &                       + TTT(iPnt,nOrd*(iQM-1)+2)*ExtPot( 8,iAt)
     &                       + TTT(iPnt,nOrd*(iQM-1)+3)*ExtPot( 6,iAt)
     &                       + TTT(iPnt,nOrd*(iQM-1)+4)*ExtPot(10,iAt)
               DB(iPnt,3,iQM)= DB(iPnt,3,iQM)
     &                       + TTT(iPnt,nOrd*(iQM-1)+2)*ExtPot( 9,iAt)
     &                       + TTT(iPnt,nOrd*(iQM-1)+3)*ExtPot(10,iAt)
     &                       + TTT(iPnt,nOrd*(iQM-1)+4)*ExtPot( 7,iAt)
            End If
            jQM = 0
            Do jAt = 1, natom
               If (iWork(ipIsMM+jAt-1) .ne. 0) Goto 2
               jQM = jQM + 1
               Do iOrd = 1, nOrd
                  iMlt = nOrd*(jQM-1) + iOrd
                  DB(iPnt,1,iQM)= DB(iPnt,1,iQM)
     &                          + DTTT(iMlt,iPnt,1,iQM)*ExtPot(iOrd,jAt)
                  DB(iPnt,2,iQM)= DB(iPnt,2,iQM)
     &                          + DTTT(iMlt,iPnt,2,iQM)*ExtPot(iOrd,jAt)
                  DB(iPnt,3,iQM)= DB(iPnt,3,iQM)
     &                          + DTTT(iMlt,iPnt,3,iQM)*ExtPot(iOrd,jAt)
               End Do
2              Continue
            End Do
1           Continue
         End Do
      End Do
*
*     Some printing for debug
*
      If (iPL.ge.4) Then
         Do iQM = 1, nAtQM
            Write(6,*) 'dB/dq_i for i = ',iQM
            Do jPnt = 1, nGrdPt
               Write(6,1234) jPnt,(DB(jPnt,iXYZ,iQM),iXYZ=1,3)
 1234       Format(I6,3D13.6)
            End Do
         End Do
      End If
      Return
      End
