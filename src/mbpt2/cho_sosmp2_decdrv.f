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
* Copyright (C) 2007, Francesco Aquilante                              *
************************************************************************
      SubRoutine Cho_SOSmp2_DecDrv(irc,DelOrig,Diag)
C
C     Francesco Aquilante, May 2007.
C
C     Purpose: decompose M(ai,bj) = (ai|bj)^2 for use in
C              SOS-MP2 approach.
C
C     DelOrig: flag for deleting files with original vectors after
C              decomposition completes.
C
#include "implicit.fh"
      External Cho_SOSmp2_Col, ChoMP2_Vec
      Integer  irc
      Logical  DelOrig
      Real*8   Diag(*)
#include "cholesky.fh"
#include "chomp2_cfg.fh"
#include "chomp2.fh"
#include "chomp2_dec.fh"
#include "WrkSpc.fh"

      Character*6  ThisNm
      Character*17 SecNam
      Parameter (SecNam = 'Cho_SOSmp2_DecDrv', ThisNm = 'DecDrv')

      Logical Restart, Failed
      Parameter (Restart = .false.)

      Integer iClos(2)

C     Initializations.
C     ----------------

      irc = 0

      lErrStat = 3
      Call GetMem('ErrStat','Allo','Real',ipErrStat,lErrStat)
      If (Verbose) Then
         nBin = 18
         Call GetMem('Bin','Allo','Real',ipBin,nBin)
      Else
         ipBin = -999999
         nBin  = 0
      End If

      Do iSym = 1,nSym
         nMP2Vec(iSym) = 0
         InCore(iSym)  = .false.
      End Do
      NowSym    = -999999
      ip_OldVec = -999999
      l_OldVec  = 0

      If (DelOrig) Then
         iClos(1) = 3  ! signals close and delete original vectors
      Else
         iClos(1) = 2  ! signals close and keep original vectors
      End If
      iClos(2) = 2     ! signals close and keep result vectors


C     Print.
C     ------

      If (Verbose) Then
         Write(6,*)
         Call Cho_Head('Cholesky Decomposition of  M(ai,bj) = (ai|bj)^2'
     &                 //' for SOS-MP2','=',80,6)
         Write(6,'(/,1X,A)') 'Configuration of decomposition:'
         Write(6,'(1X,A,1P,D15.6)') 'Threshold: ',ThrMP2
         Write(6,'(1X,A,1P,D15.6)') 'Span     : ',SpanMP2
         If (ChkDecoMP2) Then
            Write(6,'(1X,A)') 'Full decomposition check activated.'
         End If
      End If

C     Start symmetry loop.
C     --------------------

      kOffD = 1
      Do iSym = 1,nSym

         nDim = nT1am(iSym)
         If (nDim.gt.0 .and. NumCho(iSym).gt.0) Then

            If (Verbose .and. nBin.gt.0) Then
               Work(ipBin) = 1.0D2
               Do iBin = ipBin+1,ipBin+nBin-1
                  Work(iBin) = Work(iBin-1)*1.0D-1
               End Do
               Write(6,'(//,1X,A,I2,A,I9)')
     &         '>>> Cholesky decomposing symmetry block ',iSym,
     &         ', dimension: ',nDim
               Write(6,'(/,1X,A)') 'Analysis of initial diagonal:'
               Call Cho_AnaSize(Diag(kOffD),nDim,Work(ipBin),nBin,6)
            End If

C           Open files.
C           -----------

            Do iTyp = 1,2
               Call ChoMP2_OpenF(1,iTyp,iSym)
            End Do

C           Setup decomposition.
C           --------------------

            NowSym = iSym

            If (MxQualMP2 .ne. MxQual_Def) Then ! user-defined
               MxQual = min(max(MxQualMP2,1),nDim)
            Else ! default
               If (nDim .gt. 10) Then
                  MxQual = max(min(nDim/10,MxQualMP2),1)
               Else
                  MxQual = max(min(nDim,MxQualMP2),1)
               End If
            End If
#if !defined (_I8_)
            lTstBuf = (nDim+MxQual)*MxQual
            lTstQua = nDim*(MxQual+1)
            Do While ((lTstBuf.lt.0.or.lTstQua.lt.0) .and. MxQual.gt.0)
               MxQual  = MxQual - 1
               lTstBuf = (nDim+MxQual)*MxQual
               lTstQua = nDim*(MxQual+1)
            End Do
            If (MxQual .lt. 1) Then
               Write(6,*) SecNam,': MxQual causes integer overflow!'
               Write(6,*) SecNam,': parameters:'
               Write(6,*) 'Symmetry block: ',iSym
               Write(6,*) 'Dimension     : ',nDim
               Write(6,*) 'MxQual        : ',MxQual
               irc = -99
               Go To 1 ! exit
            End If
#endif

            lQual   = nDim*(MxQual + 1)
            liQual  = MxQual
            liPivot = nDim
            Call GetMem('Qual','Allo','Real',ipQual,lQual)
            Call GetMem('iQual','Allo','Inte',ipiQual,liQual)
            Call GetMem('iPivot','Allo','Inte',ipiPivot,liPivot)

            Call GetMem('GetMax','Max ','Real',ipB,lB)
            lBuf = min((nDim+MxQual)*MxQual,lB)
            Left = lB - lBuf
            nInC = Left/nDim
            If (nInC .ge. NumCho(iSym)) Then
               InCore(iSym) = .true.
               l_OldVec = nDim*NumCho(iSym)
               Call GetMem('OldVec','Allo','Real',ip_OldVec,l_OldVec)
               iOpt = 2
               lTot = l_OldVec
               iAdr = 1
               Call ddaFile(lUnit_F(iSym,1),iOpt,
     &                      Work(ip_OldVec),lTot,iAdr)
            End If
            Call GetMem('GetMx2','Max ','Real',ipB,lBuf)
            Call GetMem('DecBuf','Allo','Real',ipBuf,lBuf)

C           Decompose this symmetry block.
C           ------------------------------

            Thr  = ThrMP2
            Span = SpanMP2
            Call ChoDec(Cho_SOSmp2_Col,ChoMP2_Vec,Restart,Thr,Span,
     &                  MxQual,Diag(kOffD),Work(ipQual),Work(ipBuf),
     &                  iWork(ipiPivot),iWork(ipiQual),nDim,lBuf,
     &                  Work(ipErrStat),nMP2Vec(iSym),irc)
            If (irc .ne. 0) Then
               Write(6,*) SecNam,': ChoDec returned ',irc,
     &                           '   Symmetry block: ',iSym
               Go To 1 ! exit...
            End If
            XMn = Work(ipErrStat)
            XMx = Work(ipErrStat+1)
            RMS = Work(ipErrStat+2)
            If (Verbose) Then
               Write(6,'(/,1X,A)')
     &         '- decomposition completed!'
               Write(6,'(1X,A,I9,A,I9,A)')
     &         'Number of vectors needed: ',nMP2Vec(iSym),
     &         ' (number of AO vectors: ',NumCho(iSym),')'
               Write(6,'(1X,A)')
     &         'Error statistics for (ai|ai)^2 [min,max,rms]:'
               Write(6,'(1X,1P,3(D15.6,1X))') XMn,XMx,RMS
            End If
            Failed = abs(Xmn).gt.Thr .or. abs(XMx).gt.thr .or.
     &               RMS.gt.Thr
            If (Failed) Then
               If (.not. Verbose) Then
                  Write(6,'(1X,A)')
     &            'Error statistics for (ai|ai)^2 [min,max,rms]:'
                  Write(6,'(1X,1P,3(D15.6,1X))') XMn,XMx,RMS
               End If
               Write(6,*) SecNam,': (ai|bj)^2 decomposition failed!'
               irc = -9999
               Go To 1 ! exit
            End If

C           If requested, check decomposition.
C           ----------------------------------

            If (ChkDecoMP2) Then
               Write(6,*)
               Write(6,*)SecNam,': Checking M(ai,bj)=(ai|bj)^2 CD.'
               Write(6,*)'Symmetry block: ',iSym
               Write(6,*)'Threshold, Span, MxQual: ',Thr,Span,MxQual
               Write(6,*)'Error statistics for (ai|ai)^2 [min,max,rms]:'
               Write(6,*) (Work(ipErrStat+i),i=0,2)
               Call Cho_SOSmp2_DecChk(irc,iSym,Work(ipQual),nDim,MxQual,
     &                            Work(ipBuf),lBuf,Work(ipErrStat))
               If (irc .ne. 0) Then
                Write(6,*) SecNam,': ChoMP2_DecChk returned ',irc,
     &                            '   Symmetry block: ',iSym
                Call ChoMP2_Quit(SecNam,'SOS-MP2 decomposition failed!',
     &                             ' ')
               Else
                  XMn = Work(ipErrStat)
                  XMx = Work(ipErrStat+1)
                  RMS = Work(ipErrStat+2)
                  Failed = Failed .or. abs(Xmn).gt.Thr .or.
     &                     abs(XMx).gt.Thr .or. RMS.gt.Thr
                  Write(6,*)
     &            'Error statistics for (ai|bj)^2 [min,max,rms]:'
                  Write(6,*) XMn,XMx,RMS
                  If (Failed) Then
                     Write(6,*) '==> DECOMPOSITION FAILURE <=='
                     irc = -9999
                     Go To 1 ! exit
                  Else
                     Write(6,*) '==> DECOMPOSITION SUCCESS <=='
                  End If
                  Call xFlush(6)
               End If
            End If

C           Free memory.
C           ------------

            Call GetMem('DecBuf','Free','Real',ipBuf,lBuf)
            If (InCore(iSym)) Then
               Call GetMem('OldVec','Free','Real',ip_OldVec,l_OldVec)
               ip_OldVec = -999999
               l_OldVec  = 0
            End If
            Call GetMem('iPivot','Free','Inte',ipiPivot,liPivot)
            Call GetMem('iQual','Free','Inte',ipiQual,liQual)
            Call GetMem('Qual','Free','Real',ipQual,lQual)

C           Close (possibly deleting original) files.
C           -----------------------------------------

            Do iTyp = 1,2
               Call ChoMP2_OpenF(iClos(iTyp),iTyp,iSym)
            End Do

C           Update pointer to diagonal block.
C           ---------------------------------

            kOffD = kOffD + nT1am(iSym)

         Else

            If (Verbose) Then
               Write(6,'(//,1X,A,I2,A)')
     &         '>>> Symmetry block',iSym,' is empty!'
            End If

         End If

      End Do

    1 If (irc .ne. 0) Then ! make sure files are closed before exit
         Do iSym = 1,nSym
            Do iTyp = 1,2
               Call ChoMP2_OpenF(2,iTyp,iSym)
            End Do
         End Do
      End If
      Call GetMem('Flush','Flush','Real',ipErrStat,lErrStat)
      Call GetMem('ErrStat','Free','Real',ipErrStat,lErrStat)
      End
