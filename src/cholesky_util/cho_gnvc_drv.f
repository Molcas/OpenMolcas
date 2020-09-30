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
      SubRoutine Cho_GnVc_Drv(irc,Diag)
C
C     Purpose: generate vectors from a known "map" and diagonal.
C              First reduced set must have been set up, which is
C              reasonable, since it is naturally done along with the
C              diagonal.
C
#include "implicit.fh"
      Real*8 Diag(*)
#include "cholesky.fh"
#include "choprint.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      Character*12 SecNam
      Parameter (SecNam = 'Cho_GnVc_Drv')

      Character*2  Unt
      Character*26 String

      Integer LastRed(8), nScrV(8)
      Integer ip_mapRS2RS, l_mapRS2RS
      Common / GnVcMp / ip_mapRS2RS(8), l_mapRS2RS(8)

      Parameter (N2 = InfVec_N2)

      InfVec(i,j,k)=iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
      nVecRS(i,j)=iWork(ip_nVecRS-1+nSym*(j-1)+i)
      iVecRS(i,j)=iWork(ip_iVecRS-1+nSym*(j-1)+i)
      mapRS2RS(i,j)=iWork(ip_mapRS2RS(i)-1+j)

#if defined (_DEBUG_)
#endif

C     Start timing.
C     -------------

      Call Cho_Timer(tCPU1,tWall1)

C     Set return code.
C     ----------------

      irc = 0

C     Make a dummy allocation.
C     ------------------------

      l_GnVcDum = 1
      Call Cho_Mem('GnVcDum','Allo','Real',ip_GnVcDum,l_GnVcDum)

C     Copy first reduced set to current reduced set stored at location
C     2.
C     ----------------------------------------------------------------

      Call Cho_X_RSCopy(irc,1,2)
      If (irc .ne. 0) Then
         Write(Lupri,*) SecNam,': Cho_X_RSCopy returned ',irc
         Go To 100 ! exit
      End If

C     Allocate memory for shell pair list.
C     ------------------------------------

      l_ListSP = nnShl
      Call Cho_Mem('ListSP','Allo','Inte',ip_ListSP,l_ListSP)

C     Set number of integral passes (= number of reduced sets).
C     Set ID of last reduced set in each symmetry.
C     ---------------------------------------------------------

      nPass = 0
      Do iSym = 1,nSym
         If (NumCho(iSym) .lt. 1) Then
            LastRed(iSym) = 0
         Else
            LastRed(iSym) = InfVec(NumCho(iSym),2,iSym)
         End If
         nPass = max(nPass,LastRed(iSym))
      End Do
      If (nPass .lt. 1) Then
         Call Cho_Quit('nPass is non-positive in '//SecNam,102)
      Else If (nPass .ne. XnPass) Then
         Call Cho_Quit('nPass != XnPass in '//SecNam,102)
      End If

C     nVecRS(iSym,jRed): #vectors of sym. iSym, reduced set jRed.
C     iVecRS(iSym,jRed): 1st vec. of sym. iSym, reduced set jRed.
C                        (0 means no vectors!!)
C     -----------------------------------------------------------

      l_nVecRS = nSym*nPass
      l_iVecRS = l_nVecRS
      Call Cho_Mem('nVecRS','Allo','Inte',ip_nVecRS,l_nVecRS)
      Call Cho_Mem('iVecRS','Allo','Inte',ip_iVecRS,l_iVecRS)
      Call Cho_iZero(iWork(ip_nVecRS),l_nVecRS)
      Call Cho_iZero(iWork(ip_iVecRS),l_iVecRS)
      Do iSym = 1,nSym
         nTotVec = 0
         Do jRed = 1,LastRed(iSym)
            Call Cho_X_nVecRS(jRed,iSym,iVec1,nVec)
            If (nVec.lt.0 .or. iVec1.lt.0) Then
               Call Cho_Quit(SecNam
     &                       //'Cho_X_nVecRS returned negative nVec',
     &                       103)
            Else
               iWork(ip_nVecRS-1+nSym*(jRed-1)+iSym) = nVec
               iWork(ip_iVecRS-1+nSym*(jRed-1)+iSym) = iVec1
            End If
            nTotVec = nTotVec + nVecRS(iSym,jRed)
         End Do
         If (nTotVec .ne. NumCho(iSym)) Then
            Call Cho_Quit('Initialization error in '//SecNam,102)
         End If
      End Do

C     Allocate RS-to-RS mapping array.
C     --------------------------------

      Do iSym = 1,nSym
         l_mapRS2RS(iSym) = nnBstR(iSym,1)
         Call Cho_Mem('mapRS2RS','Allo','Inte',
     &                ip_mapRS2RS(iSym),l_mapRS2RS(iSym))
      End Do

C     Split remaining memory.
C     -----------------------

      Call Cho_Mem('GnVcMx','GetM','Real',ip_Wrk,l_WrkT)
      If (l_WrkT .lt. 2) Then
         Write(Lupri,*) SecNam,': max. allocatable memory is ',
     &                  l_WrkT
         Write(Lupri,*) 'Please increase available memory.'
         irc = 101
         Go To 100 ! exit
      Else
         l_Wrk = l_WrkT/2
      End If

#if defined (_DEBUG_)
C     Debug: force batching.
C     ----------------------

      lNdMx = 0
      Do iPass = 1,nPass
         lNeed = nnBstR(1,2)*nVecRS(1,iPass)
         Do iSym = 2,nSym
            lNeed = lNeed + nnBstR(iSym,2)*nVecRS(iSym,iPass)
         End Do
         lNdMx = max(lNdMx,lNeed)
      End Do
      l_Wrk = min(l_Wrk,lNdMx)
#endif

C     Reinitialize vector counters.
C     -----------------------------

      Call Cho_iZero(NumCho,nSym)
      NumChT = 0

C     Start batch loop over integral passes.
C     --------------------------------------

      iPass = 0
      nBatch= 0
      Do While (iPass .lt. nPass)

         If (iPrint .ge. INF_PASS) Call Cho_Timer(TlTot1,WlTot1)

C        Update batch counter.
C        ---------------------

         nBatch = nBatch + 1
         iPass1 = iPass + 1
         If (nBatch .gt. nPass) Then
            Write(Lupri,*) SecNam,': batch counter exceeds ',
     &      'pass counter: ',nBatch,' > ',nPass
            irc = 103
            Go To 100 ! exit
         End If

C        Set up this batch based on current reduced set.
C        NumPass: the number of passes treated in this batch.
C        --------------------------------------------------------

         jRed    = iPass
         NumPass = 0
         l_Int   = 0
         lThis   = 0
         Do While (jRed .lt. nPass)
            jRed  = jRed + 1
            lThis = nnBstR(1,2)*nVecRS(1,jRed)
            Do iSym = 2,nSym
               lThis = lThis + nnBstR(iSym,2)*nVecRS(iSym,jRed)
            End Do
            l_Int = l_Int + lThis
            If (l_Int .gt. l_Wrk) Then
               l_Int = l_Int - lThis ! reset memory need
               jRed  = nPass ! break loop
            Else
               NumPass = NumPass + 1
            End If
         End Do
         If (NumPass .lt. 1) Then
            Write(Lupri,*) SecNam,': insufficient memory for batch ',
     &                     nBatch
            Write(Lupri,*) SecNam,': at least  ',lThis,
     &                     ' double precision words needed.'
            Write(Lupri,*) SecNam,': available ',l_Wrk,
     &                     ' double precision words.'
            irc = 101
            Go To 100 ! exit
         End If

C        Print.
C        ------

         If (iPrint .ge. INF_PASS) Then
            Write(String,'(A19,I7)') 'Integral Pass Batch',nBatch
            Call Cho_Head(String,'*',80,Lupri)
            Write(Lupri,'(A,I5,A,I5,A,I5,A)')
     &      'Integral passes treated:',iPass1,' to',iPass+NumPass,
     &      ' of',nPass,' passes.'
            Call Cho_Word2Byte(l_WrkT,8,dl_WrkT,Unt)
            Write(Lupri,'(A,I10,A,F10.3,A,A)')
     &      'Total memory available           : ',l_WrkT,
     &      ' 8-byte words; ',dl_WrkT,' ',Unt
            Call Cho_Word2Byte(l_Int,8,dl_Int,Unt)
            Write(Lupri,'(A,I10,A,F10.3,A,A)')
     &      'Memory used for integrals/vectors: ',l_Int,
     &      ' 8-byte words; ',dl_Int,' ',Unt
            Call Cho_iZero(nScrV,nSym)
            Do i = iPass1,iPass+NumPass
               Do iSym = 1,nSym
                  nScrV(iSym) = nScrV(iSym) + nVecRS(iSym,i)
               End Do
            End Do
            Write(Lupri,'(A,8I8)')
     &      '#vec. gener.  : ',(nScrV(i),i=1,nSym)
            Call Cho_Flush(Lupri)
         End If

C        Allocate memory for integral columns and initialize.
C        ----------------------------------------------------

         Call Cho_Mem('GnVc.Int','Allo','Real',ip_Int,l_Int)
         Call Cho_dZero(Work(ip_Int),l_Int)

C        Set up map from first reduced set to reduced set iPass1.
C        --------------------------------------------------------

         irc = 0
         Call Cho_X_RSCopy(irc,1,3)
         If (irc .ne. 0) Then
            Call Cho_Quit(SecNam
     &                    //': non-zero return code from Cho_X_RSCopy',
     &                    104)
         End If
         Do iSym = 1,nSym
            Call Cho_RS2RS(iWork(ip_mapRS2RS(iSym)),l_mapRS2RS(iSym),
     &                     3,2,iPass1,iSym)
         End Do

C        Set up "qualified column" index arrays.
C        nQual(iSym): #qualifieds in symmetry iSym
C        iQuAB(iAB,iSym): addr of qualified iAB, sym. iSym in curr.
C                         reduced set.
C        ----------------------------------------------------------

         Call Cho_iZero(nQual,nSym)
         iPass2 = iPass1 + NumPass - 1
         Do jPass = iPass1,iPass2
            Do iSym = 1,nSym
               kOff0 = ip_iQuAB - 1 + MaxQual*(iSym-1)
               iV1 = iVecRS(iSym,jPass)
               iV2 = iV1 + nVecRS(iSym,jPass) - 1
               Do iV = iV1,iV2
                  iAB = InfVec(iV,1,iSym) ! addr in 1st red. set
                  jAB = iAB - iiBstR(iSym,1)
                  kAB = mapRS2RS(iSym,jAB) ! addr in curr. red. set
#if defined (_DEBUG_)
                  If (kAB.lt.1 .or. kAB.gt.nnBstR(iSym,2)) Then
                     Write(Lupri,*) SecNam,': illegal kAB = ',kAB
                     Write(Lupri,*) 'Vector, symmetry, pass: ',
     &                              IV,iSym,jPass
                     Write(Lupri,*) 'Allowed range: 1 - ',nnBstR(iSym,2)
                     Call Cho_Quit('Index error in '//SecNam,104)
                  End If
#endif
                  nQual(iSym) = nQual(iSym) + 1
                  iWork(kOff0+nQual(iSym)) = iiBstR(iSym,2) + kAB
               End Do
            End Do
         End Do

         iOff_Col(1) = 0
         NumInt = nnBstR(1,2)*nQual(1)
         Do iSym = 2,nSym
            iOff_Col(iSym) = NumInt
            NumInt = NumInt + nnBstR(iSym,2)*nQual(iSym)
         End Do
         If (l_Int .ne. NumInt) Then
            Call Cho_Quit('Integral memory error in '//SecNam,101)
         End If

C        Compute all integrals needed for NumPass passes.
C        ------------------------------------------------

         If (iPrint .ge. INF_PASS) Call Cho_Timer(TlInt1,WlInt1)
         NumSP = 0
         Call Cho_GnVc_GetInt(Work(ip_Int),l_Int,
     &                        iWork(ip_nVecRS),iWork(ip_iVecRS),
     &                        iWork(ip_ListSP),
     &                        nSym,nPass,nnShl,iPass1,NumPass,NumSP)
         If (NumSP .lt. 1) Then
            Call Cho_Quit('No shell pairs calculated!',104)
         End If
         If (iPrint .ge. INF_PASS) Call Cho_Timer(TlInt2,WlInt2)

C        Generate vectors.
C        -----------------

         If (iPrint .ge. INF_PASS) Call Cho_Timer(TlDec1,WlDec1)
         Call Cho_GnVc_GenVec(Diag,Work(ip_Int),l_Int,
     &                        iWork(ip_nVecRS),iWork(ip_iVecRS),
     &                        nSym,nPass,iPass1,NumPass)
         If (iPrint .ge. INF_PASS) Call Cho_Timer(TlDec2,WlDec2)

C        Deallocate memory.
C        ------------------

         Call Cho_Mem('GnVc.Int','Free','Real',ip_Int,l_Int)

C        Print timing for this batch.
C        ----------------------------

         If (iPrint .ge. INF_PASS) Then
            TlInt = TlInt2 - TlInt1
            WlInt = WlInt2 - WlInt1
            TlDec = TlDec2 - TlDec1
            WlDec = WlDec2 - WlDec1
            Call Cho_Timer(TlTot2,WlTot2)
            TlTot = TlTot2 - TlTot1
            WlTot = WlTot2 - WlTot1
            Write(Lupri,'(/,A,I7,A)')
     &      'Overall timings for integral pass batch',nBatch,
     &      ' (CPU/Wall in seconds):'
            Write(Lupri,'(A,F12.2,1X,F12.2)')
     &      'Integrals    : ',TlInt,WlInt
            Write(Lupri,'(A,F12.2,1X,F12.2)')
     &      'Decomposition: ',TlDec,WlDec
            Write(Lupri,'(A,F12.2,1X,F12.2)')
     &      'Total        : ',TlTot,WlTot
            Write(Lupri,'(A,I7,A,I7,A)')
     &      'Integral passes treated:',iPass1,' to',iPass1-1+NumPass
            Call Cho_Flush(Lupri)
         End If

C        Update pass counter.
C        --------------------

         iPass = iPass + NumPass

      End Do ! integral pass loop

C     Exit after deallocating memory (by flushing).
C     ---------------------------------------------

  100 Continue
      Did_DecDrv = .true.
      Call Cho_Mem('GnVcDum','Flush','Real',ip_GnVcDum,l_GnVcDum)
      Call Cho_Timer(tCPU2,tWall2)
      tDecDrv(1) = tDecDrv(1) + tCPU2  - tCPU1
      tDecDrv(2) = tDecDrv(2) + tWall2 - tWall1
#if defined (_DEBUG_)
#endif
      Return

      End
