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
*               Francesco Aquilante                                    *
************************************************************************
*  Cho_X_CheckDiag
*
*> @brief
*>   Check diagonal
*> @author Thomas Bondo Pedersen
*> @modified_by F. Aquilante (add ::OneCenter_ChkDiag)
*> @modified_by T.B. Pedersen (If ``(Cho_1Center)``: \p Err only contains 1-center errors on exit)
*>
*> @details
*> This routine reads and analyzes (histogram and statistics)
*> the exact integral diagonal, computes and analyzes the
*> diagonal from Cholesky vectors,
*> and the difference between the two (exact minus Cholesky).
*>
*> The statistics printed are: minimum value, maximum
*> value, mean value, mean absolute value, variance (wrt mean
*> value), and standard deviation (wrt mean value).
*>
*> On exit:
*>
*> - \p Err(1) = min error
*> - \p Err(2) = max error
*> - \p Err(3) = average error
*> - \p Err(4) = RMS error
*>
*> Return code is ``0`` if successful execution. If \p irc is non-zero,
*> the contents or \p Err are ill-defined.
*> Results will only be printed to output if \c iPrint is ``-5`` or
*> greater (\c iPrint stored in choprint.inc).
*>
*> @param[out] irc Return code
*> @param[out] Err min, max, average, and RMS error
************************************************************************
      SubRoutine Cho_X_CheckDiag(irc,Err)
      Implicit None
      Integer irc
      Real*8  Err(4)
#include "cholesky.fh"
#include "choprint.fh"
#include "WrkSpc.fh"

      Character*15 SecNam
      Parameter (SecNam = 'Cho_X_CheckDiag')

      Integer ip_XD, l_XD
      Integer ip_CD, l_CD
      Integer ip_Bin, l_Bin
      Integer ip_Stat, l_Stat

      Integer iPrThr
      Parameter (iPrThr=-5)

      Real*8   dDot_
      external ddot_

      Integer i
      Real*8  Stat
      Stat(i)=Work(ip_Stat-1+i)

C     Set return code.
C     ----------------

      irc = 0
      If (nnBstRT(1) .lt. 1) Then
         Call Cho_dZero(Err,4)
         Return
      End If

C     Allocations.
C     ------------

      l_XD = nnBstRT(1)
      l_CD = nnBstRT(1)
      l_Bin = 16
      l_Stat = 7
      Call GetMem('ExactDiag','Allo','Real',ip_XD,l_XD)
      Call GetMem('ChoDiag','Allo','Real',ip_CD,l_CD)
      Call GetMem('ChoBin','Allo','Real',ip_Bin,l_Bin)
      Call GetMem('Stat','Allo','Real',ip_Stat,l_Stat)

C     Set bins for histograms.
C     ------------------------

      Work(ip_Bin) = 1.0d0
      Do i = 1,l_Bin-1
         Work(ip_Bin+i) = Work(ip_Bin-1+i)*1.0d-1
      End Do

C     Read exact diagonal.
C     --------------------

      Call Cho_IODiag(Work(ip_XD),2)

C     Print histogram of exact diagonal and get statistics.
C     -----------------------------------------------------

      If (iPrint.ge.iPrThr) Then
         Call Cho_Head('Analysis of Exact Integral Diagonal','=',80,6)
         Call Cho_AnaSize(Work(ip_XD),l_XD,Work(ip_Bin),l_Bin,6)
         Call Statistics(Work(ip_XD),l_XD,Work(ip_Stat),1,2,3,4,5,6,7)
         Call Cho_PrtSt(Work(ip_XD),l_XD,Work(ip_Stat))
      End If

C     Calculate Cholesky diagonal.
C     ----------------------------

      Call Cho_X_CalcChoDiag(irc,Work(ip_CD))
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': Cho_X_CalcChoDiag returned ',irc
         Go To 1 ! return after dealloc
      End If

C     Print histogram of Cholesky diagonal and get statistics.
C     --------------------------------------------------------

      If (iPrint.ge.iPrThr) Then
         Call Cho_Head('Analysis of Cholesky Integral Diagonal','=',80,
     &                 6)
         Call Cho_AnaSize(Work(ip_CD),l_CD,Work(ip_Bin),l_Bin,6)
         Call Statistics(Work(ip_CD),l_CD,Work(ip_Stat),1,2,3,4,5,6,7)
         Call Cho_PrtSt(Work(ip_CD),l_CD,Work(ip_Stat))
      End If

C     Subtract Cholesky diagonal from exact diagonal.
C     -----------------------------------------------

      Call dAXPY_(nnBstRT(1),-1.0d0,Work(ip_CD),1,Work(ip_XD),1)

C     Print histogram of difference array and get statistics.
C     -------------------------------------------------------

      If (iPrint.ge.iPrThr) Then
         Call Cho_Head('Analysis of Difference (Exact-Cholesky)','=',80,
     &                 6)
         Call Cho_AnaSize(Work(ip_XD),l_XD,Work(ip_Bin),l_Bin,6)
      End If
      Call Statistics(Work(ip_XD),l_XD,Work(ip_Stat),1,2,3,4,5,6,7)
      If (iPrint.ge.iPrThr) Then
         Call Cho_PrtSt(Work(ip_XD),l_XD,Work(ip_Stat))
      End If

C     Set Err array.
C     --------------

      Err(1) = Stat(3)
      Err(2) = Stat(4)
      Err(3) = Stat(1)
      Err(4) = sqrt(dDot_(nnBstRT(1),Work(ip_XD),1,
     &                              Work(ip_XD),1)/dble(nnBstRT(1)))

      If (iPrint.ge.iPrThr) Then
         Write(6,'(/,1X,A,1P,D15.6)')
     &   'Minimum error   : ',Err(1)
         Write(6,'(1X,A,1P,D15.6)')
     &   'Maximum error   : ',Err(2)
         Write(6,'(1X,A,1P,D15.6)')
     &   'Average error   : ',Err(3)
         Write(6,'(1X,A,1P,D15.6)')
     &   'RMS error       : ',Err(4)
      End If

C     Error analysis for the 1-center diagonals only.
C     If this is a one-center calculation, use statistics from 1-center
C     diagonals only as elements of Err array.
C     -----------------------------------------------------------------

      If (nSym.eq.1) Then
         Call OneCenter_ChkDiag(Work(ip_XD),l_XD,Work(ip_Stat),
     &                          iPrint.ge.iPrThr)
         If (Cho_1Center) Then
            Err(1) = Stat(3)
            Err(2) = Stat(4)
            Err(3) = Stat(1)
            Err(4) = sqrt(dDot_(nnBstRT(1),Work(ip_XD),1,
     &                    Work(ip_XD),1)/dble(nnBstRT(1)))
         End If
      End If

C     Deallocations.
C     --------------

    1 Continue
      Call GetMem('Stat','Free','Real',ip_Stat,l_Stat)
      Call GetMem('ChoBin','Free','Real',ip_Bin,l_Bin)
      Call GetMem('ChoDiag','Free','Real',ip_CD,l_CD)
      Call GetMem('ExactDiag','Free','Real',ip_XD,l_XD)

      End
      SubRoutine Cho_PrtSt(Vec,lVec,Stat)
      Implicit None
      Integer lVec
      Real*8  Vec(lVec)
      Real*8  Stat(7)

      Real*8   dDot_
      external ddot_

      Write(6,'(/,1X,A,I15)')
     & 'No. of elements: ',lVec
      Write(6,'(1X,A,1P,D15.6)')
     & 'Frobenius norm : ',sqrt(dDot_(lVec,Vec,1,Vec,1))
      Write(6,'(1X,A,1P,D15.6)')
     & 'Minimum value  : ',Stat(3)
      Write(6,'(1X,A,1P,D15.6)')
     & 'Maximum value  : ',Stat(4)
      Write(6,'(1X,A,1P,D15.6)')
     & 'Mean value     : ',Stat(1)
      Write(6,'(1X,A,1P,D15.6)')
     & 'Mean abs. value: ',Stat(2)
      Write(6,'(1X,A,1P,D15.6)')
     & 'Max. abs. value: ',Stat(5)
      Write(6,'(1X,A,1P,D15.6)')
     & 'Biased variance: ',Stat(6)
      Write(6,'(1X,A,1P,D15.6,A)')
     & 'Standard dev.  : ',Stat(7),' (unbiased variance)'

      End
      Subroutine OneCenter_ChkDiag(Diag,l_D,Stat,DoPrint)
      use ChoArr, only: iRS2F
      Implicit Real*8 (a-h,o-z)
      Real*8 Diag(l_D), Stat(7)
      Logical DoPrint
#include "Molcas.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "choptr.fh"
      Character*(LENIN8)  Name(maxbfn)
      Character*(LENIN)  ctmp1, ctmp2
      Real*8 Err(4)

      Call Get_cArray('Unique Basis Names',Name,LENIN8*nBasT)

      Do krs=1,nnBstRT(1)
         ia = iRS2F(1,krs)
         ctmp1=Name(ia)(1:LENIN)
         ib = iRS2F(2,krs)
         ctmp2=Name(ib)(1:LENIN)
         If (ctmp1.ne.ctmp2) Diag(krs)=0.0d0
      End Do

      If (DoPrint) Then
         Call Cho_Head('Analysis of Difference (1-Center only)','=',80,
     &                 6)
      End If
      Call Statistics(Diag,l_D,Stat,1,2,3,4,5,6,7)
      If (DoPrint) Then
         Call Cho_PrtSt(Diag,l_D,Stat)
      End If

C     Set Err array.
C     --------------

      Err(1) = Stat(3)
      Err(2) = Stat(4)
      Err(3) = Stat(1)
      Err(4) = sqrt(dDot_(nnBstRT(1),Diag,1,
     &                              Diag,1)/dble(nnBstRT(1)))

      If (DoPrint) Then
         Write(6,'(/,1X,A,1P,D15.6)')
     &   'Minimum error   : ',Err(1)
         Write(6,'(1X,A,1P,D15.6)')
     &   'Maximum error   : ',Err(2)
         Write(6,'(1X,A,1P,D15.6)')
     &   'Average error   : ',Err(3)
         Write(6,'(1X,A,1P,D15.6)')
     &   'RMS error       : ',Err(4)
      End If

      Return
      End
