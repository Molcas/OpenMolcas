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
*  Cho_X_CalculateGMat
*
*> @brief
*>   Calculate Cholesky \f$ G \f$ matrix (metric matrix)
*> @author Thomas Bondo Pedersen
*>
*> @details
*> Calculate the metric matrix from Cholesky vectors (i.e. exact).
*>
*> \f[ G_{IJ} = \sum_{K=1}^{\min(I,J)} L_{IK} L_{JK} \f]
*>
*> The matrix is symmetric and stored on disk (file ``AVECXX``) in
*> triangular storage. The file is opened and closed in this
*> routine using routines ::DAName_MF_WA and ::DAClos, respectively
*> (i.e. \c FileName is opened as a word addressable multifile).
*> The calculation failed if \p irc is different from zero on exit.
*>
*> @note
*> This routine should *not* be used with DF.
*>
*> @param[out] irc Return code
************************************************************************
      SubRoutine Cho_X_CalculateGMat(irc)
      Implicit Real*8 (A-H,O-Z)
      Character*(6) FileName

#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

#if defined (_DEBUG_)
      real*8 ddot_
      external ddot_
#endif

      Logical isDF

      Parameter (N2 = InfVec_N2)
      InfVcT(i,j,k)=iWork(ip_InfVec_T-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
      InfVec(i,j,k)=iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
      NVT(i)=iWork(ip_NVT-1+i)
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
      iRS2RS(i)=iWork(ip_iRS2RS-1+i)

C     Set return code.
C     ----------------

      irc = 0

C     Refuse to perform calculation for DF.
C     -------------------------------------

      Call DecideOnDF(isDF)
      If (isDF) Then
         irc = -1
         Return
      End If


C     Scratch location in index arrays.
C     ---------------------------------

      iLoc = 3 ! do NOT change (used implicitly by reading routine)

C     Create memory pointer for flushing memory.
C     ------------------------------------------

      l_Flush = 1
      Call GetMem('FLUSH','Allo','Inte',ip_Flush,l_Flush)

C     Get pointer to InfVec array for all vectors (needed for parallel
C     runs) and the total number of vectors.
C     ----------------------------------------------------------------

      l_NVT = nSym
      Call GetMem('NVT','Allo','Inte',ip_NVT,l_NVT)
      Call Cho_CGM_InfVec(ip_InfVec_T,iWork(ip_NVT),l_NVT)

C     Copy rs1 to location 2.
C     -----------------------

      Call Cho_X_RSCopy(irc,1,2)
      If (irc .ne. 0) Then
         irc = 1
         Go To 1 ! exit after deallocation
      End If

C     Calculate triangular G matrix.
C     G(IJ)=sum_K L(I,K)*L(J,K).
C     ------------------------------

      iRedC = -1
      Do iSym = 1,nSym

C        Open file.
C        ----------

         Write (FileName,'(A4,I2.2)') 'AVEC',iSym-1
         lUnit = 7
         Call DAName_MF_WA(lUnit,FileName)
         iDisk=0

         l_iRS2RS = nnBstR(iSym,1)
         l_G = NVT(iSym)*(NVT(iSym)+1)/2
         Call GetMem('RS-TO-RS','Allo','Inte',ip_iRS2RS,l_iRS2RS)
         Call GetMem('G','Allo','Real',ip_G,l_G)
         Call GetMem('MX','Max ','Real',ip_Wrk,l_Wrk)
         Call GetMem('Wrk','Allo','Real',ip_Wrk,l_Wrk)
         Call fZero(Work(ip_G),l_G)
         Call iZero(iWork(ip_iRS2RS),l_iRS2RS)
         idRS2RS = -2
         kOffG = ip_G - 1
         KK1 = 1
         Do While (KK1 .le. NumCho(iSym))
            nVRead = 0
            mUsed = 0
            Call Cho_X_VecRd(Work(ip_Wrk),l_Wrk,KK1,NumCho(iSym),iSym,
     &                       nVRead,iRedC,mUsed)
            If (nVRead .lt. 1) Then
               irc = 2
               Go To 1 ! exit after deallocation
            End If
            kOffV = ip_Wrk - 1
            Do KKK = 0,nVRead-1
               KK = KK1 + KKK
               iRed = InfVec(KK,2,iSym)
               If (iRedC .ne. iRed) Then
                  Call Cho_X_SetRed(irc,iLoc,iRed)
                  If (irc .ne. 0) Then
                     irc = 3
                     Go To 1 ! exit after deallocation
                  End If
                  iRedC = iRed
               End If
               If (idRS2RS .ne. iRedC) Then
                  Call Cho_RS2RS(iWork(ip_iRS2RS),l_iRS2RS,
     &                           2,iLoc,iRedC,iSym)
                  idRS2RS = iRedC
               End If
               K = InfVec(KK,5,iSym)
               Do J = K,NVT(iSym)
                  iJ = iRS2RS(InfVcT(J,1,iSym)-iiBstR(iSym,1))
                  V_J = Work(kOffV+iJ)
                  Do I = K,J
                     iI = iRS2RS(InfVcT(I,1,iSym)-iiBstR(iSym,1))
                     kG_IJ = kOffG + iTri(I,J)
                     Work(kG_IJ) = Work(kG_IJ) + Work(kOffV+iI)*V_J
                  End Do
               End Do
               kOffV = kOffV + nnBstR(iSym,iLoc)
            End Do
            KK1 = KK1 + nVRead
         End Do
         Call Cho_GADGOP(Work(ip_G),l_G,'+')
         iOpt = 1
         Call DDAFile(lUnit,iOpt,Work(ip_G),l_G,iDisk)
#if defined (_DEBUG_)
         Call TriPrt('G-matix',' ',Work(ip_G),NVT(iSym))
         Write(6,'(A,I2,A,1P,D16.7)')
     &   'G matrix, sym.',iSym,': Norm = ',
     &   sqrt(dDot_(l_G,Work(ip_G),1,Work(ip_G),1))
#endif
         Call GetMem('Wrk','Free','Real',ip_Wrk,l_Wrk)
         Call GetMem('G','Free','Real',ip_G,l_G)
         Call GetMem('RS-TO-RS','Free','Inte',ip_iRS2RS,l_iRS2RS)

C        Close file
C        ----------
         Call DAClos(lUnit)
      End Do

C     Flush memory, close file, and exit.
C     -----------------------------------

    1 Continue
      Call GetMem('FLUSH','Flush','Inte',ip_Flush,l_Flush)
      Call GetMem('FLUSH','Free','Inte',ip_Flush,l_Flush)

      End
C
C**********************************************************************C
C
      SubRoutine Cho_CGM_InfVec(ip_InfVec_T,NVT,n)
      Implicit None
      Integer ip_InfVec_T
      Integer n
      Integer NVT(n)

      Call Cho_X_GetIP_InfVec(ip_InfVec_T)
      Call Cho_X_GetTotV(NVT,n)

      End
