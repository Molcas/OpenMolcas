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
*               2020,2021, Roland Lindh                                *
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
      use ChoSwp, only: InfVec
      Implicit Real*8 (A-H,O-Z)
      Character*(6) FileName

#include "real.fh"
#include "cholesky.fh"
#include "stdalloc.fh"

#if defined (_DEBUGPRINT_)
      real*8 ddot_
      external ddot_
#endif

      Logical isDF
      Integer, Pointer:: InfVcT(:,:,:)
      Integer, Allocatable:: NVT(:), iRS2RS(:)
      Real*8, Allocatable:: Wrk(:), G(:)
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      SubRoutine Cho_CGM_InfVec(InfVcT,NVT,n)
      Integer, Pointer:: InfVcT(:,:,:)
      Integer :: n, NVT(n)
      End SubRoutine Cho_CGM_InfVec
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

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

C     Get pointer to InfVec array for all vectors (needed for parallel
C     runs) and the total number of vectors.
C     ----------------------------------------------------------------

      Call mma_allocate(NVT,nSym,Label='NVT')
      Call Cho_CGM_InfVec(InfVcT,NVT,SIZE(NVT))

C     Copy rs1 to location 2.
C     -----------------------

      Call Cho_X_RSCopy(irc,1,2)
      If (irc .ne. 0) Then
         irc = 1
         Go To 1 ! exit
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
         Call mma_allocate(iRS2RS,l_iRS2RS,Label="iRS2RS")
         iRS2RS(:)=0
         l_G = NVT(iSym)*(NVT(iSym)+1)/2
         Call mma_allocate(G,l_G,Label='G')
         Call mma_MaxDBLE(l_Wrk)
         Call mma_allocate(Wrk,l_Wrk,Label='Wrk')
         Wrk(:)=Zero
         G(:)=Zero

         idRS2RS = -2
         KK1 = 1
         Do While (KK1 .le. NumCho(iSym))
            nVRead = 0
            mUsed = 0
            Call Cho_X_VecRd(Wrk,SIZE(Wrk),KK1,NumCho(iSym),iSym,
     &                       nVRead,iRedC,mUsed)
            If (nVRead .lt. 1) Then
               irc = 2
               Go To 1 ! exit after deallocation
            End If
            kOffV = 0
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
                  Call Cho_RS2RS(iRS2RS,SIZE(iRS2RS),2,iLoc,iRedC,iSym)
                  idRS2RS = iRedC
               End If
               K = InfVec(KK,5,iSym)
               Do J = K,NVT(iSym)
                  iJ = iRS2RS(InfVcT(J,1,iSym)-iiBstR(iSym,1))
                  V_J = Wrk(kOffV+iJ)
                  Do I = K,J
                     iI = iRS2RS(InfVcT(I,1,iSym)-iiBstR(iSym,1))
                     kG_IJ = iTri(I,J)
                     G(kG_IJ) = G(kG_IJ) + Wrk(kOffV+iI)*V_J
                  End Do
               End Do
               kOffV = kOffV + nnBstR(iSym,iLoc)
            End Do
            KK1 = KK1 + nVRead
         End Do
         Call Cho_GADGOP(G,SIZE(G),'+')
         iOpt = 1
         Call DDAFile(lUnit,iOpt,G,SIZE(G),iDisk)
#if defined (_DEBUGPRINT_)
         Call TriPrt('G-matix',' ',G,NVT(iSym))
         Write(6,'(A,I2,A,1P,D16.7)')
     &   'G matrix, sym.',iSym,': Norm = ',
     &   sqrt(dDot_(SIZE(G),G,1,G,1))
#endif
         Call mma_deallocate(Wrk)
         Call mma_deallocate(G)
         Call mma_deallocate(iRS2RS)

C        Close file
C        ----------
         Call DAClos(lUnit)
      End Do

    1 Continue
      Call mma_deallocate(NVT)

      End
C
C**********************************************************************C
C
      SubRoutine Cho_CGM_InfVec(InfVcT,NVT,n)
      Implicit None
      Integer, Pointer:: InfVcT(:,:,:)
      Integer n
      Integer NVT(n)
*                                                                      *
************************************************************************
*                                                                      *
      Interface
      Subroutine Cho_X_GetIP_InfVec(InfVcT)
      Integer, Pointer:: InfVct(:,:,:)
      End Subroutine Cho_X_GetIP_InfVec
      End Interface
*                                                                      *
************************************************************************
*                                                                      *
      Call Cho_X_GetIP_InfVec(InfVcT)
      Call Cho_X_GetTotV(NVT,n)

      End
