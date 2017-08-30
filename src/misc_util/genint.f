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
* Copyright (C) 2004, Francesco Aquilante                              *
************************************************************************
*  GEN_INT
*
*> @brief
*>   Generates integrals from Cholesky vectors
*> @author F. Aquilante, May 2004
*> @modified_by F. Aquilante, Sep. 2004
*>
*> @note
*> The transposition ``L(ab,J)`` &rarr; ``L(ba,J)`` of the vectors
*> ``(syma.ne.symb)`` is necessary because the calling routine
*> requires the integrals in the order \f$ (sr|qp) \f$ which
*> is reversed compared to the order of the symmetries
*> given as input arguments.
*>
*> @param[out] rc
*> @param[in]  iSymp
*> @param[in]  iSymq
*> @param[in]  iSymr
*> @param[in]  iSyms
*> @param[in]  ipq1
*> @param[in]  numpq
*> @param[out] Xint
************************************************************************
      SUBROUTINE GEN_INT(rc,iSymp,iSymq,iSymr,iSyms,ipq1,numpq,Xint)
************************************************************************
*
*   Modified  September 2004
*   Reason:
*   the transposition L(ab,J) --> L(ba,J) of the vectors
*   (syma.ne.symb) is necessary because the calling routine
*   requires the integrals in the order (sr|qp) which
*   is reversed compared to the order of the symmetries
*   given as input arguments
*
************************************************************************

      Implicit Real*8 (a-h,o-z)
      INTEGER   rc
      INTEGER   iSymp,iSymq,iSymr,iSyms
      INTEGER   pq,pq2,numpq,pq1_save
      Real*8    Xint(*)

#include "RdOrd.fh"

#include "WrkSpc.fh"
#include "TwoRc.fh"

C *************************************
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
C *************************************

      jSym = MulD2h(iSymp,iSymq)

      If (NumCho(jSym) .lt. 1) Return

C --- save the value of pq1 because it belongs to a Common block
      pq1_save = pq1
      pq1      = ipq1

      pq2 = pq1 + numpq - 1
      If (iSymp .eq. iSymq) Then
         Npq = nBas(iSymp)*(nBas(iSymp) + 1)/2
      Else
         Npq = nBas(iSymp)*nBas(iSymq)
      End If

      If (iSymr .eq. iSyms) Then
         Nrs = nBas(iSymr)*(nBas(iSymr) + 1)/2
      Else
         Nrs = nBas(iSymr)*nBas(iSyms)
      End If

C --- Set up the batch procedure
C ------------------------------
      Call GetMem('Maxmem','MAX ','REAL',KDUM,LWORK)

C------- Memory management ----------------------------
      if (iSymp.ne.iSymr) then
              mNeed = 2*Max(Npq,Nrs) + Nrs
      else
              mNeed = 2*Npq
      endif
C------------------------------------------------------
      If (mNeed .gt. 0) Then
         nVec = Min(LWORK/mNeed,NumCho(jSym))
      Else
C         ***QUIT*** bad initialization
         WRITE(6,*) 'Gen_Int: bad initialization'
         rc=99
         CALL QTrace()
         CALL Abend()
         nVec = -9999  ! dummy assignment - avoid compiler warnings
      End If
      If (nVec .gt. 0) Then
         nBatch = (NumCho(jSym) - 1)/nVec + 1
      Else
C         ***QUIT*** insufficient memory
         WRITE(6,*) 'Gen_Int: Insufficient memory for batch'
         WRITE(6,*) 'LWORK= ',LWORK
         WRITE(6,*) 'mNeed= ',mNeed
         WRITE(6,*) 'NumCho= ',NumCho(jsym)
         WRITE(6,*) 'jsym= ',jsym
         rc = rcRD05
         CALL QTrace()
         CALL Abend()
         nBatch = -9999  ! dummy assignment
      End If

C  --- Start the batch procedure for reading the vectors and computing
C  the integrals ---

c      Call FZero(Xint,numpq*Nrs)
      Do i=1,numpq*Nrs
         Xint(i) = ZERO
      End Do

      DO iBatch = 1,nBatch

         If (iBatch .eq. nBatch) Then
            NumV = NumCho(jSym) - nVec*(nBatch - 1)
         Else
            NumV = nVec
         End If

         iVec1 = nVec*(iBatch - 1) + 1

         if (iSymp.ne.iSymr) then

            LenMem1 = Max(Npq,Nrs)*NumV
            LenMem2 = Max(Npq,Nrs)*NumV
            LenMem3 = Nrs*NumV
C --- Allocate memory for reading the vectors and do the transposition
            Call GetMem('MemC1','ALLO','REAL',kVec1,LenMem1)
            Call GetMem('MemC2','ALLO','REAL',kVec2,LenMem2)
            Call GetMem('MemC3','ALLO','REAL',kVec3,LenMem3)

         else

            LenMem1 = Npq*NumV  ! equal to LenMem2
C --- Allocate memory for reading the vectors and do the transposition
            Call GetMem('MemC1','ALLO','REAL',kVec1,LenMem1)
            Call GetMem('MemC2','ALLO','REAL',kVec2,LenMem1)

         endif

C !--- Copying out (and transpose) the elements of the 1st vector ---!
C ---------------- L(pq,J) ---> L(qp,J)  -----------------------------
C --------------------------------------------------------------------
         if (iSymp.ne.iSymq) then  ! transposition needed
            Call RdChoVec(Work(kVec1),Npq,NumV,iVec1,LuCVec(1))
            koff1=0
            koff2=0
            do jvec=1,NumV
             do jq=1,nBas(iSymq)
               do jp=1,nBas(iSymp)

                  koff1 = kVec1 + Npq*(jvec-1) + nBas(iSymp)*(jq-1)
     &                          + (jp-1)
                  koff2 = kVec2 + Npq*(jvec-1) + nBas(iSymq)*(jp-1)
     &                          + (jq-1)
                  work(koff2) = work(koff1)

               end do
             end do
            end do
         else  ! no need to transpose "diagonal" vectors
            Call RdChoVec(Work(kVec2),Npq,NumV,iVec1,LuCVec(1))
         endif
         kWqp=kVec2
C --------------------------------------------------------------------

      IF (numpq.eq.Npq) THEN

       If (iSymp.ne.iSymr) THEN !need to read the 2nd vector also
C --------------------------------------------------------------------
         if (iSymr.ne.iSyms) then
            Call RdChoVec(Work(kVec1),Nrs,NumV,iVec1,LuCVec(2))
            koff1=0
            koff2=0
            do jvec=1,NumV
             do js=1,nBas(iSyms)
               do jr=1,nBas(iSymr)

                  koff1 = kVec1 + Nrs*(jvec-1) + nBas(iSymr)*(js-1)
     &                          + (jr-1)
                  koff2 = kVec3 + Nrs*(jvec-1) + nBas(iSyms)*(jr-1)
     &                          + (js-1)
                  work(koff2) = work(koff1)

               end do
             end do
            end do
         else
            Call RdChoVec(Work(kVec3),Nrs,NumV,iVec1,LuCVec(2))
         endif
         kWsr = kVec3
C --------------------------------------------------------------------

C --- Computing the integrals (II|JJ)
C -----------------------------------
C --- (sr|{qp}) <- (sr|{qp}) + sum_I L(sr,#I)* L({qp},#I)
C==========================================================
         CALL DGEMM_('N','T',Nrs,numpq,NumV,
     &              ONE,Work(kWsr),Nrs,
     &              WORK(kWqp),numpq,ONE,Xint,Nrs)


       ELSE ! isymp = isymr   (Npq=Nrs)

C --- Computing integrals of the type (II|II) and (IJ|IJ)

         CALL DGEMM_('N','T',Nrs,numpq,NumV,
     &              ONE,Work(kWqp),Nrs,
     &              WORK(kWqp),numpq,ONE,Xint,Nrs)


       END If


      ELSE  ! numpq.ne.Npq

C !--- Copying out the elements of the 1st vector ---!
C ----------------------------------------------------
            Do J = 1,NumV
               Do jpq = 1,numpq
                  pq = pq1 + jpq - 1
C --- Address of the matrix element (pq,J) in the full matrix
                  kOff1 = kVec2 + Npq*(J - 1) + (pq - 1)
C --- Address of the matrix element (pq,J) in the sub-block matrix
                  kOff2 = kVec1 + numpq*(J - 1) + (jpq - 1)
C --- Copy out the elements of the sub-block matrix if not the full
C --- matrix
                  Work(kOff2) = Work(kOff1)
               End Do
            End Do
            kXqp = kVec1
            kWsr = kVec2

       If (iSymp.ne.iSymr) Then
C --------------------------------------------------------------------
         if (iSymr.ne.iSyms) then  !   L(rs,J) ---> L(sr,J)
            Call RdChoVec(Work(kVec2),Nrs,NumV,iVec1,LuCVec(2))
            koff1=0
            koff2=0
            do jvec=1,NumV
             do js=1,nBas(iSyms)
               do jr=1,nBas(iSymr)

                  koff1 = kVec2 + Nrs*(jvec-1) + nBas(iSymr)*(js-1)
     &                          + (jr-1)
                  koff2 = kVec3 + Nrs*(jvec-1) + nBas(iSyms)*(jr-1)
     &                          + (js-1)
                  work(koff2) = work(koff1)

               end do
             end do
            end do
         else
            Call RdChoVec(Work(kVec3),Nrs,NumV,iVec1,LuCVec(2))
         endif
         kWsr = kVec3
C --------------------------------------------------------------------
       Endif

C --- Computing the integrals
C -----------------------------------
C --- (sr|{qp}) <- (sr|{qp}) + sum_I L(sr,#I)* L({qp},#I)
C==========================================================
         CALL DGEMM_('N','T',Nrs,numpq,NumV,
     &              ONE,Work(kWsr),Nrs,
     &              WORK(kXqp),numpq,ONE,Xint,Nrs)


      ENDIF


C --- Free the memory
       if (iSymp.ne.iSymr) then
          Call GetMem('MemC3','FREE','REAL',kVec3,LenMem3)
          Call GetMem('MemC2','FREE','REAL',kVec2,LenMem2)
          Call GetMem('MemC1','FREE','REAL',kVec1,LenMem1)
       else
          Call GetMem('MemC2','FREE','REAL',kVec2,LenMem1)
          Call GetMem('MemC1','FREE','REAL',kVec1,LenMem1)
       endif

      END DO  ! end of the batch procedure

      rc  = rc0000
      pq1 = pq1_save

      Return
      END

**************************************************************
