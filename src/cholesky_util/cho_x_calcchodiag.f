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
* Copyright (C) Francesco Aquilante                                    *
*               Thomas Bondo Pedersen                                  *
************************************************************************
      SUBROUTINE Cho_X_CalcChoDiag(rc,Diag)
************************************************************
*
*   <DOC>
*     <Name>Cho\_X\_CalcChoDiag</Name>
*     <Syntax>Call Cho\_X\_CalcChoDiag(rc,Diag)</Syntax>
*     <Arguments>
*       \Argument{rc}{Return code}{Integer}{out}
*       \Argument{Diag}{Array containing diagonal on exit}{Real*8}{out}
*     </Arguments>
*     <Purpose>Calculate integral diagonal from Cholesky vectors
*     </Purpose>
*     <Dependencies>Cho\_X\_Init must have been called</Dependencies>
*     <Author>Francesco Aquilante</Author>
*     <Modified_by>Thomas Bondo Pedersen</Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*         This routine calculates the integral diagonal from Cholesky
*         vectors,
*
*         (ab|ab) = sum(J)  (Lab,J)**2   (a,b:  AO-indices)
*
*         The diagonal calculation is parallelized.
*         The diagonal is returned in first reduced set storage and must
*         be allocated before calling this routine.
*         Return code is 0 if successful execution.
*     </Description>
*    </DOC>
**********************************************************************

**********************************************************************
*  Author : F. Aquilante
*
**********************************************************************
C
C  Computes diagonal elements of the ERI matrix from Cholesky vectors
C          and returns them in Diag (1st red set storage)
C
C      (ab|ab) = sum_J  (Lab,J)^2
C
C      a,b:  AO-index
C
**********************************************************************

      Implicit Real*8 (a-h,o-z)

      Integer   rc
      Real*8    Diag(*)
      Character*17 SECNAM
      Parameter (SECNAM = 'Cho_X_CalcChoDiag')

#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      parameter ( N2 = InfVec_N2 )

************************************************************************
      InfVec(i,j,k) = iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
******
      IndRed(i,k) = iWork(ip_IndRed-1+nnBstrT(1)*(k-1)+i)
******
      nDimRS(i,j) = iWork(ip_nDimRS-1+nSym*(j-1)+i)
************************************************************************

#if defined (_DEBUG_)
      Call QEnter(SECNAM)
#endif

      Call fZero(Diag,nnBstRT(1))


      IREDC= -1  ! unknwn reduced set

      iLoc = 3 ! use scratch location in reduced index arrays

C *************** BIG LOOP OVER VECTORS SYMMETRY *******************
      DO jSym=1,nSym

         If (NumCho(jSym) .lt. 1) GOTO 1000

         JRED1 = InfVec(1,2,jSym)  ! red set of the 1st vec
         JRED2 = InfVec(NumCho(jSym),2,jSym) !red set of the last vec
         Do JRED=JRED1,JRED2

            CALL Cho_X_nVecRS(JRED,JSYM,iVrs,nVrs)

            If (nVrs.eq.0) GOTO 999  ! no vectors in that (jred,jsym)

            if (nVrs.lt.0) then
               Write(6,*)SECNAM//': Cho_X_nVecRS returned nVrs<0. STOP!'
               rc = 77
               Return
            endif

            Call Cho_X_SetRed(irc,iLoc,JRED) !set index arrays at iLoc
            if(irc.ne.0)then
              Write(6,*)SECNAM//'cho_X_setred non-zero return code.',
     &                         '  rc= ',irc
              rc = irc
              Return
            endif

            IREDC=JRED

            nRS = nDimRS(JSYM,JRED)

            Call GetMem('MaxM','Max','Real',KDUM,LWORK)

            nVec  = Min(LWORK/Max(nRS,1),nVrs)

            If (nVec.lt.1) Then
               WRITE(6,*) SECNAM//': Insufficient memory for batch'
               WRITE(6,*) ' LWORK= ',LWORK
               WRITE(6,*) ' jsym= ',jsym
               WRITE(6,*) ' min. mem. need for reading= ',nRS
               rc = 33
               Return
               nBatch = -9999  ! dummy assignment
            End If

            LREAD = nRS*nVec

            Call GetMem('rsL','Allo','Real',ipLrs,LREAD)

C --- BATCH over the vectors ----------------------------

            nBatch = (nVrs-1)/nVec + 1

            DO iBatch=1,nBatch

               If (iBatch.eq.nBatch) Then
                  JNUM = nVrs - nVec*(nBatch-1)
               else
                  JNUM = nVec
               endif

               JVEC = nVec*(iBatch-1) + iVrs
               IVEC2 = JVEC - 1 + JNUM

               CALL CHO_VECRD(Work(ipLrs),LREAD,JVEC,IVEC2,JSYM,
     &                        NUMV,IREDC,MUSED)

               If (NUMV.le.0 .or.NUMV.ne.JNUM ) then
                  Call GetMem('rsL','Free','Real',ipLrs,LREAD)
                  rc=77
                  Return
               End If


C ---------------------------------------------------------------------
C --- Compute the diagonals :   D(ab) = D(ab) + sum_J (Lab,J)^2
C
C --- Stored in the 1st reduced set

               Do krs=1,nRS

                  mrs = iiBstR(JSYM,iLoc) + krs
                  jrs = IndRed(mrs,iLoc) ! address in 1st red set

                  Do jvc=1,JNUM

                     ipL = ipLrs + nRS*(jvc-1)
                     Diag(jrs) = Diag(jrs)
     &                         + Work(ipL+krs-1)**2
                  End Do

               End Do

C --------------------------------------------------------------------
C --------------------------------------------------------------------

            END DO  ! end batch loop

C --- free memory
            Call GetMem('rsL','Free','Real',ipLrs,LREAD)

999         Continue


         END DO   ! loop over red sets

1000     CONTINUE

      END DO   !loop over JSYM


      Call Cho_GAdGOp(Diag(1),NNBSTRT(1),'+')

      rc  = 0

#if defined (_DEBUG_)
      CAll QExit(SECNAM)
#endif

      Return
      END

**************************************************************
**************************************************************
