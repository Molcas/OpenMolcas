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
      SubRoutine Cho_GnVc_GenVec(Diag,xInt,lInt,nVecRS,iVecRS,
     &                           mSym,mPass,iPass1,NumPass)
C
C     Purpose: generate Cholesky vectors from raw integral columns.
C
      use ChoArr, only: nDimRS
      use ChoSwp, only: InfVec, IndRed
#include "implicit.fh"
      Real*8  Diag(*), xInt(lInt)
      Integer nVecRS(mSym,mPass), iVecRS(mSym,mPass)
#include "cholesky.fh"
#include "choprint.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Character(LEN=15), Parameter:: SecNam = 'Cho_GnVc_GenVec'

      Integer NumCho_OLD(8), iOff1(8), iOff2(8)

      Real*8, Allocatable:: Wrk(:), VecTmp(:)

      Integer ip_mapRS2RS, l_mapRS2RS
      Common / GnVcMp / ip_mapRS2RS(8), l_mapRS2RS(8)

      mapRS2RS(i,j)=iWork(ip_mapRS2RS(i)-1+j)

C     Check input.
C     ------------

      If (NumPass .lt. 1) return

      If (mSym .ne. nSym) Then
         Call Cho_Quit('Input error [1] in '//SecNam,103)
      End If

      If (iPass1 .lt. 1) Then
         Call Cho_Quit('Input error [2] in '//SecNam,103)
      End If

      iPass2 = iPass1 + NumPass - 1
      If (iPass2 .gt. mPass) Then
         Call Cho_Quit('Input error [3] in '//SecNam,103)
      End If

      nPass = XnPass
      If (mPass .ne. nPass) Then
         Call Cho_Quit('Input error [4] in '//SecNam,103)
      End If

      NumVec = 0
      Do iPass = iPass1,iPass2
         Do iSym = 1,nSym
            NumVec = NumVec + nVecRS(iSym,iPass)
         End Do
      End Do
      If (NumVec .lt. 1) return ! exit

C     Subtract previous vectors.
C     --------------------------

      Call mma_maxDBLE(l_Wrk)
      Call mma_allocate(Wrk,l_Wrk,Label='Wrk')
      Do iSym = 1,nSym
         kOff = iOff_Col(iSym) + 1
         Call Cho_Subtr(xInt(kOff),Wrk,SIZE(Wrk),iSym)
      End Do
      Call mma_deallocate(Wrk)

C     Initialize vector generation.
C     -----------------------------

      Do iSym = 1,nSym
         iOff1(iSym) = iOff_Col(iSym) + 1
         iOff2(iSym) = iOff_Col(iSym) + 1
      End Do

      l_VecTmp = 0
      Do iPass = iPass1,iPass2
         Do iSym = 1,nSym
            Need = nnBstR(iSym,2)*nVecRS(iSym,iPass)
            l_VecTmp = max(l_VecTmp,Need)
         End Do
      End Do
      MxSubtr = 0
      Do iPass = iPass1,iPass2
         Do iSym = 1,nSym
            nAB = 0
            Do jPass = iPass+1,iPass2
               nAB = nAB + nVecRS(iSym,jPass)
            End Do
            Need = nAB*nVecRS(iSym,iPass)
            MxSubtr = max(MxSubtr,Need)
         End Do
      End Do
      l_VecTmp = max(l_VecTmp,MxSubtr)
      Call mma_allocate(VecTmp,l_VecTmp,Label='VecTmp')

C     Copy reduced set iPass1 to location 3.
C     --------------------------------------

      irc = 0
      Call Cho_X_RSCopy(irc,2,3)
      If (irc .ne. 0) Then
         Write(Lupri,*) SecNam,': Cho_X_RSCopy returned ',irc
         Call Cho_Quit('Error termination in '//SecNam,104)
      End If

C     Decomposition pass loop.
C     ------------------------

      Do iPass = iPass1,iPass2

C        Print header.
C        -------------

         LenLin = 0 ! to avoid compiler warnings
         If (iPrint .ge. INF_PROGRESS) Then
            Call Cho_Head(SecNam//
     &                    ': Generation of Vectors from Map','=',
     &                    80,Lupri)
            Write(Lupri,'(/,A,I5)')
     &      'Integral pass number',iPass
            Write(Lupri,'(A,8I8)')
     &      '#Cholesky vec.: ',(NumCho(iSym),iSym=1,nSym)
            Write(Lupri,'(A,8I8)')
     &      '#qualified    : ',(nVecRS(iSym,iPass),iSym=1,nSym)
            Write(Lupri,'(A,8I8)')
     &      'Current  dim. : ',(nnBstR(iSym,3),iSym=1,nSym)
            Write(Lupri,'(A,8I8)')
     &      'Original dim. : ',(nnBstR(iSym,1),iSym=1,nSym)
            Write(Lupri,'(/,A,/,A,A)')
     &      '           #Vectors             Treated Diagonal',
     &      'Sym.     Sym.     Total     Index     Before      After',
     &      '   Conv. Neg.   New Max'
            LenLin = 79
            Write(Lupri,'(80A)') ('-',i=1,LenLin)
            Call Cho_Flush(Lupri)
            Call iCopy(nSym,NumCho,1,NumCho_OLD,1)
         Else If (iPrint .ge. INF_PASS) Then
            Write(Lupri,'(/,A,I5)')
     &      'Integral pass number',iPass
            Write(LUPRI,'(A,8I8)')
     &      '#Cholesky vec.: ',(NumCho(iSym),iSym=1,nSym)
            Write(LUPRI,'(A,8I8)')
     &      '#qualified    : ',(nVecRS(iSym,iPass),iSym=1,nSym)
            Call Cho_Flush(Lupri)
            Call iCopy(nSym,NumCho,1,NumCho_OLD,1)
         End If

C        Zero entries in integral matrix that are not part of this
C        reduced set.
C        ---------------------------------------------------------

         If (iPass .gt. iPass1) Then
            Do iSym = 1,nSym
               Do iV = 1,nVecRS(iSym,iPass)
                  lTot = nnBstR(iSym,2)
                  VecTmp(1:lTot)=0.0D0
                  lOff0 = iOff1(iSym) + nnBstR(iSym,2)*(iV-1) - 1
                  Do lAB = 1,nnBstR(iSym,3)
                     jAB = IndRed(iiBstR(iSym,3)+lAB,3) - iiBstR(iSym,1)
                     kAB = mapRS2RS(iSym,jAB)
                     VecTmp(kAB) = xInt(lOff0+kAB)
                  End Do
                  Call dCopy_(lTot,VecTmp,1,xInt(lOff0+1),1)
               End Do
            End Do
         End If

C        Write reduced set info for this pass (index arrays are stored
C        at location 3).
C        -------------------------------------------------------------

         Call Cho_P_PutRed(iPass,3)

C        Start symmetry loop.
C        --------------------

         Do iSym = 1,nSym

            If (nVecRS(iSym,iPass) .lt. 1) Go To 100 ! cycle sym. loop

C           Generate vectors for this pass and symmetry.
C           --------------------------------------------

            Do iV = 1,nVecRS(iSym,iPass)

               kOff0 = iOff1(iSym) + nnBstR(iSym,2)*(iV-1) - 1

               iVec = iVecRS(iSym,iPass) + iV - 1
               iAB  = InfVec(iVec,1,iSym) ! addr in 1st red. set

               XC = Diag(iAB)
               If (abs(Diag(iAB)) .gt. 1.0d-14) Then ! TODO/FIXME
                  Fac = 1.0d0/sqrt(abs(Diag(iAB)))
               Else
                  Fac = 1.0d7
               End If
               kOff = kOff0 + 1
               Call dScal_(nnBstR(iSym,2),Fac,xInt(kOff),1)

               Do i = 1,nnBstR(iSym,2)
                  ii = iiBstR(iSym,2) + i
                  jj = IndRed(ii,2)
                  If (Diag(jj) .eq. 0.0d0) Then
                     xInt(kOff0+i) = 0.0d0
                  End If
               End Do

               Do i = 1,nnBstR(iSym,2)
                  ii = iiBstR(iSym,2) + i
                  jj = IndRed(ii,2)
                  Diag(jj) = Diag(jj) - xInt(kOff0+i)*xInt(kOff0+i)
               End Do

               olDiag    = Diag(iAB)
               Diag(iAB) = 0.0d0
               Call Cho_ChkDia(Diag,iSym,xMin,xMax,xM,nNegT,nNeg,nConv)
               nNZTot = nNZTot + nNeg

               kOff1 = kOff0 + 1
               Do jV = iV+1,nVecRS(iSym,iPass)
                  jVec = iVecRS(iSym,iPass) + jV - 1
                  jAB  = InfVec(jVec,1,iSym)
                  kOff2 = iOff1(iSym) + nnBstR(iSym,2)*(jV-1)
                  Fac   = -xInt(kOff0+mapRS2RS(iSym,jAB-iiBstR(iSym,1)))
                  Call dAXPY_(nnBstR(iSym,2),Fac,xInt(kOff1),1,
     &                                          xInt(kOff2),1)
               End Do

               Call Cho_SetVecInf(iVec,iSym,iAB,iPass,3)

               If (iPrint .ge. INF_PROGRESS) Then
                  iVecT = NumChT + iV
              Write(Lupri,'(I3,3(1X,I9),2(1X,D11.3),2(1X,I4),1X,D11.3)')
     &            iSym,iVec,iVecT,iAB,XC,olDiag,nConv,nNeg,xM
               End If

            End Do

C           Subtract contributions to later vectors.
C           ----------------------------------------

            nAB = nQual(iSym)
            Do jPass = iPass1,iPass
               nAB = nAB - nVecRS(iSym,jPass)
            End Do
            If (nAB .gt. 0) Then
               ip_Scr = 1
               iP = iPass
               jVec0 = -1
               Do While (iP.lt.iPass2 .and. jVec0.lt.0)
                  iP = iP +1
                  jVec0 = iVecRS(iSym,iP) - 1
               End Do
               If (jVec0 .lt. 0) Then ! should never happen
                  Call Cho_Quit('jVec0 < 0 in '//SecNam,103)
               End If
               Do iV = 1,nVecRS(iSym,iPass)
                  kOff1 = ip_Scr + nAB*(iV-1) - 1
                  kOff2 = iOff1(iSym) + nnBstR(iSym,2)*(iV-1) - 1
                  Do iAB = 1,nAB
                     jVec = jVec0 + iAB
                     jAB  = InfVec(jVec,1,iSym)
                     kAB  = mapRS2RS(iSym,jAB-iiBstR(iSym,1))
                     VecTmp(kOff1+iAB) = xInt(kOff2+kAB)
                  End Do
               End Do
               kOff1 = iOff1(iSym)
               kOff2 = iOff1(iSym)
     &               + nnBstR(iSym,2)*nVecRS(iSym,iPass)
               Call DGEMM_('N','T',
     &                    nnBstR(iSym,2),nAB,nVecRS(iSym,iPass),
     &                    -1.0d0,xInt(kOff1),nnBstR(iSym,2),
     &                           VecTmp,nAB,
     &                     1.0d0,xInt(kOff2),nnBstR(iSym,2))
            End If

C           Reorder vectors to appropriate reduced set.
C           Skipped for iPass1, as they are already in correct storage.
C           -----------------------------------------------------------

            If (iPass .gt. iPass1) Then
               lTot = nnBstR(iSym,2)*nVecRS(iSym,iPass)
               Call dCopy_(lTot,xInt(iOff1(iSym)),1,VecTmp,1)
               Do iV = 1,nVecRS(iSym,iPass)
                  kOff0 = iOff2(iSym) + nnBstR(iSym,3)*(iV-1) - 1
                  lOff0 = nnBstR(iSym,2)*(iV-1)
                  Do kAB = 1,nnBstR(iSym,3)
                     jAB = IndRed(iiBstR(iSym,3)+kAB,3)
                     lAB = mapRS2RS(iSym,jAB-iiBstR(iSym,1))
                     xInt(kOff0+kAB) = VecTmp(lOff0+lAB)
                  End Do
               End Do
            End If

C           Update vector counters.
C           -----------------------

            NumCho(iSym) = NumCho(iSym) + nVecRS(iSym,iPass)
            NumChT = NumChT + nVecRS(iSym,iPass)

C           Update pointer arrays.
C           iOff1: pointer to integral columns (in xInt).
C           iOff2: pointer to vectors (also in xInt).
C           ---------------------------------------------

            iOff1(iSym) = iOff1(iSym)
     &                  + nnBstR(iSym,2)*nVecRS(iSym,iPass)
            iOff2(iSym) = iOff2(iSym)
     &                  + nnBstR(iSym,3)*nVecRS(iSym,iPass)

C           Cycle point for empty symmetry.
C           -------------------------------

  100       Continue
            If (iPrint .ge. INF_PROGRESS) Call Cho_Flush(Lupri)

         End Do ! symmetry

C        Print.
C        ------

         If (iPrint .ge. INF_PROGRESS) Then
            Do iSym = 1,nSym
               NumCho_OLD(iSym) = NumCho(iSym) - NumCho_OLD(iSym)
            End Do
            Write(Lupri,'(80A)') ('-',I=1,LenLin)
            Write(Lupri,'(A,8I8)')
     &      '#vec. gener.  : ',(NumCho_OLD(iSym),iSym=1,nSym)
            Call Cho_Flush(Lupri)
         Else If (iPrint .GE. INF_PASS) Then
            Do iSym = 1,nSym
               NumCho_OLD(iSym) = NumCho(iSym) - NumCho_OLD(iSym)
            End Do
            Write(Lupri,'(A,8I8)')
     &      '#vec. gener.  : ',(NumCho_OLD(iSym),iSym=1,nSym)
            Call Cho_Flush(Lupri)
         End If

C        Analyze diagonal.
C        -----------------

         If (iPrint .ge. INF_PASS) Then
            Bin1 = 1.0d2
            Step = 1.0d-1
            nBin = 18
            Call Cho_AnaDia(Diag,Bin1,Step,nBin,.false.)
         End If

C        Set next (iPass+1) reduced set at location 2.
C        Reduced set iPass1 is now stored at location 3.
C        -----------------------------------------------

         Call Cho_SetRed(Diag)
         jPass = iPass + 1
         Call Cho_SetRSDim(nDimRS,nSym,MaxRed,jPass,2)
         If (iPrint .ge. INF_PASS) Then
            Call Cho_PrtRed(2)
            Call Cho_Flush(Lupri)
         End If

C        Swap locations so that:
C        location 2 contains reduced set iPass1 and
C        location 3 contains next (iPass+1) reduced set.
C        -----------------------------------------------

         irc = 0
         Call Cho_X_RSSwap(irc,2,3)
         If (irc .ne. 0) Then
            Write(Lupri,*) SecNam,': Cho_X_RSSwap returned ',irc
            Call Cho_Quit('Error termination in '//SecNam,104)
         End If

      End Do ! integral pass

C     Deallocate temporary vector array.
C     ----------------------------------

      Call mma_deallocate(VecTmp)

C     Write vectors to disk.
C     ----------------------

      Call Cho_Timer(C1,W1)
      Do iSym = 1,nSym
         NumVec = nVecRS(iSym,iPass1)
         Do iPass = iPass1+1,iPass2
            NumVec = NumVec + nVecRS(iSym,iPass)
         End Do
         If (NumVec .gt. 0) Then
            iPass = iPass1
            iVec1 = iVecRS(iSym,iPass)
            Do While (iVec1.lt.1 .and. iPass.lt.iPass2)
               iPass = iPass + 1
               iVec1 = iVecRS(iSym,iPass)
            End Do
            If (iVec1 .lt. 1) Then
               Call Cho_Quit('Logical error in '//SecNam,103)
            Else
               Call Cho_PutVec2(xInt(iOff_Col(iSym)+1),NumVec,iVec1,
     &                          iSym)
            End If
         End If
      End Do
      Call Cho_Timer(C2,W2)
      tDecom(1,2) = tDecom(1,2) + C2 - C1
      tDecom(2,2) = tDecom(2,2) + W2 - W1

C     Write restart files.
C     --------------------

      Call Cho_P_WrRstC(iPass2)

C     Store next (iPass2+1) reduced set at location 2.
C     ------------------------------------------------

      irc = 0
      Call Cho_X_RSCopy(irc,3,2)
      If (irc .ne. 0) Then
         Write(Lupri,*) SecNam,': Cho_X_RSCopy returned ',irc
         Call Cho_Quit('Error termination in '//SecNam,104)
      End If

      End
