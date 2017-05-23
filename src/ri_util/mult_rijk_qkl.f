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
* Copyright (C) Jonas Bostrom                                          *
************************************************************************
      Subroutine Mult_RijK_QKL(iSO,nBas_aux,nIrrep)
**************************************************************************
*     Author: Jonas Bostrom
*
*     Purpose: Computation of the DF coefficient vectors in MO-basis.
*
*     Equation:   C(il,K)  =  sum_J  R(il,J) * Q(K,J)
*
*     Input:
*            iSO : alpha (iSO=1) or beta (iSO=2) orbitals
*            nBas_aux : number of aux bsfs in each irrep.
*            nIrrep : number of irreps.
**************************************************************************
      Implicit Real*8 (a-h,o-z)
      Integer nBas_Aux(1:nIrrep), nVec(1:nIrrep)
      Character  Fname*6, Fname2*6, Name_Q*6
      Character*50 CFmt
      Character*13 SECNAM
      Parameter (SECNAM = 'MULT_RIJK_QKL')
      Logical timings
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"
#include "exterm.fh"
#include "pso.fh"
*#define _DEBUG_
*#define _CD_TIMING_
#ifdef _CD_TIMING_
#include "temptime.fh"
#endif
#include "para_info.fh"
*
      parameter ( N2 = InfVec_N2 )
      COMMON  /CHOTIME /timings
*
*************************
*     Define some indeces
      MulD2h(i,j) = iEOR(i-1,j-1) + 1
      InfVec(i,j,k) = iWork(ip_InfVec-1+MaxVec*N2*(k-1)+MaxVec*(j-1)+i)
*************************

      CALL CWTime(TotCPU1,TotWall1)

      Do i = 1, nIrrep
         nVec(i) = NumCho(i)
      End Do
      Call GAIGOP(nVec,nIrrep,'+')

      nTotCho = 0
      MaxCho = 0
      MaxLocCho = 0
      Do jSym = 1, nIrrep
         nTotCho = nTotCho + nVec(jSym)
         MaxCho = Max(MaxCho,nVec(jSym))
         MaxLocCho= Max(MaxLocCho,NumCho(jSym))
         Do iSym = 1, nIrrep
            iAdrCVec(jSym,iSym,iSO) = 0
         End Do
      End Do


*     Loop over the first cholesky symmetry
*
      kCount=0
      Do jSym = 1, nIrrep
*
***      Check so the symmetry contains vectors
*-------------------------------------------------
         NumCV = NumCho(jSym)
         NumAux = nBas_Aux(jSym)
         If(jSym.eq.1) NumAux = NumAux - 1
*     Save the number of auxiliary basis functions to be
*     accessed later
         NumAuxVec(jSym) = NumAux


         Call GAIGOP_SCAL(NumCV,'max')
         If(NumCV .lt. 1) goto 1000
*
         nTotFIorb = 0
         MaxMOprod = 0
         MaxMOprodR = 0
         Do iSym = 1, nIrrep
            kSym = MulD2h(JSym,iSym)

            nTotFIorb=nTotFIorb + nIJ1(iSym,kSym,iSO)
            MaxMOprod = Max(MaxMOprod,nIJ1(iSym,kSym,iSO))
            MaxMOProdR = Max(MaxMOprodR,nIJR(iSym,kSym,iSO))
         End Do
*
         Call GetMem('MemChk','Max','Real',iDum,MemMax)
         nJvec1 = (MemMax-NumAux*MaxMOprod)/(MaxMOprod + NumAux)
         If(nJvec1.lt.1) Then
            Write(6,*) 'Too little memory in:',SECNAM
            Call Abend
         End If
         nJvec1 = min(nJvec1,NumCho(jSym))
         nJbat = NumCho(jSym)/nJvec1
         iRest = mod(NumCho(jSym),nJvec1)
         If(iRest.ne.0) Then
            nJbat = nJbat+1
            nJvecLast = iRest
         Else
            nJvecLast = nJvec1
         End If
*
         If(Is_Real_Par()) Then
            iFirstCho = InfVec(1,5,jSym)
         Else
            iFirstCho = 1
         End If
*
         l_QVector = nJVec1*NumAux
         l_RVector = MaxMOprod*nJVec1
         l_CVector = MaxMOprodR*NumAux
*
         Call GetMem('Q_Vector','Allo','Real',ip_QVector,l_QVector)
         Call GetMem('R_Vector','Allo','Real',ip_RVector,l_RVector)
         Call GetMem('C_Vector','Allo','Real',ip_CVector,l_CVector)
*

         iSeed2 = 8
         LuCVec = IsFreeUnit(iSeed2)
         If (iSO.eq.1) Then
            Write(Fname2,'(A4,I1,I1)') 'CVEA',jSym
         ElseIf (iSO.eq.2) Then
            Write(Fname2,'(A4,I1,I1)') 'CVEB',jSym
         EndIf
         Call DANAME_MF_WA(LuCVec,Fname2)
         iAdrC = 0
*
***   Get Q Vectors from Disk
*----------------------------------
      Do iSym = 1, nIrrep
         lSym = MulD2h(iSym,jSym)
*
         If(nIJ1(iSym,lSym,iSO).lt.1) Go To 2000
         Call FZero(Work(ip_CVector),l_CVector)
         Do iJBat = 1, nJBat
            If(iJBat.eq.nJBat) Then
               njVec = nJVecLast
            Else
               nJvec = nJvec1
            End If
*
            iSeed=55+jSym-1
            Lu_Q=IsFreeUnit(iSeed)
            Write(Name_Q,'(A4,I2.2)') 'QVEC',jSym-1
            Call DaName_MF_WA(Lu_Q,Name_Q)
            l_Q = nJvec*NumAux
            iAdrQ=(iFirstCho-1)*NumAux + (iJBat-1)*nJVec*NumAux
            Call dDaFile(Lu_Q,2,Work(ip_Qvector),l_Q,iAdrQ)

#ifdef _DEBUG_
            Call RecPrt('Q-vectors',' ',Work(ip_QVector),
     &                  nJVec,NumAux)
#endif
*
*
            iSeed=7
            LuRVec = IsFreeUnit(iSeed)
            If (iSO.eq.1) Then
               Write(Fname,'(A4,I1,I1)') 'CHTA',iSym,lSym
            ElseIf (iSO.eq.2) Then
               Write(Fname,'(A4,I1,I1)') 'CHTB',iSym,lSym
            EndIf
            Call DANAME_MF_WA(LuRVec,Fname)
*
***         Loop over all cholesky vectors on all nodes
*------------------------------------------------------

*
*                 Get R-Vectors from disk
*-------------------------------------------

            iAdrR = nIJ1(iSym,lSym,iSO)*nJVec1*(iJBat-1)
            l_RVec = nJvec * nIJ1(iSym,lSym,iSO)
            Call dDaFile(LuRVec,2,Work(ip_RVector),
     &                   l_RVec,iAdrR)

            Call dGemm_('N','T',nIJ1(iSym,lSym,iSO),NumAux,nJVec,1.0d0,
     &                 Work(ip_RVector),nIJ1(iSym,lSym,iSO),
     &                 Work(ip_QVector),NumAux,0.0d0,
     &                 Work(ip_CVector),nIJ1(iSym,lSym,iSO))
         End Do



#ifdef _DEBUG_
         Write (6,*) 'jSym=',jSym
         Call RecPrt('R-Vectors',' ',Work(ip_RVector),
     &               nIJ1(iSym,lSym,iSO),NumAux)
         Call RecPrt('C-Vectors',' ',Work(ip_CVector),
     &               nIJ1(iSym,lSym,iSO),NumAux)
#endif
         If ((.not.lSA).and.(iSym.eq.lSym)) Then
            Do iAux = 1, NumAux
               index = -1
               Do i = 1, nChOrb(iSym-1,iSO)
                  index = index+i
                  index2 = index + (iAux-1)*nIJ1(iSym,lSym,iSO)
                  Work(ip_Cvector+index2) =
     &                 Work(ip_Cvector+index2)/sqrt(2.0d0)
               End Do
            End Do
         End If

         Call DaClos(Lu_Q)
         Call DACLOS(LuRVec)
*
 2000    Continue
*
         Call GADGOP(Work(ip_CVector),l_CVector,'+')
         iAdrCVec(jSym,iSym,iSO) = iAdrC
         Call dDaFile(LuCVec,1,Work(ip_CVector),
     &                nIJ1(iSym,lSym,iSO)*NumAux,iAdrC)
         If(nIJ1(iSym,lSym,iSO) .lt. nIJR(iSym,lSym,iSO)) Then
            iAdrC = iAdrC +
     &          (nIJR(iSym,lSym,iSO)-nIJ1(iSym,lSym,iSO))*NumAux
         End If

         End Do !iSym
*
         Call GetMem('C_Vector','Free','Real',ip_CVector,l_CVector)
         Call GetMem('R_Vector','Free','Real',ip_RVector,l_RVector)
         Call GetMem('Q_Vector','Free','Real',ip_QVector,l_QVector)
*
         Call DACLOS(LuCVec)
 1000    Continue

      End Do  ! jSym
*
      CALL CWTime(TotCPU2,TotWall2)
      TotCPU = TotCPU2 - TotCPU1
      TotWall= TotWall2 - TotWall1
#ifdef _CD_TIMING_
      rMult_CPU = TOTCPU
      rMult_Wall = TOTWALL
#endif
*
      If(timings)then

      CFmt='(2x,A)'
      Write(6,*)
      Write(6,CFmt)'Cholesky Gradients timing from '//SECNAM
      Write(6,CFmt)'----------------------------------------'
      Write(6,*)
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'                                CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
         Write(6,'(2x,A26,2f10.2)')'TOTAL                            '
     &                           //'         ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)

      endif

      Return
      End
