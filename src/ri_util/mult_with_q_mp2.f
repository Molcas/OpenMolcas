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
      Subroutine Mult_with_Q_MP2(nBas_aux,nBas,nIrrep)
**************************************************************************
*     Author: Jonas Bostrom
*
*     Purpose: Multiply MP2 A~_sep and B~_sep with inverse cholesky factors.
*
**************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"
#include "exterm.fh"
#include "chomp2g_alaska.fh"
*
      Integer nBas_Aux(1:nIrrep), nBas(1:nIrrep)
      Integer nLRb(8)
      Character*6 Name_Q, Name
      Integer Lu_B(4), Lu_A(2) , iAdrA_in(8), iAdrA_Out(8)
*
      Character*15 SECNAM
      Parameter (SECNAM = 'Mult_with_Q_MP2')
      Logical timings
*
      COMMON  /CHOTIME /timings
*                                                                      *
************************************************************************
*                                                                      *
      CALL CWTime(TotCPU1,TotWall1)
*
      Do i = 1,2
         iSeed=7
         Lu_A(i) = IsFreeUnit(iSeed)
         Write(Name,'(A5,I1)') 'AMP2V', i
         Call DaName_MF_WA(Lu_A(i),Name)
      End Do
*
      Do i = 1, 4
         iSeed = 7
         Lu_B(i) = IsFreeUnit(iSeed)
         Write(Name,'(A5,I1)') 'BMP2V', i
         Call DaName_MF_WA(Lu_B(i),Name)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      iA_in = 1
      iA_Out = 1
      Do iSym = 1, nSym
         NumCV = NumCho(iSym)
         NumAux = nBas_Aux(iSym)-1
         nLR = 0
         Do jSym = 1, nSym
            kSym = iEor(iSym-1,jSym-1)+1
            nLR = nLR + nBas(jSym)*nBas(kSym)
         End Do
         nLRb(iSym) = nLR
         iAdrA_in(iSym)  = iA_in
         iA_in = iA_in + NumCV*NumCV
         iAdrA_out(iSym) = iA_out
         iA_out= iA_out+ NumAux*NumAux
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      MaxValue = 200
      Do iSym = 1, nSym
*
         nBas2 = nLRb(iSym)
         NumCV = NumCho(iSym)
         NumAux = nBas_Aux(iSym)-1
         nAuxVe = NumAux
*
*     Get Q-vectors from disk
*     -----------------------
*
      l_Q = NumCV*NumAux
      Call GetMem('Q_Vector','Allo','Real',ip_Q,l_Q)
*
      iSeed=7
      Lu_Q = IsFreeUnit(iSeed)
      Write(Name_Q,'(A4,I2.2)') 'QVEC', iSym-1
      Call DaName_MF_WA(Lu_Q,Name_Q)
*
      iOpt = 2
      iAdrQ=0
      Call dDaFile(Lu_Q,iOpt,Work(ip_Q),l_Q,iAdrQ)
*
*     Get MP2 A-tilde vectors from disk
*     ---------------------------------
*
      l_A_t  = NumCV*NumCV
      l_A    = NumAux*NumAux
      l_A_ht = NumAux*NumCV
      Call GetMem('A_Tilde_vec','Allo','Real',ip_A_t,l_A_t)
      Call GetMem('A_vec','Allo','Real',ip_A,l_A)
      Call GetMem('A_HalfT_vec','Allo','Real',ip_A_ht,l_A_ht)
*
      Do iType = 1,2
*
         iOpt = 2
         iAdrA = iAdrA_in(iSym)
         Call dDaFile(Lu_A(iType),iOpt,Work(ip_A_t),l_A_t,iAdrA)
#ifdef _DEBUG_
         Write(6,*) 'Q-vectors'
         Do i = 1, l_Q
            Write(6,*) Work(ip_Q+i-1)
         End Do

         Write(6,*) 'A-vectors'
         Do i = 1, l_A
            Write(6,*) Work(ip_A_t+i-1)
         End Do
#endif
*
*     Make first halftransformation to cholesky-base
*     ----------------------------------------------
*
         Call dGemm_('N','N', NumAux, NumCV, NumCV,
     &              1.0d0, Work(ip_Q),NumAux,
     &                     Work(ip_A_t), NumCV,
     &              0.0d0, Work(ip_A_ht), NumAux)
*
         Call dGemm_('N','T', NumAux, NumAux, NumCV,
     &              1.0d0, Work(ip_A_ht), NumAux,
     &                     Work(ip_Q),NumAux,
     &              0.0d0, Work(ip_A), NumAux)
*
*     Put transformed A-vectors back on disk
*
         iOpt = 1
         iAdrA = iAdrA_out(iSym)
         Call dDaFile(Lu_A(iType),iOpt,Work(ip_A),l_A,iAdrA)
*
      End Do
*
      Call GetMem('A_Tilde_vec','Free','Real',ip_A_t,l_A_t)
      Call GetMem('A_vec','Free','Real',ip_A,l_A)
      Call GetMem('A_HalfT_vec','Free','Real',ip_A_ht,l_A_ht)
*                                                                      *
************************************************************************
*                                                                      *
         Call GetMem('MaxMem','Max','Real',iDum,MaxMem)
         MaxMem=9*(MaxMem/10)
         Call GetMem('MaxMem','Allo','Real',ip_B_t,MaxMem)
*
         nVec = MaxMem / ( 2*nLRb(iSym) )
         nVec = min(Max(NumCV,NumAux),nVec)
         If(nVec .lt. 1) Then
            Call ChoMP2_Quit(SecNam,'nVec is non-positive','[1]')
         End If
*
         l_B   = nLRb(iSym)*nVec
         l_B_t = nLRb(iSym)*nVec
         ip_B = ip_B_t + l_B_t
*
      Do iType = 1,2
*
*     The B-vectors should be read one batch at the time
*     --------------------------------------------------
*
         Do kVec = 1, NumAux, nVec
            NumVecK = Min(nVec,NumAux-(kVec-1))
*
            Do jVec = 1, NumCV, nVec
               NumVecJ = Min(nVec,NumCV - (jVec-1))
*
               iOpt = 2
               lTot = NumVecJ*nLRb(iSym)
               iAdr = 1 + nLRb(iSym)*(jVec-1)
               Call dDaFile(Lu_B(iType),iOpt,Work(ip_B_t),lTot,iAdr)
*
               Fac = 0.0D0
               If (jVec.ne.1) Fac = 1.0D0
               iOffQ1 = kVec-1 + NumAux*(jVec-1)
               Call dGemm_('N','T',nBas2, NumVecK, NumVecJ,
     &                    1.0d0,Work(ip_B_t), nBas2,
     &                          Work(ip_Q+iOffQ1),NumAux,
     &                    Fac,  Work(ip_B),nBas2)
            End Do
*
            iOpt = 1
            lTot = NumVecK*nBas2
            iAdr = 1 + nBas2*(kVec-1)
            Call dDaFile(Lu_B(iType+2),iOpt,Work(ip_B),lTot,iAdr)
*
         End Do
*
      End Do
*
         Call GetMem('MaxMem','Free','Real',ip_B_t,MaxMem)
         Call GetMem('Q_Vector','Free','Real',ip_Q,l_Q)
*
         Call DaClos(Lu_Q)
*
      End Do ! iSym
*
      Do i = 1,2
         Call DaClos(Lu_A(i))
      End Do
      Do i = 1, 4
         Call DaClos(Lu_B(i))
      End Do
*
      Return
      End
