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
      Subroutine Mult_with_Q_CASPT2(nBas_aux,nBas,nIrrep,SubAux)
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
      Logical SubAux
      Integer nLRb(8)
      Character*6 Name_Q
C     Integer Lu_B(4), Lu_A(2) , iAdrA_in(8), iAdrA_Out(8)
*
      Character*15 SECNAM
      Parameter (SECNAM = 'Mult_with_Q_MP2')
      Logical timings,is_error
*
      COMMON  /CHOTIME /timings
*
      Character*4096 RealName
      Integer, Parameter :: LuCMOPT2=61, !! The A-vector
     *                      LuGAMMA=60   !! The B-vector
*                                                                      *
************************************************************************
*                                                                      *
      CALL CWTime(TotCPU0,TotWall0)
*                                                                      *
************************************************************************
*                                                                      *
      iA_in = 1
      iA_Out = 1
      Do iSym = 1, nSym
         NumCV = NumCho(iSym)
         NumAux = nBas_Aux(iSym)-1
         If (SubAux) NumAux = nBas_Aux(iSym)-1
         nLR = 0
         Do jSym = 1, nSym
            kSym = iEor(iSym-1,jSym-1)+1
            nLR = nLR + nBas(jSym)*nBas(kSym)
         End Do
         nLRb(iSym) = nLR
C        iAdrA_in(iSym)  = iA_in
         iA_in = iA_in + NumCV*NumCV
C        iAdrA_out(iSym) = iA_out
         iA_out= iA_out+ NumAux*NumAux
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      MaxValue = 200
      Do iSym = 1, nSym
*
         nBas2 = nLRb(iSym)
C        write (*,*) "nBas2 = ", nBas2
         NumCV = NumCho(iSym)
         NumAux = nBas_Aux(iSym)-1
         If (SubAux) NumAux = nBas_Aux(iSym)-1
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
      Call PrgmTranslate('CMOPT2',RealName,lRealName)
C     call molcas_Open(LuCMOPT2,RealName(1:lRealName))
C     Open (Unit=LuCMOPT2,
C    *      File=RealName(1:lRealName),
C    *      Status='OLD',
C    *      Form='UNFORMATTED')
      Call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),
     &                      'DIRECT','UNFORMATTED',
     &                      iost,.FALSE.,
     &                      1,'OLD',is_error)
C     If (is_error) then
C       write (6,'(x,"Maybe, you did not add GRAD or GRDT ",
C    *               "keyword in &CASPT2?")')
C       write (6,'(x,"Please add either one, if this is single-point ",
C    *               "gradient calculation.")')
C       write (6,'(x,"Otherwise, something is wrong...")')
C       call abend()
C     end if
      Do i = 1, l_A_t
        Read (LuCMOPT2,END=100) Work(ip_A_t+i-1)
      End Do
*
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
     &           1.0d0, Work(ip_Q),NumAux,
     &                  Work(ip_A_t), NumCV,
     &           0.0d0, Work(ip_A_ht), NumAux)
*
      Call dGemm_('N','T', NumAux, NumAux, NumCV,
     &           1.0d0, Work(ip_A_ht), NumAux,
     &                  Work(ip_Q),NumAux,
     &           0.0d0, Work(ip_A), NumAux)
*
*     Put transformed A-vectors back on disk
*
      Rewind (LuCMOPT2)
      Do i = 1, l_A
        Write (LuCMOPT2) Work(ip_A+i-1)
      End Do
*
      Close (LuCMOPT2)
C     write (*,*) "l_A_t = ", l_A_t
C     write (*,*) "l_A   = ", l_A
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
C     write (*,*) "nvec in mult = ", nvec
      If(nVec .lt. 1) Then
         Call ChoMP2_Quit(SecNam,'nVec is non-positive','[1]')
      End If
*
      l_B   = nLRb(iSym)*nVec
      l_B_t = nLRb(iSym)*nVec
      ip_B = ip_B_t + l_B_t
*
      Call PrgmTranslate('GAMMA',RealName,lRealName)
C     call molcas_Open(LuGAMMA,RealName(1:lRealName))
C     Open (Unit=LuGAMMA,
C    *      File=RealName(1:lRealName),
C    *      Status='OLD',
C    *      Form='UNFORMATTED',
C    *      Access='DIRECT',
C    *      Recl=nBas2*8)
        Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                        'DIRECT','UNFORMATTED',
     &                        iost,.TRUE.,
     &                        nBas2*8,'OLD',is_error)
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
            Do lVec = 1, NumCV
              Read (Unit=LuGAMMA,Rec=lVec)
     *       (Work(ip_B_t+i-1+nBas2*(lVec-1)),i=1,nBas2)
            End Do
*
            Fac = 0.0D0
            If (jVec.ne.1) Fac = 1.0D0
            iOffQ1 = kVec-1 + NumAux*(jVec-1)
            Call dGemm_('N','T',nBas2, NumVecK, NumVecJ,
     &                 1.0d0,Work(ip_B_t), nBas2,
     &                       Work(ip_Q+iOffQ1),NumAux,
     &                 Fac,  Work(ip_B),nBas2)
         End Do
*
         Do lVec = 1, NumAux
           Write (Unit=LuGAMMA,Rec=lVec)
     *       (Work(ip_B+i-1+nBas2*(lVec-1)),i=1,nBas2)
         End Do
      End Do
*
      Call GetMem('MaxMem','Free','Real',ip_B_t,MaxMem)
      Call GetMem('Q_Vector','Free','Real',ip_Q,l_Q)
*
      Call DaClos(Lu_Q)
      Close (LuGAMMA)
*
      End Do ! iSym
*
      CALL CWTime(TotCPU1,TotWall1)
C     write (6,*) 'CPU/Wall Time for mult_with_q_caspt2:',
C    *  totcpu1-totcpu0, totwall1-totwall0
      Return
C
  100 continue
      write (6,'(x,"Maybe, you did not add GRAD or GRDT ",
     *             "keyword in &CASPT2?")')
      write (6,'(x,"Please add either one, if this is single-point ",
     *             "gradient calculation.")')
      write (6,'(x,"Otherwise, something is wrong...")')
      call abend()
C
      End
