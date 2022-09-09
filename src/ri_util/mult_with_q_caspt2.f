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
#include "WrkSpc.fh"
#include "exterm.fh"
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
      ! Integer LUCMOPT2  !! The A-vector
      Integer LUGAMMA   !! The B-vector
      Integer LUAPT2
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
      Do iSym = 1, nSym
*
         nBas2 = nLRb(iSym)
C        write (*,*) "nBas2 = ", nBas2
         NumCV = NumCho(iSym)
         NumAux = nBas_Aux(iSym)-1
         If (SubAux) NumAux = nBas_Aux(iSym)-1
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
    !   LUCMOPT2 = 61
    !   Call PrgmTranslate('CMOPT2',RealName,lRealName)
    !   Call MOLCAS_Open_Ext2(LuCMOPT2,RealName(1:lRealName),
    !  &                      'DIRECT','UNFORMATTED',
    !  &                      iost,.FALSE.,
    !  &                      1,'OLD',is_error)

    !   Read (LuCMOPT2,END=100) Work(ip_A_t:ip_A_t+l_A_t-1)

      ! Read A_PT2 from LUAPT2
      LUAPT2 = 77
      call daname_mf_wa(LUAPT2, 'A_PT2')
      id = 0
      call ddafile(LUAPT2, 2, work(ip_A_t), l_A_t, id)

      !! Symmetrized A_PT2
      Do i = 1, NumCV
         Do j = 1, i
           aaa = Work(ip_A_t+i-1+NumCV*(j-1))
     *         + Work(ip_A_t+j-1+NumCV*(i-1))
           aaa = aaa*0.5d+00
           Work(ip_A_t+i-1+NumCV*(j-1)) = aaa
           Work(ip_A_t+j-1+NumCV*(i-1)) = aaa
         End Do
      End Do
*
#ifdef _DEBUGPRINT_
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
      ! Rewind (LuCMOPT2)
      ! Write (LuCMOPT2) Work(ip_A:ip_A+l_A-1)
      ! Close (LuCMOPT2)

      ! write A_PT2 to LUAPT2
      id = 0
      call ddafile(LUAPT2, 1, work(ip_A), l_A, id)

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
      !! Applicable to C1 only
      nBasTri = nBas(1)*(nBas(1)+1)/2
      nVec = MaxMem / ( 2*nBasTri+1 ) ! MaxMem / ( 2*nLRb(iSym)+1 )
      nVec = min(Max(NumCV,NumAux),nVec)
      If(nVec .lt. 1) Then
         Call ChoMP2_Quit(SecNam,'nVec is non-positive','[1]')
      End If
*
      l_B_t = nBasTri*nVec ! nLRb(iSym)*nVec
      ip_B = ip_B_t + l_B_t
      ip_B2= ip_B   + l_B_t

      LUGAMMA = 60
      Call PrgmTranslate('GAMMA',RealName,lRealName)
      Call MOLCAS_Open_Ext2(LuGamma,RealName(1:lRealName),
     &                     'DIRECT','UNFORMATTED',
     &                      iost,.TRUE.,
     &                      nBas2*8,'OLD',is_error)
*
      LuGamma2 = 62
      Call PrgmTranslate('GAMMA2',RealName,lRealName)
      Call MOLCAS_Open_Ext2(LuGamma2,RealName(1:lRealName),
     &                      'DIRECT','UNFORMATTED',
     &                      iost,.TRUE.,
     &                      NumAux*8,'REPLACE',is_error)
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
          Do lVec = 1, NumVecJ
            Read (Unit=LuGAMMA,Rec=jVec+lVec-1)Work(ip_B2:ip_B2+nBas2-1)
            !! symmetrize (mu nu | P)
            !! only the lower triangle part is used
            nseq = 0
            Do i = 1, nBas(1)
               Do j = 1, i
                 aaa = Work(ip_B2+i-1+nBas(1)*(j-1))
     *               + Work(ip_B2+j-1+nBas(1)*(i-1))
                 Work(ip_B_t+nseq+nBasTri*(lVec-1)) = aaa
                 nseq = nseq + 1
               End Do
            End Do
          End Do
*
          Fac = 0.0D0
          If (jVec.ne.1) Fac = 1.0D0
          iOffQ1 = kVec-1 + NumAux*(jVec-1)
          Call dGemm_('N','T',nBasTri, NumVecK, NumVecJ,
     &               1.0d0,Work(ip_B_t), nBasTri,
     &                     Work(ip_Q+iOffQ1),NumAux,
     &               Fac,  Work(ip_B),nBasTri)
        End Do
*
        !! Because scaling with 0.5 is omitted when symmetrized
        Call DScal_(nBasTri*NumVecK,0.5D+00,Work(ip_B),1)
*
        !! (mu nu | P) --> (P | mu nu)
        If (Max(NumCV,NumAux).eq.nVec) Then
          nseq = 0
          Do lVec = 1, nBas(1)
            Do jVec = 1, lVec
              nseq = nseq + 1
            Call DCopy_(NumAux,Work(ip_B+nseq-1),nBasTri,Work(ip_B_t),1)
             Write (Unit=LuGAMMA2,Rec=nseq) Work(ip_B_t:ip_B_t+NumAux-1)
            End Do
          End Do
        Else
          nseq = 0
          Do lVec = 1, nBas(1)
            Do jVec = 1 ,lVec
              nseq = nseq + 1
              If (kVec.ne.1) Read (Unit=LuGAMMA2,Rec=nseq)
     *          Work(ip_B_t:ip_B_t+kVec-2)
              Call DCopy_(NumVecK,Work(ip_B+nseq-1),nBasTri,
     *                            Work(ip_B_t+kVec-1),1)
              Write (Unit=LuGAMMA2,Rec=nseq)
     *          (Work(ip_B_t:ip_B_t+kVec+NumVecK-2))
            End Do
          End Do
        End If
      End Do
*
      Close (LuGAMMA2)
*
      Call GetMem('MaxMem','Free','Real',ip_B_t,MaxMem)
      Call GetMem('Q_Vector','Free','Real',ip_Q,l_Q)
*
      Call DaClos(Lu_Q)
      call daclos(LUAPT2)
      Close (LuGAMMA)
*
      End Do ! iSym
*
      CALL CWTime(TotCPU1,TotWall1)
C     write (6,*) 'CPU/Wall Time for mult_with_q_caspt2:',
C    *  totcpu1-totcpu0, totwall1-totwall0
      Return
C
  ! 100 continue
    !   write (6,'(x,"Maybe, you did not add GRAD or GRDT ",
    !  *             "keyword in &CASPT2?")')
    !   write (6,'(x,"Please add either one, if this is single-point ",
    !  *             "gradient calculation.")')
    !   write (6,'(x,"Otherwise, something is wrong...")')
    !   call abend()
C
      End
