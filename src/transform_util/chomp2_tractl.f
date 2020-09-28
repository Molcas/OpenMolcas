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
* Copyright (C) 2004,2005, Giovanni Ghigo                              *
************************************************************************
*  ChoMP2_TraCtl
*
*> @brief
*>   Driver for the generation of the two-electrons integrals file (``MOLINT``)
*>   from the AO-based Cholesky Full Vectors for MBPT2 program
*> @author Giovanni Ghigo
*>
*> @details
*> This routine is similar to ::Cho_TraCtl routine but only the
*> transformed vectors type ``C`` (TCVC) are generated. See ::Cho_TraCtl
*> and related routines for more details.
*>
*> @note
*> The number of frozen and deleted MO used in the post-SCF
*> must be written in the RunFile in \c nFroPT and \c nDelPT arrays.
*>
*> @param[in] LUINTM Unit number of two-electrons integrals file (``MOLINT``)
*> @param[in] CMO    MO coefficients
*> @param[in] NCMO   Total number of MO coefficients
************************************************************************
      Subroutine ChoMP2_TraCtl(LUINTM,CMO,NCMO)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden                                    *
*           October 2004                                               *
* Modified for Cholesky-MP2 May 2005                                   *
*----------------------------------------------------------------------*
* This is the routine for the transformation from AO basis to MO basis *
* of the Cholesky Full Vectors & Two-electrons integral generation.    *
* Reads the Cholesky Full Vectors (CHFV), generate new Transformed     *
* Cholesky Full Vectors (TCVx) then generate the Two-electrons inte-   *
* gral file (MOLINT).                                                  *
*  <p,q|k,l>  where p,q: All MO; k,l: Occupied                         *
* THIS CODE IS ONLY FOR MBPT2 AND IT IS SPLIT FROM THE GENERAL CODE    *
* BUT IT IS STILL INTEGRATED AND DEPENDENT ON THE GENERAL CODE         *
************************************************************************
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
#include "rasdim.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "cho_tra.fh"
      Dimension CMO(NCMO)
      Character*4 CHNm
      Character*6 CHName
      Parameter (CHNm='CHFV')
      Logical Found

      MulD2h(i,j) = iEor(i-1,j-1)+1

      IfTest=.False.

      Call Timing(CPU0,CPE,TIO0,TIOE)

***   INIZIALIZATION   *************************************************

* --- Define what has to be calculated. ---
*      DoTCVA flag for the generation of TCVA
*      DoFull flag for the generation of TCVF
*      DoCoul flag for the generation of coulomb integrals

*     MBPT2
        DoTCVA = .False.
        DoFull = .False.
        DoCoul = .False.

* --- Get Informations from RunFile. ---
*      The following informations must be passed to the Cholesky
*      transformation routine through RunFile. COMMON blocks could
*      not be used due to several conflicts.
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      Call Get_iArray('nFroPT',nFro,nSym)
      Call Get_iArray('nDelPT',nDel,nSym)
      Call Get_iArray('nIsh',nIsh,nSym)
*PAM07      Call Get_iArray('nAsh',nAsh,nSym)
* Replaced by:
      Do i=1,nSym
        nAsh(i)=0
      End Do
      Call qpg_iArray('nAsh',Found,nData)
      If(Found .and. nData.eq.nSym) Then
        Call Get_iArray('nAsh',nAsh,nSym)
      End If
* End of replacement
      Call Get_iArray('NumCho',NumCho,nSym)
      nBasT=0
      Do i=1,nSym
        nBasT   = nBasT + nBas(i)
        nOrb(i) = nBas(i) - nFro(i) - nDel(i)
        nIsh(i) = nIsh(i) - nFro(i)
        nOsh(i) = nIsh(i) + nAsh(i)
        nSsh(i) = nOrb(i) - nOsh(i)
      EndDo

* --- Inizialize information arrays. ---

*     The TCVx existing flag and the Memory Allocation & Length array:
      Do iType=1,MxTCVx
        Do iSym=1,MaxSym
          Do jSym=1,MaxSym
            TCVXist(iType,iSym,jSym)=.False. ! TCVx existing flag.
            iMemTCVX(iType,iSym,jSym,1)=0 ! Memory Address and
            iMemTCVX(iType,iSym,jSym,2)=0 ! Length in Work(TCVx)
          EndDo
        EndDo
      EndDo

*     The Address Field for MOLINT:
      LenIAD2M=3*36*36
      Do i=1,36*36
       IAD2M(1,i)=0
       IAD2M(2,i)=0
       IAD2M(3,i)=0
      EndDo
      iAddrIAD2M=0
      Call iDaFile(LUINTM,1,IAD2M,LenIAD2M,iAddrIAD2M)

*     The Timing:
      CPU_Tra=0.0d0
      TIO_Tra=0.0d0
      CPU_Gen=0.0d0
      TIO_Gen=0.0d0

* *** START LOOP iSymL on TOTAL SYMMETRY of L (iSym * jSym)   **********
      DO iSymL=1,nSym

*       Re-Inizialize the TCVx & iMemTCVX
        Do iSym=1,MaxSym
          Do jSym=1,MaxSym
            TCVXist(3,iSym,jSym)=.False. ! TCVx existing flag.
            iMemTCVX(3,iSym,jSym,1)=0 ! Memory Address and
            iMemTCVX(3,iSym,jSym,2)=0 ! Length in Work(TCVx)
          EndDo
        EndDo
        Call Mem_Est(iSymL,nVec,nFVec)
        If (nVec.GT.0 .and. nFVec.GT.0) then
          nBatch =(NumCho(iSymL)-1)/nVec + 1
        else
          Write(6,*)
          Write(6,*) ' ************************************'
          Write(6,*) ' *  Insufficient memory for batch ! *'
          Write(6,*) ' ************************************'
          Write(6,*)
          Call XFlush(6)
          Call QTrace()
          Call Abend()
        EndIf

*  ---  START LOOP iBatch  ---------------------------------------------
        DO iBatch=1,nBatch
          If (iBatch.EQ.nBatch) then
            NumV = NumCho(iSymL) -nVec * (nBatch-1)
          else
            NumV = nVec
          EndIf
          nFBatch = (NumV-1) / nFVec + 1

*   ---   Start Transformation of Cholesy Vectors  CHFV -> TCVx
          Call Timing(CPU1,CPE,TIO1,TIOE)

*         Start Loop on CHFV-iSym-jSym
          Do iSym=1,nSym
          lUCHFV = -1
           If(nBas(iSym).GT.0) then
            Do jSym=1,iSym
            lUCHFV = -1
             If(nBas(jSym).GT.0 .and. MulD2h(iSym,jSym).EQ.iSymL) then
              lUCHFV=7
              iStrtVec_AB = nVec * (iBatch-1) + 1
              Write(CHName,'(A4,I1,I1)') CHNm,iSym,jSym
              Call dAName_MF_WA(lUCHFV,CHName)
              If (iSym.EQ.jSym) then
               Call ChoMP2_TraS(iSymL, iSym,jSym, NumV, CMO,NCMO,
     &                      lUCHFV, iStrtVec_AB, nFVec,nFBatch)
              Else
               Call ChoMP2_TraA(iSymL, iSym,jSym ,NumV, CMO,NCMO,
     &                      lUCHFV, iStrtVec_AB, nFVec,nFBatch)
              EndIf
              Call dAClos(lUCHFV)
             EndIf
            EndDo
           EndIf
          EndDo
*         End Loop on CHFV-iSym-jSym

          Call Timing(CPU2,CPE,TIO2,TIOE)
          CPU_Tra=CPU_Tra+CPU2-CPU1
          TIO_Tra=TIO_Tra+TIO2-TIO1
*   ---   End Transformation of Cholesky Vectors  CHFV -> TCVx

*   ---   Start Generation of Integrals files  TCVx -> MOLINT

*         Start Loop on I, J, A, B Symmetries
          nSymP=(nSym**2+nSym)/2
          Do iSymI = 1, nSym
            Do iSymJ = 1, iSymI
              Do iSymA = 1, nSym
                Do iSymB = 1, nSym
                  iSymAI = MulD2h(iSymA,iSymI)
                  iSymBJ = MulD2h(iSymB,iSymJ)

                  If (iSymAI.EQ.iSymL .and. iSymBJ.EQ.iSymL)
     &              Call ChoMP2_TwoEl(iBatch,nBatch,NumV,LUINTM,
     &                        iAddrIAD2M,iSymI,iSymJ,iSymA,iSymB, iSymL)

                EndDo
              EndDo
            EndDo
          EndDo
*         End Loop on I, J, A, B Symmetries

          Do iSym=1,MaxSym
           Do jSym=1,MaxSym
            If(iMemTCVX(3,iSym,jSym,1).GT.0) then
             iAddr=iMemTCVX(3,iSym,jSym,1)
             iLen =iMemTCVX(3,iSym,jSym,2)
             Call GetMem('iAddr','Free','Real',iAddr,iLen)
             iMemTCVX(3,iSym,jSym,1)=0
             iMemTCVX(3,iSym,jSym,2)=0
            EndIf
           EndDo
          EndDo
          Call Timing(CPU3,CPE,TIO3,TIOE)
          CPU_Gen=CPU_Gen+CPU3-CPU2
          TIO_Gen=TIO_Gen+TIO3-TIO2
*   ---   End Generation of Two-Electrons Integrals files

        ENDDO
*  ---  END LOOP iBatch  -----------------------------------------------

      ENDDO
* *** END MAIN LOOP iSymL on TOTAL SYMMETRY of L (iSym * jSym)   *******

      iAddrIAD2M=0
      Call iDaFile(LUINTM,1,IAD2M,LenIAD2M,iAddrIAD2M)

      Write(6,*)'TIMING INFORMATION:   CPU(s)   %CPU   Elapsed(s)'
      Write(6,'(A,F9.2,1X,F6.1,1X,F12.2)') ' Transformation     ',
     & CPU_Tra, 1.0d2*(CPU_Tra+5.0d-5)/(TIO_Tra+5.0d-5), TIO_Tra
      Write(6,'(A,F9.2,1X,F6.1,1X,F12.2)') ' Generation         ',
     & CPU_Gen, 1.0d2*(CPU_Gen+5.0d-5)/(TIO_Gen+5.0d-5), TIO_Gen
      Call Timing(CPU4,CPE,TIO4,TIOE)
      CPU_Tot=CPU4-CPU0
      TIO_Tot=TIO4-TIO0
      Write(6,'(A,F9.2,1X,F6.1,1X,F12.2)') ' TOTAL              ',
     & CPU_Tot, 1.0d2*(CPU_Tot+5.0d-5)/(TIO_Tot+5.0d-5), TIO_Tot
      Call XFlush(6)

      Return
      End
