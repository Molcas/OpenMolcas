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
* Copyright (C) 2004, Giovanni Ghigo                                   *
************************************************************************
      Subroutine Cho_TraCtl(iTraType,LUINTM,CMO,NCMO,DoExch2)
************************************************************************
* Author :  Giovanni Ghigo                                             *
*           Lund University, Sweden                                    *
*           October 2004                                               *
*----------------------------------------------------------------------*
* This is the routine for the transformation from AO basis to MO basis *
* of the Cholesky Full Vectors & Two-electrons integral generation.    *
* Reads the Cholesky Full Vectors (CHFV), generate new Transformed     *
* Cholesky Full Vectors (TCVx) then generate the Two-electrons inte-   *
* gral file (MOLINT).                                                  *
*  <p,q|k,l>  where p,q: All MO; k,l: Occupied                         *
* If DoFull=.True also <p,q|r,s> where p,q,r,s: All MO (for CC)        *
************************************************************************
*
*   <DOC>
*     <Name>Cho\_TraCtl</Name>
*     <Syntax>Call Cho\_TraCtl(iType,LUINTM,CMO,NCMO)
*     </Syntax>
*     <Arguments>
*      \Argument{iType}{Caller program (see below)}{Integer}{in}
*      \Argument{LUINTM}{Unit number of Two-electrons integrals file
*      (MOLINT)}{Integer}{in}
*      \Argument{CMO}{MO coefficients}{Array Real*8}{in}
*      \Argument{NCMO}{Total number of MO coefficients}{Integer}{in}
*      \Argument{DoExch2}{Flag for the generation of Exch-2 integrals}
*      {logical}{in}
*     </Arguments>
*     <Purpose>
*      Driver for the generation of the Two-electrons integrals file
*      (MOLINT) from the AO-based Cholesky Full Vectors.\\
*      Called by TraCtl\_Drv.
*     </Purpose>
*     <Dependencies>
*      The number of frozen and deleted MO used in the post-CASSFC/SCF
*      must be written in the RunFile in nFroPT and nDelPT arrays.\\
*      Routines called: Mem\_Est, Cho\_TraA, Cho\_TraS, Cho\_TwoEl.
*     </Dependencies>
*     <Author>
*      G. Ghigo
*     </Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*      All programs that need the generation of the Two-electrons
*      integrals file in MO basis must tell to the Cholesky routine
*      who they are through iType:
*      \begin{itemize}
*       \item 1 MBPT2
*       \item 2 CASPT2
*       \item 3 MCLR
*       \item 4 CC
*      \end{itemize}
*      The AO-based Cholesky Full Vectors (CHFV) are first transformed
*      in a new set of MO-based Transformed Cholesky Full Vectors
*      (TCVx) then the Two-electrons integrals files (MOLINT) is
*      generated. The integrals in the file (MOLINT) are differentely
*      ordered but addresses can be read from IAD2M as usual.\\
*      The main loop runs over all the symmetries iSymL and for each
*      CHFV$_{pq}$ (p, q are AO indices) the following MO-based
*      Transformed Cholesky Full Vectors (TCVx) can be generated:
*      \begin{itemize}
*       \item  (1) TCVA: L$_{ij}$, L$_{ji}$ if Sym(p).NE.Sym(q)
*       \item  (2) TCVB: L$_{tj}$, L$_{ui}$ if Sym(p).NE.Sym(q)
*       \item  (3) TCVC: L$_{aj}$, L$_{bi}$ if Sym(p).NE.Sym(q)
*       \item  (4) TCVD: L$_{tu}$, L$_{ut}$ if Sym(p).NE.Sym(q)
*       \item  (5) TCVE: L$_{au}$, L$_{bt}$ if Sym(p).NE.Sym(q)
*       \item  (6) TCVF: L$_{ab}$  Only if DoFull=.True. (Not
*                                                      implemented yet!)
*      \end{itemize}
*      MO Indices  i,j: Inactive;   t,u: Active;   a,b: Secondary\\
*      Which TCVx have to be generated is defined in the logical array
*      TCVXist by the routine Mem\_Est and the memory pointer and length
*      are contained in the array iMemTCVX.\\
*      Note that when Sym(p).NE.Sym(q) the order of symmetry is
*      exchanged while the first index is always the MO we excite into
*      and the second index is always the MO we excite out.\\
*      TCV type A and D, L$_{ji}$, and L$_{ut}$, are not really generated
*      because they are just the transposed L$_{ij}$ and L$_{tu}$.
*      However, they are stored in memory.\\
*      The main loop for the generation of the integrals (MOLINT) is
*      based on the total symmetry of the MO-based TCVx and for each
*      symmetry, a batch procedure is performed (nBatch). Within this
*      batch cycle, the TCVx are generated by Cho\_TraS and Cho\_TraA
*      depending if Sym(p).EQ.Sym(q). In both these routines, a second
*      inner-batch procedure (nFBatch) is established reading the Full
*      vectors and tranforming them in the TCVx.
*      Both the number of vectors
*      transformed in the two batch procedure are defined in Mem\_Est
*      (nVec for nBatch and nFVec for nFBatch).\\
*      The MO-based two-electron integrals generation is performed by
*      Cho\_TwoEl. The integrals are (p,k|q,l) where p,q: All MO; k,l:
*      Occupied.
*      For MBPT2 integrals where p,q are occupied are not generated and
*      suitable routines called by ChoMP2\_TraCtl are used.
*      For CC all p,q,k,l are all MO (not implemented yet).
*     </Description>
*    </DOC>
*
******************************************************************
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
      Logical DoExch2
      Logical Found

      MulD2h(i,j) = iEor(i-1,j-1)+1

CGG   ------------------------------------------------------------------
      IfTest=.False.
c      IfTest=.True.
c      DoExch2=.True.
CGG   ------------------------------------------------------------------

      Call QENTER('Cho_TraCtl')
      Call Timing(CPU0,CPE,TIO0,TIOE)

***   INIZIALIZATION   *************************************************

      CALL CWTIME(TCR1,TWR1)
      Call Cho_X_init(irc,0.0)
      If (irc.ne.0) Then
        write(6,*) ' In Cho_TraCtl: Cho_X_Init returned non-zero'//
     &             ' rc = ',irc
        Call Abend()
      EndIf
      Call Cho_X_ReoVec(irc) ! get (if not there) CD vects in full stor
      If (irc.ne.0) Then
        write(6,*) ' In Cho_TraCtl: Cho_X_ReoVec returned non-zero'//
     &             ' rc = ',irc
        Call Abend()
      EndIf
      Call Cho_X_final(irc)
      CALL CWTIME(TCR2,TWR2)
      tcpu_reo=(TCR2-TCR1)
      twal_reo=(TWR2-TWR1)
      write(6,*) ' Reordering of the Cholesky vectors to full storage. '
      write(6,*) ' Elapsed time for the reordering section: ',tcpu_reo
      write(6,*) ' CPU time for the reordering section: ',tcpu_reo
      write(6,*)

* --- Define what has to be calculated. ---
*      DoExc2 flag for the generation of Exch-2 integrals
*      DoTCVA flag for the generation of TCVA
*      DoFull flag for the generation of TCVF
*      DoCoul flag for the generation of coulomb integrals

      DoExc2=DoExch2

*     MBPT2
      If (iTraType.EQ.1) then
        DoTCVA = .False.
        DoFull = .False.
        DoCoul = .False.
      EndIf

*     CASPT2
      If (iTraType.EQ.2) then
        DoTCVA = .True.
        DoFull = .False.
        DoCoul = .True.
      EndIf

*     MCLR
      If (iTraType.EQ.3) then
        DoTCVA = .True.
        DoFull = .False.
        DoCoul = .True.
      EndIf

*     CC
      If (iTraType.EQ.4) then
        DoTCVA = .True.
        DoFull = .True.
        DoCoul = .True.
      EndIf

*     If Tra is not recognized it is assumed is CASPT2 type.
      If (iTraType.GE.5) then
        DoTCVA = .True.
        DoFull = .False.
        DoCoul = .True.
      EndIf

* --- Get Informations from RunFile. ---
*      The following informations must be passed to the Cholesky
*      transformation routine through RunFile. COMMON blocks could
*      not be used due to several conflicts.
      Call Get_iScalar('nSym',nSym)
      Call Get_iArray('nBas',nBas,nSym)
      Call Get_iArray('nFroPT',nFro,nSym)
      Call Get_iArray('nDelPT',nDel,nSym)
      Call Get_iArray('nIsh',nIsh,nSym)
* PAM 2007      Call Get_iArray('nAsh',nAsh,nSym)
* Replaced by:
      Do i=1,nSym
        nAsh(i)=0
      End Do
      Call qpg_iArray('nAsh',Found,nData)
      If(Found .and. nData.eq.nSym) Then
        Call Get_iArray('nAsh',nAsh,nSym)
      End If
* End of replacement
      Call Get_NVnode(NumCho)
      nBasT=0
      Do i=1,nSym
        nBasT   = nBasT + nBas(i)
        nOrb(i) = nBas(i) - nFro(i) - nDel(i)
        nOsh(i) = nIsh(i) + nAsh(i)
        nSsh(i) = nOrb(i) - nOsh(i)
      EndDo

* --- Copy data to common ERI. ---
      NSYMZ=NSYM
      LUINTMZ=LUINTM
      DO I=1,NSYM
        NORBZ(I)=NORB(I)
        NOSHZ(I)=NOSH(I)
      END DO

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

CGG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*)
      Write(6,'(A,8I5)')'           Symmetries :',(i,i=1,nSym)
      Write(6,*)
      Write(6,'(A,8I5)')'               Frozen :',(nFro(i),i=1,nSym)
      Write(6,'(A,8I5)')'         Inactive (I) :',(nIsh(i),i=1,nSym)
      Write(6,'(A,8I5)')'           Active (A) :',(nAsh(i),i=1,nSym)
      Write(6,'(A,8I5)')'        Secondary (S) :',(nSsh(i),i=1,nSym)
      Write(6,'(A,8I5)')'              Deleted :',(nDel(i),i=1,nSym)
      Write(6,*)
      Write(6,'(A,8I5)')'      Total correlated:',(nOrb(i),i=1,nSym)
      Write(6,*)
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------


* *** START LOOP iSymL on TOTAL SYMMETRY of L (iSym * jSym)   **********
      DO iSymL=1,nSym

*       Re-Inizialize the TCVx & iMemTCVX
        Do iType=1,MxTCVx
          Do iSym=1,MaxSym
            Do jSym=1,MaxSym
              TCVXist(iType,iSym,jSym)=.False. ! TCVx existing flag.
              iMemTCVX(iType,iSym,jSym,1)=0 ! Memory Address and
              iMemTCVX(iType,iSym,jSym,2)=0 ! Length in Work(TCVx)
            EndDo
          EndDo
        EndDo
        Call Mem_Est(iSymL,nVec,nFVec)
CGG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*)
      Write(6,*)' TCVx generated in Symmetry',iSymL
      Do i=1,nSym
      Do j=1,nSym
      Do k=1,6
      If (TCVXist(k,i,j)) Write(6,*)' Type=',k,' Symmetries:',i,j
      EndDo
      EndDo
      EndDo
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------
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
CGG   ------------------------------------------------------------------
c      If(IfTest) then
c      Write(6,*)
c      Write(6,*)' iBatch=',iBatch,' of',nBatch,' - NumV=',NumV,
c     &                                             ' - nFBatch=',nFBatch
c      Call XFlush(6)
c      EndIf
CGG   ------------------------------------------------------------------

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
CGG   ------------------------------------------------------------------
c      If(IfTest) then
c      Write(6,*) ' - Open ',CHName,' unit=',lUCHFV,' Vect=',iStrtVec_AB
c      Call XFlush(6)
c      EndIf
CGG   ------------------------------------------------------------------

              If (iSym.EQ.jSym) then
               Call Cho_TraS(iSymL, iSym,jSym, NumV, CMO,NCMO,
     &                      lUCHFV, iStrtVec_AB, nFVec,nFBatch)
              Else
               Call Cho_TraA(iSymL, iSym,jSym ,NumV, CMO,NCMO,
     &                      lUCHFV, iStrtVec_AB, nFVec,nFBatch)
              EndIf

              Call dAClos(lUCHFV)
CGG   ------------------------------------------------------------------
c      If(IfTest) then
c      Write(6,*) ' - Closed ',CHName
c      Call XFlush(6)
c      EndIf
CGG   ------------------------------------------------------------------
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
CGG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*) ' - Generation of Integrals:'
      Call XFlush(6)
      EndIf
CGG   ------------------------------------------------------------------

*         Start Loop on I, J, A, B Symmetries
          nSymP=(nSym**2+nSym)/2
          Do iSymI = 1, nSym
            Do iSymJ = 1, iSymI
              Do iSymA = 1, nSym
                Do iSymB = 1, nSym
                  iSymAI = MulD2h(iSymA,iSymI)
                  iSymBJ = MulD2h(iSymB,iSymJ)

                  If (iSymAI.EQ.iSymL .and. iSymBJ.EQ.iSymL)
     &              Call Cho_TwoEl(iBatch,nBatch,NumV,LUINTM,iAddrIAD2M,
     &                                 iSymI,iSymJ,iSymA,iSymB, iSymL)

                EndDo
              EndDo
            EndDo
          EndDo
*         End Loop on I, J, A, B Symmetries

          Do iType=1,MxTCVx
           Do iSym=1,MaxSym
            Do jSym=1,MaxSym
             If(iMemTCVX(iType,iSym,jSym,1).GT.0) then
              iAddr=iMemTCVX(iType,iSym,jSym,1)
              iLen =iMemTCVX(iType,iSym,jSym,2)
              Call GetMem('iAddr','Free','Real',iAddr,iLen)
              iMemTCVX(iType,iSym,jSym,1)=0
              iMemTCVX(iType,iSym,jSym,2)=0
             EndIf
            EndDo
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
     & CPU_Tra, 1.0d2*CPU_Tra/Max(1.0d0,TIO_Tra), TIO_Tra
      Write(6,'(A,F9.2,1X,F6.1,1X,F12.2)') ' Generation         ',
     & CPU_Gen, 1.0d2*CPU_Gen/Max(1.0d0,TIO_Gen), TIO_Gen
      Call Timing(CPU4,CPE,TIO4,TIOE)
      CPU_Tot=CPU4-CPU0
      TIO_Tot=TIO4-TIO0
      Write(6,'(A,F9.2,1X,F6.1,1X,F12.2)') ' TOTAL              ',
     & CPU_Tot, 1.0d2*CPU_Tot/Max(1.0d0,TIO_Tot), TIO_Tot
      write(6,*)
      Call XFlush(6)
CGG   ------------------------------------------------------------------
      If(IfTest) then
c      IPRX=0    ! Do not print Coulomb nor Exchange Integrals
      IPRX=1    ! Print only Coulomb Integrals
c      IPRX=2    ! Print all Integrals
c      IPRX=3    ! Print only Exchange Integrals
       Call RdInt2(IPRX,DoTCVA)
      EndIf
CGG   ------------------------------------------------------------------

      Call put_tra_comm(IAD2M,NSYMZ,NORBZ,NOSHZ,LUINTMZ)

      Call QEXIT('Cho_TraCtl')
      Return
      End

************************************************************************
      Subroutine put_tra_comm(IBD2M,NSYMX,NORBX,NOSHX,LUINTMX)

      Implicit Real*8 (a-h,o-z)
      Integer IBD2M(3,36*36),NSYMX,NORBX(8),NOSHX(8),LUINTMX
#include "intgrl.fh"

      Do j=1,36*36
         Do i=1,3
            IAD2M(i,j) = IBD2M(i,j)
         End Do
      End Do
      NSYMZ = NSYMX
      Do i=1,8
         NORBZ(i) = NORBX(i)
         NOSHZ(i) = NOSHX(i)
      End Do
      LUINTMZ = LUINTMX

      Return
      End

************************************************************************
      Subroutine Get_NVnode(NVEC)

      Integer NVEC(*)
#include "cholesky.fh"

      Do i=1,nSym
         NVEC(i) =  NumCho(i)
      End Do

      Return
      End
