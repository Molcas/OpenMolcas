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
************************************************************************
*  Cho_X_Init
*
*> @brief
*>   Initialize Cholesky vector information for external use.
*> @author Thomas Bondo Pedersen
*>
*> @details
*> This routine reads and processes the information
*> stored on the runfile/restart files by the Cholesky
*> decomposition utility.
*> This routine is also used for setting up the environment in
*> density fitting (DF or RI) runs.
*> Index arrays are allocated and initialized.
*> All information is stored in the include files
*> cholesky.fh, choorb.fh, choarr.f90, and choswp.f90.
*>
*> \p BufFrac is the fraction of total available memory that will be
*> allocated as Cholesky vector buffer. For example, \p BufFrac = ``0.35``
*> implies that 35% of the total available memory (after the
*> allocations of this routine) will be allocated for Cholesky
*> vectors. The vectors will be read into the buffer as part of
*> the initialization. The reading routines can then be used as
*> usual; the buffer is automatically taken care of by the reading
*> routines. If \p BufFrac is less than or equal to zero, no buffer
*> will be used.
*>
*> Return codes:
*>
*> - \p irc = ``-2``: Local DF not implemented here
*> - \p irc = ``-1``: Cholesky flag not found on runfile
*> - \p irc =  ``0``: initialization success
*> - \p irc =  ``1``: runfile info corrupted
*> - \p irc =  ``2``: restart file info corrupted
*> - \p irc =  ``3``: inconsistent include file(s) detected (typically an internal error/bug)
*> - \p irc =  ``4``: error in parallel setup
*>
*> @note
*> The two-electron repulsion integrals must have been decomposed by Seward.
*>
*> @param[out] irc     Return code
*> @param[in]  BufFrac Fraction of memory to be used as buffer
************************************************************************
      Subroutine Cho_X_Init(irc,BufFrac)
      use ChoArr, only: iSOShl, iBasSh, nBasSh, nBstSh, iSP2F, iShlSO,
     &                  iRS2F, nDimRS
      use ChoSwp, only: nnBstRSh, nnBstRSh_Hidden
      use ChoSwp, only: iiBstRSh, iiBstRSh_Hidden
      use ChoSwp, only:   IndRSh,   IndRSh_Hidden
      use ChoSwp, only:   IndRed,   IndRed_Hidden
#include "implicit.fh"
#include "choorb.fh"
#include "cholesky.fh"
#include "choptr2.fh"
#include "chosp.fh"
#include "choini.fh"
#include "choprint.fh"
#include "chobkm.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Character*10 SecNam
      Parameter (SecNam = 'Cho_X_Init')

#if defined (_DEBUGPRINT_)
      Character*2 Unt
#endif

      Logical DidCholesky, FirstCall
      Logical isDF, isLocalDF, DoDummy

      Integer ChoIsIni, ip, l

      Integer  Cho_iSumElm
      External Cho_iSumElm

      Save FirstCall
      Data FirstCall /.true./

C     Register entry.
C     ---------------

#if defined (_DEBUGPRINT_)
      Call GetMem('CXI_MX1','Max ','Real',ip_Max,l_Max)
      Call Cho_Word2Byte(l_Max,8,Byte,Unt)
      Write(6,*) '>>>>> Available memory on entry to ',SecNam,': ',
     &           l_Max,' = ',Byte,Unt
#endif

C     Check that this is a Cholesky run.
C     ----------------------------------

      Call DecideOnCholesky(DidCholesky)
      If (.not. DidCholesky) Go To 100

C     Check if already initialized.
C     -----------------------------

      If (FirstCall) Then ! it cannot be already done
         FirstCall = .false.
      Else ! might be already done
         Call Get_iScalar('ChoIni',ChoIsIni)
         If (ChoIsIni .eq. ChoIniCheck) Then ! already done
            irc = 0
            Go To 1
         End If
      End If

C     Check if this is density fitting (DF).
C     --------------------------------------

      Call DecideOnDF(isDF)
      If (isDF) Then
         Call DecideOnLocalDF(isLocalDF)
         If (isLocalDF) Go To 99
      End If

C     Define all entries in include files choorb.fh, cholesky.fh,
C     and choprint.fh.
C     -------------------------------------------------------------

      Call Cho_X_SetInc(irc)
      If (irc .ne. 0) Then
         Go To 103  ! include file inconsistency detected
      End If

C     Set parallel info (picked up from para_info).
C     -------------------------------------------------

      CHO_FAKE_PAR = .False.
      Call Cho_ParConf(CHO_FAKE_PAR)

C     Define entries in choptr2.fh.
C     ------------------------------

      Call Cho_SetPtr2()

C     Set run mode to "external".
C     ---------------------------

      Run_Mode = Run_External

C     Set output unit used by the decomposition core routines.
C     --------------------------------------------------------

      LuPri = 6

C     Set print level to -5 (ensuring that Cho_X_Checkdiag will
C     print information, if called). All other routines will be
C     silent.
C     ---------------------------------------------------------

      iPrint=-5

C     Get nSym: the number of irreps.
C     -------------------------------

      Call Get_iScalar('nSym',nSym)
      If (nSym.lt.1 .or. nSym.gt.8) Then
         Write(6,*) SecNam,': nSym out of bounds: ',nSym
         Go To 101
      End If

C     Get Cho_AdrVec: addressing of vector files (1=WA, 2=DA).
C     --------------------------------------------------------

      Call Get_iScalar('ChoVec Address',Cho_AdrVec)

C     Open files with red. set and vector info.
C     -----------------------------------------

      Call Cho_UnIni()
      Call Cho_OpenVR(1,2)

C     Set vector I/O model etc.
C     -------------------------

      Cho_IOVec = 3
      N1_VecRd  = 2
      N2_VecRd  = 3
      nSys_Call = 0
      nDGM_Call = 0

C     Get info (derived) from the runfile.
C     nBas  : #basis functions in each irrep
C     iSOShl: shell index for each basis function (SO)
C     NumCho: #Cholesky vectors in each irrep
C     MaxVec: max. element in NumCho (used to allocate InfVec)
C     --------------------------------------------------------

      Call Get_iArray('nBas',nBas,nSym)
      iBas(1) = 0
      nBasT   = nBas(1)
      Do iSym = 2,nSym
         iBas(iSym) = nBasT
         nBasT = nBasT + nBas(iSym)
      End Do
      If (nBasT .lt. 1) Then
         Write(6,*) SecNam,': nBasT out of bounds: ',nBasT
         Go To 101
      End If
      Call mma_allocate(iSOShl,nBasT,Label='iSOShl')
      Call Get_iArray('ISOSHL',iSOShl,nBasT)

      Call Get_iArray('NumCho',NumCho,nSym)
      NumChT = Cho_iSumElm(NumCho,nSym)
      MaxVec = NumCho(1)
      Do iSym = 2,nSym
         MaxVec = max(MaxVec,NumCho(iSym))
      End Do

C     Read Cholesky restart file. Allocate InfRed and InfVec.
C     nShell: #shells
C     nnShl_Tot : total #shell pairs
C     nnShl : #shell pairs contributing in diagonal
C     MaxRed: #reduced sets (used to allocate InfRed)
C     InfRed: InfRed(i) is the disk address of reduced set i
C     InfVec: InfVec(i,1,iSym) is the parent index of vector i in irrep
C                              iSym in first reduced set
C             InfVec(i,2,iSym) is the reduced set of this vector
C             InfVec(i,3,iSym) is the disk address for reading this
C                              vector
C             InfVec(i,4,iSym) is the WA disk address of this vector.
C             InfVec(i,5,iSym) in a parallel run is the global index
C                              of the i-th vector in iSym for this node
C                              (in serial: InfVec(i,5,iSym) = i)
C     Note: InfVec(i,5,iSym) is treated in a quite dirty way here:
C       for DF:
C           not defined in serial, thus it will be defined by the call
C           to Cho_X_DefineInfVec_5.
C           defined in parallel (simply read from disk and not modified
C           by Cho_X_DefineInfVec_5).
C       for Cholesky:
C           not defined in serial, thus it will be defined by the call
C           to Cho_X_DefineInfVec_5.
C           not defined in parallel, thus it will be defined by the call
C           to Cho_X_DefineInfVec_5. It is later modified by
C           Cho_X_Init_Par.
C     ------------------------------------------------------------------

      ierr = 0
      Call Cho_X_RdRst(ierr)
      If (ierr .ne. 0) Go To 102
      nnShl_Tot = nShell*(nShell+1)/2
      Call Cho_X_DefineInfVec_5(isDF)

C     nnShl_SP makes it possible to use function Cho_F2SP.
C     (Stored in chosp.fh)
C     ----------------------------------------------------

      nnShl_SP = nnShl

C     Allocate and initialize index arrays.
C     -------------------------------------

      call mma_allocate(iiBstRSh_Hidden,nSym,nnShl,3,
     &                  Label='iiBstRSh_Hidden')
      iiBstRSh => iiBstRSh_Hidden
      call mma_allocate(nnBstRSh_Hidden,nSym,nnShl,3,
     &                  Label='nnBstRSh_Hidden')
      nnBstRSh => nnBstRSh_Hidden
      Call Cho_RstD_GetInd1()
      mmBstRT = nnBstRT(1)

      Call mma_allocate(IndRed_Hidden,nnBstRT(1),3,
     &                  Label='IndRed_Hidden')
      IndRed => IndRed_Hidden
      Call mma_allocate(IndRSh_Hidden,nnBstRT(1),Label='IndRSh_Hidden')
      IndRSh => IndRSh_Hidden
      Call Cho_RstD_GetInd2()

      Call mma_allocate(iSP2F,nnShl,Label='iSP2F')
      Call Cho_RstD_GetInd3(iSP2F,SIZE(iSP2F))

C     Allocate and read bookmarks (if available on runfile).
C     ------------------------------------------------------

      If (isDF) Then
         ip_BkmVec=0
         l_BkmVec=0
         nRow_BkmVec=0
         nCol_BkmVec=0
         ip_BkmThr=0
         l_BkmThr=0
         nRow_BkmThr=0
         nCol_BkmThr=0
      Else
         l=4
         Call GetMem('BkmDim','Allo','Inte',ip,l)
         Call Get_iArray('Cholesky BkmDim',iWork(ip),l)
         nRow_BkmVec=iWork(ip)
         nCol_BkmVec=iWork(ip+1)
         nRow_BkmThr=iWork(ip+2)
         nCol_BkmThr=iWork(ip+3)
         Call GetMem('BkmDim','Free','Inte',ip,l)
         If (nRow_BkmVec.gt.0 .and. nCol_BkmVec.gt.0 .and.
     &       nRow_BkmThr.gt.0 .and. nCol_BkmThr.gt.0) Then
            l_BkmVec=nRow_BkmVec*nCol_BkmVec
            Call GetMem('BkmVec','Allo','Inte',ip_BkmVec,l_BkmVec)
            Call Get_iArray('Cholesky BkmVec',iWork(ip_BkmVec),l_BkmVec)
            l_BkmThr=nRow_BkmThr*nCol_BkmThr
            Call GetMem('BkmVec','Allo','Real',ip_BkmThr,l_BkmThr)
            Call Get_dArray('Cholesky BkmThr',Work(ip_BkmThr),l_BkmThr)
         Else
            ip_BkmVec=0
            l_BkmVec=0
            nRow_BkmVec=0
            nCol_BkmVec=0
            ip_BkmThr=0
            l_BkmThr=0
            nRow_BkmThr=0
            nCol_BkmThr=0
         End If
      End If

C     mySP is set because a few core routines may use it.
C     After the decomposition is done, it must be a trivial mapping and
C     the user (programmer) should not worry about it at all.
C     -----------------------------------------------------------------

      l_mySP = nnShl
      Call GetMem('mySP','Allo','Inte',ip_mySP,l_mySP)
      Do ijShl = 1,nnShl
         iWork(ip_mySP-1+ijShl) = ijShl
      End Do

C     Copy reduced set 1 to location 2.
C     ---------------------------------

      Call Cho_RSCopy(iiBstRSh,nnBstRSh,
     &                IndRed,1,2,nSym,nnShl,mmBstRT,3)

C     Get dimensions of reduced sets.
C     -------------------------------

      Call mma_allocate(nDimRS,nSym,MaxRed,Label='nDimRS')
      Call iCopy(nSym,nnBstR(1,1),1,nDimRS,1)
      iLoc = 3
      Do iRed = 2,MaxRed
         Call Cho_GetRed(iRed,iLoc,.false.)
         Call Cho_SetRedInd(iiBstRSh,nnBstRSh,nSym,nnShl,iLoc)
         Call iCopy(nSym,nnBstR(1,iLoc),1,nDimRS(:,iRed),1)
      End Do

C     Copy reduced set 1 to location 3.
C     ---------------------------------

      Call Cho_RSCopy(iiBstRSh,nnBstRSh,
     &                IndRed,1,3,nSym,nnShl,mmBstRT,3)

C     Derive:
C     nBasSh: #basis functions in sym. block of shell.
C     nBstSh: #basis functions in each shell.
C     Mx2Sh : max. shell pair dimension.
C     iShlSO: index of SO within its shell.
C     ------------------------------------------------

      Call mma_allocate(iBasSh,nSym, nShell,Label='iBasSh')
      Call mma_allocate(nBasSh,nSym, nShell,Label='nBasSh')
      Call mma_allocate(nBstSh,nShell,Label='nBstSh')
      Call Cho_SetSh(iBasSh,nBasSh,nBstSh,
     &               iBas,nBas,iSOShl,nSym,nShell,nBasT)

      MxOrSh = nBstSh(1)
      Do iShl = 2,nShell
         MxOrSh = max(MxOrSh,nBstSh(iShl))
      End Do

      Mx2Sh = 0
      Do ijShl = 1,nnShl
         Call Cho_InvPck(iSP2F(ijShl),iShl,jShl,.True.)
         If (iShl .eq. jShl) Then
            Numij = nBstSh(iShl)*(nBstSh(iShl)+1)/2
         Else
            Numij = nBstSh(iShl)*nBstSh(jShl)
         End If
         Mx2Sh = max(Mx2Sh,Numij)
      End Do

      Call mma_allocate(iShlSO,nBasT,Label='iShlSO')
      Call Cho_SetSh2(iShlSO,iSOShl,nBstSh,nBasT,nShell)

C     Allocate and compute mapping RS1->Full.
C     ---------------------------------------

      Call mma_allocate(iRS2F,2,nnBstRT(1),Label='iRS2F')
      Call Cho_RSTOF(iRS2F,2,nnBstRT(1),1)

C     Allocate integer scratch array used for Cholesky vector I/O.
C     ------------------------------------------------------------

      DoDummy = .not.(Cho_IOVec.eq.1 .or. Cho_IOVec.eq.2 .or.
     &                Cho_IOVec.eq.3 .or. Cho_IOVec.eq.4)
      Call Cho_Allo_iScr(DoDummy)

C     Setup for parallel runs.
C     ------------------------

      Call Cho_X_Init_Par(irc,isDF)
      If (irc .ne. 0) Then
         Go To 104
      End If

#if defined (_DEBUGPRINT_)
C     Debug: test bookmarks.
C     Note that 1C-CD flag must be available on runfile
C     (make sure _DEBUGPRINT_ is defined also in Cho_Final().
C     --------------------------------------------------
      If (l_BkmVec.gt.0 .and. l_BkmThr.gt.0) Then
         Call Get_iScalar('1C-CD',is1CCD)
         Call Cho_TestBookmark(irc,.True.,is1CCD.eq.1)
         If (irc.ne.0) Call Cho_Quit('Bookmark test failed!',104)
      End If
#endif

C     Allocate and initialize (i.e. read vectors) vector buffer.
C     ----------------------------------------------------------

      Frac = min(max(BufFrac,0.0d0),1.0d0)
      Call Cho_VecBuf_Init(Frac,nnBstR(1,1)) ! allocate
      Call Cho_VecBuf_Ini2() ! read

C     Normal exit point.
C     ------------------

      ChoIsIni = ChoIniCheck
      Call Put_iScalar('ChoIni',ChoIsIni)
      irc = 0
      Go To 1

C     Error branches.
C     ===============

   99 Continue ! Local DF not implemented
C TODO/FIXME: compute Cholesky vectors from LDF coefficients here?
         irc=-2
         Write(6,'(//,A,A,//)')
     &   SecNam,': Local DF not implemented!'
      Go To 1

  100 Continue ! Cholesky flag not found on runfile
         irc = -1
         Write(6,'(//,A,A,//)')
     &   SecNam,': two-electron integrals not Cholesky decomposed!'
      Go To 1

  101 Continue ! Bad info obtained from runfile
         irc = 1
         Write(6,'(//,A,A,//)')
     &   SecNam,': WARNING: error reading runfile!'
      Go To 1

  102 Continue ! restart info corrupted
         irc = 2
         Write(6,'(//,A,A)')
     &   SecNam,': WARNING: error reading restart info!'
         Write(6,'(A,A,I6,//)')
     &   SecNam,': return code from read:',ierr
      Go To 1

  103 Continue ! include file inconsistency detected
         irc = 3
         Write(6,'(//,A,A,//)')
     &   SecNam,': WARNING: include file inconsistency detected!'
      Go To 1

  104 Continue ! error in parallel setup
         irc = 4
         Write(6,'(//,A,A,//)')
     &   SecNam,': WARNING: error in parallel setup!'
      Go To 1

C     Return.
C     =======

    1 Continue
#if defined (_DEBUGPRINT_)
      Call GetMem('CXI_MX2','Max ','Real',ip_Max,l_Max)
      Call Cho_Word2Byte(l_Max,8,Byte,Unt)
      Write(6,*) '>>>>> Available memory on exit from ',SecNam,': ',
     &           l_Max,' = ',Byte,Unt
      Call xFlush(6)
#endif
      End


      SubRoutine Cho_X_DefineInfVec_5(isDF)
C
C     Purpose: Trivial definition of location 5 of InfVec:
C              InfVec(i,5,iSym) = i
C              The routine does nothing in case of parallel DF.
C
      Use Para_Info, Only: Is_Real_Par
      use ChoSwp, only: InfVec
      Implicit None
      Logical isDF
#include "cholesky.fh"

      Integer iSym, i
      Logical doDefine

C Define in case of
C 1) serial Cholesky
C 2) serial DF
C 3) parallel Cholesky
C Do NOT define for parallel DF.
      doDefine = .not.Is_Real_Par() .or.
     &           (Is_Real_Par() .and. .not.isDF)
      If (doDefine) Then
         Do iSym = 1,nSym
            Do i = 1,NumCho(iSym)
               InfVec(i,5,iSym) = i
            End Do
         End Do
      End If

      End
