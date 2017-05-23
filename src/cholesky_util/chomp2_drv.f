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
* Copyright (C) 2004, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_Drv(irc,EMP2,CMO,EOcc,EVir)
C
C     Thomas Bondo Pedersen, October 2004.
C
C     Purpose: driver for computing the MP2 energy correction EMP2
C              using Cholesky decomposed two-electron integrals.
C              Input must have been processed and MO coefficients
C              and orbital energies must be passed as arguments.
C
C     Notes:
C
C       - all MO Cholesky vector files generated here are deleted before
C         exit, except for error terminations (i.e. no cleanup actions
C         are taken!)
C
#include "implicit.fh"
      Dimension CMO(*), EOcc(*), EVir(*)
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2g.fh"
#include "chomp2_cfg.fh"
#include "choorb.fh"
#include "WrkSpc.fh"

      Character*3  ThisNm
      Character*10 SecNam
      Parameter (SecNam = 'ChoMP2_Drv', ThisNm = 'Drx')

      Parameter (Chk_Mem_ChoMP2 = 0.123456789D0, Tol = 1.0D-15)
      Parameter (iFmt = 0)

      Logical DoSort, Delete, Delete_def
      Parameter (Delete_def = .true.)

#if defined (_DEBUG_)
      Verbose = .true.
#endif
      If (Verbose) Then
         Call CWTime(CPUTot1,WallTot1)
      End If

C     Initializations.
C     ----------------

      Call qEnter(ThisNm)
      irc = 0

      EMP2 = 0.0d0
      If(DoDens) Then
         EMP2_dens = 0.0d0
      End If

      If (Verbose) Then
         Call CWTime(CPUIni1,WallIni1)
      End If

      l_Dum = 1
      Call GetMem('Dummy','Allo','Real',ip_Dum,l_Dum)
      Work(ip_Dum) = Chk_Mem_ChoMP2

      FracMem = 0.0d0 ! no buffer allocated
      Call Cho_X_Init(irc,FracMem)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': Cho_X_Init returned ',irc
         Call ChoMP2_Quit(SecNam,'Cholesky initialization error',' ')
      End If

      Call ChoMP2_Setup(irc)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': ChoMP2_Setup returned ',irc
         Go To 1  ! exit
      End If

      If(DoDens) Then
C        Write(6,*) 'Run ChoMP2g_setup'
         Call ChoMP2g_Setup(irc, EOcc, EVir)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2g_Setup returned ',irc
            Go To 1             ! exit
         End If
      End If

      If (Verbose) Then
         Call ChoMP2_Setup_Prt(irc)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2_Setup_Prt returned ',irc
            Go To 1  ! exit
         End If
         Call CWTime(CPUIni2,WallIni2)
         Call Cho_PrtTim('Cholesky MP2 initialization',CPUIni2,CPUIni1,
     &                   WallIni2,WallIni1,iFmt)
      End If

C     Transform Cholesky vectors directly from reduced set to MO
C     representation. Result vectors are stored on disk.
C     If decomposition of (ai|bj) is requested, compute also the
C     (ai|ai) diagonal here.
C     ----------------------------------------------------------

      If (Verbose) Then
         Call CWTime(CPUTra1,WallTra1)
      End If
      If (DecoMP2) Then
         lDiag = nT1am(1)
         Do iSym = 2,nSym
            lDiag = lDiag + nT1am(iSym)
         End Do
      Else If(DoDens) Then
         lDiag = nMoMo(1,6)
         Do iSym = 2, nSym
            lDiag = lDiag + nMoMo(iSym,6)
         End Do
      Else
         lDiag = 1
      End If
      Call GetMem('Diag','Allo','Real',ipDiag,lDiag)
*
      If(.not.DoDens) Then
         Call ChoMP2_TraDrv(irc,CMO,Work(ipDiag),DecoMP2)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2_TraDrv returned ',irc
            Go To 1             ! exit
         End If
         If (Verbose) Then
            Call CWTime(CPUTra2,WallTra2)
            Call Cho_PrtTim('Cholesky MP2 transformation',CPUTra2,
     &                       CPUTra1,WallTra2,WallTra1,iFmt)
         End If
      Else If(DoDens) Then
         Call ChoMP2g_TraDrv(irc,CMO,Work(ipDiag),DecoMP2)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2g_TraDrv returned ',irc
            Go To 1             ! exit
         End If
         If (Verbose) Then
            Call CWTime(CPUTra2,WallTra2)
            Call Cho_PrtTim('Cholesky MP2 transformation',CPUTra2,
     &                       CPUTra1,WallTra2,WallTra1,iFmt)
         End If
      End If

C     Finalize Cholesky info (to release memory).
C     Retain essential info: LuPri, nSym, and NumCho(*).
C     --------------------------------------------------

      nSym_Sav = nSym
      Call iCopy(nSym,NumCho,1,nMP2Vec,1)

      Call Cho_X_Final(irc)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': Cho_X_Final returned ',irc
         Go To 1 ! exit
      End If

      LuPri = 6
      nSym  = nSym_Sav
      Call iCopy(nSym,nMP2Vec,1,NumCho,1)

C     Decompose (ai|bj) integrals, if requested.
C     Set number of vectors to be used in energy calculation.
C     -------------------------------------------------------

      If (DecoMP2) Then
         If (Verbose) Then
            Call CWTime(CPUDec1,WallDec1)
         End If
         Delete = Delete_def ! delete transf. vector files after dec.
         Call ChoMP2_DecDrv(irc,Delete,Work(ipDiag),'Integrals')
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2_DecDrv returned ',irc
            Call ChoMP2_Quit(SecNam,'MP2 decomposition failed!',' ')
         End If
         If (Verbose) Then
            Call CWTime(CPUDec2,WallDec2)
            Call Cho_PrtTim('Cholesky MP2 decomposition',
     &                      CPUDec2,CPUDec1,
     &                      WallDec2,WallDec1,iFmt)
         End If
      Else If(DoDens) Then
         If (Verbose) Then
            Call CWTime(CPUDec1,WallDec1)
         End If
         Call ChoMP2g_AmpDiag(irc,ipDiag,EOcc,EVir)
         Delete = .false.    ! do not delete transf. vectors.
         Call ChoMP2_DecDrv(irc,Delete,Work(ipDiag),'Amplitudes')
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2_DecDrv returned ',irc
            Call ChoMP2_Quit(SecNam,'MP2 decomposition failed!',' ')
         End If
         If (Verbose) Then
            Call CWTime(CPUDec2,WallDec2)
            Call Cho_PrtTim('Cholesky MP2 decomposition',
     &                      CPUDec2,CPUDec1,
     &                      WallDec2,WallDec1,iFmt)
         End If
      Else
         Call iCopy(nSym,NumCho,1,nMP2Vec,1)
      End If
      Call GetMem('Diag','Free','Real',ipDiag,lDiag)

C     Presort Cholesky vectors if needed.
C     -----------------------------------

      DoSort = nBatch .gt. 1
      If (DoSort.and. (.not.DoDens)) Then
         If (Verbose) Then
            Call CWTime(CPUSrt1,WallSrt1)
         End If
         Delete = Delete_def
         Call ChoMP2_SrtDrv(irc,Delete)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2_SrtDrv returned ',irc
            If (Delete) Then ! full vectors not available
               Call ChoMP2_Quit(SecNam,'MP2 presort failed!',' ')
            Else
               Write(6,*) SecNam,
     &                    ': trying to use full vectors instead...'
            End If
            DoSort = .false.
         End If
         If (Verbose) Then
            Call CWTime(CPUSrt2,WallSrt2)
            Call Cho_PrtTim('Cholesky MP2 presort',
     &                      CPUSrt2,CPUSrt1,
     &                      WallSrt2,WallSrt1,iFmt)
         End If
      End If

C     FNO section: MP2 pseudodensity
C     ------------------------------

      If (DoFNO.and.(.not.DoDens)) Then
         If (Verbose) Then
            Call CWTime(CPUDab1,WallDab1)
         End If
         Delete = Delete_def
         Call ChoMP2_FNO(irc,Work(ip_Dab),Work(ip_Dii),
     &                       EOcc,EVir,DoSort,Delete)
         call dscal_(l_Dii,-1.0d0,Work(ip_Dii),1)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2_FNO returned ',irc
            Go To 1 ! exit
         End If
         If (Verbose) Then
            Call CWTime(CPUDab2,WallDab2)
            Call Cho_PrtTim('Cholesky MP2 FNO section ',
     &                      CPUDab2,CPUDab1,
     &                      WallDab2,WallDab1,iFmt)
         End If
         Go To 1 ! exit
      End If

*     Compute MP2 Density and energy correction.
*     ------------------------------------------
      If (DoDens) Then
         If (Verbose) Then
            Call CWTime(CPUEnr1,WallEnr1)
         End If
         Delete = .false.
         Call ChoMP2g_DensDrv(irc,Work(ip_EOccu),Work(ip_EVirt),
     &                       Work(ip_EFroz),CMO)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2g_DensDrv returned ',irc
            Go To 1             ! exit
         End If
      End If

*     Compute some matrices for Mp2-gradients
*     ---------------------------------------
      If(DoGrdt) Then
         If (Verbose) Then
            Call CWTime(CPUEnr1,WallEnr1)
         End If
         Call ChoMP2g_GradSetup(irc,CMO)
         If(irc .ne. 0) Then
            Write(6,*) SecNam, ':ChoMP2g_GradSetup returned ', irc
            Go To 1 ! exit
         End If
         If (Verbose) Then
            Call CWTime(CPUEnr2,WallEnr2)
            Call Cho_PrtTim('Cholesky Grad setup',
     &                      CPUEnr2,CPUEnr1,
     &                      WallEnr2,WallEnr1,iFmt)
         End If
      End If


C     Compute MP2 energy correction.
C     ------------------------------

      If (Verbose) Then
         Call CWTime(CPUEnr1,WallEnr1)
      End If
      Delete = Delete_def
      If (Laplace .and. SOS_MP2) Then
         Call ChoLSOSMP2_Energy(irc,EMP2,EOcc,EVir,DoSort,Delete)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoLSOSMP2_Energy returned ',irc
            Go To 1 ! exit
         End If
      Else
         Call ChoMP2_Energy(irc,EMP2,EOcc,EVir,DoSort,Delete)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': ChoMP2_Energy returned ',irc
            Go To 1 ! exit
         End If
      End If
      If (Verbose) Then
         Call CWTime(CPUEnr2,WallEnr2)
         Call Cho_PrtTim('Cholesky MP2 energy',
     &                   CPUEnr2,CPUEnr1,
     &                   WallEnr2,WallEnr1,iFmt)
      End If

C     Exit.
C     -----

    1 Diff = abs(Work(ip_Dum)-Chk_Mem_ChoMP2)
      If (Diff .gt. Tol) Then
         Write(6,*) SecNam,': Memory Boundary Error!'
         If (irc .eq. 0) irc = -9999
      End If
      If (Verbose) Then
         Call CWTime(CPUTot2,WallTot2)
         Call Cho_PrtTim('Cholesky MP2',CPUTot2,CPUTot1,
     &                   WallTot2,WallTot1,iFmt)
      End If
      Call GetMem('Flush','Flush','Real',ip_Dum,l_Dum)
      Call GetMem('Dummy','Free','Real',ip_Dum,l_Dum)
      Call qExit(ThisNm)
      Return
      End
