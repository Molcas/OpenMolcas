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
* Copyright (C) 2007, Francesco Aquilante                              *
************************************************************************
      SubRoutine Cho_SOSmp2_Drv(irc,EMP2,CMO,EOcc,EVir)
C
C     Francesco Aquilante, May 2007.
C
C     Purpose: driver for computing the Scaled Opposite-Spin (SOS)
C              MP2 energy correction EMP2
C              using Cholesky (or RI) representation for the
C              two-electron integrals.
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
#include "chomp2_cfg.fh"
#include "WrkSpc.fh"

      Character*3  ThisNm
      Character*14 SecNam
      Parameter (SecNam = 'Cho_SOSmp2_Drv', ThisNm = 'Drv')

      Parameter (Chk_Mem_ChoMP2 = 0.123456789D0, Tol = 1.0D-15)
      Parameter (iFmt = 0)

      Logical Delete, Delete_def
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

      Call Cho_SOSmp2_Setup(irc)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': Cho_SOSmp2_Setup returned ',irc
         Go To 1  ! exit
      End If

      If (Verbose) Then
         Call Cho_SOSmp2_Setup_Prt(irc)
         If (irc .ne. 0) Then
            Write(6,*) SecNam,': Cho_SOSmp2_Setup_Prt returned ',irc
            Go To 1  ! exit
         End If
         Call CWTime(CPUIni2,WallIni2)
         Call Cho_PrtTim('Cholesky SOS-MP2 initialization',CPUIni2,
     &                   CPUIni1,WallIni2,WallIni1,iFmt)
      End If

C     Transform Cholesky vectors directly from reduced set to MO
C     representation. Result vectors are stored on disk.
C     Compute also the (ai|ai)^2 diagonal here.
C     ----------------------------------------------------------

      If (Verbose) Then
         Call CWTime(CPUTra1,WallTra1)
      End If
      lDiag = nT1am(1)
      Do iSym = 2,nSym
         lDiag = lDiag + nT1am(iSym)
      End Do
      Call GetMem('Diag','Allo','Real',ipDiag,lDiag)

      Call ChoMP2_TraDrv(irc,CMO,Work(ipDiag),.true.)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': ChoMP2_TraDrv returned ',irc
         Go To 1  ! exit
      End If
C
C     Squaring each diagonal element
C     ------------------------------
      Do ia=0,lDiag-1
         Work(ipDiag+ia)=Work(ipDiag+ia)**2
      End Do
      If (set_cd_thr) ThrMP2=ddot_(lDiag,[1.0d0],0,
     &                                  Work(ipDiag),1)/(5.0d0*lDiag)

      If (Verbose) Then
         Call CWTime(CPUTra2,WallTra2)
         Call Cho_PrtTim('Cholesky MP2 transformation',CPUTra2,CPUTra1,
     &                   WallTra2,WallTra1,iFmt)
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

C     Decompose M(ai,bj) = (ai|bj)^2 .
C     Set number of vectors to be used in energy calculation.
C     -------------------------------------------------------

      If (Verbose) Then
         Call CWTime(CPUDec1,WallDec1)
      End If
      Delete = Delete_def ! delete transf. vector files after dec.
      Call Cho_SOSmp2_DecDrv(irc,Delete,Work(ipDiag))
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': Cho_SOSmp2_DecDrv returned ',irc
         Call ChoMP2_Quit(SecNam,'SOS-MP2 decomposition failed!',' ')
      End If
      If (Verbose) Then
         Call CWTime(CPUDec2,WallDec2)
         Call Cho_PrtTim('Cholesky SOS-MP2 decomposition',
     &                   CPUDec2,CPUDec1,
     &                   WallDec2,WallDec1,iFmt)
      End If
      Call GetMem('Diag','Free','Real',ipDiag,lDiag)

C     Compute SOS-MP2 energy correction.
C     ----------------------------------

      If (Verbose) Then
         Call CWTime(CPUEnr1,WallEnr1)
      End If
      Delete = Delete_def
      Call Cho_SOSmp2_Energy(irc,EMP2,EOcc,EVir,Delete)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': Cho_SOSmp2_Energy returned ',irc
         Go To 1 ! exit
      End If
      If (Verbose) Then
         Call CWTime(CPUEnr2,WallEnr2)
         Call Cho_PrtTim('Cholesky SOS-MP2 energy',
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
         Call Cho_PrtTim('Cholesky SOS-MP2',CPUTot2,CPUTot1,
     &                   WallTot2,WallTot1,iFmt)
      End If
      Call GetMem('Flush','Flush','Real',ip_Dum,l_Dum)
      Call GetMem('Dummy','Free','Real',ip_Dum,l_Dum)
      Call qExit(ThisNm)
      Return
      End
