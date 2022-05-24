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
* Copyright (C) 2008, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine ChoMP2_O4_Drv(irc,EMP2,CMO,EOcc,EVir)
C
C     Thomas Bondo Pedersen, Jan. 2008.
C     - based on ChoMP2_Drv() by T. B. Pedersen.
C
C     Purpose: driver for computing the MP2 energy correction EMP2
C              using Cholesky decomposed two-electron integrals
C              in a quartic scaling fashion.
C              Input must have been processed and MO coefficients
C              and orbital energies must be passed as arguments.
C
C     Notes:
C
C       - all MO Cholesky vector files generated here are deleted before
C         exit, except for error terminations (i.e. no cleanup actions
C         are taken!)
C       - the amplitude decomposition does NOT use EOcc and EVir; these
C         are accessed through pointers ip_EOc and ip_EVir stored in
C         chomp2_dec.fh
C
#include "implicit.fh"
      Dimension CMO(*), EOcc(*), EVir(*)
#include "cholesky.fh"
#include "chomp2.fh"
#include "chomp2_cfg.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"

      Character(LEN=6), Parameter:: ThisNm = 'O4_Drv'
      Character(LEN=13), Parameter:: SecNam = 'ChoMP2_O4_Drv'

      Real*8, Parameter:: Chk_Mem_ChoMP2 = 0.123456789D0, Tol = 1.0D-15
      Integer, Parameter:: iFmt = 0

      Logical Delete, DoAmpDiag
      Logical, Parameter:: Delete_def = .true.

      Integer a, ai
      Integer lU_AO(8)

      Character(LEN=3) BaseName_AO
      Real*8, Allocatable:: Check(:)

      MulD2h(k,l)=iEor(k-1,l-1)+1

#if defined (_DEBUGPRINT_)
      Verbose = .true.
#endif
      If (Verbose) Then
         Call CWTime(CPUTot1,WallTot1)
      End If

C     Initializations.
C     ----------------

      irc = 0

      EMP2 = 0.0d0

      If (Verbose) Then
         Call CWTime(CPUIni1,WallIni1)
      End If

      Call mma_allocate(Check,1,Label='Check')
      Check(1) = Chk_Mem_ChoMP2

      FracMem = 0.0d0 ! no buffer allocated
      Call Cho_X_Init(irc,FracMem)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': Cho_X_Init returned ',irc
         Call ChoMP2_Quit(SecNam,'Cholesky initialization error',' ')
      End If

C-TBP:
C Frankie,
C The setup is still the same here (i.e. batching etc. is included)
C I don't use the batching info at all, though, so it should be save
C to remove it - unless you need it, of course.

      Call ChoMP2_Setup(irc)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': ChoMP2_Setup returned ',irc
         Go To 1  ! exit
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
C     Compute also amplitude diagonal here.
C     ----------------------------------------------------------

      If (Verbose) Then
         Call CWTime(CPUTra1,WallTra1)
      End If

      lDiag = nT1am(1)
      Do iSym = 2,nSym
         lDiag = lDiag + nT1am(iSym)
      End Do

      Call GetMem('Diag','Allo','Real',ipDiag,lDiag)

      Call ChoMP2_TraDrv(irc,CMO,Work(ipDiag),.True.)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': ChoMP2_TraDrv returned ',irc
         Go To 1  ! exit
      End If
      kD0 = ipDiag - 1
      Do iSym = 1,nSym
         Do iSymi = 1,nSym
            iSyma = MulD2h(iSymi,iSym)
            kD1 = kD0 + iT1Am(iSyma,iSymi)
            Do i = 1,nOcc(iSymi)
               kD2 = kD1 + nVir(iSyma)*(i-1)
               Ei  = EOcc(iOcc(iSymi)+i)
               Do a = 1,nVir(iSyma)
                  ai = kD2 + a
                  DE = 2.0d0*(EVir(iVir(iSyma)+a)-Ei)
                  Work(ai) = Work(ai)/DE
               End Do
            End Do
         End Do
         kD0 = kD0 + nT1Am(iSym)
      End Do

      If (Verbose) Then
         Call CWTime(CPUTra2,WallTra2)
         Call Cho_PrtTim('Cholesky MP2 transformation',CPUTra2,CPUTra1,
     &                   WallTra2,WallTra1,iFmt)
      End If

C     Decompose MP2 amplitudes (times -1).
C     ------------------------------------

      If (Verbose) Then
         Call CWTime(CPUDec1,WallDec1)
      End If

C-TBP:
C Frankie,
C I just modified the decomposition slightly so that it treats
C amplitudes, too. You should be aware, though, that result vectors
C are always written on the same files, so if you do another
C decomposition of, say, integrals (or squared integrals), then
C the vector files will be overwritten!!
C The number of vectors is always written to nMP2Vec(iSym) in
C chomp2.fh - this is overwritten too, if you do another CD!!

      Delete = Delete_def ! delete transf. vector files after dec.
      Call ChoMP2_DecDrv(irc,Delete,Work(ipDiag),'Amplitudes')
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': ChoMP2_DecDrv returned ',irc
         Call ChoMP2_Quit(SecNam,'MP2 decomposition failed!',
     &                    ' ')
      End If
      Call GetMem('Diag','Free','Real',ipDiag,lDiag)

      If (Verbose) Then
         Call CWTime(CPUDec2,WallDec2)
         Call Cho_PrtTim('Cholesky MP2 decomposition',
     &                   CPUDec2,CPUDec1,
     &                   WallDec2,WallDec1,iFmt)
      End If

C     Backtransform amplitude vectors to AO basis.
C     Calculate also backtransformed amplitude diagonal.
C     --------------------------------------------------

      If (Verbose) Then
         Call CWTime(CPUBT1,WallBT1)
      End If

      iTyp = 2 ! type of MO vectors (i.e. amp. vectors)
      Delete = Delete_def ! delete amp. vectors after backtransf.
      BaseName_AO = 'AAO' ! basename for multifiles of AO vectors
      DoAmpDiag = .True. ! calculate backtransf. amp. diagonal
      lDiag = 0
      Do iSym = 1,nSym
         Do iSymb = 1,nSym
            iSyma = MulD2h(iSymb,iSym)
            lDiag = lDiag + nBas(iSyma)*nBas(iSymb)
         End Do
      End Do

      Call GetMem('AODiag','Allo','Real',ipDiag,lDiag)
      Call ChoMP2_VectorMO2AO(iTyp,Delete,BaseName_AO,CMO,DoAmpDiag,
     &                        Work(ipDiag),lDiag,lU_AO,irc)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': ChoMP2_VectorMO2AO returned ',irc
         Call ChoMP2_Quit(SecNam,
     &                'MP2 amplitude vector backtransformation failed!',
     &                ' ')
      End If
      Call GetMem('Diag','Free','Real',ipDiag,lDiag)

      If (Verbose) Then
         Call CWTime(CPUBT2,WallBT2)
         Call Cho_PrtTim('Cholesky MP2 backtransformation',
     &                   CPUBT2,CPUBT1,
     &                   WallBT2,WallBT1,iFmt)
      End If

C-TBP:
C Frankie,
C You can now read the AO vectors from the units lU_AO(iSym)
C using ddaFile(). The files are word-addressable so that it
C is possible to read from the file as if addressing an array.
C The number of vectors is nMP2Vec(iSym) stored in chomp2.fh.
C To save memory, you may want to finalize Cholesky info before
C continuing with the backtransformed vectors - but remember
C that all Cholesky information (from the AO integral CD) is
C then lost (f.ex. NumCho(iSym) becomes useless).

C     Finalize Cholesky info.
C     -----------------------

      Call Cho_X_Final(irc)
      If (irc .ne. 0) Then
         Write(6,*) SecNam,': Cho_X_Final returned ',irc
         Go To 1 ! exit
      End If

C     Close and delete files containing backtransformed amplitude
C     vectors.
C     -----------------------------------------------------------

      Do iSym = 1,nSym
         Call daEras(lU_AO(iSym))
      End Do

C     Exit.
C     -----

    1 Continue
      Diff = abs(Check(1)-Chk_Mem_ChoMP2)
      If (Diff .gt. Tol) Then
         Write(6,*) SecNam,': Memory Boundary Error!'
         If (irc .eq. 0) irc = -9999
      End If
      If (Verbose) Then
         Call CWTime(CPUTot2,WallTot2)
         Call Cho_PrtTim('Cholesky MP2',CPUTot2,CPUTot1,
     &                   WallTot2,WallTot1,iFmt)
      End If

      Call mma_deallocate(Check)
      End
