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
* Copyright (C) 2007,2008, Thomas Bondo Pedersen                       *
************************************************************************
      SubRoutine ChoMP2_VectorMO2AO(iTyp,Delete,BaseName_AO,CMO,DoDiag,
     &                              Diag,lDiag,lU_AO,irc)
C
C     Thomas Bondo Pedersen, Dec. 2007 - Jan. 2008.
C
C     Purpose: backtransform vectors from MO to AO basis.
C
C     Input:
C        iTyp......... specifies type of vectors (iTyp is used to open
C                      MO vector files through ChoMP2_OpenF()).
C        Delete....... Flag specifiyng whether MO vector files are to
C                      be deleted before exiting this routine.
C        BaseName_AO.. Base name (3 characters) for the files containing
C                      AO vectors. Symmetry index will be appended!
C        CMO.......... MO coefficient array
C        DoDiag....... if .True., calculate AO diagonal elements as
C                      D(ab) = sum_J L(J,ab)*L(J,ab)
C     Output:
C        Diag......... Contains the diagonal if requested (flag DoDiag).
C        lDiag........ Dimension of Diag
C        lU_AO........ Array containing units of open AO vector files.
C        irc.......... return code.
C                      = 0 if successful, non-zero otherwise
C                      (must be checked by caller).
C
C***********************************************************************
C!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C***********************************************************************
C*** CHOLESKY INFORMATION MUST BE INITIALIZED WHEN CALLING THIS ROUTINE.
C    --> I.e. Cho_X_Init() must have been called.
C    --> It is actually sufficient that choorb.fh + nSym in
C        cholesky.fh are available (which they are when Cho_X_Init()
C        has been called).
C*** CHOLESKY MP2 INFORMATION MUST BE INITIALIZED WHEN CALLING THIS
C    ROUTINE.
C    --> I.e. ChoMP2_Setup() must have been called.
C*** AO VECTORS ARE STORED IN LOWER TRIANGULAR [M(I,J), I.GE.J] FORMAT.
C    --> I.e. vectors are stored as L(J,ab) where a>=b
C*** DIAGONAL IS STORED IN LOWER TRIANGULAR FORMAT.
C    --> I.e. diagonal is stored as D(ab) where a>=b (same as vectors).
C***********************************************************************
C!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C***********************************************************************
C
      Implicit None
      Integer     iTyp
      Logical     Delete
      Character*3 BaseName_AO
      Real*8      CMO(*)
      Logical     DoDiag
      Integer     lDiag
      Real*8      Diag(lDiag)
      Integer     lU_AO(*)
      Integer     irc
#include "cholesky.fh"
#include "choorb.fh"
#include "chomp2.fh"
#include "WrkSpc.fh"

      Character*18 SecNam
      Character*11 ThisNm
      Parameter (SecNam = 'ChoMP2_VectorMO2AO')
      Parameter (ThisNm = 'VectorMO2AO')

      Character*4 FullName_AO

      Integer iSym, iSyma, iSymb, iCount, iOpen, iClose
      Integer ip_COcc, ip_CVir
      Integer l_COcc,  l_CVir

      Logical Debug
#if defined (_DEBUG_)
      Parameter (Debug = .True.)
#else
      Parameter (Debug = .False.)
#endif

      Integer MulD2h, k, l
      MulD2h(k,l)=iEOr(k-1,l-1)+1

C     Initializations.
C     ----------------

      irc = 0
      Do iSym = 1,nSym
         lU_AO(iSym) = -999999
      End Do
      If (DoDiag) Then
         iCount = 0
         Do iSym = 1,nSym
            Do iSymb = 1,nSym
               iSyma = MulD2h(iSymb,iSym)
               iCount = iCount + nBas(iSyma)*nBas(iSymb)
            End Do
         End Do
         If (iCount .ne. lDiag) Then
            Write(6,*) SecNam,': WARNING: ',
     &                 'inconsistent diagonal allocation!'
            If (iCount .gt. lDiag) Then
               Write(6,*) '   - insufficient memory, will return now...'
               irc = 1
               Return
            Else
               Write(6,*) '   - sufficient memory, going to continue...'
            End If
         End If
      End If

C     Reorder CMO. This also removes frozen orbitals.
C     -----------------------------------------------

      l_COcc = nT1AOT(1)
      l_CVir = nAOVir(1)
      Call GetMem('COcc','Allo','Real',ip_COcc,l_COcc)
      Call GetMem('CVir','Allo','Real',ip_CVir,l_CVir)
      Call ChoMP2_MOReOrd(CMO,Work(ip_COcc),Work(ip_CVir))

C     Backtransform.
C     --------------

      Call ChoMP2_BackTra(iTyp,Work(ip_COcc),Work(ip_CVir),
     &                    BaseName_AO,DoDiag,Diag)

C     Open AO vector files (i.e. get units to return).
C     ------------------------------------------------

      Do iSym = 1,nSym
         Write(FullName_AO,'(A3,I1)') BaseName_AO,iSym
         lU_AO(iSym) = 7
         Call daName_MF_WA(lU_AO(iSym),FullName_AO)
      End Do

C     Debug: check backtransformation.
C     --------------------------------

      If (Debug) Then
         Call ChoMP2_CheckBackTra(iTyp,Work(ip_COcc),Work(ip_CVir),
     &                            lU_AO)
      End If

C     Delete MO files if requested.
C     -----------------------------

      If (Delete) Then
         iOpen = 1
         iClose = 3
         Do iSym = 1,nSym
            Call ChoMP2_OpenF(iOpen,iTyp,iSym)
            Call ChoMP2_OpenF(iClose,iTyp,iSym)
         End Do
      End If

C     Deallocate and exit.
C     --------------------

      Call GetMem('CVir','Free','Real',ip_CVir,l_CVir)
      Call GetMem('COcc','Free','Real',ip_COcc,l_COcc)
      End
