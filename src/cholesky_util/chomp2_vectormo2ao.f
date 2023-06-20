!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2007,2008, Thomas Bondo Pedersen                       *
!***********************************************************************
      SubRoutine ChoMP2_VectorMO2AO(iTyp,Delete,BaseName_AO,CMO,DoDiag,
     &                              Diag,lDiag,lU_AO,irc)
!
!     Thomas Bondo Pedersen, Dec. 2007 - Jan. 2008.
!
!     Purpose: backtransform vectors from MO to AO basis.
!
!     Input:
!        iTyp......... specifies type of vectors (iTyp is used to open
!                      MO vector files through ChoMP2_OpenF()).
!        Delete....... Flag specifiyng whether MO vector files are to
!                      be deleted before exiting this routine.
!        BaseName_AO.. Base name (3 characters) for the files containing
!                      AO vectors. Symmetry index will be appended!
!        CMO.......... MO coefficient array
!        DoDiag....... if .True., calculate AO diagonal elements as
!                      D(ab) = sum_J L(J,ab)*L(J,ab)
!     Output:
!        Diag......... Contains the diagonal if requested (flag DoDiag).
!        lDiag........ Dimension of Diag
!        lU_AO........ Array containing units of open AO vector files.
!        irc.......... return code.
!                      = 0 if successful, non-zero otherwise
!                      (must be checked by caller).
!
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!*** CHOLESKY INFORMATION MUST BE INITIALIZED WHEN CALLING THIS ROUTINE.
!    --> I.e. Cho_X_Init() must have been called.
!    --> It is actually sufficient that choorb.fh + nSym in
!        cholesky.fh are available (which they are when Cho_X_Init()
!        has been called).
!*** CHOLESKY MP2 INFORMATION MUST BE INITIALIZED WHEN CALLING THIS
!    ROUTINE.
!    --> I.e. ChoMP2_Setup() must have been called.
!*** AO VECTORS ARE STORED IN LOWER TRIANGULAR [M(I,J), I.GE.J] FORMAT.
!    --> I.e. vectors are stored as L(J,ab) where a>=b
!*** DIAGONAL IS STORED IN LOWER TRIANGULAR FORMAT.
!    --> I.e. diagonal is stored as D(ab) where a>=b (same as vectors).
!***********************************************************************
!!!!!!!!!!!!!!!!!!!!!!!!!! IMPORTANT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************************
!
      use stdalloc
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

      Character(LEN=18), Parameter:: SecNam = 'ChoMP2_VectorMO2AO'
      Character(LEN=11), Parameter:: ThisNm = 'VectorMO2AO'

      Character(LEN=4) FullName_AO

      Integer iSym, iSyma, iSymb, iCount, iOpen, iClose

#if defined (_DEBUGPRINT_)
      Logical, Parameter:: Debug = .True.
#else
      Logical, Parameter:: Debug = .False.
#endif

      Real*8, Allocatable:: COcc(:), CVir(:)

      Integer MulD2h, k, l
      MulD2h(k,l)=iEOr(k-1,l-1)+1

!     Initializations.
!     ----------------

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

!     Reorder CMO. This also removes frozen orbitals.
!     -----------------------------------------------

      Call mma_allocate(COcc,nT1AOT(1),Label='COcc')
      Call mma_allocate(CVir,nAOVir(1),Label='CVir')
      Call ChoMP2_MOReOrd(CMO,COcc,CVir)

!     Backtransform.
!     --------------

      Call ChoMP2_BackTra(iTyp,COcc,CVir,BaseName_AO,DoDiag,Diag)

!     Open AO vector files (i.e. get units to return).
!     ------------------------------------------------

      Do iSym = 1,nSym
         Write(FullName_AO,'(A3,I1)') BaseName_AO,iSym
         lU_AO(iSym) = 7
         Call daName_MF_WA(lU_AO(iSym),FullName_AO)
      End Do

!     Debug: check backtransformation.
!     --------------------------------

      If (Debug) Then
         Call ChoMP2_CheckBackTra(iTyp,COcc,CVir,lU_AO)
      End If

!     Delete MO files if requested.
!     -----------------------------

      If (Delete) Then
         iOpen = 1
         iClose = 3
         Do iSym = 1,nSym
            Call ChoMP2_OpenF(iOpen,iTyp,iSym)
            Call ChoMP2_OpenF(iClose,iTyp,iSym)
         End Do
      End If

!     Deallocate and exit.
!     --------------------

      Call mma_deallocate(CVir)
      Call mma_deallocate(COcc)
      End
