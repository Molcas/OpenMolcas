!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      SubRoutine Cho_X_Init_Par(irc,isDF)
!
!     Purpose: setup for parallel Cholesky/DF.
!
      Implicit None
      Integer irc
      Logical isDF

      If (isDF) Then
         Call Cho_X_Init_Par_DF(irc)
      Else
         Call Cho_X_Init_Par_GenBak()
         Call Cho_X_Init_Par_Cho(irc)
      End If

      End
      SubRoutine Cho_X_Init_Par_DF(irc)
!
!     Purpose: setup for parallel DF.
!
#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: MyRank, nProcs, Is_Real_Par
#endif
      Implicit None
      Integer irc

      Character(LEN=17), Parameter:: SecNam = 'Cho_X_Init_Par_DF'

#if defined (_DEBUGPRINT_)
      Logical, Parameter:: LocDbg = .True.
#else
      Logical, Parameter:: LocDbg = .False.
#endif

#if defined (_MOLCAS_MPP_)
#include "cholesky.fh"
      Integer nV(8)
      Integer iSym
      Logical isSerial

      irc = 0

!     Return if serial.
!     -----------------

      isSerial = nProcs.eq.1 .or. .not.Is_Real_Par()
      If (isSerial) Then
         If (LocDbg) Then
            Write(6,*) SecNam,': serial run, nothing to do...'
            Write(6,*) '#nodes: ',nProcs,'  myRank: ',myRank
         End If
         Return
      Else
         If (LocDbg) Then
            Write(6,*) SecNam,': parallel run...'
            Write(6,*) '#nodes: ',nProcs,'  myRank: ',myRank
         End If
      End If

!     Reset number of vectors to the number on this node as stored on
!     the runfile.
!     ---------------------------------------------------------------

      Call iCopy(nSym,NumCho,1,nV,1)
      Call Get_iArray('nVec_RI',NumCho,nSym)
      NumChT = NumCho(1)
      Do iSym = 2,nSym
         NumChT = NumChT + NumCho(iSym)
      End Do

!     Debug print.
!     ------------

      If (LocDbg) Then
         Write(6,*)
         Write(6,*) 'Output from ',SecNam,':'
         Write(6,*) 'NumCho before: ',(nV(iSym),iSym=1,nSym)
         Write(6,*) 'NumCho after : ',(NumCho(iSym),iSym=1,nSym)
      End If

#else

      irc = 0
      If (LocDbg) Then
         Write(6,*) SecNam,': serial run, nothing to do...'
      End If

#endif

      End
      SubRoutine Cho_X_Init_Par_Cho(irc)
!
!     Purpose: setup for parallel Cholesky.
!
#if defined (_MOLCAS_MPP_)
      Use Para_Info, Only: MyRank, nProcs, Is_Real_Par
      use ChoSwp, only: InfVec
      use stdalloc
#endif
      Implicit None
      Integer irc

      Character(LEN=18), Parameter:: SecNam = 'Cho_X_Init_Par_Cho'

#if defined (_DEBUGPRINT_)
      Logical, Parameter:: LocDbg = .True.
#else
      Logical, Parameter:: LocDbg = .False.
#endif

#if defined (_MOLCAS_MPP_)
#include "cholesky.fh"

      Integer nV(8)
      Integer iSym, i, j
      Logical isSerial


      Integer, Allocatable:: IDV(:), myInfV(:)

      irc = 0

!     Return if serial.
!     -----------------

      isSerial = nProcs.eq.1 .or. .not.Is_Real_Par()
      If (isSerial) Then
         If (LocDbg) Then
            Write(6,*) SecNam,': serial run, nothing to do...'
            Write(6,*) '#nodes: ',nProcs,'  myRank: ',myRank
         End If
         Return
      Else
         If (LocDbg) Then
            Write(6,*) SecNam,': parallel run...'
            Write(6,*) '#nodes: ',nProcs,'  myRank: ',myRank
         End If
      End If

!     Reset vector info to fit vectors stored on this node.
!     -----------------------------------------------------

      Do iSym = 1,nSym
         nV(iSym) = 0
         If (NumCho(iSym) .gt. 0) Then
            Call mma_allocate(IDV,NumCho(iSym),Label='IDV')
            Call Cho_Distrib_Vec(1,NumCho(iSym),IDV,nV(iSym))
            If (nV(iSym) .gt. 0) Then
               Call mma_allocate(myInfV,nV(iSym),Label='myInfV')
               Do j = 1,SIZE(InfVec,2)
                  If (j .ne. 3) Then
                     Do i = 1,nV(iSym)
                        myInfV(i) = InfVec(IDV(i),j,iSym)
                     End Do
                     InfVec(:,j,iSym) = myInfV(:)
                  End If
               End Do
               Call mma_deallocate(myInfV)
            End If
            Call mma_deallocate(IDV)
         End If
      End Do

!     Reset number of vectors.
!     ------------------------

      Call iSwap(nSym,NumCho,1,nV,1)
      NumChT = NumCho(1)
      Do iSym = 2,nSym
         NumChT = NumChT + NumCho(iSym)
      End Do

!     Debug print.
!     ------------

      If (LocDbg) Then
         Write(6,*)
         Write(6,*) 'Output from ',SecNam,':'
         Write(6,*) 'NumCho before: ',(nV(iSym),iSym=1,nSym)
         Write(6,*) 'NumCho after : ',(NumCho(iSym),iSym=1,nSym)
      End If

#else

      irc = 0
      If (LocDbg) Then
         Write(6,*) SecNam,': serial run, nothing to do...'
      End If

#endif

      End
      SubRoutine Cho_X_Init_Par_GenBak()
      Use Para_Info, Only: Is_Real_Par
      use ChoSwp, only: InfVec, InfVec_Bak
      use ChPari
      use stdalloc, only: mma_allocate
      Implicit None
#include "cholesky.fh"

      NumCho_Bak(:)=0
      If (Is_Real_Par()) Then
         Call mma_allocate(InfVec_Bak,SIZE(InfVec,1),SIZE(InfVec,2),    &
     &                     SIZE(InfVec,3),Label='InfVec_Bak')
         InfVec_Bak(:,:,:)=InfVec(:,:,:)
         NumCho_Bak(1:nSym)=NumCho(1:nSym)
      End If

      End
