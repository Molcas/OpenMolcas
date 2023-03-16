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
! Copyright (C) 1992, Per-Olof Widmark                                 *
!               1992, Markus P. Fuelscher                              *
!               1992, Piotr Borowski                                   *
!               2003, Valera Veryazov                                  *
!               1998,2022, Roland Lindh                                *
!***********************************************************************
      SubRoutine SOrb(LuOrb,SIntTh,iTerm)
!***********************************************************************
!                                                                      *
!     purpose: Get starting orbitals from:                             *
!              -1) default choice                                      *
!               0) diagonalizaton of the core                          *
!               1) via intermediate calculation of HF AOs              *
!               2) input orbitals                                      *
!               3) input density matrix                                *
!                                                                      *
!***********************************************************************
#ifdef _HDF5_
      Use mh5, Only: mh5_close_file
      use InfSCF, only: IsHDF5, FileOrb_ID
#endif
      use InfSCF, only: DoCholesky, InVec, nBB, KSDFT, iUHF, One_Grid,
     &                  SCF_FileOrb, StVec, Scrmbl,
     &                  nBas, nOrb, nSym, ScrFac, nBB, nBT, nnB
      Use SCF_Arrays, only: CMO, TrM, FockAO, Ovrlp, EOrb, OccNo,
     &                      OneHam
      use Files
      Implicit None
!
      Real*8 SIntTh
      Integer iTerm, LuOrb, nD
      Integer IsUHF, iD, nData
      Character FName*512, KSDFT_save*80
      Logical FstItr
      Logical found
!
      nD = iUHF + 1
      CALL DecideonCholesky(DoCholesky)
!-------- Cholesky and NDDO are incompatible
      IF (DoCholesky.and.InVec.eq.1) THEN
         call WarningMessage(1,
     &    ' In SORB: Cholesky and NDDO not implemented !!!; '//
     &    ' NDDO option ignored')
         InVec=-1
      ENDIF
!
!---- Is default clause choosen?
!
      If(InVec.eq.-1) Then
         Call qpg_darray('SCF orbitals',found,ndata)
         If(found .and. nData.eq.nBB) Then
            Call qpg_darray('OrbE',found,ndata)
            If(found) Then
               InVec=8
            End If
         End If
      End If
      If(InVec.eq.-1) Then
         Call qpg_darray('Guessorb',found,ndata)
         If(found .and. nData.eq.nBB) Then
            Call qpg_darray('Guessorb energies',found,ndata)
            If(found) Then
               InVec=9
            End If
         End If
      End If
      If(InVec.eq.-1) Then
         InVec=0
      End If
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Has the user selected a method?
!
100   Continue
!                                                                      *
!***********************************************************************
!                                                                      *
      Select Case (InVec)
!                                                                      *
!***********************************************************************
!                                                                      *
      Case (0)

!-------- Diagonalize core
          Call Start0()
!                                                                      *
!***********************************************************************
!                                                                      *
      Case (1)
!
!-------- HF AO orbitals as intermediate step...
!
!------- NDDO, always none-DFT
!
         Call SwiOpt(.False.,OneHam,Ovrlp,nBT,CMO,nBB,nD)
         Call Start0()
         InVec=0
         Call SOrbCHk(OneHam,FockAO,nBT,nD)
         KSDFT_save=KSDFT
         KSDFT='SCF'
         Call WrInp_SCF(SIntTh)
         FstItr=.True.
         Call WfCtl_SCF(iTerm,'NDDO      ',FstItr,SIntTh)
         KSDFT=KSDFT_save
         Call Free_TLists
         If (iTerm.ne.0) Call Quit(iTerm)
         Write(6,*)
         Write(6,'(A)') 'Generation of NDDO vectors completed!'
         Write(6,*)
         Write(6,*) '2nd step: optimizing HF MOs...'
         Write(6,*) '------------------------------'
         Call SwiOpt(.TRUE.,OneHam,Ovrlp,nBT,CMO,nBB,nD)
!------- Reset to to start from the current MO set
         Call Init_SCF()
         InVec=5
!IFG: I presume the arguments after LuOut in these two calls are correct,
!     they were missing!
         If(iUHF.eq.0) Then
            FName='SCFORB'
            Call Start2(FName,LuOut,CMO,nBB,nD,Ovrlp,nBT,
     &               EOrb,OccNo,nnB)
         Else
            FName='UHFORB'
            Call Start2(FName,LuOut,CMO,nBB,nD,Ovrlp,nBT,
     &               EOrb,OccNo,nnB)
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
      Case (2)
!                                                                      *
!-------- Read INPORB
         One_Grid=.True.
         FName=SCF_FileOrb
         Call Start2(FName,LuOrb,CMO,nBB,nD,Ovrlp,nBT,
     &               EOrb,OccNo,nnB)
!                                                                      *
!***********************************************************************
!                                                                      *
      Case (3)

!-------- Read COMOLD
         One_Grid=.True.
         Call Start3(CMO,TrM,nBB,nD,OneHam,Ovrlp,nBT)
!                                                                      *
!***********************************************************************
!                                                                      *
      Case (6)

         write(6,*)
         write(6,*) '     Constrained SCF calculation '
         write(6,*)
         StVec='Constrained orbitals'
         One_Grid=.True.
         FName=SCF_FileOrb
         Call Chk_Vec_UHF(FName,LuOrb,isUHF)
         If (IsUHF.eq.1) Then
            InVec=2
            Go To 100
         EndIf
         Call Start6(FName,LuOrb,CMO,nBB,nD,EOrb,OccNo,nnB)
!                                                                      *
!***********************************************************************
!                                                                      *
      Case (8)

         StVec='Detected old SCF orbitals'
         One_Grid=.True.
         Call start0y(CMO,nBB,nD,EOrb,nnB)
!                                                                      *
!***********************************************************************
!                                                                      *
      Case (9)

         StVec='Detected guessorb starting orbitals'
!        One_Grid=.True.
         Call start0x(CMO,nBB,nD,EOrb,nnB)
!                                                                      *
!***********************************************************************
!                                                                      *
      Case Default

         Write (6,*) 'Illegal inVec value:',InVec
         Call Abend()
!                                                                      *
!***********************************************************************
!                                                                      *
      End Select
!                                                                      *
!***********************************************************************
!                                                                      *
      If (Scrmbl) Then
         Do iD = 1, nD
            Call Scram(CMO(1,iD),nSym,nBas,nOrb,ScrFac)
         End Do
      End If

      Call SOrbCHk(OneHam,FockAO,nBT,nD)
#ifdef _HDF5_
      If (isHDF5) Call mh5_close_file(fileorb_id)
#endif

      End subroutine SOrb
