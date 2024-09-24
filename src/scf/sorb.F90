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
      use InfSCF, only: DoCholesky, InVec, nBB, KSDFT, nD, One_Grid, SCF_FileOrb, StVec, &
                        Scrmbl, nBas, nOrb, nSym, ScrFac, nBB, nBT, nnB
      Use SCF_Arrays, only: CMO, TrM, FockAO, Ovrlp, EOrb, OccNo, OneHam
      use Files, only: LuOut

      Implicit None
!
      Real*8 SIntTh
      Integer iTerm, LuOrb, IsUHF, iD, nData
      Character FName*512, KSDFT_save*80
      Logical FstItr
      Logical found
!
      CALL DecideonCholesky(DoCholesky)
!-------- Cholesky and NDDO are incompatible
      IF (DoCholesky.and.InVec.eq.1) THEN
         call WarningMessage(1,' In SORB: Cholesky and NDDO not implemented !!!;  NDDO option ignored')
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
         If(nD==1) Then
            FName='SCFORB'
            Call Start2(FName,LuOut,CMO,nBB,nD,Ovrlp,nBT,EOrb,OccNo,nnB)
         Else
            FName='UHFORB'
            Call Start2(FName,LuOut,CMO,nBB,nD,Ovrlp,nBT,EOrb,OccNo,nnB)
         End If
!                                                                      *
!***********************************************************************
!                                                                      *
      Case (2)
!                                                                      *
!-------- Read INPORB
         One_Grid=.True.
         FName=SCF_FileOrb
         Call Start2(FName,LuOrb,CMO,nBB,nD,Ovrlp,nBT,EOrb,OccNo,nnB)
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

      Contains

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
!***********************************************************************
      SubRoutine Start0()
!***********************************************************************
!                                                                      *
!     purpose: Get starting orbitals from diagonalization of the core. *
!              The result is stored in CMO. Those orbitals are         *
!              then optimized in SubRoutine WfCtl. A set of orthogonal,*
!              symmetry adapted AOs is stored in TrM.                  *
!                                                                      *
!***********************************************************************
      use InfSCF, only: nBB, nBO, nBT, nOcc, nnB, nD
      use SCF_Arrays, only: CMO, TrM, OneHam, Ovrlp, EOrb
      Implicit None

      Integer iD
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
!---- Form transformation matrix
      Call TrGen(TrM(1,1),nBB,Ovrlp,OneHam,nBT)
      if(nD.eq.2) Call dCopy_(nBB,TrM(1,1),1,TrM(1,2),1)
!
!---- Diagonalize core
      Do iD = 1, nD
         Call DCore(OneHam,nBT,CMO(1,iD),TrM(1,iD),nBO, EOrb(1,iD),nnB,nOcc(1,iD),Ovrlp)
      End Do
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
!
      Return
      End SubRoutine Start0

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
!***********************************************************************
      SubRoutine Start0x(CMO,mBB,nD,EOr,mmB)
!***********************************************************************
!                                                                      *
! This routine reads start orbitals generated by guessorb.             *
!                                                                      *
!***********************************************************************
      use InfSCF, only: nSym, nDel, nOrb, nBas
      Implicit None
!
      Integer mBB, nD, mmB
      Real*8 CMO(mBB,nD), EOr(mmB,nD)

      Logical Found
      Integer nData, iD, iRC, iSym, nSum
!
!----------------------------------------------------------------------*
! Get start orbitals.                                                  *
!----------------------------------------------------------------------*
!
      iRC=0
      Call qpg_darray('Guessorb',Found,ndata)
      If (Found) Then
         If (nData.ne.mBB) Then
            Write (6,*) 'Start0x: nData.ne.mBB'
            Write (6,*) '         nData=',nData
            Write (6,*) '         mBB  =',mBB
            Call Abend()
         End If
         Call get_darray('Guessorb',CMO(1,1),ndata)
      Else
         iRC=-1
      End If
      If (iRC.ne.0) Then
         Write (6,*) 'Start0x: no orbitals found!'
         Call Abend()
      End If
!
      iRC=0
      Call qpg_darray('Guessorb energies',found,ndata)
      If (Found) Then
         If (nData.ne.mmB) Then
            Write (6,*) 'Start0x: nData.ne.mmB'
            Write (6,*) '         nData=',nData
            Write (6,*) '         mmB  =',mmB
            Call Abend()
         End If
         Call get_darray('Guessorb energies',EOr(1,1),ndata)
      Else
         iRC=-1
      End If
      If (iRC.ne.0) Then
         Write (6,*) 'Start0x: no energies found!'
         Call Abend()
      End If
!
      If (nD.eq.2) Then
         Call dCopy_(mBB,CMO(1,1),1,CMO(1,2),1)
         Call dCopy_(mmB,EOr(1,1),1,EOr(1,2),1)
      End If
!
      Call qpg_iarray('nDel_go',Found,ndata)
      nSum=0
      If (Found) Then
         Call Get_iArray('nDel_go',nDel,ndata)
         Call Put_iArray('nDel',nDel,ndata)
         Do iSym=1,nSym
            nSum=nSum+nDel(iSym)
         End Do
      End If
!
      If (nSum.gt.0) Then
         Do iSym=1,nSym
            nOrb(iSym)=nBas(iSym)-nDel(iSym)
         End Do
         Do iD = 1, nD
            Call TrimCMO(CMO(1,iD),CMO(1,iD),nSym,nBas,nOrb)
            Call TrimEor(EOr(1,iD),EOr(1,iD),nSym,nBas,nOrb)
         End Do
      End If
!
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
      Return
      End SubRoutine Start0x

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
!***********************************************************************
      SubRoutine Start0y(CMO,mBB,nD,EOr,mmB)
!***********************************************************************
!                                                                      *
! This routine reads old SCf orbitals as start orbitals.               *
!                                                                      *
!***********************************************************************
      use InfSCF, only: nSym, nDel, nOrb, nBas
      Implicit None
!
      Integer mBB, nD, mmB
      Real*8 CMO(mBB,nD), EOr(mmB,nD)

      Logical found
      Integer nData, iD, iSym, nSum
!
!----------------------------------------------------------------------*
! Get start orbitals.                                                  *
!----------------------------------------------------------------------*
!
      Call qpg_darray('SCF orbitals',found,ndata)
      If (Found) Then
         Call get_darray('SCF orbitals',CMO(1,1),ndata)
      End If
      Call qpg_darray('OrbE',found,ndata)
      If (Found) Then
         Call get_darray('OrbE',EOr(1,1),ndata)
      End If
!
      If (nD.eq.2) Then
!
         Call DCopy_(mBB,CMO(1,1),1,CMO(1,2),1)
         Call DCopy_(mmB,EOr(1,1),1,EOr(1,2),1)
!
         Call qpg_darray('SCF orbitals_ab',found,ndata)
         If(Found) Then
            Call get_darray('SCF orbitals_ab',CMO(1,2),ndata)
         End If
         Call qpg_darray('OrbE_ab',found,ndata)
         If(found) Then
            Call get_darray('OrbE_ab',EOr(1,2),ndata)
         End If
      End If
!
      Call qpg_iarray('nDel',Found,ndata)
      nSum=0
      If(Found) Then
         Call Get_iArray('nDel',nDel,ndata)
         Do iSym=1,nSym
            nSum=nSum+nDel(iSym)
         End Do
      End If
      If (nSum.gt.0) Then
         Do iSym=1,nSym
            nOrb(iSym)=nBas(iSym)-nDel(iSym)
         End Do
         Do iD = 1, nD
            Call TrimCMO(CMO(1,iD),CMO(1,iD),nSym,nBas,nOrb)
            Call TrimEor(Eor(1,iD),Eor(1,iD),nSym,nBas,nOrb)
         End Do
      End If
!----------------------------------------------------------------------*
! Done                                                                 *
!----------------------------------------------------------------------*
      Return
      End SubRoutine Start0y





      SubRoutine SOrbChk(OneHam,Fock,mBT,nD)
      Implicit None
!
      Integer mBT, nD
      Real*8 OneHam(mBT), Fock(mBT,nD)

      Integer iD
      Real*8 Whatever
!
      Do iD = 1, nD
!----    Check orthonormality of start orbitals
         Call ChkOrt(iD,Whatever)
!
!----    Form the first Fock matrix
         Call DCopy_(mBT,OneHam,1,Fock(1,iD),1)
      End Do
!
      Return
      End SubRoutine SOrbChk
      End subroutine SOrb
