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
!               2022, Roland Lindh                                     *
!***********************************************************************
      SubRoutine MODens()
!***********************************************************************
!                                                                      *
!     purpose: Compute density matrix in molecular orbital basis       *
!                                                                      *
!***********************************************************************
      use InfSCF, only: MaxBas, MaxOrb, DMOMax, nSym, TEEE, nDens, TimFld, MaxBXO, nBas, nOcc, nOrb, nD
#ifdef _DEBUGPRINT_
      Use InfSCF, only: nBO, nBT
#endif
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: Zero, One
      Use SCF_Arrays, only: Dens, CMO, Ovrlp
      Implicit None
!
      Integer iD, jD, iT, iOvl, iSym, iiBO, iiBT, i, j
      Real*8, Dimension(:), Allocatable:: DnsS, OvlS, DMoO, Aux1, Aux2
      Real*8 CPU1, CPU2, Tim1, Tim2, Tim3
!
!----------------------------------------------------------------------*
!     Start                                                            *
!----------------------------------------------------------------------*
!
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
!
!---- Allocate memory for squared density matrix
      Call mma_allocate(DnsS,MaxBas**2,Label='DnsS')
!
!---- Allocate memory for squared overlap matrix
      Call mma_allocate(OvlS,MaxBas**2,Label='OvlS')
!
!---- Allocate memory for density matrix in MO basis
      Call mma_allocate(DMoO,MaxOrb*(MaxOrb+1)/2,Label='DMoO')
!
!---- Allocate memory for auxiliary matrix
      Call mma_allocate(Aux1,MaxBxO,Label='Aux1')
!
!---- Allocate memory for auxiliary matrix
      Call mma_allocate(Aux2,MaxBxO,Label='Aux2')
!
!
      DMOMax = Zero
      Do jD = 1, nD
!
#ifdef _DEBUGPRINT_
         Call NrmClc(Dens(1,jD,nDens),nBT,'MoDens','D in AO   ')
         Call NrmClc(Ovrlp         ,nBT,'MoDens','Overlap   ')
         Call NrmClc(CMO(1,jD)    ,nBO,'MoDens','CMOs      ')
         Write (6,*) 'nOcc=',(nOcc(i,jD),i=1,nSym)
!        Write (6,'(F16.8)') DXot(MaxBxO,CMO(1,jD),1,CMO(1,jD),1)
#endif
         it   = 1
         id   = 1
         iOvl = 1
         Do iSym = 1, nSym
!
            iiBO = nBas(iSym)*nOrb(iSym)
            iiBT = nBas(iSym)*(nBas(iSym) + 1)/2
!
            If (nOcc(iSym,jD)>0 .or. (Teee.and.nBas(iSym)>0)) Then
               Call DSq(Dens(id,jD,nDens),DnsS,1,nBas(iSym),nBas(iSym))
               Call Square(Ovrlp(iOvl),OvlS,1,nBas(iSym),nBas(iSym))
!
               Call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nBas(iSym),      &
                           One,OvlS,nBas(iSym),                           &
                                 CMO(it,jD),nBas(iSym),                   &
                           Zero,Aux1,nBas(iSym))
               Call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nBas(iSym),      &
                           One,DnsS,nBas(iSym),                           &
                                 Aux1,nBas(iSym),                         &
                           Zero,Aux2,nBas(iSym))
               Call DGEMM_('N','N',nBas(iSym),nOrb(iSym),nBas(iSym),      &
                           One,OvlS,nBas(iSym),                           &
                                 Aux2,nBas(iSym),                         &
                           Zero,Aux1,nBas(iSym))
               Call DGEMM_Tri('T','N',nOrb(iSym),nOrb(iSym),nBas(iSym),   &
                              One,CMO(it,jD),nBas(iSym),                  &
                                  Aux1,nBas(iSym),                        &
                              Zero,DMoO,nOrb(iSym))
!
!              Call TriPrt('D(mo)','(8F12.6)',DMoO,nOrb(iSym))
               Do i = nOcc(iSym,jD)+1 , nOrb(iSym)
                  Do j = 1, nOcc(iSym,jD)
                     DMOMax = Max(DMOMax,DBLE(nD)*Abs(DMoO(i*(i-1)/2+j)))
                  End Do
               End Do
            End If
!
            it   = it   + iiBO
            id   = id   + iiBT
            iOvl = iOvl + iiBT
         End Do
!
      End Do
!
!---- Deallocate memory
      Call mma_deallocate(Aux2)
      Call mma_deallocate(Aux1)
      Call mma_deallocate(DMoO)
      Call mma_deallocate(OvlS)
      Call mma_deallocate(DnsS)
!
#ifdef _DEBUGPRINT_
      Write(6,*)' DMOMax in MODens',DMOMax
#endif
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld(14) = TimFld(14) + (Cpu2 - Cpu1)
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
!
      Return
      End SubRoutine MODens
