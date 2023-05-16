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
! Copyright (C) 1994, Martin Schuetz                                   *
!               2017,2022, Roland Lindh                                *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine TraFck(canorb,FOVMax)
!***********************************************************************
!                                                                      *
!     purpose: Transform Fock Matrix to get Orbital energies           *
!                                                                      *
!     input:                                                           *
!       canorb  : Boolean: TRUE, if canonical orbs are desired,        *
!                          FALSE otherwise                             *
!                                                                      *
!     output:                                                          *
!       FOVMax  : Max Element of occ/virt block in Fock Matrix,        *
!                 transformed into MO basis                            *
!                                                                      *
!***********************************************************************
      use SpinAV, only: Do_SpinAV
      use InfSCF, only: MaxBas, nBO, nBT, nnFr, nSym, FckAuf, nBas, nFro, nConstr, nOcc, nOrb, TimFld
      use SCF_Arrays, only: FockAO, Ovrlp, CMO, EOrb
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
!
      Integer nD
      Real*8 FOVMax
      Logical CanOrb
!
      Real*8, Dimension(:), Allocatable:: FckM, FckS, HlfF, EigV, Ctmp, Scratch, CMOOld, Scrt, COvrlp
      Real*8 Fia, Dummy
      Real*8 Cpu1,Cpu2,Tim1,Tim2,Tim3
      Real*8 Tmp, Tmp0, Tmp1
      Integer ioFckM,iCMO,jEOr,iptr,iptr2,nOrbmF,nOccmF, nVrt,ia, iD, iSym, iDiag, iDum, iErr, iiCMO,     &
              iiEigV, iiScratch, iOff, jjEOr, kk, kOcc, kOff, n2Sort, n2Zero, nFound
      Integer, External:: iDaMax_
      Real*8, External:: DDot_
!
      nD   =Size(FockAO,2)

      Call Timing(Cpu1,Tim1,Tim2,Tim3)
#ifdef _DEBUGPRINT_
      Call NrmClc(FockAO,SIZE(FockAO),'TraFck','FockAO')
#endif
!---- allocate memory for modified Fock matrix
      Call mma_allocate(FckM,nBT,Label='FckM')
!---- allocate memory for squared Fock matrix
      Call mma_allocate(FckS,MaxBas**2,Label='FckS')
!
      FOVMax=Zero
!                                                                      *
!***********************************************************************
!                                                                      *
      Do iD = 1, nD
!                                                                      *
!***********************************************************************
!                                                                      *
!----    modify Fock matrix
         FckM(:)=FockAO(:,iD)
         If (nnFr.gt.0) then
           Call ModFck(FckM,Ovrlp,nBT,CMO(1,iD),nBO,nOcc(1,iD))
         endif
!
         ioFckM=1
         iCMO=1
         jEOr=1
!                                                                      *
!***********************************************************************
!                                                                      *
         Do iSym = 1,nSym
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
            Write (6,*) 'iD=',iD
            Call RecPrt('TraFck: Old CMO',' ',CMO(iCMO,iD),nBas(iSym),nOrb(iSym))
            jCMO=iCMO
#endif

            nOrbmF=nOrb(iSym)   -nFro(iSym)
            nOccmF=nOcc(iSym,iD)-nFro(iSym)
            nVrt  =nOrb(iSym)   -nOcc(iSym,iD)

            iCMO=iCMO+nBas(iSym)*nFro(iSym) ! pointer to CMO

            If (nOrbmF.gt.0) Then

!              allocate memory for half-transformed Fock matrix
               Call mma_allocate(HlfF,nBas(iSym)*nOrbmF,Label='HlfF')
!--------      transform Fock matrix into new MO space
               Call Square(FckM(ioFckM),FckS,1,nBas(iSym),nBas(iSym))
               Call DGEMM_('N','N',nBas(iSym),nOrbmF,nBas(iSym),    &
                           One,FckS,nBas(iSym),                     &
                           CMO(iCMO,iD),nBas(iSym),                 &
                           Zero,HlfF,nBas(iSym))
               Call DGEMM_Tri('T','N',nOrbmF,nOrbmF,nBas(iSym),     &
                              One,CMO(iCMO,iD),nBas(iSym),          &
                              HlfF,nBas(iSym),                      &
                              Zero,FckS,nOrbmF)
#ifdef _DEBUGPRINT_
               Write(6,*) 'transformed Fck in trafck:'
               Call TriPrt(' ',' ',FckS,nOrbmF)
#endif
!              dispose memory of half-transformed Fock matrix
               Call mma_deallocate(HlfF)

               If ((nOccmF.gt.0).AND.(nVrt.gt.0)) Then
                  iptr=1+nOccmF*(nOccmF+1)/2
!
!                 get max Fock Matrix Element in OV block...
!
                  Do ia=1,nVrt
                     Fia=abs(FckS(iptr+IDAMAX_(nOccmF,FckS(iptr),1)-1))
                     FOVMax=Max(Fia,FOVMax)
                     iptr=iptr+nOccmF+ia
                  End Do
               End If
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!--------      eventually diagonalize occ/occ & virt/virt Block
!              separately to form canonical orbitals
!
               If (canorb) Call Mk_CanOrb()

            End If
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!------    Update pointers
           iCMO=iCMO+nOrbmF*nBas(iSym)
           ioFckM=ioFckM+nBas(iSym)*(nBas(iSym)+1)/2
           jEOr=jEOr+nOrbmF
#ifdef _DEBUGPRINT_
           Call RecPrt('TraFck: New CMO',' ',CMO(jCMO,iD),nBas(iSym),nOrb(iSym))
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
         End Do ! iSym
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGCHECK__
!------  Check orthogonality
         If (canorb) Call ChkOrt(iD,Whatever)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
      End Do ! iD
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Deallocate memory
      Call mma_deallocate(FckS)
      Call mma_deallocate(FckM)
!
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld(13) = TimFld(13) + (Cpu2 - Cpu1)

      Return
      Contains

!                                                                      *
!***********************************************************************
!                                                                      *
      Subroutine Mk_CanOrb()
      Integer :: ii, ia, i, j, iOcc, jOcc
!                                                                      *
!***********************************************************************
!                                                                      *
!     allocate space for Eigenvectors of occ/occ & virt/virt block
      Call mma_allocate(EigV,Max(nOccmF,nVrt)**2,Label='EigV')
!     allocate space for temporary CMO (for MatMul)
      Call mma_allocate(CTmp,Max(nOccmF,nVrt)*nBas(iSym),Label='CTmpX')
      If (.NOT.FckAuf) Then
         Call mma_allocate(CMOOld,nOccmF*nBas(iSym),Label='CMOOld')
         Call DCopy_(nOccmF*nBas(iSym),CMO(iCMO,iD),1,CMOOld,1)
         Call mma_allocate(Scrt,nBas(iSym)**2,Label='Scrt')
         Call mma_allocate(COvrlp,nBas(iSym)*nOccmF,Label='COvrlp')
      End If
!
      If (nOccmF.gt.0) Then
!
!--------diagonalize occ/occ first find the proper pointer to EOr
!
         jEOr=jEOr+nFro(iSym)
         Call mma_allocate(Scratch,nOccmF**2,Label='Scratch')
         Dummy=Zero
         iDum=0
         nOccmF=nOccmF-nConstr(iSym)
!
         Call Diag_Driver('V','A','L',nOccmF,FckS,                &
                          Scratch,nOccmF,Dummy,Dummy,iDum,        &
                          iDum,EOrb(jEOr,iD),EigV,nOccmF,         &
     &                    1,0,'J',nFound,iErr)
!
         If (nConstr(iSym).gt.0) Then
            Call FZero(Scratch,(nOccmF+nConstr(iSym))**2)
            Do j=1,nOccmF
               iiEigV=1+nOccmF*(j-1)
               iiScratch=1+(nOccmF+nConstr(iSym))*(j-1)
               call dcopy_(nOccmF,EigV(iiEigV),1,Scratch(iiScratch),1)
            End Do
            Do j=nOccmF+1,nOccmF+nConstr(iSym)
               iiScratch=1+(nOccmF+nConstr(iSym))*(j-1)+j-1
               Scratch(iiScratch)=One
            End Do
            call dcopy_((nOccmF+nConstr(iSym))**2,Scratch,1,EigV,1)
         End If
         nOccmF=nOccmF+nConstr(iSym)
         Call mma_deallocate(Scratch)
         n2zero=nOccmF
         If (Do_SpinAV) n2zero=n2zero+nConstr(iSym)
         Call dCopy_(n2zero*(n2zero+1)/2,[Zero],0,FckS,1)
!
         iDiag = 0
         Do i = 1, n2zero
            iDiag = iDiag + i
            FckS(iDiag) = EOrb(jEOr+i-1,iD)
         End Do
!
!------- Rotate MOs to diagonalize occ/occ block
         Do ii=0,nOccmF-1
            call dcopy_(nBas(iSym),CMO(iCMO+nBas(iSym)*ii,iD),1,CTmp(1+nBas(iSym)*ii),1)
         End Do
         Call DGEMM_('N','N',nBas(iSym),nOccmF,nOccmF,      &
                     One,Ctmp,nBas(iSym),                   &
                         EigV,nOccmF,                       &
                     Zero,CMO(iCMO,iD),nBas(iSym))
!
!        Fix standard phase pf the orbitals
!
         Do i = 1, nOccmF
            call VecPhase(CMO(iCMO+(i-1)*nBas(iSym),iD),nBas(iSym))
         End Do
!
!        Order the occupied orbitals by maximum overlap with
!        the old MOs.
!
         If (.NOT.FckAuf) Then
!
            Call FZero(COvrlp,nOccmF*nBas(iSym))
            Call Square(Ovrlp(ioFckM),Scrt,1,nBas(iSym),nBas(iSym))
            Call DGEMM_('T','N', nOccmF,nBas(iSym),nBas(iSym),     &
                        One,CMO(iCMO,iD),nBas(iSym),               &
                            Scrt,nBas(iSym),                       &
                        Zero,COvrlp,nOccmF)
!
            Do iOcc = 1, nOccmF-1  !  Loop over the new MOs
!
               iOff = (iOcc-1)*nBas(iSym) + iCMO
               kOcc=0
               Tmp0=Zero
               Do jOcc = 1, nOccmF !  Loop over the new MOs
                  Tmp1=Abs(DDot_(nBas(iSym),COvrlp(jOcc),nOccmF,CMO(iOff,iD),1))
                  If (Tmp1.gt.Tmp0) Then
                     Tmp0=Tmp1
                     kOcc=jOcc
                  End If
               End Do
!
               If (iOcc.ne.kOcc) Then
                  ii = iOcc + jEOr - 1
                  kk = kOcc + jEOr - 1
                  tmp = EOrb(ii,iD)
                  EOrb(ii,iD) = EOrb(kk,iD)
                  EOrb(kk,iD) = tmp
                  kOff = (kOcc-1)*nBas(iSym) + iCMO
                  Call DSwap_(nBas(iSym),CMO(iOff,iD),1,CMO(kOff,iD),1)
               End If
            End Do
!
            iDiag = 0
            Do i = 1, n2zero
               iDiag = iDiag + i
               FckS(iDiag) = EOrb(jEOr+i-1,iD)
            End Do
!
         End If
      End If
!
      If (nVrt.gt.0) Then
!
!--------now diagonalize virt/virt block
!        setup virt/virt block in triangular Fock Matrix
!
         iptr=1+nOccmF*(nOccmF+3)/2
         iptr2=1
         Do ia=1,nVrt
            call dcopy_(ia,FckS(iptr),1,FckS(iptr2),1)
            iptr=iptr+nOccmF+ia
            iptr2=iptr2+ia
         End Do
         Call mma_allocate(Scratch,nVrt**2,Label='Scratch')
         If (Do_SpinAV) Then
            nVrt=nVrt-nConstr(iSym)
            jEOr = jEOr + nConstr(iSym)
            iptr=1+nOccmF*(nOccmF+3)/2
            Do ia=1,nConstr(iSym)
               iptr=iptr+nOccmF+ia
               FckS(iptr)=-0.666d3*dble(1000-ia)
            End Do
         End If
         Dummy=Zero
         iDum=0
         Call Diag_Driver('V','A','L',nVrt,FckS,             &
                          Scratch,nVrt,Dummy,Dummy,iDum,     &
                          iDum,EOrb(jEOr+nOccmF,iD),EigV,    &
                          nVrt,1,0,'J',nFound,iErr)
         If (Do_SpinAV) Then
            nVrt=nVrt+nConstr(iSym)
            jEOr = jEOr - nConstr(iSym)
            Call FZero(Scratch,nVrt**2)
            Do j=1,nVrt-nConstr(iSym)
               iiEigV=1+(nVrt-nConstr(iSym))*(j-1)
               iiScratch=1+nVrt*(nConstr(iSym)+j-1)
               call dcopy_(nVrt-nConstr(iSym),EigV(iiEigV),1,Scratch(iiScratch),1)
            End Do
            Do j=1,nConstr(iSym)
               iiScratch=1+nVrt*(j-1)+j-1
               Scratch(iiScratch)=One
            End Do
            call dcopy_(nVrt**2,Scratch,1,EigV,1)
         End If
         Call mma_deallocate(Scratch)
!--------rotate MOs to diagonalize virt/virt block
         iptr=iCMO+nOccmF*nBas(iSym)
         Do ia=0,nVrt-1
            call dcopy_(nBas(iSym),CMO(iptr+nBas(iSym)*ia,iD),1,CTmp(1+nBas(iSym)*ia),1)
         End Do
         Call DGEMM_('N','N',nBas(iSym),nVrt,nVrt,       &
                     One, Ctmp,nBas(iSym),               &
                          EigV,nVrt,                     &
                     Zero,CMO(iptr,iD),nBas(iSym))
      End If
!
      If (nConstr(iSym).gt.0) Then
!--------Sort non-wavelet eigenvalues/eigenvectors
         n2sort=nOccmF-nConstr(iSym)
         If (FckAuf) Call SortEig(EOrb(jEOr,iD),CMO(iCMO,iD),n2sort,nBas(iSym),1,.true.)
         jjEOr=jEOr+nOccmF+nConstr(iSym)
         iiCMO=iCMO+nBas(iSym)*(nOccmF+nConstr(iSym))
         n2sort=nVrt-nConstr(iSym)
         If (FckAuf) Call SortEig(EOrb(jjEOr,iD),CMO(iiCMO,iD),n2sort,nBas(iSym),1,.true.)
      Else
!--------Sort all eigenvalues and eigenvectors
         If (FckAuf) Call SortEig(EOrb(jEOr,iD),CMO(iCMO,iD),nOrbmF,nBas(iSym),1,.true.)
      End If
!     dispose memory
      If (.NOT.FckAuf) Then
         Call mma_deallocate(COvrlp)
         Call mma_deallocate(Scrt)
         Call mma_deallocate(CMOOld)
      End If
      Call mma_deallocate(CTmp)
      Call mma_deallocate(EigV)

      End Subroutine Mk_CanOrb

      End Subroutine TraFck
