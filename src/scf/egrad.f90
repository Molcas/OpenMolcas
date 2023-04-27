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
!               2017,2022, Roland Lindh                                *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine EGrad(O,S,nOTSD,C,nC,G,nG,nD,iOpt)
!***********************************************************************
!                                                                      *
!     purpose: This routine calculates the gradient of the SCF energy  *
!              with respect to the rotation parameters.                *
!              Grad(E) = C(t)(FDS-SDF)C                                *
!                                                                      *
!                                                                      *
!     input:                                                           *
!       O       : one-electron hamiltonian of length nOTSD             *
!       S       : overlap in AO basis of length nOTSD                  *
!       C       : matrix transforming to the set of orthonormal        *
!                 (and spherical, if needed) functions of length nC    *
!                                                                      *
!     output:                                                          *
!       G       : gradient of the SCF energy with respect to the       *
!                 rotation parameters of length nG                     *
!***********************************************************************
      Use Orb_Type, only: OrbType
      use InfSCF, only: MaxBas, nBO, nBT, nnFr, nSym, nBas, nOrb, nFro, nOcc, MapDns, iDisk
      use SCF_Arrays, only: Dens, TwoHam, Vxc
      use Constants, only: Zero, One, Two
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nOTSD, nD, nC, nG, iOpt
      Real*8 O(nOTSD),S(nOTSD),C(nC,nD), G(nG,nD)
!
      Integer i, j, k, l, iD, ig, ih, ij, iOff, it, nOr, nOrbmF, nBs, iSym, jDT
      Real*8, Dimension(:,:), Allocatable:: FckM
      Real*8, Dimension(:), Allocatable:: Aux1, Aux2, Aux3
      Real*8, Dimension(:,:), Allocatable, Target::AuxD, AuxT, AuxV
      Real*8, Dimension(:,:), Pointer:: D, T, V
!
!----------------------------------------------------------------------*
!     Start
#ifdef _DEBUGPRINT_
      Write (6,*) 'EGrad: input arrays'
      Write (6,*) '==================================================='
      Call NrmClc(O,nOTSD   ,'EGrad','O')
      Call NrmClc(S,nOTSD   ,'EGrad','S')
      Call NrmClc(D,nOTSD*nD,'EGrad','D')
      Call NrmClc(T,nOTSD*nD,'EGrad','T')
      Call NrmClc(V,nOTSD*nD,'EGrad','V')
      Call NrmClc(C,nC   *nD,'EGrad','C')
!     Do iD = 1, nD
!        Write (*,*) 'OrbType(:,iD)', OrbType(:,iD)
!     End Do ! iD
      Write (6,*) '==================================================='
      Write (6,*)
#endif
      jDT=MapDns(iOpt)
      If (jDT<0) Then
         Call mma_allocate(AuxD,nOTSD,nD,Label='AuxD')
         Call mma_allocate(AuxT,nOTSD,nD,Label='AuxT')
         Call mma_allocate(AuxV,nOTSD,nD,Label='AuxV')

         Call RWDTG(-jDT,AuxD,nOTSD*nD,'R','DENS  ',iDisk,SIZE(iDisk,1))
         Call RWDTG(-jDT,AuxT,nOTSD*nD,'R','TWOHAM',iDisk,SIZE(iDisk,1))
         Call RWDTG(-jDT,AuxV,nOTSD*nD,'R','dVxcdR',iDisk,SIZE(iDisk,1))

         D(1:nOTSD,1:nD) => AuxD(:,:)
         T(1:nOTSD,1:nD) => AuxT(:,:)
         V(1:nOTSD,1:nD) => AuxV(:,:)
      Else
         D(1:nOTSD,1:nD) =>   Dens(:,:,jDT)
         T(1:nOTSD,1:nD) => TwoHam(:,:,jDT)
         V(1:nOTSD,1:nD) =>    Vxc(:,:,jDT)
      End If
!----------------------------------------------------------------------*
!
!---- Allocate memory for modified fock matrix
      Call mma_allocate(FckM,nBT,nD,Label='FckM')
      FckM(:,:)=Zero
      G(:,:)=Zero
!
!---- Allocate memory for auxiliary matrices
      Call mma_allocate(Aux1,MaxBas**2,Label='Aux1')
      Call mma_allocate(Aux2,MaxBas**2,Label='Aux2')
      Call mma_allocate(Aux3,MaxBas**2,Label='Aux3')
!
      Do iD = 1, nD
!
         FckM(:,iD) = O(:) + T(:,iD)
#ifdef _DEBUGPRINT_
         Write (6,*) 'iD=',iD
         Call NrmClc(FckM(1,iD),nBT,'EGrad','FckM')
#endif
         If (nnFr.gt.0) Call ModFck(FckM(:,iD),S,nBT,C(:,iD),nBO,nOcc(:,iD))
!
         FckM(:,iD) = FckM(:,iD) + V(:,iD)
#ifdef _DEBUGPRINT_
         Call NrmClc(FckM(1,iD),nBT,'EGrad','FckM')
#endif
!
         iOff = 0
         ij = 1
         it = 1
         ig = 1
         Do iSym = 1, nSym
            nBs = nBas(iSym)
            nOr = nOrb(iSym)
            nOrbmF = nOrb(iSym)-nFro(iSym)
!
            If (nOrb(iSym).gt.0) Then
!
!----------    Square Fock matrix and perform C(T)F
               Aux2(:)=Zero
               Call Square(FckM(ij:,iD),Aux2,1,nBs,nBs)
#ifdef _DEBUGPRINT_
         Write (6,*) 'iSym=',iSym
         Call NrmClc(Aux2,nBs*nBs,'EGrad','Aux2')
#endif
               Aux1(:)=Zero
               Call DGEMM_('T','N',               &
                           nOr,nBs,nBs,           &
                           One,C(it,iD),nBs,      &
                               Aux2,nBs,          &
                          Zero,Aux1,nOr)
#ifdef _DEBUGPRINT_
         Call NrmClc(Aux1,nOr*nBs,'EGrad','Aux1')
#endif
!
!----------    Square density matrix and perform C(T)FD
               Aux2(:)=Zero
               Call DSq(D(ij:,iD),Aux2,1,nBs,nBs)
#ifdef _DEBUGPRINT_
         Call NrmClc(Aux2,nBs*nBs,'EGrad','Aux2')
#endif
               Aux3(:)=Zero
               Call DGEMM_('N','N',               &
                           nOr,nBs,nBs,           &
                           One,Aux1,nOr,          &
                               Aux2,nBs,          &
                          Zero,Aux3,nOr)
#ifdef _DEBUGPRINT_
         Call NrmClc(Aux3,nOr*nBs,'EGrad','Aux3')
#endif
!
!----------    Square overlap matrix and perform C(T)FDS
               Aux2(:)=Zero
               Call Square(S(ij:),Aux2,1,nBs,nBs)
#ifdef _DEBUGPRINT_
         Call NrmClc(Aux2,nBs*nBs,'EGrad','Aux2')
#endif
               Aux1(:)=Zero
               Call DGEMM_('N','N',               &
                           nOr,nBs,nBs,           &
                           One,Aux3,nOr,          &
                                 Aux2,nBs,        &
                           Zero,Aux1,nOr)
#ifdef _DEBUGPRINT_
         Call NrmClc(Aux1,nOr*nBs,'EGrad','Aux1')
#endif
!----------    C(T)FDSC
               Aux2(:)=Zero
               Call DGEMM_('N','N',               &
                           nOr,nOr,nBs,           &
                           One,Aux1,nOr,          &
                                 C(it,iD),nBs,    &
                           Zero,Aux2,nOr)
#ifdef _DEBUGPRINT_
         Call NrmClc(Aux2,nOr*nOr,'EGrad','Aux2')
#endif
!
               Call Asym(Aux2,G(ig,iD),nOr)
#ifdef _DEBUGPRINT_
      Write (6,*)
      Call NrmClc(G,nG   *nD,'EGrad','G')
      Write (6,*)
#endif
!
!              At this point enforce that the gradient is exactly zero
!              for elements corresponding to orbitals of different
!              fermion types.
!
               Do i = 1, nOr
                  If (i.le.nFro(iSym)) Then
                     k=-1
                  Else
                     k=OrbType(iOff+i-nFro(iSym),iD)
                  End If
!
                  Do j = 1, nOr
                     If (j.le.nFro(iSym)) Then
                        l=-1
                     Else
                        l=OrbType(iOff+j-nFro(iSym),iD)
                     End If
!
                     ih = ig + (i-1)*nOr + j - 1
                     If (k/=l) G(ih,iD)=Zero
!
                  End Do
               End Do
!
            End If
            ij = ij + nBs*(nBs + 1)/2
            it = it + nBs*nOr
            ig = ig + nOr*nOr
            iOff = iOff + nOrbmF
!
         End Do ! iSym
!
      End Do ! iD
!
!---- Deallocate memory
      Call mma_deallocate(Aux3)
      Call mma_deallocate(Aux2)
      Call mma_deallocate(Aux1)
      Call mma_deallocate(FckM)
!
      G(:,:)=Two*G(:,:)
!
#ifdef _DEBUGPRINT_
      Block
         Real*8 GMax
         Integer  ::i_Max=0, j_Max=0
         GMax=Zero
         Do i = 1, nD
           Do j = 1, nG
             If (GMax<Abs(G(j,i))) Then
                 GMax=Abs(G(j,i))
                 i_Max=i
                 j_Max=j
               End If
           End Do
         End Do
         Write (6,*) 'GMax,i_Max,j_Max=',GMax,i_Max, j_Max
      End Block
      Call NrmClc(G,nG   *nD,'EGrad','G')
#endif
      If (jDT<0) Then
         Call mma_deallocate(AuxD)
         Call mma_deallocate(AuxT)
         Call mma_deallocate(AuxV)
      End If
      Nullify(D,T,V)
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
!
      Contains
      SubRoutine Asym(H,A,n)
      Implicit None

      Integer n, i, j
      Real*8 H(n,n), A(n,n)

      Do i = 1, n
         Do j = 1, i
            A(i,j) = H(i,j) - H(j,i)
         End Do
      End Do
      End Subroutine ASym

      End Subroutine EGrad
