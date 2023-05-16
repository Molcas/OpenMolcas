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
! Copyright (C) 2022, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine Mk_FockMO(O,S,nOTSD,C,nC,FockMO,nFockMO,nD,iOpt)
      Use Orb_Type, only: OrbType
      use InfSCF, only: MaxBas, nBO, nBT, nnFr, nSym, nBas, nOrb, nFro, nOcc, MapDns, iDisk
      use SCF_Arrays, only: TwoHam, Vxc
      use Constants, only: Zero, One
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
!
      Integer nOTSD, nD, nC, nFockMO, iOpt
      Real*8 O(nOTSD),S(nOTSD),C(nC,nD), FockMO(nFockMO,nD)
!
      Integer i, j, k, l, iD, ig, ih, ij, iOff, it, nOr, nOrbmF, nBs, iSym, jDT
      Real*8, Dimension(:,:), Allocatable:: FckM
      Real*8, Dimension(:), Allocatable:: Aux1, Aux2
      Real*8, Dimension(:,:), Allocatable, Target:: AuxT, AuxV
      Real*8, Dimension(:,:), Pointer:: T, V
!
!----------------------------------------------------------------------*
!     Start
      jDT=MapDns(iOpt)
      If (jDT<0) Then
         Call mma_allocate(AuxT,nOTSD,nD,Label='AuxT')
         Call mma_allocate(AuxV,nOTSD,nD,Label='AuxV')

         Call RWDTG(-jDT,AuxT,nOTSD*nD,'R','TWOHAM',iDisk,SIZE(iDisk,1))
         Call RWDTG(-jDT,AuxV,nOTSD*nD,'R','dVxcdR',iDisk,SIZE(iDisk,1))

         T(1:nOTSD,1:nD) => AuxT(:,:)
         V(1:nOTSD,1:nD) => AuxV(:,:)
      Else
         T(1:nOTSD,1:nD) => TwoHam(:,:,jDT)
         V(1:nOTSD,1:nD) =>    Vxc(:,:,jDT)
      End If
#ifdef _DEBUGPRINT_
      Write (6,*) 'EGrad: input arrays'
      Write (6,*) '==================================================='
      Call NrmClc(O,nOTSD   ,'EGrad','O')
      Call NrmClc(T,nOTSD*nD,'EGrad','T')
      Call NrmClc(V,nOTSD*nD,'EGrad','V')
      Call NrmClc(C,nC   *nD,'EGrad','C')
!     Do iD = 1, nD
!        Write (*,*) 'OrbType(:,iD)', OrbType(:,iD)
!     End Do ! iD
      Write (6,*) '==================================================='
      Write (6,*)
#endif
!----------------------------------------------------------------------*
!
!---- Allocate memory for modified fock matrix
      Call mma_allocate(FckM,nBT,nD,Label='FckM')
      FckM(:,:)=Zero
      FockMO(:,:)=Zero
!
!---- Allocate memory for auxiliary matrices
      Call mma_allocate(Aux1,MaxBas**2,Label='Aux1')
      Call mma_allocate(Aux2,MaxBas**2,Label='Aux2')
!
      Do iD = 1, nD
!
         FckM(:,iD) = O(:) + T(:,iD)
#ifdef _DEBUGPRINT_
         Write (6,*) 'iD=',iD
         Call NrmClc(FckM(1,iD),nBT,'Mk_FockMO','FckM')
#endif
         If (nnFr.gt.0) Call ModFck(FckM(1,iD),S,nBT,C(1,iD),nBO,nOcc(1,1))
!
         FckM(:,iD) = FckM(:,iD) + V(:,iD)
#ifdef _DEBUGPRINT_
         Call NrmClc(FckM(1,iD),nBT,'Mk_FockMO','FckM')
#endif
!
         iOff = 0
         ij = 1
         it = 1
         iG = 1
         Do iSym = 1, nSym
            nBs = nBas(iSym)
            nOr = nOrb(iSym)
            nOrbmF = nOrb(iSym)-nFro(iSym)
!
            If (nOrb(iSym).gt.0) Then
!
!----------    Square Fock matrix and perform C(T)F
               Aux2(:)=Zero
               Call Square(FckM(ij,iD),Aux2,1,nBs,nBs)
#ifdef _DEBUGPRINT_
         Write (6,*) 'iSym=',iSym
         Call NrmClc(Aux2,nBs*nBs,'Mk_FockMO','Aux2')
#endif
               Aux1(:)=Zero
               Call DGEMM_('T','N',nOr,nBs,nBs,        &
                           One,C(it,iD),nBs,           &
                               Aux2,nBs,               &
                           Zero,Aux1,nOr)
#ifdef _DEBUGPRINT_
         Call NrmClc(Aux1,nOr*nBs,'Mk_FockMO','Aux1')
#endif
!----------    C(T)FDSC
               Aux2(:)=Zero
               Call DGEMM_('N','N',nOr,nOr,nBs,        &
                           One,Aux1,nOr,               &
                               C(it,iD),nBs,           &
                           Zero,FockMO(iG,iD),nOr)
#ifdef _DEBUGPRINT_
         Call NrmClc(FockMO(iG,iD),nOr*nOr,'Mk_FockMO','FockMO')
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
                     ih = iG + (i-1)*nOr + j - 1
                     If (k/=l) FockMO(ih,iD)=Zero
!
                  End Do
               End Do
!
            End If
            ij = ij + nBs*(nBs + 1)/2
            it = it + nBs*nOr
            iG = iG + nOr*nOr
            iOff = iOff + nOrbmF
!
         End Do ! iSym
!
      End Do ! iD
!
!---- Deallocate memory
      Call mma_deallocate(Aux2)
      Call mma_deallocate(Aux1)
      Call mma_deallocate(FckM)
!
      If (jDT<0) Then
         Call mma_deallocate(AuxT)
         Call mma_deallocate(AuxV)
      End If
      Nullify(T,V)
#ifdef _DEBUGPRINT_
      Call RecPrt('Mk_FockMO: FockMO',' ',FockMO,1,Size(FockMO))
#endif
!
!----------------------------------------------------------------------*
!     Exit                                                             *
!----------------------------------------------------------------------*
!
      End Subroutine Mk_FockMO
