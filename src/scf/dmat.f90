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
!               2022, Roland Lindh                                     *
!***********************************************************************
!#define _DEBUGPRINT_
      SubRoutine DMat(XCf,nXCf,nD)
!***********************************************************************
!                                                                      *
! Purpose: Compute aufbau density matrix                               *
!                                                                      *
!***********************************************************************
      Use Interfaces_SCF, Only: dOne_SCF, MinDns
      use InfSCF, only: DDnOFF, DNorm, iDisk, InVec, iPsLst, Iter, MiniDN, nBas
      use InfSCF, only: nBT, nDens, nDsk, nFrz, nIterP, nMem, nOrb, nSym, MapDns, nIter
      use stdalloc, only: mma_allocate, mma_deallocate
      use Constants, only: One
      use SCF_Arrays, only: Dens, TwoHam, Vxc, CMO, OccNo
!   Dens    : density matrix - vector containing some (NumDT) last     *
!             (optimized) density matrix differences - (nDT,nD,NumDT)  *
!   TwoHam  : two-el. part of the Fock matrix - vector containing      *
!             corresponding 2-el. contributions - (nDT,nD,NumDT)       *
!   Vxc     : Vxc     part of the Fock matrix - vector containing      *
!             corresponding 2-el. contributions - (nDT,nD,NumDT)       *
!   CMO     : molecular orbitals of length nCMO                        *
!   OccNo   : occupation numbers of length lthO                        *
      Implicit None
      Integer nXCF, nD
      Real*8  XCf(nXCf,nD)
!----------------------------------------------------------------------*
! Local variables                                                      *
!----------------------------------------------------------------------*
      Logical alpha_density
      Real*8, Dimension(:), Allocatable:: Aux
      Integer iFrom, iOnDsk, iD
      Real*8, External:: DDot_
      Integer nCMO
!----------------------------------------------------------------------*
! Start                                                                *
!----------------------------------------------------------------------*
      nCMO=Size(CMO,1)
!     Call Timing(Cpu1,Tim1,Tim2,Tim3)
!
!     Form proper MapDns vector
!
      If (MapDns(iter).eq.0) Then   ! Position not defined
!
!        Update MapDns and eventually write earlier Dens, TwoHam,
!        and Vxc matrices to disk.
!
         nDsk=Max(0,iter-nMem) ! is there too many densities?
!
         If (nDsk.eq.0) Then    ! keep the array in memory
!
            MapDns(iter)=iter
!
         Else
!
            iFrom=MapDns(iter-nMem) ! get index of array to dump
!
            MapDns(iter)=iFrom      ! assign to new array
!
            If (iter-nMem.eq.1) Then
               MapDns(iter-nMem)=-1 !  Initiate
            Else
               MapDns(iter-nMem)=MapDns(iter-nMem-1)-1
            End If
!
!           Dump the vectors
!
            iOnDsk=-MapDns(iter-nMem)
            Call RWDTG(iOnDsk,Dens(1,1,iFrom),  nBT*nD,'W','DENS  ',iDisk,SIZE(iDisk,1))
            Call RWDTG(iOnDsk,TwoHam(1,1,iFrom),nBT*nD,'W','TWOHAM',iDisk,SIZE(iDisk,1))
            Call RWDTG(iOnDsk,Vxc   (1,1,iFrom),nBT*nD,'W','dVxcdR',iDisk,SIZE(iDisk,1))
         End If
!
      End If
!
! Check if MapDns is correct
!
      iPsLst=MapDns(iter)
      If (iPsLst.le.0) Then
         Write (6,*) 'DMat: iPsLst.le.0'
         Write (6,*) 'iPsLst=',iPsLst
         Call Abend()
      End If
!
!
! Form i-th density matrix in the position iPsLst
!
      If (InVec.eq.3 .and. iter.eq.1) Then
!
! First density matrix is actually in the first position
! (read from RUNFILE) on the first iteration.
!
      Else
!
!        Using the CMOs generate the new density in position iPsLst
!
         alpha_density=.True.
         Do iD = 1, nD
            Call DOne_SCF(nSym,nBas,nOrb,nFrz,CMO(1,iD),nCMO,OccNo(1,iD),Dens(1,iD,iPsLst),alpha_density)
            alpha_density=.False.
         End Do ! iD
      End If
!
      Do iD = 1, nD
         Call ChkTrD(nSym,nBas,nOrb,OccNo(:,iD),Size(OccNo,1),Dens(:,iD,iPsLst),Size(Dens,1))
      End Do ! iD
!
#ifdef _DEBUGPRINT_
      Call NrmClc(Dens(1,1,iPsLst),nBT*nD,'DMat  ','D Iter    ')
#endif
!
!     Put the actual densities on the run file.
!
!     Note, this is from where the DFT part of the DrvXV code gets the
!     total density when it computes the DFT contributions to the
!     total energy and the Fock matrix.
!
      Call DensAB(nBT,iPsLst,nD,Dens)
!
! Form density difference (normal or minimized). Notice, that
! for minimized differences i-th density is copied to position
! nDens and kept there till the end of iteration.
! For normal differences position nDens is occupied by the density
! of the previous iteration, until
! (i)  OptClc is called (DIIS-only case), or
! (ii) after Diis (QNR case)
! after that, the actual density is stored at that memory location
!
      If(MiniDn.and. Max(0,nIter(nIterP)-1).gt.0) Then
!
!        Minimized density option
!
         Call DCopy_(nBT*nD,Dens(1,1,iPsLst),1,Dens(1,1,nDens),1)
         If (iter.gt.1) Call MinDns(Dens,nBT,nDens,XCf,nXCf,nD)
!
      Else If (.not.DDnOFF) Then
!
!        Do the density difference, D(iPsLst)-D(nDens)=D(k+1)-D(k)
!
         Call mma_allocate(Aux,nBT*nD,Label='Aux')
         Call DCopy_(nBT*nD,Dens(1,1,iPsLst),1,Aux,1)
         Call Daxpy_(nBT*nD,-One,Dens(1,1,nDens ),1,Dens(1,1,iPsLst),1)
         Call DCopy_(nBT*nD,Aux,1,Dens(1,1,nDens),1)
         Call mma_deallocate(Aux)
!
      Else
!
         call DCopy_(nBT*nD,Dens(1,1,iPsLst),1,Dens(1,1,nDens ),1)
!
      End If
!
      DNorm=DBLE(nD)*DDot_(nBT*nD,Dens(1,1,iPsLst),1,Dens(1,1,iPsLst),1)
!
#ifdef _DEBUGPRINT_
      Write (6,*) 'DNorm=',DNorm
      Call NrmClc(Dens(1,1,iPsLst),nBT*nD,'DMat  ','D iPsLst  ')
      Call NrmClc(Dens(1,1,nDens), nBT*nD,'DMat  ','D nDens   ')
#endif
!
!      Call Timing(Cpu2,Tim1,Tim2,Tim3)
!      TimFld( 4) = TimFld( 4) + (Cpu2 - Cpu1)
!----------------------------------------------------------------------*
! Exit                                                                 *
!----------------------------------------------------------------------*
      Return
      End SubRoutine DMat
