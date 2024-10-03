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
subroutine TraFck(canorb,FOVMax)
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

use Index_Functions, only: nTri_Elem
use SpinAV, only: Do_SpinAV
use InfSCF, only: CMO, EOrb, FckAuf, FockAO, MaxBas, nBas, nBO, nBT, nConstr, nFro, nnFr, nOcc, nOrb, nSym, Ovrlp, TimFld
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp), intent(out) :: FOVMax
logical(kind=iwp), intent(in) :: CanOrb
integer(kind=iwp) :: ia, iCMO, iD, iDiag, iDum, iErr, iiCMO, iiEigV, iiScratch, ioFckM, iOff, iptr, iptr2, iSym, jEOr, jjEOr, kk, &
                     kOcc, kOff, n2Sort, n2Zero, nD, nFound, nOccmF, nOrbmF, nVrt
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: jCMO
#endif
real(kind=wp) :: Cpu1, Cpu2, Dummy, Fia, Tim1, Tim2, Tim3, Tmp, Tmp0, Tmp1
real(kind=wp), allocatable :: CMOOld(:), COvrlp(:), Ctmp(:), EigV(:), FckM(:), FckS(:), HlfF(:), Scratch(:), Scrt(:)
integer(kind=iwp), external :: iDaMax_
real(kind=wp), external :: DDot_

nD = size(FockAO,2)

call Timing(Cpu1,Tim1,Tim2,Tim3)
#ifdef _DEBUGPRINT_
call NrmClc(FockAO,size(FockAO),'TraFck','FockAO')
#endif
! allocate memory for modified Fock matrix
call mma_allocate(FckM,nBT,Label='FckM')
! allocate memory for squared Fock matrix
call mma_allocate(FckS,MaxBas**2,Label='FckS')

FOVMax = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
do iD=1,nD
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! modify Fock matrix
  FckM(:) = FockAO(:,iD)
  if (nnFr > 0) call ModFck(FckM,Ovrlp,nBT,CMO(1,iD),nBO,nOcc(1,iD))

  ioFckM = 1
  iCMO = 1
  jEOr = 1
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  do iSym=1,nSym
    !                                                                  *
    !*******************************************************************
    !                                                                  *
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iD=',iD
    call RecPrt('TraFck: Old CMO',' ',CMO(iCMO,iD),nBas(iSym),nOrb(iSym))
    jCMO = iCMO
#   endif

    nOrbmF = nOrb(iSym)-nFro(iSym)
    nOccmF = nOcc(iSym,iD)-nFro(iSym)
    nVrt = nOrb(iSym)-nOcc(iSym,iD)

    iCMO = iCMO+nBas(iSym)*nFro(iSym) ! pointer to CMO

    if (nOrbmF > 0) then

      ! allocate memory for half-transformed Fock matrix
      call mma_allocate(HlfF,nBas(iSym)*nOrbmF,Label='HlfF')
      ! transform Fock matrix into new MO space
      call Square(FckM(ioFckM),FckS,1,nBas(iSym),nBas(iSym))
      call DGEMM_('N','N',nBas(iSym),nOrbmF,nBas(iSym), &
                  One,FckS,nBas(iSym), &
                  CMO(iCMO,iD),nBas(iSym), &
                  Zero,HlfF,nBas(iSym))
      call DGEMM_Tri('T','N',nOrbmF,nOrbmF,nBas(iSym), &
                     One,CMO(iCMO,iD),nBas(iSym), &
                     HlfF,nBas(iSym), &
                     Zero,FckS,nOrbmF)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'transformed Fck in trafck:'
      call TriPrt(' ',' ',FckS,nOrbmF)
#     endif
      ! dispose memory of half-transformed Fock matrix
      call mma_deallocate(HlfF)

      if ((nOccmF > 0) .and. (nVrt > 0)) then
        iptr = 1+nTri_Elem(nOccmF)

        ! get max Fock Matrix Element in OV block...

        do ia=1,nVrt
          Fia = abs(FckS(iptr+IDAMAX_(nOccmF,FckS(iptr),1)-1))
          FOVMax = max(Fia,FOVMax)
          iptr = iptr+nOccmF+ia
        end do
      end if
      !                                                                *
      !----------------------------------------------------------------*
      !                                                                *
      ! eventually diagonalize occ/occ & virt/virt Block
      ! separately to form canonical orbitals

      if (canorb) call Mk_CanOrb()

    end if
    !                                                                  *
    !------------------------------------------------------------------*
    !                                                                  *
    ! Update pointers
    iCMO = iCMO+nOrbmF*nBas(iSym)
    ioFckM = ioFckM+nTri_Elem(nBas(iSym))
    jEOr = jEOr+nOrbmF
#   ifdef _DEBUGPRINT_
    call RecPrt('TraFck: New CMO',' ',CMO(jCMO,iD),nBas(iSym),nOrb(iSym))
#   endif
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  end do ! iSym
  !                                                                    *
  !*********************************************************************
  !                                                                    *
# ifdef _DEBUGCHECK__
  ! Check orthogonality
  if (canorb) call ChkOrt(iD,Whatever)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end do ! iD
!                                                                      *
!***********************************************************************
!                                                                      *
! Deallocate memory
call mma_deallocate(FckS)
call mma_deallocate(FckM)

call Timing(Cpu2,Tim1,Tim2,Tim3)
TimFld(13) = TimFld(13)+(Cpu2-Cpu1)

return

contains

subroutine Mk_CanOrb()

  integer(kind=iwp) :: i, ia, ii, iOcc, j, jOcc

  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! allocate space for Eigenvectors of occ/occ & virt/virt block
  call mma_allocate(EigV,max(nOccmF,nVrt)**2,Label='EigV')
  ! allocate space for temporary CMO (for MatMul)
  call mma_allocate(CTmp,max(nOccmF,nVrt)*nBas(iSym),Label='CTmpX')
  if (.not. FckAuf) then
    call mma_allocate(CMOOld,nOccmF*nBas(iSym),Label='CMOOld')
    CMOOld(:) = CMO(iCMO:iCMO+nOccmF*nBas(iSym)-1,iD)
    call mma_allocate(Scrt,nBas(iSym)**2,Label='Scrt')
    call mma_allocate(COvrlp,nBas(iSym)*nOccmF,Label='COvrlp')
  end if

  if (nOccmF > 0) then

    ! diagonalize occ/occ first find the proper pointer to EOr

    jEOr = jEOr+nFro(iSym)
    call mma_allocate(Scratch,nOccmF**2,Label='Scratch')
    Dummy = Zero
    iDum = 0
    nOccmF = nOccmF-nConstr(iSym)

    call Diag_Driver('V','A','L',nOccmF,FckS,Scratch,nOccmF,Dummy,Dummy,iDum,iDum,EOrb(jEOr,iD),EigV,nOccmF,1,0,'J',nFound,iErr)

    if (nConstr(iSym) > 0) then
      Scratch(1:(nOccmF+nConstr(iSym))**2) = Zero
      do j=1,nOccmF
        iiEigV = 1+nOccmF*(j-1)
        iiScratch = 1+(nOccmF+nConstr(iSym))*(j-1)
        Scratch(iiScratch:iiScratch+nOccmF-1) = EigV(iiEigV:iiEigV+nOccmF-1)
      end do
      do j=nOccmF+1,nOccmF+nConstr(iSym)
        iiScratch = 1+(nOccmF+nConstr(iSym))*(j-1)+j-1
        Scratch(iiScratch) = One
      end do
      EigV(1:(nOccmF+nConstr(iSym))**2) = Scratch(1:(nOccmF+nConstr(iSym))**2)
    end if
    nOccmF = nOccmF+nConstr(iSym)
    call mma_deallocate(Scratch)
    n2zero = nOccmF
    if (Do_SpinAV) n2zero = n2zero+nConstr(iSym)
    FckS(1:nTri_Elem(n2zero)) = Zero

    iDiag = 0
    do i=1,n2zero
      iDiag = iDiag+i
      FckS(iDiag) = EOrb(jEOr+i-1,iD)
    end do

    ! Rotate MOs to diagonalize occ/occ block
    do ii=0,nOccmF-1
      CTmp(nBas(iSym)*ii+1:nBas(iSym)*(ii+1)) = CMO(iCMO+nBas(iSym)*ii:iCMO+nBas(iSym)*(ii+1)-1,iD)
    end do
    call DGEMM_('N','N',nBas(iSym),nOccmF,nOccmF, &
                One,Ctmp,nBas(iSym), &
                EigV,nOccmF, &
                Zero,CMO(iCMO,iD),nBas(iSym))

    ! Fix standard phase pf the orbitals

    do i=1,nOccmF
      call VecPhase(CMO(iCMO+(i-1)*nBas(iSym),iD),nBas(iSym))
    end do

    ! Order the occupied orbitals by maximum overlap with the old MOs.

    if (.not. FckAuf) then

      COvrlp(:) = Zero
      call Square(Ovrlp(ioFckM),Scrt,1,nBas(iSym),nBas(iSym))
      call DGEMM_('T','N',nOccmF,nBas(iSym),nBas(iSym), &
                  One,CMO(iCMO,iD),nBas(iSym), &
                  Scrt,nBas(iSym), &
                  Zero,COvrlp,nOccmF)

      do iOcc=1,nOccmF-1  !  Loop over the new MOs

        iOff = (iOcc-1)*nBas(iSym)+iCMO
        kOcc = 0
        Tmp0 = Zero
        do jOcc=1,nOccmF !  Loop over the new MOs
          Tmp1 = abs(DDot_(nBas(iSym),COvrlp(jOcc),nOccmF,CMO(iOff,iD),1))
          if (Tmp1 > Tmp0) then
            Tmp0 = Tmp1
            kOcc = jOcc
          end if
        end do

        if (iOcc /= kOcc) then
          ii = iOcc+jEOr-1
          kk = kOcc+jEOr-1
          tmp = EOrb(ii,iD)
          EOrb(ii,iD) = EOrb(kk,iD)
          EOrb(kk,iD) = tmp
          kOff = (kOcc-1)*nBas(iSym)+iCMO
          call DSwap_(nBas(iSym),CMO(iOff,iD),1,CMO(kOff,iD),1)
        end if
      end do

      iDiag = 0
      do i=1,n2zero
        iDiag = iDiag+i
        FckS(iDiag) = EOrb(jEOr+i-1,iD)
      end do

    end if
  end if

  if (nVrt > 0) then

    ! now diagonalize virt/virt block
    ! setup virt/virt block in triangular Fock Matrix

    iptr = 1+nTri_Elem(nOccmF)+nOccmF
    iptr2 = 1
    do ia=1,nVrt
      FckS(iptr2:iptr2+ia-1) = FckS(iptr:iptr+ia-1)
      iptr = iptr+nOccmF+ia
      iptr2 = iptr2+ia
    end do
    call mma_allocate(Scratch,nVrt**2,Label='Scratch')
    if (Do_SpinAV) then
      nVrt = nVrt-nConstr(iSym)
      jEOr = jEOr+nConstr(iSym)
      iptr = 1+nTri_Elem(nOccmF)+nOccmF
      do ia=1,nConstr(iSym)
        iptr = iptr+nOccmF+ia
        FckS(iptr) = -0.666e3_wp*real(1000-ia,kind=wp)
      end do
    end if
    Dummy = Zero
    iDum = 0
    call Diag_Driver('V','A','L',nVrt,FckS,Scratch,nVrt,Dummy,Dummy,iDum,iDum,EOrb(jEOr+nOccmF,iD),EigV,nVrt,1,0,'J',nFound,iErr)
    if (Do_SpinAV) then
      nVrt = nVrt+nConstr(iSym)
      jEOr = jEOr-nConstr(iSym)
      Scratch(1:nVrt**2) = Zero
      do j=1,nVrt-nConstr(iSym)
        iiEigV = 1+(nVrt-nConstr(iSym))*(j-1)
        iiScratch = 1+nVrt*(nConstr(iSym)+j-1)
        Scratch(iiScratch:iiScratch+nVrt-nConstr(iSym)-1) = EigV(iiEigV:iiEigV+nVrt-nConstr(iSym)-1)
      end do
      do j=1,nConstr(iSym)
        iiScratch = 1+nVrt*(j-1)+j-1
        Scratch(iiScratch) = One
      end do
      EigV(1:nVrt**2) = Scratch(1:nVrt**2)
    end if
    call mma_deallocate(Scratch)
    ! rotate MOs to diagonalize virt/virt block
    iptr = iCMO+nOccmF*nBas(iSym)
    do ia=0,nVrt-1
      CTmp(nBas(iSym)*ia+1:nBas(iSym)*(ia+1)) = CMO(iptr+nBas(iSym)*ia:iptr+nBas(iSym)*(ia+1)-1,iD)
    end do
    call DGEMM_('N','N',nBas(iSym),nVrt,nVrt, &
                One,Ctmp,nBas(iSym), &
                EigV,nVrt, &
                Zero,CMO(iptr,iD),nBas(iSym))
  end if

  if (nConstr(iSym) > 0) then
    ! Sort non-wavelet eigenvalues/eigenvectors
    n2sort = nOccmF-nConstr(iSym)
    if (FckAuf) call SortEig(EOrb(jEOr,iD),CMO(iCMO,iD),n2sort,nBas(iSym),1,.true.)
    jjEOr = jEOr+nOccmF+nConstr(iSym)
    iiCMO = iCMO+nBas(iSym)*(nOccmF+nConstr(iSym))
    n2sort = nVrt-nConstr(iSym)
    if (FckAuf) call SortEig(EOrb(jjEOr,iD),CMO(iiCMO,iD),n2sort,nBas(iSym),1,.true.)
  else
    ! Sort all eigenvalues and eigenvectors
    if (FckAuf) call SortEig(EOrb(jEOr,iD),CMO(iCMO,iD),nOrbmF,nBas(iSym),1,.true.)
  end if
  ! dispose memory
  if (.not. FckAuf) then
    call mma_deallocate(COvrlp)
    call mma_deallocate(Scrt)
    call mma_deallocate(CMOOld)
  end if
  call mma_deallocate(CTmp)
  call mma_deallocate(EigV)

end subroutine Mk_CanOrb

end subroutine TraFck
