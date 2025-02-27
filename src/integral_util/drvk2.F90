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
! Copyright (C) 1990,1991,1993, Roland Lindh                           *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Drvk2(DoFock,DoGrad)
!***********************************************************************
!                                                                      *
!  Object: to precompute all pair entites as zeta, kappa, P and the    *
!          integral prescreening vector to be used with the Schwarz    *
!          inequality.                                                 *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             March '90                                                *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry,            *
!             University of Lund, SWEDEN.                              *
!             June '91, modified for k2 loop.                          *
!             Modified for direct SCF, January '93                     *
!***********************************************************************

use Index_Functions, only: iTri, nTri_Elem1, nTri3_Elem1
use setup, only: mSkal
use iSD_data, only: iSD, nSD
use k2_structure, only: Indk2, k2_Processed, k2Data
use k2_arrays, only: BraKet, Create_BraKet, Destroy_BraKet, DoGrad_, ipOffD, Sew_Scr
use Dens_Stuff, only: ipDij, mDCRij, mDij
use Basis_Info, only: Shells
use Symmetry_Info, only: iOper, nIrrep
use Gateway_global, only: force_part_c
use Sizes_of_Seward, only: S
use UnixInfo, only: ProgName
use stdalloc, only: mma_allocate, mma_deallocate, mma_maxDBLE
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Gateway_Info, only: lSchw
use Definitions, only: u6
#endif

implicit none
logical(kind=iwp), intent(in) :: DoFock, DoGrad
integer(kind=iwp) :: iAng, iAng4(4), iBas, iCmp, iDCRR(0:7), ijCmp, ijInc, ijS, ik2, ipMem1, ipMem2, iPrim, iPrimi, iPrimSave(4), &
                     iS, iSD4(0:nSD,4), iShell, iShll, jAng, jBas, jCmp, jPrim, jPrimj, jS, jShell, jShll, la_, mabMax_, mabMin_, &
                     Mem1, Mem2, MemMax, MemPrm, MemTmp, mk2, mScree, nBasi, nBasj, nDCRR, ne_, nHm, nHrrMtrx, nScree, nSO, nZeta
real(kind=wp) :: Coor(3,4)
logical(kind=iwp) :: force_part_save, ReOrder, Rls
character(len=8) :: Method
real(kind=wp), allocatable :: HRRMtrx(:,:), Knew(:), Lnew(:), Pnew(:), Qnew(:), Scr(:,:)

!                                                                      *
!***********************************************************************
!                                                                      *
if (k2_processed) return

ReOrder = .false.
if (index(ProgName,'scf') /= 0) then
  call Get_cArray('Relax Method',Method,8)
  ReOrder = Method == 'KS-DFT'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
DoGrad_ = DoGrad
la_ = S%iAngMx
mabMin_ = nTri3_Elem1(max(la_,la_)-1)
mabMax_ = nTri3_Elem1(la_+la_)-1
ne_ = (mabMax_-mabMin_+1)
nHrrMtrx = ne_*nTri_Elem1(la_)*nTri_Elem1(la_)
call mma_allocate(HRRMtrx,nHRRMtrx,2,Label='HRRMatrix')
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate memory for k2 data. Observe that the call to Allok2
! can be done elsewhere and then this call will simply result in
! a return.

call Allok2()
!                                                                      *
!***********************************************************************
!                                                                      *
nScree = 0
mScree = 0
mk2 = 0

! Allocate memory for zeta, kappa, and P.
call Create_BraKet(S%m2Max,S%m2Max)
!                                                                      *
!***********************************************************************
!                                                                      *
MemTmp = 0
do iAng=0,S%iAngMx
  MemTmp = max(MemTmp,(S%MaxPrm(iAng)*nTri_Elem1(iAng))**2)
end do
call mma_allocate(Scr,MemTmp,3,Label='Scr')
call mma_allocate(Knew,S%m2Max,Label='Knew')
call mma_allocate(Lnew,S%m2Max,Label='Lnew')
call mma_allocate(Pnew,S%m2Max*3,Label='Pnew')
call mma_allocate(Qnew,S%m2Max*3,Label='Qnew')
!                                                                      *
!***********************************************************************
!                                                                      *
if (allocated(Sew_Scr)) then
  Rls = .false.
  MemMax = size(Sew_Scr)
  !write(u6,*) 'Drvk2: Memory already allocated:',MemMax
else
  Rls = .true.
  call mma_maxDBLE(MemMax)
  if (MemMax > 8000) MemMax = MemMax-8000
  call mma_allocate(Sew_Scr,MemMax,Label='Sew_Scr')
  !write(u6,*) 'Drvk2: Memory allocated:',MemMax
end if
ipMem1 = 1
!                                                                      *
!***********************************************************************
!                                                                      *
! Canonical double loop over shells. This includes both valence and
! auxiliary basis functions

do iS=1,mSkal
  iSD4(:,1) = iSD(:,iS)
  iSD4(:,3) = iSD(:,iS)

  iShll = iSD4(0,1)

  ! In case of auxiliary basis sets we want the iS index to point at the dummay shell.
  if (Shells(iShll)%Aux .and. (iS /= mSkal)) cycle

  iBas = iSD4(3,1)
  iPrim = iSD4(5,1)
  iShell = iSD4(11,1)

  if (ReOrder) call OrdExpD2C(iPrim,Shells(iShll)%Exp,iBas,Shells(iShll)%pCff)

  do jS=1,iS

    call Gen_iSD4(iS,jS,iS,jS,iSD,nSD,iSD4)

    jShll = iSD4(0,2)

    ! In case the first shell is the dummy auxiliary basis shell make sure that
    ! the second shell, jS, also is a auxiliary basis shell.
    if (Shells(iShll)%Aux .and. (.not. Shells(jShll)%Aux)) cycle
    ! Make sure that the second shell never is the dummy auxiliary basis shell.
    if (Shells(jShll)%Aux .and. (jS == mSkal)) cycle

    jBas = iSD4(3,2)
    jPrim = iSD4(5,2)
    jShell = iSD4(11,2)

    call Coor_Setup(iSD4,nSD,Coor)

    iPrimi = iPrim
    jPrimj = jPrim

    nBasi = iBas
    nBasj = jBas

    ! Fake shell 3 and 4
    iSD4(5,3:4) = 1
    iSD4(3,3:4) = 1

    nZeta = iPrim*jPrim

    call ConMax(BraKet%Eta(:),iPrim,jPrim,Shells(iShll)%pCff,nBasi,Shells(jShll)%pCff,nBasj)

#   ifdef _DEBUGPRINT_
    call RecPrt(' Sym. Dist. Centers',' ',Coor,3,4)
#   endif

    ijS = iTri(iShell,jShell)
    if (DoFock) then
      ipDij = ipOffD(1,ijS)
      mDCRij = ipOffD(2,ijS)
      mDij = ipOffD(3,ijS)
    else
      ipDij = -1
      mDCRij = 1
      mDij = 1
    end if

    nSO = 1

    ! Compute memory request for the primitives, i.e. how much
    ! memory is needed up to the transfer equation.

    iAng4(:) = iSD4(1,:)
    call MemRys(iAng4,MemPrm)

    ! Decide on the partitioning of the shells based on
    ! on the available memory and the requested memory

    ! Now do a dirty trick to avoid splitting of the first
    ! contracted index. Move all over on the second index.

    iPrimSave(:) = iSD4(5,:) ! Store away original setting

    iSD4(3,1) = 1
    iSD4(3,2) = nZeta

    iSD4(5,1) = 1
    iSD4(5,2) = nZeta

    force_part_save = force_part_c
    force_part_c = .false.

    call PSOAO0(nSO,MemPrm,MemMax,ipMem1,ipMem2,Mem1,Mem2,.false.,nSD,iSD4)

    force_part_c = force_part_save
    ijInc = min(iSD4(4,2),iSD4(6,2))

    ! restore correct index
    iSD4(5,:) = iPrimSave(:)

    iPrimi = iSD4(5,1)
    jPrimj = iSD4(5,2)

#   ifdef _DEBUGPRINT_
    write(u6,*) ' ************** Memory partitioning **************'
    write(u6,*) ' ipMem1=',ipMem1
    write(u6,*) ' ipMem2=',ipMem2
    write(u6,*) ' Mem1=',Mem1
    write(u6,*) ' Mem2=',Mem2
    write(u6,*) ' ***********************************************'
#   endif

    ! Find the Double Coset Representatives for center A and B.

    iDCRR(0:nIrrep-1) = iOper(0:nIrrep-1)
    nDCRR = nIrrep

    ! Compute all pair entities (zeta, kappa, P, and [nm|nm],
    ! total of six types) for all possible unique pairs of
    ! centers generated for the symmetry unique centers A and B.

    iAng = iSD4(1,1)
    jAng = iSD4(1,2)
    iCmp = iSD4(2,1)
    jCmp = iSD4(2,2)
    nHm = iCmp*jCmp*(nTri3_Elem1(iAng+jAng)-nTri3_Elem1(max(iAng,jAng)-1))
    nHm = nHm*nIrrep
    ijCmp = nTri_Elem1(iAng)*nTri_Elem1(jAng)
    if (.not. DoGrad_) ijCmp = 0

    ik2 = Indk2(3,ijS)

    call k2Loop(Coor,iDCRR,nDCRR,k2data(:,ik2),Shells(iShll)%Exp,iPrimi,Shells(jShll)%Exp,jPrimj,BraKet%xA(:),BraKet%xB(:), &
                Shells(iShll)%pCff,nBasi,Shells(jShll)%pCff,nBasj,BraKet%Zeta(:),BraKet%ZInv(:),BraKet%KappaAB(:),BraKet%P(:,:), &
                BraKet%IndZet(:),nZeta,ijInc,BraKet%Eta(:),Sew_Scr(ipMem2),Mem2,nScree,mScree,ijCmp,DoFock,Scr,MemTmp,Knew,Lnew, &
                Pnew,Qnew,S%m2Max,DoGrad,HrrMtrx,nHrrMtrx,nSD,iSD4)

    Indk2(2,ijS) = nDCRR
    mk2 = mk2+nDCRR

  end do
end do

call Destroy_Braket()
!                                                                      *
!***********************************************************************
!                                                                      *
if (Rls) call mma_deallocate(Sew_Scr)
call mma_deallocate(Qnew)
call mma_deallocate(Pnew)
call mma_deallocate(Lnew)
call mma_deallocate(Knew)
call mma_deallocate(Scr)
call mma_deallocate(HRRMtrx)
!                                                                      *
!***********************************************************************
!                                                                      *
!rScree = One-(One*mScree)/(One*nScree)
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) ' *** The k2 entities has been precomputed ***'
write(u6,'(I7,A)') mk2,' blocks of k2 data were computed.'
if (lSchw) then
  write(u6,*) ' Prescreening based on primitive integrals.'
else
  write(u6,*) ' Prescreening based on radial overlap.'
end if
!write(u6,'(1X,A,F7.5)') 'Pair screening ratio:',rScree
write(u6,*)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *

k2_processed = .true.

end subroutine Drvk2
