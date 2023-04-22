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
! Copyright (C) 1990, Roland Lindh                                     *
!               1995, Anders Bernhardsson                              *
!***********************************************************************

subroutine Cnt1El2(Kernel,KrnlMm,Label,iDCnt,iDCar,loper,rHrmt,DiffOp,Lab_Dsk,iadd,isym,kcar,nordop)
!***********************************************************************
!                                                                      *
! Object: to compute the one-electron integrals. The method employed at*
!         this point is not necessarily the fastest. However, the total*
!         time for the computation of integrals will depend on the time*
!         spent in computing the two-electron integrals.               *
!         The memory at this point is assumed to be large enough to do *
!         the computation in core.                                     *
!         The data is structured with respect to four indices, two (my *
!         ny or i j) refer to primitives or basis functions and two (a *
!         b) refer to the components of the cartesian or spherical     *
!         harmonic gaussians.                                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January '90                                              *
!             Rewritten for gradients needed in hessian calculations   *
!             and general operator treatment                           *
!             May '95 By:                                              *
!             Anders Bernhardsson , Dept. of Theoretical Chemistry,    *
!             University  of Lund, SWEDEN.                             *
!***********************************************************************

use McKinley_global, only: nFck, sIrrep
use mck_interface, only: mck_mem, oneel_mck_kernel
use Index_Functions, only: nTri_Elem, nTri_Elem1
use Real_Spherical, only: ipSph, RSph
use iSD_data, only: iSD
use Basis_Info, only: dbsc, MolWgh, nBas, Shells
use Center_Info, only: dc
use Symmetry_Info, only: iOper, nIrrep
use Sizes_of_Seward, only: S
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
procedure(oneel_mck_kernel) :: Kernel
procedure(mck_mem) :: KrnlMm
character(len=8), intent(in) :: Label, Lab_Dsk
integer(kind=iwp), intent(in) :: iDCnt, iDCar, iadd, isym, kcar, nordop
integer(kind=iwp), intent(out) :: loper
real(kind=wp), intent(in) :: rHrmt
logical(kind=iwp), intent(in) :: DiffOp
#include "Molcas.fh"
#include "disp.fh"
integer(kind=iwp) :: iAng, iAO, iBas, iCar, iCmp, iCnt, iCnttp, iComp, iDCRR(0:7), iDCRT(0:7), iI, iIC, iIrrep, IndGrd(0:7), iopt, &
                     ip(8), iPrim, irc, iS, iShell, iShll, iSmLbl, iSOBlk, iStabM(0:7), iStabO(0:7), iStart, iuv, jAng, jAO, jBas, &
                     jCmp, jCnt, jCnttp, jdisp, jIrrep, jPrim, jS, jShell, jShll, kk, kOper, lDCRR, LenInt, LenInt_Tot, lFinal, &
                     LmbdR, LmbdT, maxi, mdci, mdcj, MemKer, MemKrn, mSO, nDCRR, NDCRT, nDens, ndenssq, nDisp, nIC, nnIrrep, &
                     nOp(2), nOrder, nrOp, nScr1, nSkal, nSO, nStabM, nStabO
real(kind=wp) :: A(3), B(3), CCoor(3), Fact, RB(3)
logical(kind=iwp) :: DiffCnt, IfGrd(3,2), Trans(2)
character(len=8) :: LabDsk
real(kind=wp), allocatable :: Fnl(:), Integrals(:), Kappa(:), Kern(:), PCoor(:,:), Scr(:), ScrSph(:), SO(:), Zeta(:), ZI(:)
integer(kind=iwp), external :: MemSO1, NrOpr
logical(kind=iwp), external :: EQ, TF

! Compute the number of blocks from each component of the operator
! and the irreps it will span.

LabDsk = Lab_Dsk
CCoor(:) = Zero
IndGrd(0:nIrrep-1) = 0
loper = 0
nnIrrep = nIrrep
if (sIrrep) nnIrrep = 1
do iIrrep=0,nnIrrep-1
  jIrrep = nropr(ieor(ioper(iIrrep),ioper(isym)))
  nDisp = IndDsp(iDcnt,iIrrep)
  do iCar=1,3
    iComp = 2**(iCar-1)
    if (TF(iDCnt,iIrrep,iComp)) then
      ndisp = ndisp+1
      if (iDCar == icar) then
        loper = loper+2**jIrrep
        IndGrd(jIrrep) = nDisp
      end if
    end if
  end do
end do
nIC = 0
if (loper == 0) return

ip(1:nIrrep) = 0

iStart = 1
do iIrrep=0,nIrrep-1
  if (btest(loper,iIrrep)) then
    LenInt = nFck(iIrrep)
    nIc = nIC+1
    ip(NIC) = iStart
    iStart = iStart+LenInt
  end if
end do
LenInt_Tot = iStart-1
call mma_allocate(Integrals,LenInt_Tot,Label='Integrals')
Integrals(:) = Zero

call SOS(iStabO,nStabO,1)

! Auxiliary memory allocation.

!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Basis_Mode('Valence')
call Nr_Shells(nSkal)
call Setup_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *

! Double loop over shells.

do iS=1,nSkal
  iShll = iSD(0,iS)
  iAng = iSD(1,iS)
  iCmp = iSD(2,iS)
  iBas = iSD(3,iS)
  iPrim = iSD(5,iS)
  iAO = iSD(7,iS)
  mdci = iSD(10,iS)
  iShell = iSD(11,iS)
  iCnttp = iSD(13,iS)
  iCnt = iSD(14,iS)
  A(1:3) = dbsc(iCnttp)%Coor(1:3,iCnt)

  do jS=1,iS
    jShll = iSD(0,jS)
    jAng = iSD(1,jS)
    jCmp = iSD(2,jS)
    jBas = iSD(3,jS)
    jPrim = iSD(5,jS)
    jAO = iSD(7,jS)
    mdcj = iSD(10,jS)
    jShell = iSD(11,jS)
    jCnttp = iSD(13,jS)
    jCnt = iSD(14,jS)
    B(1:3) = dbsc(jCnttp)%Coor(1:3,jCnt)

    ! Call kernel routine to get memory requirement. Observe, however
    ! that kernels which will use the HRR will allocate that
    ! memory internally.

    maxi = S%maxPrm(iAng)*S%maxprm(jang)
    call mma_allocate(Zeta,maxi,Label='Zeta')
    call mma_allocate(ZI,maxi,Label='ZI')
    call mma_allocate(Kappa,maxi,Label='Kappa')
    call mma_allocate(PCoor,maxi,3,Label='PCoor')
    call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)

    ! Memory requirements for contraction and symmetry
    ! adaptation of derivatives.

    lFinal = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)*nIrrep

    MemKrn = max(MemKer*Maxi,lFinal)
    call mma_allocate(Kern,MemKrn,Label='Kern')

    ! Save some memory and use Scrt area for transformation

    ! Allocate memory for the final integrals all in the primitive basis.

    call mma_allocate(Fnl,lFinal,Label='Fnl')

    ! Scratch area for the transformation to spherical gaussians

    nScr1 = S%MaxBas(iAng)*S%MaxBas(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)*nIC
    call mma_allocate(ScrSph,nScr1,Label='ScfSph')

    ! At this point we can compute Zeta.
    ! This is now computed in the ij or ji order.

    call ZXia(Zeta,ZI,iPrim,jPrim,Shells(iShll)%Exp,Shells(jShll)%Exp)

    DiffCnt = (mdci == iDCnt) .or. (mdcj == iDCnt)
    if (DiffCnt .or. DiffOp) then
      IfGrd(:,:) = .false.
      trans(:) = .false.
      if (mdci == iDCnt) then
        IfGrd(idCar,1) = .true.
      end if
      if (mdcj == iDCnt) then
        IfGrd(idCar,2) = .true.
      end if

      if (IfGrd(iDCar,1) .and. IfGrd(iDCar,2) .and. (.not. DiffOp)) then
        IfGrd(iDCar,2) = .false.
        Trans(2) = .true.
      end if
      if (Label == 'CONNECTI') Trans(2) = .false.

      ! Allocate memory for SO integrals that will be generated by
      ! this batch of AO integrals.

      nSO = 0
      do iIrrep=0,nIrrep-1
        if (btest(loper,iIrrep)) then
          iSmLbl = 2**iIrrep
          nSO = nSO+MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
        end if
      end do
      !if (iPrint >= 29) write(u6,*) ' nSO=',nSO
      if (nSO /= 0) then
        call mma_allocate(SO,iBas*jBas*nSO,Label='SO')
        SO(:) = Zero

        ! Find the DCR for A and B

        call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)

        ! Find the stabilizer for A and B

        call Inter(dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)

        call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

        ! Compute normalization factor

        iuv = dc(mdci)%nStab*dc(mdcj)%nStab
        Fact = real(iuv*nStabO,kind=wp)/real(nIrrep**2*LmbdT,kind=wp)
        if (MolWgh == 1) then
          Fact = Fact*real(nIrrep,kind=wp)**2/real(iuv,kind=wp)
        else if (MolWgh == 2) then
          Fact = sqrt(real(iuv,kind=wp))*real(nStabO,kind=wp)/real(nIrrep*LmbdT,kind=wp)
        end if

        ! Loops over symmetry operations acting on the basis.

        nOp(1) = NrOpr(0)
        if (jBas < -999999) write(u6,*) 'gcc overoptimization',nDCRR
        do lDCRR=0,nDCRR-1
          call OA(iDCRR(lDCRR),B,RB)
          nOp(2) = NrOpr(iDCRR(lDCRR))
          if ((Label /= 'CONNECTI') .and. EQ(A,RB) .and. (.not. DiffOp)) cycle

          ! Compute kappa and P.

          call Setup1(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,A,RB,Kappa,PCoor,ZI)

          ! Compute AO integrals.
          ! for easy implementation of NA integrals.

          call Kernel(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,Zeta,Kappa,PCoor,Fnl,iPrim*jPrim,iAng,jAng,A,RB,nOrder,Kern, &
                      MemKrn,Ccoor,nOrdOp,IfGrd,IndGrd,nop,dc(mdci)%nStab,dc(mdcj)%nStab,nic,idcar,trans,kcar,isym)

          ! Transform from primitive to contracted basis functions.
          ! Order of transformation is fixed. It has been shown through
          ! testing that the index order ij,ab will give a performance
          ! that is up to 20% faster than the ab,ij index order.

          ! Transform i,jabx to jabx,I

          kk = nTri_Elem1(iAng)*nTri_Elem1(jAng)*nIC
          call DGEMM_('T','N',jPrim*kk,iBas,iPrim,One,Fnl,iPrim,Shells(iShll)%pCff,iPrim,Zero,Kern,jPrim*kk)

          ! Transform j,abxI to abxI,J

          call DGEMM_('T','N',kk*iBas,jBas,jPrim,One,Kern,jPrim,Shells(jShll)%pCff,jPrim,Zero,Fnl,kk*iBas)

          ! Transform to spherical gaussians if needed.

          kk = nTri_Elem1(iAng)*nTri_Elem1(jAng)

          if (Shells(iShll)%Transf .or. Shells(jShll)%Transf) then

            ! Result comes back as IJAB or IJAb

            call CarSph(Fnl,kk,iBas*jBas*nIC,Kern,nScr1,RSph(ipSph(iAng)),iAng,Shells(iShll)%Transf, &
                        Shells(iShll)%Prjct,RSph(ipSph(jAng)),jAng,Shells(jShll)%Transf,Shells(jShll)%Prjct,ScrSph,iCmp*jCmp)

            call DGeTmO(ScrSph,nIC,nIC,iBas*jBas*iCmp*jCmp,Kern,iBas*jBas*iCmp*jCmp)

          else

            ! Transpose abx,IJ back to IJ,abx

            call DGeTmO(Fnl,kk*nIC,kk*nIC,iBas*jBas,Kern,iBas*jBas)
          end if

          ! At this point accumulate the batch of integrals onto the
          ! final symmetry adapted integrals.

          !if (iPrint >= 99) then
          !  call RecPrt (' Accumulated SO integrals, so far...',' ',SO,iBas*jBas,nSO)
          !end if

          ! Symmetry adapt component by component

          iSOBlk = 1
          iIC = 1
          do iIrrep=0,nIrrep-1
            iSmLbl = iand(lOper,2**iIrrep)
            mSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
            if (mSO == 0) then
              do jIrrep=0,nIrrep-1
                if (btest(iSmLbl,jIrrep)) iIC = iIC+1
              end do
            else
              call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,Kern,iBas,jBas,nIC,iIC,SO(iSOBlk),mSO,nOp)
              iSOBlk = iSOBlk+mSO*iBas*jBas
            end if
          end do

        end do

        ! Multiply with factors due to projection operators

        if (Fact /= One) SO(:) = Fact*SO

        ! Scatter the SO's on to the non-zero blocks of the lower triangle.

        iSOBlk = 1
        iIC = 0
        do iIrrep=0,nIrrep-1
          if (btest(lOper,iIrrep)) then
            iSmlbl = 2**iIrrep
            iiC = iiC+1
            mSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
            if ((nfck(iirrep) /= 0) .and. (mSO /= 0)) call SOSctt(SO(iSOBlk),iBas,jBas,mSO,Integrals(ip(iIC)),nFck(iIrrep),iSmLbl, &
                                                                  iCmp,jCmp,iShell,jShell,iAO,jAO,nIC,Label,2**iIrrep,rHrmt)
            iSOBlk = iSOBlk+mSO*iBas*jBas
          end if
        end do

        call mma_deallocate(SO)
      end if
    end if
    call mma_deallocate(pCoor)
    call mma_deallocate(Kappa)
    call mma_deallocate(ZI)
    call mma_deallocate(Zeta)
    call mma_deallocate(ScrSph)
    call mma_deallocate(Fnl)
    call mma_deallocate(Kern)

  end do
end do

call Free_iSD()

! Compute properties or write integrals to disc and deallocate core.

nDens = 0
ndenssq = 0
do iI=0,nIrrep-1
  ndenssq = ndenssq+nbas(ii)**2
  nDens = nDens+nTri_Elem(nBas(iI))
end do
nrOp = 0

call mma_allocate(Scr,ndenssq,Label='Scr')
do iIrrep=0,nIrrep-1
  if (btest(loper,iIrrep)) then
    nrOp = nrOp+1
    jdisp = indgrd(iirrep)
    kOper = 2**iIrrep
    !write(u6,*) koper,isym,jdisp,iirrep
    if (iadd /= 0) then
      irc = -1
      iopt = 0
      call drdmck(irc,iOpt,LabDsk,jdisp,Scr,koper)
      if (irc /= 0) call SysAbendMsg('cnt1el2','error during read in rdmck',' ')
      Integrals(ip(nrop):ip(nrop)+nfck(iIrrep)-1) = Integrals(ip(nrop):ip(nrop)+nfck(iIrrep)-1)+Scr(1:nfck(iIrrep))
    end if
    irc = -1
    iopt = 0
    !write(u6,*) LabDsk,jdisp,koper
    call dwrmck(irc,iOpt,LabDsk,jdisp,Integrals(ip(nrop)),koper)
    if (irc /= 0) call SysAbendMsg('cnt1el2','error during write in dwrmck',' ')
  end if
end do

call mma_deallocate(Scr)
call mma_deallocate(Integrals)

return

end subroutine Cnt1El2
