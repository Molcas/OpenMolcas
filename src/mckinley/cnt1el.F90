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

subroutine Cnt1El(Kernel,KrnlMm,Label,iDCnt,iDCar,loper,rHrmt,DiffOp,dens,Lab_Dsk,iadd)
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
use mck_interface, only: grd_mck_kernel, mck_mem
use Index_Functions, only: nTri_Elem, nTri_Elem1
use Real_Spherical, only: ipSph, RSph
use iSD_data, only: iSD
use Basis_Info, only: dbsc, MolWgh, nBas, Shells
use Center_Info, only: dc
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
procedure(grd_mck_kernel) :: Kernel
procedure(mck_mem) :: KrnlMm
character(len=8), intent(in) :: Label, Lab_Dsk
integer(kind=iwp), intent(in) :: iDCnt, iDCar, iadd
integer(kind=iwp), intent(out) :: loper
real(kind=wp), intent(in) :: rHrmt, dens(*)
logical(kind=iwp), intent(in) :: DiffOp
#include "Molcas.fh"
#include "print.fh"
#include "disp.fh"
integer(kind=iwp) :: iAng, iAO, iBas, iCar, iCmp, iCnt, iCnttp, iComp, iDCRR(0:7), iDCRT(0:7), iI, iIC, iIrrep, IndGrd(0:7), iopt, &
                     ip(8), iPrim, irc, iS, iShell, iShll, iSmLbl, iSOBlk, iStabM(0:7), iStabO(0:7), iStart, iuv, jAng, jAO, jBas, &
                     jCmp, jCnt, jCnttp, jdisp, jIrrep, jPrim, jS, jShell, jShll, kk, kOper, lDCRR, LenInt, LenInt_Tot, lFinal, &
                     LmbdR, LmbdT, maxi, mdci, mdcj, mDens, MemKer, MemKrn, mSO, nDCRR, nDCRT, nDens, nDenssq, nDisp, nIC, &
                     nnIrrep, nOp(2), nOrder, nOrdOp, nrOp, nScr1, nSkal, nSO, nStabM, nStabO
real(kind=wp) :: A(3), B(3), CCoor(3), Fact, RB(3)
logical(kind=iwp) :: IfGrd(3,2), DiffCnt, Trans(2)
character(len=8) :: LabDsk
real(kind=wp), allocatable :: Fnl(:), Integrals(:), Kappa(:), Kern(:), PCoor(:,:), Scr(:), ScrSph(:), SO(:), Zeta(:), ZI(:)
integer(kind=iwp), external :: MemSO1, NrOpr
real(kind=wp), external :: DDot_
logical(kind=iwp), external :: EQ, TF

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the number of blocks from each component of the operator
! and the irreps it will span.

! Differentiated symmetry-unique center IDCNT
! Derivative wrt component IDCAR=1,2,3 (d/dx,d/dy,d/dz)
! INDDSP(IDCNT,IIRREP) is the number of displacements in
! earlier center/irrep. Thus it is an offset.

LabDsk = Lab_Dsk
nOrdOp = 0
IndGrd(0:nIrrep-1) = 0
loper = 0
#ifdef _DEBUGPRINT_
iprint = 99
#endif
nnIrrep = nIrrep
if (sIrrep) nnIrrep = 1
do iIrrep=0,nnIrrep-1
  nDisp = IndDsp(iDcnt,iIrrep)
  ! First set NDISP=ordering number of this displacement.
  ! Then loop over directions d/dx,d/dy,d/dz
  do iCar=1,3
    iComp = 2**(iCar-1)
    if (TF(iDCnt,iIrrep,iComp)) then
      ndisp = ndisp+1
      ! NDISP is now the ordering number of this displacement.
      if (iDCar == icar) then
        loper = loper+2**iIrrep
        IndGrd(iIrrep) = nDisp
      end if
    end if
  end do
end do
nIC = 0
if (loper == 0) return
! For the displacement represented by this symmetry-unique
! center IDCNT and this component IDCAR, the differentiation
! operator has components with irreps that have been marked
! with '1' in LOPER, regarded as a flag array.
! INDGRD(IIRREP) will be zero, except for those irreps, and
! will then contain the ordering number of the displacement.

! Allocate one integral array for each of these irreps.
! The address is kept in array IP().
nIC = 0
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

! Obtain ISTABO, the stabilizer of the totally symmetric irrep(!)
! Note: 3rd parameter is bit-packed set of irreps
! so '1' contains only irrep nr 0.
! But then ISTABO will be the whole group!? and NSTABO=NIRREP?!
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

    ! Memory requirements for contraction and Symmetry
    ! adaption of derivatives.

    lFinal = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)*nIrrep

    MemKrn = max(MemKer*Maxi,lFinal)
    call mma_Allocate(Kern,MemKrn,Label='Kern')

    ! Save some memory and use Scrt area for transformation

    ! Allocate memory for the final integrals all in the primitive basis.

    call mma_allocate(Fnl,lFinal,Label='Fnl')

    ! Scratch area for the transformation to spherical gaussians

    nScr1 = S%MaxBas(iAng)*S%MaxBas(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)*nIC
    call mma_allocate(ScrSph,nScr1,Label='ScrSph')

    ! At this point we can compute Zeta.
    ! This is now computed in the ij or ji order.

    call ZXia(Zeta,ZI,iPrim,jPrim,Shells(iShll)%Exp,Shells(jShll)%Exp)

    DiffCnt = (mdci == iDCnt) .or. (mdcj == iDCnt)
    if (DiffCnt .or. DiffOp) then
      IfGrd(:,:) = .false.
      ! trans(iCnt) is true means there will be a sign shift in the SYMADO
      ! routine for the contribution to the integral from the
      ! differentiation wrt center iCnt
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
      if (Label == 'OVRGRDA') Trans(2) = .false.

      ! Allocate memory for SO integrals that will be generated by
      ! this batch of AO integrals.

      nSO = 0
      do iIrrep=0,nIrrep-1
        if (btest(loper,iIrrep)) then
          iSmLbl = 2**iIrrep
          nSO = nSO+MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
        end if
      end do
#     ifdef _DEBUGPRINT_
      if (iPrint >= 29) write(u6,*) ' nSO=',nSO
#     endif
      if (nSO /= 0) then
        call mma_Allocate(SO,nSO*iBas*jBas,Label='SO')
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

          Fnl(:) = Zero
          call Kernel(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,Zeta,ZI,Kappa,PCoor,Fnl,iPrim*jPrim,iAng,jAng,A,RB,nOrder, &
                      Kern,MemKrn,Ccoor,nOrdOp,IfGrd,IndGrd,nop,loper,dc(mdci)%nStab,dc(mdcj)%nStab,nic,idcar,idcnt,iStabM,nStabM, &
                      trans,nIrrep)

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

            call CarSph(Fnl,kk,iBas*jBas*nIC,Kern,nScr1,RSph(ipSph(iAng)),iAng,Shells(iShll)%Transf,Shells(iShll)%Prjct, &
                        RSph(ipSph(jAng)),jAng,Shells(jShll)%Transf,Shells(jShll)%Prjct,ScrSph,iCmp*jCmp)

            call DGeTmO(ScrSph,nIC,nIC,iBas*jBas*iCmp*jCmp,Kern,iBas*jBas*iCmp*jCmp)

          else

            ! Transpose abx,IJ back to IJ,abx

            call DGeTmO(Fnl,kk*nIC,kk*nIC,iBas*jBas,Kern,iBas*jBas)
          end if

          ! At this point accumulate the batch of integrals onto the
          ! final symmetry adapted integrals.

#         ifdef _DEBUGPRINT_
          if (iPrint >= 99) then
            call RecPrt(' Accumulated SO integrals, so far...',' ',SO,iBas*jBas,nSO)
          end if
#         endif

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
        iiC = 0
        do iIrrep=0,nIrrep-1
          if (btest(lOper,iIrrep)) then
            iSmlbl = 2**iIrrep
            iiC = iiC+1
            mSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
            if ((nfck(iIrrep) /= 0) .and. (mSO /= 0)) call SOSctt(SO(iSOBlk),iBas,jBas,mSO,Integrals(ip(iIC)),nFck(iIrrep),iSmLbl, &
                                                                  iCmp,jCmp,iShell,jShell,iAO,jAO,nIC,Label,2**iIrrep,rHrmt)
            iSOBlk = iSOBlk+mSO*iBas*jBas
          end if
        end do

        call mma_deallocate(SO)
      end if
    end if
    call mma_deallocate(ScrSph)
    call mma_deallocate(Kern)
    call mma_deallocate(Fnl)
    call mma_deallocate(PCoor)
    call mma_deallocate(Kappa)
    call mma_deallocate(ZI)
    call mma_deallocate(Zeta)
  end do
end do
call Free_iSD()

! Compute properties or write integrals to disc and deallocate core.

nDens = 0
nDenssq = 0
do iI=0,nIrrep-1
  nDenssq = nDenssq+nBas(ii)**2+nBas(ii)
  nDens = nDens+nTri_Elem(nBas(iI))
end do
nrOp = 0

call mma_allocate(Scr,2*nDenssq,Label='Scr')
do iIrrep=0,nIrrep-1
  if (btest(loper,iIrrep)) then
    nrOp = nrOp+1
    jdisp = indgrd(iIrrep)
    kOper = 2**iIrrep
    if (show .and. (iIrrep == 0)) then
      write(u6,*) Label,': ',ddot_(nDens,Dens,1,Integrals(ip(nrop)),1)
      write(u6,*) 'oper: ',ddot_(nDens,Integrals(ip(nrop)),1,Integrals(ip(nrop)),1)
      write(u6,*) 'Dens: ',ddot_(nDens,Dens,1,Dens,1)
    else if (show) then
      mDens = nFck(iIrrep)
      write(u6,*) Label
      write(u6,'(A,G20.10)') 'oper: ',ddot_(mDens,Integrals(ip(nrop)),1,Integrals(ip(nrop)),1)
    end if

    if (iadd /= 0) then
      irc = -1
      iopt = 0
      call dRdMck(irc,iOpt,LabDsk,jdisp,Scr,koper)
      if (irc /= 0) call SysAbendMsg('cnt1el','error during read in rdmck',' ')
      Integrals(ip(nrop):ip(nrop)+nfck(iIrrep)-1) = Integrals(ip(nrop):ip(nrop)+nfck(iIrrep)-1)+Scr(1:nfck(iIrrep))
    end if
    irc = -1
    iopt = 0
#   ifdef _DEBUGPRINT_
    write(u6,'(2A,2I8)') 'LabDsk,jdisp,koper',LabDsk,jdisp,koper
#   endif
    call dWrMck(irc,iOpt,LabDsk,jdisp,Integrals(ip(nrop)),koper)
    if (irc /= 0) call SysAbendMsg('cnt1el','error during write in dwrmck',' ')
  end if
end do
call mma_deallocate(Scr)
call mma_deallocate(Integrals)

return

end subroutine Cnt1El
