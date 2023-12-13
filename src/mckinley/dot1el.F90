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
! Copyright (C) 1990,1991, Roland Lindh                                *
!               1990, IBM                                              *
!               1994, Anders Bernhardsson                              *
!***********************************************************************

subroutine Dot1El(Kernel,KrnlMm,Hess,nHess,DiffOp,CCoor,FD,nFD,lOper,nComp)
!***********************************************************************
!                                                                      *
! Object: to compute gradients of the one electron integrals.          *
!         The memory at this point is assumed to be large enough to do *
!         the computation in core.                                     *
!         The data is structured with respect to four indices, two (my *
!         ny or i j) refer to primitives or basis functions and two (a *
!         b) refer to the components of the cartesian or spherical     *
!         harmonic gaussians.                                          *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January '90                                              *
!             Anders Bernhardsson  Dec '94                             *
!                                                                      *
!             Modified for general kernel routines January '91         *
!             Modified for nonsymmetrical operators February '91       *
!             Modified for gradients October '91                       *
!             Modified for Hermite-Gauss quadrature November '90       *
!             Modified for Rys quadrature November '90                 *
!             Modified for multipole moments November '90              *
!***********************************************************************

use McKinley_global, only: sIrrep
use mck_interface, only: hss_kernel, mck_mem
use Index_Functions, only: iTri, nTri_Elem, nTri_Elem1
use Real_Spherical, only: ipSph, RSph
use iSD_data, only: iSD
use Basis_Info, only: dbsc, MolWgh, Shells
use Center_Info, only: dc
use Symmetry_Info, only: iOper, nIrrep
use Sizes_of_Seward, only: S
use Disp, only: IndDsp
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
procedure(hss_kernel) :: Kernel
procedure(mck_mem) :: KrnlMm
integer(kind=iwp), intent(in) :: nHess, nFD, nComp, lOper(nComp)
real(kind=wp), intent(out) :: Hess(nHess)
logical(kind=iwp), intent(in) :: DiffOp
real(kind=wp), intent(in) :: CCoor(3,nComp), FD(nFD)
integer(kind=iwp) :: i, iAng, iAO, iAtom, iBas, iCar, iCmp, iCnt, iCnttp, iCoM(0:7,0:7), iComp, iComp1, iComp2, iDCRR(0:7), &
                     iDCRT(0:7), ielem, iIrrep, ijS, IndGrd(0:2,0:1,0:7), IndHss(0:1,0:2,0:1,0:2,0:7), iPrim, iS, iShell, iShll, &
                     iSmLbl, iStabM(0:7), iStabO(0:7), iStop, iTmp(0:7), iuv, j, jAng, jAO, jAtom, jBas, jCar, jCmp, jCnt, jCnttp, &
                     jPrim, jS, jShell, jShll, kk, lDCRR, lFinal, llOper, LmbdR, LmbdT, mdci, mdcj, MemKer, MemKrn, nDAO, nDCRR, &
                     nDCRT, nDisp1, nDisp2, nMax, nnIrrep, nOp(2), nOrder, nOrdOp, nScrt1, nScrt2, nSkal, nSO, nStabM, nStabO, &
                     nTasks
real(kind=wp) :: A(3), B(3), FactNd, RB(3)
logical(kind=iwp) :: AeqB, Chck, ifgrd(0:2,0:1), IfHss(0:1,0:2,0:1,0:2)
!character(len=3) :: ChOper(0:7)
real(kind=wp), allocatable :: DAO(:), DSO(:), DSOpr(:), Fnl(:), Kappa(:), Kern(:), PCoor(:,:), Scrt1(:), Scrt2(:), Zeta(:), ZI(:)
integer(kind=iwp), external :: MemSO1, n2Tri, NrOpr
logical(kind=iwp), external :: EQ, TF, TstFnc

Hess(:) = Zero

! Auxiliary memory allocation.

call mma_allocate(Zeta,S%m2Max,Label='Zeta')
call mma_allocate(ZI,S%m2Max,Label='ZI')
call mma_allocate(Kappa,S%m2Max,Label='Kappa')
call mma_allocate(PCoor,S%m2Max,3,Label='PCoor')
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

nTasks = nTri_Elem(nSkal)
iS = 0
jS = 0
do ijS=1,nTasks
  jS = jS+1
  if (jS > iS) then
    iS = jS
    jS = 1
  end if

  !do iS=1,nSkal
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

  !  do jS=1,iS
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
  !write(u6,*) 'iShll,iAng,iCmp,iBas,iPrim,iAO,ixyz,mdci,iShell'
  !write(u6,*) (iSD(i,iS),i=0,11)
  !write(u6,*) 'jShll,jAng,jCmp,jBas,jPrim,jAO,jxyz,mdcj,jShell'
  !write(u6,*) (iSD(i,jS),i=0,11)

  ! Call kernel routine to get memory requirement.

  nOrdOp = 0  ! not used in this implementation
  call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)
  MemKrn = MemKer*S%m2Max
  call mma_allocate(Kern,MemKrn,Label='Kern')

  ! Allocate memory for the final integrals, all in the primitive basis.

  lFinal = 21*S%MaxPrm(iAng)*S%MaxPrm(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)
  call mma_allocate(Fnl,lFinal,Label='Fnl')

  ! Scratch area for contraction step

  nScrt1 = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)
  call mma_allocate(Scrt1,nScrt1,Label='Scrt1')

  ! Scratch area for the transformation to spherical gaussians

  nScrt2 = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)
  call mma_allocate(Scrt2,nScrt2,Label='Scrt2')

  nDAO = iPrim*jPrim*nTri_Elem1(iAng)*nTri_Elem1(jAng)
  call mma_allocate(DAO,nDAO,Label='DAO')

  ! At this point we can compute Zeta.

  call ZXia(Zeta,ZI,iPrim,jPrim,Shells(iShll)%Exp,Shells(jShll)%Exp)

  AeqB = iS == jS

  ! Find the DCR for A and B

  call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
  if (DiffOp .or. (nDCRR /= 1) .or. (.not. EQ(A,B))) then
    !if (iPrint >= 49) write(u6,'(10A)') ' {R}=(',(ChOper(iDCRR(i)),i=0,nDCRR-1),')'

    ! Find the stabilizer for A and B

    call Inter(dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)

    ! Generate all possible (left) CoSet
    ! To the stabilizer of A and B

    do j=0,nStabM-1
      iCoM(0:nIrrep-1,j) = ieor(iOper(0:nIrrep-1),iStabM(j))
    end do
    ! Order the Coset so we will have the unique ones first
    nMax = 1
    loop1: do j=1,nIrrep-1
      ! Check uniqueness
      do i=0,nMax-1
        do ielem=0,nStabM-1
          if (iCoM(i,1) == iCoM(j,ielem)) cycle loop1
        end do
      end do
      ! Move unique CoSet
      nMax = nMax+1
      iTmp(0:nStabM-1) = iCoM(nMax-1,0:nStabM-1)
      iCoM(nMax-1,0:nStabM-1) = iCoM(j,0:nStabM-1)
      iCoM(j,0:nStabM-1) = iTmp(0:nStabM-1)
      if (nMax == nIrrep/nStabM) exit loop1
    end do loop1
    IfHss(:,:,:,:) = .false.
    do iAtom=0,1
      do iCar=0,2
        do jAtom=0,iAtom
          if (iAtom == jAtom) then
            iStop = iCar
          else
            iStop = 2
          end if
          do jCar=0,iStop
            iComp1 = 2**iCar
            iComp2 = 2**jCar
            iComp = ieor(iComp1,iComp2)
            Chck = TstFnc(iCoM,0,iComp,nStabM)
            if (Chck) IfHss(iAtom,iCar,jAtom,jCar) = .true.
          end do
        end do
      end do
    end do
    IndHss(:,:,:,:,0:nirrep-1) = 0
    IndGrd(:,:,0:nirrep-1) = 0

    ! Determine which displacement in all IR's, each center is associated with.

    nnIrrep = nIrrep
    if (sIrrep) nnIrrep = 1

    do iIrrep=0,nnIrrep-1
      nDisp1 = IndDsp(mdci,iIrrep)
      nDisp2 = IndDsp(mdcj,iIrrep)
      do iCar=0,2
        iComp = 2**iCar
        if (TF(mdci,iIrrep,iComp)) then
          nDisp1 = nDisp1+1
          IndGrd(iCar,0,iIrrep) = nDisp1
          if (iIrrep == 0) IfGrd(iCar,0) = .true.
        else
          if (iIrrep == 0) IfGrd(iCar,0) = .false.
          IndGrd(iCar,0,iIrrep) = 0
        end if
        iComp = 2**iCar
        if (TF(mdcj,iIrrep,iComp)) then
          nDisp2 = nDisp2+1
          IndGrd(iCar,1,iIrrep) = nDisp2
          if (iIrrep == 0) IfGrd(iCar,1) = .true.
        else
          IndGrd(iCar,1,iIrrep) = 0
          if (iIrrep == 0) IfGrd(iCar,1) = .true.
        end if
      end do
    end do

    ! Determine index for each 2nd derivative

    do iIrrep=0,nIrrep-1
      do iAtom=0,1
        do iCar=0,2
          do jAtom=0,iAtom
            if (iAtom == jAtom) then
              istop = iCar
            else
              iStop = 2
            end if
            do jCar=0,istop
              if ((IndGrd(iCar,iAtom,iIrrep) > 0) .and. (IndGrd(jCar,jAtom,iIrrep) > 0)) then
                IndHss(iAtom,iCar,jAtom,jCar,iIrrep) = itri(IndGrd(iCar,iAtom,iIrrep),IndGrd(jCar,jAtom,iIrrep))
              else
                IndHss(iAtom,iCar,jAtom,jCar,iIrrep) = 0
              end if
            end do
          end do
        end do
      end do
    end do
    if (.not. DiffOp) then
      iAtom = 1
      do iCar=0,2
        do jAtom=0,1
          if (iAtom == jAtom) then
            iStop = iCar
          else
            iStop = 2
          end if
          do jCar=0,iStop
            if (IfHss(0,iCar,0,jCar) .or. IfHss(0,jCar,0,iCar)) then
              IfHss(iAtom,iCar,jAtom,jCar) = .false.
              IndHss(iAtom,iCar,jAtom,jCar,0:nIrrep-1) = -IndHss(iAtom,iCar,jAtom,jCar,0:nIrrep-1)
            end if
          end do
        end do
      end do

    end if

    ! Allocate memory for the elements of the Fock or 1st order
    ! denisty matrix which are associated with the current shell pair.

    iSmLbl = 1
    nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
    if (nSO /= 0) then
      call mma_allocate(DSOpr,nSO*iPrim*jPrim,Label='DSOpr')
      DSOpr(:) = Zero
      call mma_allocate(DSO,nSO*iPrim*jPrim,Label='DSO')
      DSO(:) = Zero

      ! Gather the elements from 1st order density / Fock matrix.

      call SOGthr(DSO,iBas,jBas,nSO,FD,n2Tri(iSmLbl),iSmLbl,iCmp,jCmp,iShell,jShell,AeqB,iAO,jAO)

      ! Project the Fock/1st order density matrix in AO
      ! basis on to the primitive basis.

      !if (iPrint >= 99) then
      !  call RecPrt(' Left side contraction',' ',Shells(iShll)%pCff,iPrim,iBas)
      !  call RecPrt(' Right side contraction',' ',Shells(jShll)%pCff,jPrim,jBas)
      !end if

      ! Transform IJ,AB to J,ABi
      call DGEMM_('T','T',jBas*nSO,iPrim,iBas,One,DSO,iBas,Shells(iShll)%pCff,iPrim,Zero,DSOpr,jBas*nSO)
      ! Transform J,ABi to AB,ij
      call DGEMM_('T','T',nSO*iPrim,jPrim,jBas,One,DSOpr,jBas,Shells(jShll)%pCff,jPrim,Zero,DSO,nSO*iPrim)
      ! Transpose to ij,AB
      call DGeTmO(DSO,nSO,nSO,iPrim*jPrim,DSOpr,iPrim*jPrim)
      call mma_deallocate(DSO)

      !if (iPrint >= 99) call RecPrt(' Decontracted 1st order density/Fock matrix',' ',DSOpr,iPrim*jPrim,nSO)

      ! Loops over symmetry operations.

      nOp(1) = NrOpr(0)
      if (jBas < -999999) write(u6,*) 'gcc overoptimization',nDCRR
      do lDCRR=0,nDCRR-1
        call OA(iDCRR(lDCRR),B,RB)
        nOp(2) = NrOpr(iDCRR(lDCRR))
        if (EQ(A,RB) .and. (.not. DiffOp)) cycle

        !if (iPrint >= 49) write(u6,'(10A)') ' {M}=(',(ChOper(iStabM(i)),i=0,nStabM-1),')'
        !end if

        llOper = lOper(1)
        do iComp=2,nComp
          llOper = ior(llOper,lOper(iComp))
        end do
        call SOS(iStabO,nStabO,llOper)
        call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

        ! Compute normalization factor due the DCR symmetrization
        ! of the two basis functions and the operator.

        iuv = dc(mdci)%nStab*dc(mdcj)%nStab
        FactNd = real(iuv*nStabO,kind=wp)/real(nIrrep**2*LmbdT,kind=wp)
        if (MolWgh == 1) then
          FactNd = FactNd*real(nIrrep,kind=wp)**2/real(iuv,kind=wp)
        else if (MolWgh == 2) then
          FactNd = sqrt(real(iuv,kind=wp))*real(nStabO,kind=wp)/real(nIrrep*LmbdT,kind=wp)
        end if

        !if (iPrint >= 49) write(u6,'(A,/,2(3F6.2,2X))') ' *** Centers A, RB ***',(A(i),i=1,3),(RB(i),i=1,3)

        ! Desymmetrize the matrix with which we will contract the trace.

        call DesymD(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,DAO,iPrim,jPrim,DSOpr,nSO,nOp,FactNd)

        ! Project the spherical harmonic space onto the cartesian space.

        kk = nTri_Elem1(iAng)*nTri_Elem1(jAng)
        if (Shells(iShll)%Transf .or. Shells(jShll)%Transf) then

          ! ij,AB --> AB,ij
          call DGeTmO(DAO,iPrim*jPrim,iPrim*jPrim,iCmp*jCmp,Scrt1,iCmp*jCmp)
          ! AB,ij --> ij,ab
          call SphCar(Scrt1,iCmp*jCmp,iPrim*jPrim,Scrt2,nScrt2,RSph(ipSph(iAng)),iAng,Shells(iShll)%Transf,Shells(iShll)%Prjct, &
                      RSph(ipSph(jAng)),jAng,Shells(jShll)%Transf,Shells(jShll)%Prjct,DAO,kk)
        end if
        !if (iPrint >= 99) call RecPrt(' Decontracted FD in the cartesian space',' ',DAO,iPrim*jPrim,kk)

        ! Compute kappa and P.

        call Setup1(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,A,RB,Kappa,PCoor,ZI)

        ! Compute gradients of the primitive integrals and trace the result.

        !BS write(u6,*) 'Call the Kernel'

        call Kernel(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,Zeta,ZI,Kappa,Pcoor,Fnl,iPrim*jPrim,iAng,jAng,A,RB,nOrder, &
                    Kern,MemKer*iPrim*jPrim,Ccoor,nOrdOp,Hess,nHess,IfHss,IndHss,ifgrd,indgrd,DAO,mdci,mdcj,nOp,lOper,nComp, &
                    iStabM,nStabM,nIrrep)

#       ifdef _DEBUGPRINT_
        write(u6,*) 'Hess after Kernel call in dot1el '
        call HssPrt(Hess,nHess)
#       endif

      end do

      call mma_deallocate(DSOpr)
    end if
  end if
  call mma_deallocate(DAO)
  call mma_deallocate(Scrt2)
  call mma_deallocate(Scrt1)
  call mma_deallocate(Fnl)
  call mma_deallocate(Kern)

  !  end do
  !end do
end do

call Free_iSD()

call mma_deallocate(PCoor)
call mma_deallocate(Kappa)
call mma_deallocate(ZI)
call mma_deallocate(Zeta)

return

end subroutine Dot1El
