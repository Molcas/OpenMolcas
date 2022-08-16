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

subroutine Dot1El2(Kernel,KrnlMm,Hess,nGrad,DiffOp,CCoor,FD,nordop)
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
!             Modified for Hermite-Gauss quadrature November '90       *
!             Modified for Rys quadrature November '90                 *
!             Modified for multipole moments November '90              *
!                                                                      *
!             Modified for general kernel routines January '91         *
!             Modified for nonsymmetrical operators February '91       *
!             Modified for gradients October '91                       *
!             Modified for Hessians by AB   Dec '94                    *
!***********************************************************************

use Real_Spherical
use iSD_data
use Basis_Info
use Center_Info
use Symmetry_Info, only: nIrrep, iOper
use Sizes_of_Seward, only: S

implicit real*8(A-H,O-Z)
external Kernel, KrnlMm
#include "Molcas.fh"
#include "real.fh"
#include "stdalloc.fh"
#include "disp.fh"
#include "nsd.fh"
#include "setup.fh"
real*8 A(3), B(3), Ccoor(3), FD(*), RB(3), Hess(nGrad)
integer iDCRR(0:7), iDCRT(0:7), iStabM(0:7), iCoM(0:7,0:7), nOp(2), iStabO(0:7), IndGrd(2,3,3,0:7)
logical AeqB, EQ, DiffOp
real*8, allocatable :: Zeta(:), ZI(:), Kappa(:), PCoor(:,:), Kern(:), Scrt1(:), Scrt2(:), DAO(:), DSOpr(:), DSO(:)
logical, external :: TF
! Statement function
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

call dcopy_(nGrad,[Zero],0,Hess,1)

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

nTasks = nSkal*(nSkal+1)/2
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

  ! Call kernel routine to get memory requirement.

  call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)
  MemKrn = MemKer*S%m2Max
  call mma_allocate(Kern,MemKrn,Label='Kern')

  ! Allocate memory for the final integrals, all in the primitive basis.

  ! Scratch area for contraction step

  nScrt1 = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
  call mma_allocate(Scrt1,nScrt1,Label='Scrt1')

  ! Scratch area for the transformation to spherical gaussians

  nScrt2 = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
  call mma_allocate(Scrt2,nScrt2,Label='Scrt2')

  nDAO = iPrim*jPrim*nElem(iAng)*nElem(jAng)
  call mma_allocate(DAO,nDAO,Label='DAO')

  ! At this point we can compute Zeta.

  call ZXia(Zeta,ZI,iPrim,jPrim,Shells(iShll)%Exp,Shells(jShll)%Exp)

  AeqB = iS == jS

  ! Find the DCR for A and B

  call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
  if ((.not. DiffOp) .and. (nDCRR == 1) .and. EQ(A,B)) Go To 131

  ! Find the stabilizer for A and B

  call Inter(dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)

  ! Generate all possible (left) CoSet
  ! To the stabilizer of A and B

  do i=0,nIrrep-1
    do j=0,nStabM-1
      iCoM(i,j) = ieor(iOper(i),iStabM(j))
    end do
  end do
  !  Order the Coset so we will have the unique ones first
  nMax = 1
  do j=1,nIrrep-1
    ! Check uniqueness
    do i=0,nMax-1
      do ielem=0,nStabM-1
        if (iCoM(i,1) == iCoM(j,ielem)) Go To 435
      end do
    end do
    ! Move unique CoSet
    nMax = nMax+1
    do ielem=0,nStabM-1
      iTmp = iCoM(nMax-1,ielem)
      iCoM(nMax-1,ielem) = iCoM(j,ielem)
      iCoM(j,ielem) = iTmp
    end do
    if (nMax == nIrrep/nStabM) Go To 439
435 continue
  end do
439 continue

  ! Allocate memory for the elements of the Fock or 1st order
  ! denisty matrix which are associated with the current shell pair.

  iSmLbl = 1
  nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
  if (nSO == 0) Go To 131
  call mma_allocate(DSOpr,nSO*iPrim*jPrim,Label='DSOpr')
  DSOpr(:) = Zero
  call mma_allocate(DSO,nSO*iPrim*jPrim,Label='DSO')
  DSO(:) = Zero

  ! Gather the elements from 1st order density / Fock matrix.

  call SOGthr(DSO,iBas,jBas,nSO,FD,n2Tri(iSmLbl),iSmLbl,iCmp,jCmp,iShell,jShell,AeqB,iAO,jAO)

  ! Project the Fock/1st order density matrix in AO
  ! basis on to the primitive basis.

  ! Transform IJ,AB to J,ABi
  call DGEMM_('T','T',jBas*nSO,iPrim,iBas,1.0d0,DSO,iBas,Shells(iShll)%pCff,iPrim,0.0d0,DSOpr,jBas*nSO)
  ! Transform J,ABi to AB,ij
  call DGEMM_('T','T',nSO*iPrim,jPrim,jBas,1.0d0,DSOpr,jBas,Shells(jShll)%pCff,jPrim,0.0d0,DSO,nSO*iPrim)
  ! Transpose to ij,AB
  call DGeTmO(DSO,nSO,nSO,iPrim*jPrim,DSOpr,iPrim*jPrim)
  call mma_deallocate(DSO)

  ! Loops over symmetry operations.

  nOp(1) = NrOpr(0)
  if (jBas < -999999) write(6,*) 'gcc overoptimization',nDCRR
  do lDCRR=0,nDCRR-1
    call OA(iDCRR(lDCRR),B,RB)
    nOp(2) = NrOpr(iDCRR(lDCRR))
    if (EQ(A,RB) .and. (.not. DiffOp)) Go To 140

    lloper = 1
    call SOS(iStabO,nStabO,llOper)
    call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

    ! Compute normalization factor due the DCR symmetrization
    ! of the two basis functions and the operator.

    iuv = dc(mdci)%nStab*dc(mdcj)%nStab
    FactNd = dble(iuv*nStabO)/dble(nIrrep**2*LmbdT)
    if (MolWgh == 1) then
      FactNd = FactNd*dble(nIrrep)**2/dble(iuv)
    else if (MolWgh == 2) then
      FactNd = sqrt(dble(iuv))*nStabO/dble(nIrrep*LmbdT)
    end if

    ! Desymmetrize the matrix with which we will contract the trace.

    call DesymD(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,DAO,iPrim,jPrim,DSOpr,nSO,nOp,FactNd)

    ! Project the spherical harmonic space onto the cartesian space.

    kk = nElem(iAng)*nElem(jAng)
    if (Shells(iShll)%Transf .or. Shells(jShll)%Transf) then

      ! ij,AB --> AB,ij
      call DGeTmO(DAO,iPrim*jPrim,iPrim*jPrim,iCmp*jCmp,Scrt1,iCmp*jCmp)
      ! AB,ij --> ij,ab
      call SphCar(Scrt1,iCmp*jCmp,iPrim*jPrim,Scrt2,nScr2,RSph(ipSph(iAng)),iAng,Shells(iShll)%Transf,Shells(iShll)%Prjct, &
                  RSph(ipSph(jAng)),jAng,Shells(jShll)%Transf,Shells(jShll)%Prjct,DAO,kk)
    end if

    ! Compute kappa and P.

    call Setup1(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,A,RB,Kappa,PCoor,ZI)

    call Icopy(18*nirrep,[0],0,IndGrd,1)
    kk = 0
    do jIrrep=0,nirrep-1
      do Jcar=1,3
        iirrep = irrfnc(2**(jcar-1))
        if (iirrep == jirrep) then
          jj = 0
          do i=0,jirrep-1
            jj = ldisp(i)+jj
          end do
          nDisp = IndDsp(mdci,jIrrep)-jj
          do iCar=1,3
            iComp = 2**(iCar-1)
            if (TF(mdci,jIrrep,iComp)) then
              ndisp = ndisp+1
              IndGrd(1,icar,jcar,jIrrep) = kk+nDisp
            end if
          end do
          kk = kk+ldisp(jirrep)
        end if
      end do
    end do

    kk = 0
    do jIrrep=0,nirrep-1
      do Jcar=1,3
        iirrep = irrfnc(2**(jcar-1))
        if (iirrep == jirrep) then
          jj = 0
          do i=0,jirrep-1
            jj = ldisp(i)+jj
          end do
          nDisp = IndDsp(mdcj,jIrrep)-jj
          do iCar=1,3
            iComp = 2**(iCar-1)
            if (TF(mdcj,jIrrep,iComp)) then
              ndisp = ndisp+1
              IndGrd(2,icar,jcar,jIrrep) = kk+nDisp
            end if
          end do
          kk = kk+ldisp(jirrep)
        end if
      end do
    end do

    ! Compute gradients of the primitive integrals and trace the result.

    call Kernel(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,Zeta,ZI,Kappa,Pcoor,iPrim*jPrim,iAng,jAng,A,RB,nOrder,Kern,MemKer, &
                Ccoor,nOrdOp,Hess,indgrd,DAO,mdci,mdcj,nOp,iStabM,nStabM)

140 continue
  end do

  call mma_deallocate(DSOpr)
131 continue
  call mma_deallocate(DAO)
  call mma_deallocate(Scrt2)
  call mma_deallocate(Scrt1)
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

end subroutine Dot1El2
