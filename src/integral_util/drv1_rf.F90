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
! Copyright (C) 1990-1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Drv1_RF(FactOp,nOpr,FD,nFD,CCoor,lOper,Cavxyz,lMax)
!***********************************************************************
!                                                                      *
! Object: to compute the local multipole moment, desymmetrize the 1st  *
!         order density matrix and accumulate contributions to the     *
!         global multipole expansion.                                  *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!             January '90                                              *
!             Modified for Hermite-Gauss quadrature November '90       *
!             Modified for Rys quadrature November '90                 *
!             Modified for multipole moments November '90              *
!                                                                      *
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for general kernel routines January  91         *
!             Modified for nonsymmetrical operators February  91       *
!             Modified for gradients October  91                       *
!             Modified for reaction field calculations July  92        *
!             Modified loop structure  99                              *
!***********************************************************************

use Index_Functions, only: nTri_Elem1, nTri3_Elem1
use Real_Spherical, only: ipSph, rSph
use iSD_data, only: iSD
use Basis_Info, only: DBSC, MolWgh, Shells
use Center_Info, only: DC
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use define_af, only: Angtp
use Symmetry_Info, only: ChOper
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nOpr, nFD, lOper(nOpr), lMax
real(kind=wp), intent(in) :: FactOp(nOpr), FD(nFD), CCoor(3,nOpr)
real(kind=wp), intent(inout) :: Cavxyz(nTri3_Elem1(lMax))
integer(kind=iwp) :: iAng, iAO, iBas, iCmp, iCnt, iCnttp, iDCRR(0:7), iDCRT(0:7), iOpr, iPrim, iS, iShell, iShll, iSmLbl, &
                     iStabM(0:7), iStabO(0:7), iuv, jAng, jAO, jBas, jCmp, jCnt, jCnttp, jPrim, jS, jShell, jShll, kk, lDCRR, &
                     lDCRT, lFinal, Lmbdr, Lmbdt, mdci, mdcj, MemKer, MemKrn, nComp, nDAO, nDCRR, nDCRT, nOp(3), nOrder, nOrdOp, &
                     NrOpr, nScr1, nScr2, nSkal, nSO, nStabM, nStabO
real(kind=wp) :: A(3), B(3), C(3), FactND, RB(3), TA(3), TRB(3)
real(kind=wp), allocatable :: DAO(:), DSO(:), DSOpr(:), Fnl(:), Kappa(:), Kern(:), PCoor(:,:), Scr1(:), Scr2(:), Zeta(:), ZI(:)
integer(kind=iwp), external :: MemSO1, n2Tri
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif

! Auxiliary memory allocation.

call mma_allocate(Zeta,S%m2Max,Label='Zeta')
call mma_allocate(ZI,S%m2Max,Label='ZI')
call mma_allocate(Kappa,S%m2Max,Label='Kappa')
call mma_allocate(PCoor,S%m2Max,3,Label='PCoor')
!                                                                      *
!***********************************************************************
!                                                                      *
call Nr_Shells(nSkal)
!                                                                      *
!***********************************************************************
!                                                                      *
! Double loop over shells. These loops decide the integral type

do iS=1,nSkal
  iShll = iSD(0,iS)
  if (Shells(iShll)%Aux) exit
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

    iSmLbl = 1
    nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
    if (nSO == 0) cycle
#   ifdef _DEBUGPRINT_
    write(u6,'(A,A,A,A,A)') ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
#   endif

    ! Call kernel routine to get memory requirement.

    nOrdOp = lMax
    call RFMem(nOrder,MemKer,iAng,jAng,nOrdOp)
    !write(u6,*) nOrder,MemKer,iAng,jAng,nOrdOp
    MemKrn = MemKer*S%m2Max
    call mma_allocate(Kern,MemKrn,Label='Kern')

    ! Allocate memory for the final integrals, all in the
    ! primitive basis.

    nComp = nTri3_Elem1(lMax)
    lFinal = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)*nComp
    call mma_allocate(Fnl,lFinal,Label='Fnl')

    ! Scratch area for contraction step

    nScr1 = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)
    call mma_allocate(Scr1,nScr1,Label='Scr1')

    ! Scratch area for the transformation to spherical gaussians

    nScr2 = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nTri_Elem1(iAng)*nTri_Elem1(jAng)
    call mma_allocate(Scr2,nScr2,Label='Scr2')

    nDAO = iPrim*jPrim*nTri_Elem1(iAng)*nTri_Elem1(jAng)
    call mma_allocate(DAO,nDAO,Label='DAO')

    ! At this point we can compute Zeta.

    call ZXia(Zeta,ZI,iPrim,jPrim,Shells(iShll)%Exp,Shells(jShll)%Exp)

    ! Find the DCR for A and B

    call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
#   ifdef _DEBUGPRINT_
    write(u6,'(10A)') ' {R}=(',(ChOper(iDCRR(i)),i=0,nDCRR-1),')'
#   endif

!-----------Find the stabilizer for A and B

    call Inter(dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)

    ! Allocate memory for the elements of the Fock or 1st order
    ! denisty matrix which are associated with the current shell pair.

    call mma_allocate(DSOpr,nSO*iPrim*jPrim,Label='DSOpr')
    call mma_allocate(DSO,nSO*iPrim*jPrim,Label='DSO')

    ! Gather the elements from 1st order density / Fock matrix.

    call SOGthr(DSO,iBas,jBas,nSO,FD,n2Tri(iSmLbl),iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)

    ! Project the Fock/1st order density matrix in AO
    ! basis on to the primitive basis.

#   ifdef _DEBUGPRINT_
    call RecPrt(' Left side contraction',' ',Shells(iShll)%pCff,iPrim,iBas)
    call RecPrt(' Right side contraction',' ',Shells(jShll)%pCff,jPrim,jBas)
#   endif

    ! Transform IJ,AB to J,ABi
    call DGEMM_('T','T',jBas*nSO,iPrim,iBas,One,DSO,iBas,Shells(iShll)%pCff,iPrim,Zero,DSOpr,jBas*nSO)
    ! Transform J,ABi to AB,ij
    call DGEMM_('T','T',nSO*iPrim,jPrim,jBas,One,DSOpr,jBas,Shells(jShll)%pCff,jPrim,Zero,DSO,nSO*iPrim)
    ! Transpose to ij,AB
    call DGeTmO(DSO,nSO,nSO,iPrim*jPrim,DSOpr,iPrim*jPrim)
    call mma_deallocate(DSO)

#   ifdef _DEBUGPRINT_
    call RecPrt(' Decontracted 1st order density/Fock matrix',' ',DSOpr,iPrim*jPrim,nSO)
#   endif

    ! Loops over symmetry operations.

    do lDCRR=0,nDCRR-1
      call OA(iDCRR(lDCRR),B,RB)

      ! Loop over operators

      do iOpr=1,nOpr
        if (FactOp(iOpr) == Zero) cycle
        C(:) = Ccoor(:,iOpr)

        ! Generate stabilizer of the operator.

        call SOS(iStabO,nStabO,lOper(iOpr))

        ! Find the DCR for M and S

        call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
#       ifdef _DEBUGPRINT_
        write(u6,'(10A)') ' {M}=(',(ChOper(iStabM(i)),i=0,nStabM-1),')'
        write(u6,'(10A)') ' {O}=(',(ChOper(iStabO(i)),i=0,nStabO-1),')'
        write(u6,'(10A)') ' {T}=(',(ChOper(iDCRT(i)),i=0,nDCRT-1),')'
#       endif

        ! Compute normalization factor due the DCR symmetrization
        ! of the two basis functions and the operator.

        iuv = dc(mdci)%nStab*dc(mdcj)%nStab
        FactNd = real(iuv*nStabO,kind=wp)/real(nIrrep**2*LmbdT,kind=wp)
        if (MolWgh == 1) then
          FactNd = FactNd*real(nIrrep,kind=wp)**2/real(iuv,kind=wp)
        else if (MolWgh == 2) then
          FactNd = sqrt(real(iuv,kind=wp))*real(nStabO,kind=wp)/real(nIrrep*LmbdT,kind=wp)
        end if
        FactNd = FactNd*FactOp(iOpr)

        do lDCRT=0,nDCRT-1
          nOp(1) = NrOpr(iDCRT(lDCRT))
          nOp(2) = NrOpr(ieor(iDCRT(lDCRT),iDCRR(lDCRR)))
          nOp(3) = NrOpr(0)

          call OA(iDCRT(lDCRT),A,TA)
          call OA(iDCRT(lDCRT),RB,TRB)
#         ifdef _DEBUGPRINT_
          write(u6,'(A,/,3(3F6.2,2X))') ' *** Centers A, B, C ***',(TA(i),i=1,3),(TRB(i),i=1,3),(C(i),i=1,3)
          write(u6,*) ' nOp=',nOp
#         endif

          ! Desymmetrize the matrix with which we will contract the trace.

          call DesymD(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,DAO,iPrim,jPrim,DSOpr,nSO,nOp,FactNd)

          ! Project the spherical harmonic space onto the
          ! cartesian space.

          kk = nTri_Elem1(iAng)*nTri_Elem1(jAng)
          if (Shells(iShll)%Transf .or. Shells(jShll)%Transf) then

            ! ij,AB --> AB,ij
            call DGeTmO(DAO,iPrim*jPrim,iPrim*jPrim,iCmp*jCmp,Scr1,iCmp*jCmp)
            ! AB,ij --> ij,ab
            call SphCar(Scr1,iCmp*jCmp,iPrim*jPrim,Scr2,nScr2,RSph(ipSph(iAng)),iAng,Shells(iShll)%Transf,Shells(iShll)%Prjct, &
                        RSph(ipSph(jAng)),jAng,Shells(jShll)%Transf,Shells(jShll)%Prjct,DAO,kk)
          end if
#         ifdef _DEBUGPRINT_
          call RecPrt(' Decontracted FD in the cartesian space',' ',DAO,iPrim*jPrim,kk)
#         endif

          ! Compute kappa and P.

          call Setup1(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,TA,TRB,Kappa,PCoor,ZI)

          ! Compute primitive multipole moments.

          call RFInt(Zeta,Kappa,Pcoor,Fnl,iPrim*jPrim,nComp,iAng,jAng,TA,TRB,nOrder,Kern,MemKer,C,nOrdOp)
#         ifdef _DEBUGPRINT_
          call RecPrt(' Final Integrals',' ',Fnl,nDAO,nComp)
#         endif

          ! Trace with 1st order density matrix and accumulate
          ! to the multipole expansion around center Q.

#         ifdef _DEBUGPRINT_
          call RecPrt(' Decontracted FD in the cartesian space',' ',DAO,nDAO,1)
          call RecPrt('Cavxyz',' ',Cavxyz,1,nComp)
#         endif
          call dGeMV_('T',nDAO,nComp,-One,Fnl,nDAO,DAO,1,One,Cavxyz,1)
#         ifdef _DEBUGPRINT_
          call RecPrt('Cavxyz',' ',Cavxyz,1,nComp)
#         endif

        end do
      end do
    end do

    call mma_deallocate(DSOpr)
    call mma_deallocate(DAO)
    call mma_deallocate(Scr2)
    call mma_deallocate(Scr1)
    call mma_deallocate(Fnl)
    call mma_deallocate(Kern)
  end do
end do

call mma_deallocate(PCoor)
call mma_deallocate(Kappa)
call mma_deallocate(ZI)
call mma_deallocate(Zeta)

return

end subroutine Drv1_RF
