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
! Copyright (C) 1990-1992,1999, Roland Lindh                           *
!***********************************************************************

subroutine drv_ef_PCM(FactOp,nTs,FD,nFD,CCoor,lOper,VTessera,nOrdOp)
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

use Real_Spherical, only: ipSph, RSph
use iSD_data, only: iSD
use Basis_Info, only: dbsc, MolWgh, Shells
use Center_Info, only: dc
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use Index_Functions, only: nTri_Elem1
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nTs, nFD, lOper(nTs), nOrdOp
real(kind=wp), intent(in) :: FactOp(nTs), FD(nFD), CCoor(4,nTs)
real(kind=wp), intent(inout) :: VTessera(3,nTs)
#include "angtp.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"
integer(kind=iwp) :: i, iAng, iAO, iBas, iCmp, iCnt, iCnttp, iComp, iDCRR(0:7), iDCRT(0:7), ipFnlc, iPrim, iPrint, iRout, iS, &
                     iShell, iShll, iSmLbl, iStabM(0:7), iStabO(0:7), iTile, iuv, jAng, jAO, jBas, jCmp, jCnt, jCnttp, jPrim, jS, &
                     jShell, jShll, kk, lDCRR, lDCRT, lFinal, LmbdR, LmbdT, mdci, mdcj, MemKer, MemKrn, nComp, nDAO, nDCRR, nDCRT, &
                     niAng, njAng, nOp(3), nOrder, nScr1, nScr2, nSkal, nSO, nStabM, nStabO
real(kind=wp) :: A(3), B(3), C(3), FactNd, RB(3), TA(3), TRB(3)
logical(kind=iwp) :: AeqB
real(kind=wp), allocatable :: DAO(:), DSO(:), DSOp(:), Fnl(:), Kappa(:), Kern(:), PCoor(:), Scrt1(:), Scrt2(:), Zeta(:), ZI(:)
character(len=3), parameter :: ChOper(0:7) = ['E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz']
integer(kind=iwp), external :: MemSO1, n2Tri, NrOpr
real(kind=wp), external :: DDot_

iRout = 112
iPrint = nPrint(iRout)

! Auxiliary memory allocation.

call mma_allocate(Zeta,S%m2Max,label='Zeta')
call mma_allocate(ZI,S%m2Max,label='ZI')
call mma_allocate(Kappa,S%m2Max,label='Kappa')
call mma_allocate(PCoor,3*S%m2Max,label='PCoor')
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
  niAng = nTri_Elem1(iAng)
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
    njAng = nTri_Elem1(jAng)

    iSmLbl = 1
    nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
    if (nSO == 0) cycle
    if (iPrint >= 19) write(u6,'(A,A,A,A,A)') ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'

    ! Call kernel routine to get memory requirement.

    call EFMmP(nOrder,MemKer,iAng,jAng,nOrdOp)
    !write(u6,*) nOrder,MemKer,iAng,jAng,nOrdOp
    MemKrn = MemKer*S%m2Max
    call mma_allocate(Kern,MemKrn,label='Kernel')

    ! Allocate memory for the final integrals, all in the
    ! primitive basis.

    nComp = (nOrdOp+1)*(nOrdOp+2)/2
    lFinal = S%MaxPrm(iAng)*S%MaxPrm(jAng)*niAng*njAng*nComp
    call mma_allocate(Fnl,lFinal,label='Final')

    ! Scratch area for contraction step

    nScr1 = S%MaxPrm(iAng)*S%MaxPrm(jAng)*niAng*njAng
    call mma_allocate(Scrt1,nScr1,label='Scrtch')

    ! Scratch area for the transformation to spherical gaussians

    nScr2 = S%MaxPrm(iAng)*S%MaxPrm(jAng)*niAng*njAng
    call mma_allocate(Scrt2,nScr2,label='ScrSph')

    nDAO = iPrim*jPrim*niAng*njAng
    call mma_allocate(DAO,nDAO,label='DAO')

    ! At this point we can compute Zeta.

    call ZXia(Zeta,ZI,iPrim,jPrim,Shells(iShll)%Exp,Shells(jShll)%Exp)

    AeqB = iS == jS

    ! Find the DCR for A and B

    call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
    if (iPrint >= 49) write(u6,'(10A)') ' {R}=(',(ChOper(iDCRR(i)),i=0,nDCRR-1),')'

    ! Find the stabilizer for A and B

    call Inter(dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)

    ! Allocate memory for the elements of the Fock or 1st order
    ! denisty matrix which are associated with the current shell
    ! pair.

    call mma_allocate(DSOp,nSO*iPrim*jPrim,label='DSOpr')
    call mma_allocate(DSO,nSO*iPrim*jPrim,label='DSO')

    ! Gather the elements from 1st order density / Fock matrix.

    call SOGthr(DSO,iBas,jBas,nSO,FD,n2Tri(iSmLbl),iSmLbl,iCmp,jCmp,iShell,jShell,AeqB,iAO,jAO)

    ! Project the Fock/1st order density matrix in AO
    ! basis on to the primitive basis.

    if (iPrint >= 99) then
      call RecPrt(' Left side contraction',' ',Shells(iShll)%pCff,iPrim,iBas)
      call RecPrt(' Right side contraction',' ',Shells(jShll)%pCff,jPrim,jBas)
    end if

    ! Transform IJ,AB to J,ABi
    call DGEMM_('T','T',jBas*nSO,iPrim,iBas,One,DSO,iBas,Shells(iShll)%pCff,iPrim,Zero,DSOp,jBas*nSO)
    ! Transform J,ABi to AB,ij
    call DGEMM_('T','T',nSO*iPrim,jPrim,jBas,One,DSOp,jBas,Shells(jShll)%pCff,jPrim,Zero,DSO,nSO*iPrim)
    ! Transpose to ij,AB
    call DGeTmO(DSO,nSO,nSO,iPrim*jPrim,DSOp,iPrim*jPrim)
    call mma_deallocate(DSO)

    if (iPrint >= 99) call RecPrt(' Decontracted 1st order density/Fock matrix',' ',DSOp,iPrim*jPrim,nSO)

    ! Loops over symmetry operations.

    do lDCRR=0,nDCRR-1
      call OA(iDCRR(lDCRR),B,RB)

      ! Loop over operators

      do iTile=1,nTs
        if (FactOp(iTile) == Zero) cycle
        C = Ccoor(1:3,iTile)

        ! Generate stabilizer of the operator.

        call SOS(iStabO,nStabO,lOper(iTile))

        ! Find the DCR for M and S

        call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
        if (iPrint >= 49) then
          write(u6,'(10A)') ' {M}=(',(ChOper(iStabM(i)),i=0,nStabM-1),')'
          write(u6,'(10A)') ' {O}=(',(ChOper(iStabO(i)),i=0,nStabO-1),')'
          write(u6,'(10A)') ' {T}=(',(ChOper(iDCRT(i)),i=0,nDCRT-1),')'
        end if

        ! Compute normalization factor due the DCR symmetrization
        ! of the two basis functions and the operator.

        iuv = dc(mdci)%nStab*dc(mdcj)%nStab
        FactNd = real(iuv*nStabO,kind=wp)/real(nIrrep**2*LmbdT,kind=wp)
        if (MolWgh == 1) then
          FactNd = FactNd*real(nIrrep,kind=wp)**2/real(iuv,kind=wp)
        else if (MolWgh == 2) then
          FactNd = sqrt(real(iuv,kind=wp))*real(nStabO,kind=wp)/real(nIrrep*LmbdT,kind=wp)
        end if
        FactNd = FactNd*FactOp(iTile)

        do lDCRT=0,nDCRT-1
          nOp(1) = NrOpr(iDCRT(lDCRT))
          nOp(2) = NrOpr(ieor(iDCRT(lDCRT),iDCRR(lDCRR)))
          nOp(3) = NrOpr(0)

          call OA(iDCRT(lDCRT),A,TA)
          call OA(iDCRT(lDCRT),RB,TRB)
          if (iPrint >= 49) then
            write(u6,'(A,/,3(3F6.2,2X))') ' *** Centers A, B, C ***',(TA(i),i=1,3),(TRB(i),i=1,3),(C(i),i=1,3)
            write(u6,*) ' nOp=',nOp
          end if

          ! Desymmetrize the matrix with which we will contract the trace.

          call DesymD(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,DAO,iPrim,jPrim,DSOp,nSO,nOp,FactNd)

          ! Project the spherical harmonic space onto the cartesian space.

          kk = niAng*njAng
          if (Shells(iShll)%Transf .or. Shells(jShll)%Transf) then

            ! ij,AB --> AB,ij
            call DGeTmO(DAO,iPrim*jPrim,iPrim*jPrim,iCmp*jCmp,Scrt1,iCmp*jCmp)
            ! AB,ij --> ij,ab
            call SphCar(Scrt1,iCmp*jCmp,iPrim*jPrim,Scrt2,nScr2,RSph(ipSph(iAng)),iAng,Shells(iShll)%Transf,Shells(iShll)%Prjct, &
                        RSph(ipSph(jAng)),jAng,Shells(jShll)%Transf,Shells(jShll)%Prjct,DAO,kk)
          end if
          if (iPrint >= 99) call RecPrt(' Decontracted FD in the cartesian space',' ',DAO,iPrim*jPrim,kk)

          ! Compute kappa and P.

          call Setup1(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,TA,TRB,Kappa,PCoor,ZI)

          ! Compute the potential at a tessera.

          ! pcm_solvent
          !write(u6,*) 'Shells(iShll)%Exp(1),iPrim,Shells(jShll)%Exp(1),jPrim'
          !write(u6,*) Shells(iShll)%Exp(1),iPrim,Shells(jShll)%Exp(1),jPrim
          !write(u6,*) 'Zeta,ZI,Kappa,Pcoor'
          !write(u6,*) Zeta(1),ZI(1),Kappa(1),Pcoor(:)
          !write(u6,*) 'Fnl,iPrim*jPrim,nComp,iAng,jAng,norder'
          !write(u6,*) Fnl(1),iPrim*jPrim,nComp,iAng,jAng,norder
          ! pcm_solvent end
          call EFPrm(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,Zeta,ZI,Kappa,Pcoor,Fnl,iPrim*jPrim,nComp,iAng,jAng,TA,TRB, &
                     nOrder,Kern,MemKer,C,nOrdOp)
          if (iPrint >= 49) call RecPrt(' Final Integrals',' ',Fnl,nDAO,nComp)

          ! Trace with 1st order density matrix and accumulate
          ! to the potenital at tessera iTile

          if (iPrint >= 49) call RecPrt(' Decontracted FD in the cartesian space',' ',DAO,nDAO,1)
          ipFnlc = 1
          do iComp=1,nComp
            if (iPrint >= 49) call RecPrt('VTessera(iComp,iTile)',' ',VTessera(iComp,iTile),1,1)
            VTessera(iComp,iTile) = VTessera(iComp,iTile)+DDot_(nDAO,DAO,1,Fnl(ipFnlc),1)
            if (iPrint >= 49) call RecPrt('VTessera(iComp,iTile)',' ',VTessera(iComp,iTile),1,1)
            ipFnlc = ipFnlc+nDAO
          end do

        end do
      end do
    end do

    call mma_deallocate(Kern)
    call mma_deallocate(Fnl)
    call mma_deallocate(Scrt1)
    call mma_deallocate(Scrt1)
    call mma_deallocate(DAO)
    call mma_deallocate(DSOp)
  end do
end do

call mma_deallocate(Zeta)
call mma_deallocate(ZI)
call mma_deallocate(Kappa)
call mma_deallocate(PCoor)

return

end subroutine drv_ef_PCM
