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

subroutine V_EF_PCM(nAt,nTs,DoPot,DoFld,AtmC,Tessera,V,EF_n,EF_e)

implicit real*8(a-h,o-z)
dimension AtmC(3,nAt)
dimension V(*), EF_n(3,*), EF_e(3,*), Tessera(4,*)
logical DoPot, DoFld

! Compute potential on tesserae

if (DoPot) then
  call FZero(V,nTs)
  nOrdOp = 0
  call Mlt_PCM(nAt,nTs,nOrdOp,Tessera,AtmC,V,EF_n,EF_e)
end if

! Compute electric field on tesserae

if (DoFld) then
  call FZero(EF_n,3*nTs)
  call FZero(EF_e,3*nTs)
  nOrdOp = 1
  call Mlt_PCM(nAt,nTs,nOrdOp,Tessera,AtmC,V,EF_n,EF_e)
end if

return

end subroutine V_EF_PCM
!====
subroutine Mlt_PCM(nAt,nTs,nOrdOp,Tessera,AtmC,V,EF_n,EF_e)

implicit real*8(A-H,O-Z)
dimension V(*), EF_n(3,*), EF_e(3,*)
dimension Tessera(4,*)
dimension AtmC(3,*)
#include "real.fh"
#include "print.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
real*8 Temp(3)
real*8, allocatable :: Chrg(:)
logical Found
real*8, allocatable :: D1ao(:)

call mma_allocate(Chrg,nat)
call Get_dArray('Nuclear charge',Chrg,nAt)
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute the desired multipole on the tiles

! 1) The nuclear contribution

do iTs=1,nTs
  call EFNuc(Tessera(1,iTs),Chrg,AtmC,nAt,Temp,nOrdOp)
  if (nOrdOp == 0) then
    V(iTs) = Temp(1)
  elseif (nOrdOp == 1) then
    EF_n(1,iTs) = Temp(1)
    EF_n(2,iTs) = Temp(2)
    EF_n(3,iTs) = Temp(3)
  end if
end do

call mma_deallocate(Chrg)

! 2) The electronic contribution

! Get the total 1st order AO density matrix

call Qpg_dArray('D1ao',Found,nDens)
if (Found .and. (nDens /= 0)) then
  call mma_allocate(D1ao,nDens,Label='D1ao')
else
  write(6,*) 'Mlt_pcm: D1ao not found.'
  call Abend()
end if
call Get_D1ao(D1ao,nDens)

call Allocate_Work(ipFactOp,nTs)
call Allocate_iWork(iplOper,nTs)
call dcopy_(nTs,[One],0,Work(ipFactOp),1)
call ICopy(nTs,[255],0,iWork(iplOper),1)

call drv_ef_PCM(Work(ipFactOp),nTs,D1ao,nDens,Tessera,iWork(iplOper),EF_e,nOrdOp)
if (nOrdOp == 0) then
  do iTs=1,nTs
    V(iTs) = Ef_e(1,iTs)
  end do
end if

call Free_iWork(iplOper)
call Free_Work(ipFactOp)
call mma_deallocate(D1ao)

return

end subroutine Mlt_PCM
!====
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

use Real_Spherical
use iSD_data
use Basis_Info
use Center_Info
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
#include "angtp.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"
real*8 A(3), B(3), C(3), FD(nFD), FactOp(nTs), CCoor(4,nTs), RB(3), TRB(3), TA(3), VTessera(3,nTs)
character ChOper(0:7)*3
integer lOper(nTs), iStabO(0:7), iDCRR(0:7), iDCRT(0:7), iStabM(0:7), nOp(3)
logical AeqB
data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/
! Statement functions
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

iRout = 112
iPrint = nPrint(iRout)

! Auxiliary memory allocation.

call GetMem('Zeta','ALLO','REAL',iZeta,S%m2Max)
call GetMem('Zeta','ALLO','REAL',ipZI,S%m2Max)
call GetMem('Kappa','ALLO','REAL',iKappa,S%m2Max)
call GetMem('PCoor','ALLO','REAL',iPCoor,S%m2Max*3)
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
    if (nSO == 0) Go To 131
    if (iPrint >= 19) write(6,'(A,A,A,A,A)') ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'

    ! Call kernel routine to get memory requirement.

    call EFMmP(nOrder,MemKer,iAng,jAng,nOrdOp)
    !write(6,*) nOrder,MemKer,iAng,jAng,nOrdOp
    MemKrn = MemKer*S%m2Max
    call GetMem('Kernel','ALLO','REAL',iKern,MemKrn)

    ! Allocate memory for the final integrals, all in the
    ! primitive basis.

    nComp = (nOrdOp+1)*(nOrdOp+2)/2
    lFinal = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nElem(iAng)*nElem(jAng)*nComp
    call GetMem('Final','ALLO','REAL',ipFnl,lFinal)

    ! Scratch area for contraction step

    nScr1 = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
    call GetMem('Scrtch','ALLO','REAL',iScrt1,nScr1)

    ! Scratch area for the transformation to spherical gaussians

    nScr2 = S%MaxPrm(iAng)*S%MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
    call GetMem('ScrSph','Allo','Real',iScrt2,nScr2)

    nDAO = iPrim*jPrim*nElem(iAng)*nElem(jAng)
    call GetMem(' DAO ','Allo','Real',ipDAO,nDAO)

    ! At this point we can compute Zeta.

    call ZXia(Work(iZeta),Work(ipZI),iPrim,jPrim,Shells(iShll)%Exp,Shells(jShll)%Exp)

    AeqB = iS == jS

    ! Find the DCR for A and B

    call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
    if (iPrint >= 49) write(6,'(10A)') ' {R}=(',(ChOper(iDCRR(i)),i=0,nDCRR-1),')'

    ! Find the stabilizer for A and B

    call Inter(dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)

    ! Allocate memory for the elements of the Fock or 1st order
    ! denisty matrix which are associated with the current shell
    ! pair.

    call GetMem('DSOpr ','ALLO','REAL',ipDSOp,nSO*iPrim*jPrim)
    call GetMem('DSO ','ALLO','REAL',ipDSO,nSO*iPrim*jPrim)

    ! Gather the elements from 1st order density / Fock matrix.

    call SOGthr(Work(ipDSO),iBas,jBas,nSO,FD,n2Tri(iSmLbl),iSmLbl,iCmp,jCmp,iShell,jShell,AeqB,iAO,jAO)

    ! Project the Fock/1st order density matrix in AO
    ! basis on to the primitive basis.

    if (iPrint >= 99) then
      call RecPrt(' Left side contraction',' ',Shells(iShll)%pCff,iPrim,iBas)
      call RecPrt(' Right side contraction',' ',Shells(jShll)%pCff,jPrim,jBas)
    end if

    ! Transform IJ,AB to J,ABi
    call DGEMM_('T','T',jBas*nSO,iPrim,iBas,1.0d0,Work(ipDSO),iBas,Shells(iShll)%pCff,iPrim,0.0d0,Work(ipDSOp),jBas*nSO)
    ! Transform J,ABi to AB,ij
    call DGEMM_('T','T',nSO*iPrim,jPrim,jBas,1.0d0,Work(ipDSOp),jBas,Shells(jShll)%pCff,jPrim,0.0d0,Work(ipDSO),nSO*iPrim)
    ! Transpose to ij,AB
    call DGeTmO(Work(ipDSO),nSO,nSO,iPrim*jPrim,Work(ipDSOp),iPrim*jPrim)
    call GetMem('DSO ','Free','Real',ipDSO,nSO*iBas*jBas)

    if (iPrint >= 99) call RecPrt(' Decontracted 1st order density/Fock matrix',' ',Work(ipDSOp),iPrim*jPrim,nSO)

    ! Loops over symmetry operations.

    do lDCRR=0,nDCRR-1
      call OA(iDCRR(lDCRR),B,RB)

      ! Loop over operators

      do iTile=1,nTs
        if (FactOp(iTile) == Zero) Go To 5000
        call dcopy_(3,Ccoor(1,iTile),1,C,1)

        ! Generate stabilizer of the operator.

        call SOS(iStabO,nStabO,lOper(iTile))

        ! Find the DCR for M and S

        call DCR(LmbdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)
        if (iPrint >= 49) then
          write(6,'(10A)') ' {M}=(',(ChOper(iStabM(i)),i=0,nStabM-1),')'
          write(6,'(10A)') ' {O}=(',(ChOper(iStabO(i)),i=0,nStabO-1),')'
          write(6,'(10A)') ' {T}=(',(ChOper(iDCRT(i)),i=0,nDCRT-1),')'
        end if

        ! Compute normalization factor due the DCR symmetrization
        ! of the two basis functions and the operator.

        iuv = dc(mdci)%nStab*dc(mdcj)%nStab
        FactNd = dble(iuv*nStabO)/dble(nIrrep**2*LmbdT)
        if (MolWgh == 1) then
          FactNd = FactNd*dble(nIrrep)**2/dble(iuv)
        else if (MolWgh == 2) then
          FactNd = sqrt(dble(iuv))*dble(nStabO)/dble(nIrrep*LmbdT)
        end if
        FactNd = FactNd*FactOp(iTile)

        do lDCRT=0,nDCRT-1
          nOp(1) = NrOpr(iDCRT(lDCRT))
          nOp(2) = NrOpr(ieor(iDCRT(lDCRT),iDCRR(lDCRR)))
          nOp(3) = NrOpr(0)

          call OA(iDCRT(lDCRT),A,TA)
          call OA(iDCRT(lDCRT),RB,TRB)
          if (iPrint >= 49) then
            write(6,'(A,/,3(3F6.2,2X))') ' *** Centers A, B, C ***',(TA(i),i=1,3),(TRB(i),i=1,3),(C(i),i=1,3)
            write(6,*) ' nOp=',nOp
          end if

          ! Desymmetrize the matrix with which we will contract the trace.

          call DesymD(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,Work(ipDAO),iPrim,jPrim,Work(ipDSOp),nSO,nOp, &
                      FactNd)

          ! Project the spherical harmonic space onto the cartesian space.

          kk = nElem(iAng)*nElem(jAng)
          if (Shells(iShll)%Transf .or. Shells(jShll)%Transf) then

            ! ij,AB --> AB,ij
            call DGeTmO(Work(ipDAO),iPrim*jPrim,iPrim*jPrim,iCmp*jCmp,Work(iScrt1),iCmp*jCmp)
            ! AB,ij --> ij,ab
            call SphCar(Work(iScrt1),iCmp*jCmp,iPrim*jPrim,Work(iScrt2),nScr2,RSph(ipSph(iAng)),iAng,Shells(iShll)%Transf, &
                        Shells(iShll)%Prjct,RSph(ipSph(jAng)),jAng,Shells(jShll)%Transf,Shells(jShll)%Prjct,Work(ipDAO),kk)
          end if
          if (iPrint >= 99) call RecPrt(' Decontracted FD in the cartesian space',' ',Work(ipDAO),iPrim*jPrim,kk)

          ! Compute kappa and P.

          call Setup1(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,TA,TRB,Work(iKappa),Work(iPCoor),Work(ipZI))

          ! Compute the potential at a tessera.

          ! pcm_solvent
          !write(6,*) 'Work(iExp),iPrim,Work(jExp),jPrim'
          !write(6,*) Shells(iShll)%Exp(1),iPrim,Shells(jShll)%Exp(1),jPrim
          !write(6,*) 'Work(iZeta),Work(ipZI),Work(iKappa),Work(iPcoor)'
          !write(6,*) Work(iZeta),Work(ipZI),Work(iKappa),Work(iPcoor)
          !write(6,*) 'Work(ipFnl),iPrim*jPrim,nComp,iAng,jAng,norder'
          !write(6,*) Work(ipFnl),iPrim*jPrim,nComp,iAng,jAng,norder
          ! pcm_solvent end
          call EFPrm(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,Work(iZeta),Work(ipZI),Work(iKappa),Work(iPcoor),Work(ipFnl), &
                     iPrim*jPrim,nComp,iAng,jAng,TA,TRB,nOrder,Work(iKern),MemKer,C,nOrdOp)
          if (iPrint >= 49) call RecPrt(' Final Integrals',' ',Work(ipFnl),nDAO,nComp)

          ! Trace with 1st order density matrix and accumulate
          ! to the potenital at tessera iTile

          if (iPrint >= 49) call RecPrt(' Decontracted FD in the cartesian space',' ',Work(ipDAO),nDAO,1)
          ipFnlc = ipFnl
          do iComp=1,nComp
            if (iPrint >= 49) call RecPrt('VTessera(iComp,iTile)',' ',VTessera(iComp,iTile),1,1)
            VTessera(iComp,iTile) = VTessera(iComp,iTile)+DDot_(nDAO,Work(ipDAO),1,Work(ipFnlc),1)
            if (iPrint >= 49) call RecPrt('VTessera(iComp,iTile)',' ',VTessera(iComp,iTile),1,1)
            ipFnlc = ipFnlc+nDAO
          end do

        end do
5000    continue
      end do
    end do

    call GetMem('DSOpr ','Free','REAL',ipDSOp,nSO*iPrim*jPrim)
    call GetMem(' DAO ','Free','Real',ipDAO,iPrim*jPrim*nElem(iAng)*nElem(jAng))
    call GetMem('ScrSph','Free','Real',iScrt2,nScr2)
    call GetMem('Scrtch','Free','Real',iScrt1,nScr1)
    call GetMem('Final','Free','Real',ipFnl,lFinal)
    call GetMem('Kernel','Free','Real',iKern,MemKrn)
131 continue
  end do
end do

call GetMem('PCoor','FREE','REAL',iPCoor,S%m2Max*3)
call GetMem('Kappa','FREE','REAL',iKappa,S%m2Max)
call GetMem('Zeta','FREE','REAL',ipZI,S%m2Max)
call GetMem('Zeta','FREE','REAL',iZeta,S%m2Max)

!call GetMem('drv_ef_PCM','CHEC','REAL',iDum,iDum)

return

end subroutine drv_ef_PCM
