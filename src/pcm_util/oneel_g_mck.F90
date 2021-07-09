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
!***********************************************************************

subroutine OneEl_g_mck(Kernel,KrnlMm,Grad,nGrad,DiffOp,CCoor,FD,nFD,lOper,nComp,nOrdOp,Label)
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
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for general kernel routines January '91         *
!             Modified for nonsymmetrical operators February '91       *
!             Modified for gradients October '91                       *
!***********************************************************************

use Real_Spherical, only: ipSph, RSph
use iSD_data, only: iSD
use Basis_Info, only: dbsc, MolWgh, Shells
use Center_Info, only: dc
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
external :: Kernel, KrnlMm
integer(kind=iwp), intent(in) :: nGrad, nFD, nComp, lOper(nComp), nOrdOp
real(kind=wp), intent(out) :: Grad(nGrad)
logical(kind=iwp), intent(in) :: DiffOp
real(kind=wp), intent(in) :: CCoor(3), FD(nFD)
character(len=80), intent(in) :: Label
integer(kind=iwp) :: i, iAng, iAO, iBas, iCar, iCmp, iCnt, iCnttp, iComp, iDCRR(0:7), iDCRT(0:7), ijS, iKappa, iKern, IndGrd(3,2), &
                     iPCoor, ipDAO, ipDSO, ipDSOp, ipFnl, iPrim, iPrint, ipZI, iRout, iS, iScrt1, iScrt2, iShell, iShll, iSmLbl, &
                     iStabM(0:7), iStabO(0:7), iuv, iZeta, jAng, jAO, jBas, jCmp, jCnt, jCnttp, jPrim, jS, jShell, jShll, kk, &
                     lDCRR, lFinal, llOper, LmbdR, LmbdT, mdci, mdcj, MemKer, MemKrn, nDCRR, nDCRT, niAng, njAng, nOp(2), nOrder, &
                     nScr1, nScr2, nSkal, nSO, nStabM, nStabO, nTasks
real(kind=wp) :: A(3), B(3), FactNd, RB(3)
logical(kind=iwp) :: AeqB, EQ, IfGrad(3,3)
character(len=3), parameter :: ChOper(0:7) = ['E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz']
integer(kind=iwp), external :: MemSO1, n2Tri, NrOpr
#include "angtp.fh"
#include "Molcas.fh"
#include "WrkSpc.fh"
#include "print.fh"
#include "disp.fh"

iRout = 112
iPrint = nPrint(iRout)
call dcopy_(nGrad,[Zero],0,Grad,1)

! Auxiliary memory allocation.

call GetMem('Zeta','ALLO','REAL',iZeta,S%m2Max)
call GetMem('Zeta','ALLO','REAL',ipZI,S%m2Max)
call GetMem('Kappa','ALLO','REAL',iKappa,S%m2Max)
call GetMem('PCoor','ALLO','REAL',iPCoor,S%m2Max*3)
!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Basis_Mode('Valence')
call Nr_Shells(nSkal)
call Setup_iSD()
!                                                                      *
!***********************************************************************
!                                                                      *
! Double loop over shells. These loops decide the integral type

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
  niAng = (iAng+1)*(iAng+2)/2

  !do jS=1,iS
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
  njAng = (jAng+1)*(jAng+2)/2

  iSmLbl = 1
  nSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
  if (nSO == 0) cycle

  ! Find the DCR for A and B

  call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)
  if ((.not. DiffOp) .and. (nDCRR == 1) .and. EQ(A,B)) cycle
  if (iPrint >= 49) write(u6,'(10A)') ' {R}=(',(ChOper(iDCRR(i)),i=0,nDCRR-1),')'

  if (iPrint >= 19) write(u6,'(A,A,A,A,A)') ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'

  ! Call kernel routine to get memory requirement.

  call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)
  MemKrn = MemKer*S%m2Max
  call GetMem('Kernel','ALLO','REAL',iKern,MemKrn)

  ! Allocate memory for the final integrals, all in the primitive basis.

  lFinal = 6*S%MaxPrm(iAng)*S%MaxPrm(jAng)*niAng*njAng*nComp
  call GetMem('Final','ALLO','REAL',ipFnl,lFinal)

  ! Scratch area for contraction step

  nScr1 = S%MaxPrm(iAng)*S%MaxPrm(jAng)*niAng*njAng
  call GetMem('Scrtch','ALLO','REAL',iScrt1,nScr1)

  ! Scratch area for the transformation to spherical gaussians

  nScr2 = S%MaxPrm(iAng)*S%MaxPrm(jAng)*niAng*njAng
  call GetMem('ScrSph','Allo','Real',iScrt2,nScr2)

  call GetMem(' DAO ','Allo','Real',ipDAO,iPrim*jPrim*niAng*njAng)

  ! At this point we can compute Zeta.

  call ZXia(Work(iZeta),Work(ipZI),iPrim,jPrim,Shells(iShll)%Exp,Shells(jShll)%Exp)

  do iCar=0,2
    IndGrd(iCar+1,1) = iSD(iCar+16,iS)
    IfGrad(iCar+1,1) = iSD(iCar+16,iS) /= 0
  end do

  AeqB = iS == jS

  do iCar=0,2
    IndGrd(iCar+1,2) = iSD(iCar+16,jS)
    IfGrad(iCar+1,2) = iSD(iCar+16,jS) /= 0
  end do

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
  call DGEMM_('T','T',jBas*nSO,iPrim,iBas,One,Work(ipDSO),iBas,Shells(iShll)%pCff,iPrim,Zero,Work(ipDSOp),jBas*nSO)
  ! Transform J,ABi to AB,ij
  call DGEMM_('T','T',nSO*iPrim,jPrim,jBas,One,Work(ipDSOp),jBas,Shells(jShll)%pCff,jPrim,Zero,Work(ipDSO),nSO*iPrim)
  ! Transpose to ij,AB
  call DGeTmO(Work(ipDSO),nSO,nSO,iPrim*jPrim,Work(ipDSOp),iPrim*jPrim)
  call GetMem('DSO ','Free','Real',ipDSO,nSO*iBas*jBas)

  if (iPrint >= 99) call RecPrt(' Decontracted 1st order density/Fock matrix',' ',Work(ipDSOp),iPrim*jPrim,nSO)

  ! Loops over symmetry operations.

  nOp(1) = NrOpr(0)
  ! VV: gcc bug: one has to use this if!
  if (nDCRR >= 1) then
    do lDCRR=0,nDCRR-1
      call OA(iDCRR(lDCRR),B,RB)
      nOp(2) = NrOpr(iDCRR(lDCRR))
      if (EQ(A,RB) .and. (.not. DiffOp)) cycle
      if (.not. DiffOp) then
        ! Use the translational invariance to reduce the set of
        ! gradients to compute
        do iCar=1,3
          if (IfGrad(iCar,1) .and. IfGrad(iCar,2)) then
            IfGrad(iCar,2) = .false.
            IndGrd(iCar,2) = -IndGrd(iCar,2)
          end if
        end do
      end if

      if (iPrint >= 49) then
        write(u6,'(10A)') ' {M}=(',(ChOper(iStabM(i)),i=0,nStabM-1),')'
      end if

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
        FactNd = FactNd*real(nIrrep**2,kind=wp)/real(iuv,kind=wp)
      else if (MolWgh == 2) then
        FactNd = sqrt(real(iuv,kind=wp))*real(nStabO,kind=wp)/real(nIrrep*LmbdT,kind=wp)
      end if

      if (iPrint >= 49) then
        write(u6,'(A,/,2(3F6.2,2X))') ' *** Centers A, RB ***',(A(i),i=1,3),(RB(i),i=1,3)
      end if

      ! Desymmetrize the matrix with which we will
      ! contracte the trace.

      call DesymD(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,Work(ipDAO),iPrim,jPrim,Work(ipDSOp),nSO,nOp,FactNd)

      ! Project the spherical harmonic space onto the
      ! cartesian space.

      kk = niAng*njAng
      if (Shells(iShll)%Transf .or. Shells(jShll)%Transf) then

        ! ---ij,AB --> AB,ij
        call DGeTmO(Work(ipDAO),iPrim*jPrim,iPrim*jPrim,iCmp*jCmp,Work(iScrt1),iCmp*jCmp)
        ! ---AB,ij --> ij,ab
        call SphCar(Work(iScrt1),iCmp*jCmp,iPrim*jPrim,Work(iScrt2),nScr2,RSph(ipSph(iAng)),iAng,Shells(iShll)%Transf, &
                    Shells(iShll)%Prjct,RSph(ipSph(jAng)),jAng,Shells(jShll)%Transf,Shells(jShll)%Prjct,Work(ipDAO),kk)
      end if
      if (iPrint >= 99) call RecPrt(' Decontracted FD in the cartesian space',' ',Work(ipDAO),iPrim*jPrim,kk)

      ! Compute kappa and P.

      call Setup1(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,A,RB,Work(iKappa),Work(iPCoor),Work(ipZI))

      ! Compute gradients of the primitive integrals and
      ! trace the result.

      call Kernel(Shells(iShll)%Exp,iPrim,Shells(jShll)%Exp,jPrim,Work(iZeta),Work(ipZI),Work(iKappa),Work(iPcoor),Work(ipFnl), &
                  iPrim*jPrim,iAng,jAng,A,RB,nOrder,Work(iKern),MemKer,Ccoor,nOrdOp,Grad,nGrad,IfGrad,IndGrd,Work(ipDAO),mdci, &
                  mdcj,nOp,lOper,nComp,iStabM,nStabM)
      if (iPrint >= 49) call PrGrad_mck(' In Oneel',Grad,nGrad,ChDisp,5)

    end do
  end if
  call GetMem('DSOpr ','Free','REAL',ipDSOp,nSO*iPrim*jPrim)
  call GetMem(' DAO ','Free','Real',ipDAO,iPrim*jPrim*niAng*njAng)
  call GetMem('ScrSph','Free','Real',iScrt2,nScr2)
  call GetMem('Scrtch','Free','Real',iScrt1,nScr1)
  call GetMem('Final','Free','Real',ipFnl,lFinal)
  call GetMem('Kernel','Free','Real',iKern,MemKrn)
  !end do
  !end do
end do

call Free_iSD()
call GetMem('Kappa','FREE','REAL',iKappa,S%m2Max)
call GetMem('PCoor','FREE','REAL',iPCoor,S%m2Max*3)
call GetMem('Zeta','FREE','REAL',ipZI,S%m2Max)
call GetMem('Zeta','FREE','REAL',iZeta,S%m2Max)

if (iPrint >= 15) call PrGrad_mck(Label,Grad,nGrad,ChDisp,5)

return

end subroutine OneEl_g_mck
