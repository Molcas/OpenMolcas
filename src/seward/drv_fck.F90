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
! Copyright (C) 1990,1991,1993,1999, Roland Lindh                      *
!               1990, IBM                                              *
!***********************************************************************

subroutine Drv_Fck(Label,ip,lOper,nComp,CCoor,nOrdOp,rNuc,rHrmt,iChO,opmol,ipad,opnuc,iopadr,idirect,isyop,PtChrg,nGrid,iAddPot)

use PAM2, only: iPAMcount
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
character(len=8), intent(inout) :: Label
integer(kind=iwp), intent(in) :: nComp, lOper(nComp), nOrdOp, iChO(nComp), ipad, iopadr(*), idirect, isyop, nGrid, iAddPot
integer(kind=iwp), intent(out) :: ip(nComp)
real(kind=wp), intent(in) :: CCoor(3,nComp), rNuc(nComp), rHrmt, opmol(*), opnuc(*), PtChrg(nGrid)
#include "print.fh"
#include "warnings.h"
integer(kind=iwp) :: iadr, iComp, iIrrep, iOpt, iPrint, iRC, iRout, iSmLbl, iStabO(0:7), LenInt, LenTot, llOper, nIC, nStabO
real(kind=wp), allocatable :: Int1El(:)
integer(kind=iwp), external :: n2Tri
#include "macros.fh"
unused_var(CCoor)
unused_var(nOrdOp)
unused_var(iChO)
unused_var(opmol(1))
unused_var(opnuc(1))
unused_var(ipad)
unused_var(iopadr(1))
unused_var(idirect)
unused_var(isyop)
unused_var(PtChrg)
unused_var(iAddPot)

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 112
iPrint = nPrint(iRout)
if (iPrint >= 19) then
  write(u6,*) ' In OneEl: Label',Label
  write(u6,*) ' In OneEl: nComp'
  write(u6,'(1X,8I5)') nComp
  write(u6,*) ' In OneEl: lOper'
  write(u6,'(1X,8I5)') lOper
  write(u6,*) ' In OneEl: n2Tri'
  do iComp=1,nComp
    ip(iComp) = n2Tri(lOper(iComp))
  end do
  write(u6,'(1X,8I5)') (ip(iComp),iComp=1,nComp)
  call RecPrt(' CCoor',' ',CCoor,3,nComp)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Compute the number of blocks from each component of the operator
!     and the irreps it will span.

nIC = 0
llOper = 0
do iComp=1,nComp
  llOper = ior(llOper,lOper(iComp))
  do iIrrep=0,nIrrep-1
    if (btest(lOper(iComp),iIrrep)) nIC = nIC+1
  end do
end do
if (iPrint >= 20) write(u6,*) ' nIC =',nIC
if (nIC == 0) return
call SOS(iStabO,nStabO,llOper)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate memory for symmetry adapted one electron integrals.
! Will just store the unique elements, i.e. low triangular blocks
! and lower triangular elements in the diagonal blocks.

ip(:) = -1
LenTot = 0
do iComp=1,nComp
  LenInt = n2Tri(lOper(iComp))
  LenTot = LenTot+LenInt+4
end do
call mma_allocate(Int1El,LenTot)
ip(1) = 1
call DCopy_(LenTot,[Zero],0,Int1El(ip(1)),1)
iadr = ip(1)
do iComp=1,nComp
  LenInt = n2Tri(lOper(iComp))
  ip(icomp) = iadr
  iadr = iadr+LenInt+4
  ! Copy center of operator to work area.
  call DCopy_(3,Ccoor(1,iComp),1,Int1El(ip(iComp)+LenInt),1)
  ! Copy nuclear contribution to work area.
  Int1El(ip(iComp)+LenInt+3) = rNuc(iComp)
end do
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Compute all SO integrals for all components of the operator.

call Drv_Fck_Inner(Label,ip,Int1El,LenTot,lOper,nComp,rHrmt,iStabO,nStabO,nIC)
!                                                                      *
!***********************************************************************
!                                                                      *
!                    P O S T P R O C E S S I N G                       *
!                                                                      *
!***********************************************************************
!                                                                      *
if (iPrint >= 10) call PrMtrx(Label,lOper,nComp,ip,Int1El)
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Write integrals to disc.

do iComp=1,nComp
  iSmLbl = lOper(iComp)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  !----- Write integrals to disc

  iOpt = 0
  iRC = -1
  if (Label(1:3) == 'PAM') write(Label,'(A5,I3.3)') 'PAM  ',iPAMcount
  !write(u6,*) ' oneel *',Label,'*'

  call WrOne(iRC,iOpt,Label,iComp,Int1El(ip(iComp)),iSmLbl)

  if (Label(1:3) == 'PAM') call WrOne(iRC,iOpt,Label,1,Int1El(ip(iComp)),iSmLbl)
  iPAMcount = iPAMcount+1

  if (iRC /= 0) then
    write(u6,*) ' *** Error in subroutine ONEEL ***'
    write(u6,*) '     Abend in subroutine WrOne'
    call Quit(_RC_IO_ERROR_WRITE_)
  end if
end do  ! iComp
!                                                                      *
!***********************************************************************
!                                                                      *
!---- Deallocate memory for integral

call mma_deallocate(Int1El)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Drv_Fck

subroutine Drv_Fck_Inner(Label,ip,Int1El,LenTot,lOper,nComp,rHrmt,iStabO,nStabO,nIC)
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
!             Modified for Hermite-Gauss quadrature November '90       *
!             Modified for Rys quadrature November '90                 *
!             Modified for multipole moments November '90              *
!                                                                      *
!             Modified for general kernel routines January  91         *
!             Modified for nonsymmetrical operators February  91       *
!             Modified for better symmetry treatement October  93      *
!             Modified loop structure April 99                         *
!***********************************************************************

use iSD_data, only: iSD
use Basis_Info, only: dbsc, MolWgh, Shells
use Center_Info, only: dc
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
character(len=8), intent(in) :: Label
integer(kind=iwp), intent(in) :: nComp, ip(nComp), LenTot, lOper(nComp), iStabO(0:7), nStabO, nIC
real(kind=wp), intent(in) :: Int1El(LenTot), rHrmt
#include "angtp.fh"
#include "print.fh"
integer(kind=iwp) :: i, iAng, iAO, iB, iBas, iC, iCmp, iCnt, iCnttp, iComp, iDCRR(0:7), iDCRT(0:7), iElem, ii, iIC, iIrrep, ijB, &
                     ijC, iPrim, iPrint, iRout, iS, iShell, iShll, iSmLbl, iSOBlk, iStabM(0:7), iTo, iuv, jAng, jAO, jB, jBas, &
                     jCmp, jCnt, jCnttp, jElem, jPrim, jS, jShell, jShll, LmbdR, LambdT, lDCRR, lFinal, mdci, mdcj, mSO, nDCRR, &
                     nDCRT, nOp(2), nSkal, nSO, nStabM
real(kind=wp) :: A(3), B(3), Fact, RB(3)
real(kind=wp), allocatable :: Zeta(:), ZI(:), SO(:), Fnl(:)
character(len=*), parameter :: ChOper(0:7) = ['E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz']
integer(kind=iwp), external :: MemSO1, n2Tri, NrOpr

iRout = 112
iPrint = nPrint(iRout)
!iPrint = 99

!-----Auxiliary memory allocation.

call mma_allocate(Zeta,S%m2Max)
call mma_allocate(ZI,S%m2Max)
!                                                                      *
!***********************************************************************
!                                                                      *
call Nr_Shells(nSkal)
!                                                                      *
!***********************************************************************
!                                                                      *
!-----Double loop over shells. These loops decide the integral type

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
  do jS=iS,iS
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
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Allocate memory for SO integrals that will be generated by
    ! this batch of AO integrals.

    nSO = 0
    do iComp=1,nComp
      iSmLbl = lOper(iComp)
      nSO = nSO+MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
    end do
    if (iPrint >= 29) write(u6,*) ' nSO=',nSO
    if (nSO == 0) cycle
    call mma_allocate(SO,nSO*iBas*jBas)
    call DCopy_(nSO*iBas*jBas,[Zero],0,SO,1)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    if (iPrint >= 19) write(u6,'(A,A,A,A,A)') ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Allocate memory for the final integrals all in the
    ! primitive basis.
    iElem = (iAng+1)*(iAng+2)/2
    jElem = (jAng+1)*(jAng+2)/2
    lFinal = nIC*S%MaxPrm(iAng)*S%MaxPrm(jAng)*iElem*jElem
    call mma_allocate(Fnl,lFinal)
    call dCopy_(lFinal,[Zero],0,Fnl,1)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! At this point we can compute Zeta.
    ! This is now computed in the ij or ji order.

    call ZXia(Zeta,ZI,iPrim,jPrim,Shells(iShll)%Exp,Shells(jShll)%Exp)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Find the DCR for A and B

    call DCR(LmbdR,dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iDCRR,nDCRR)

    ! Find the stabilizer for A and B

    call Inter(dc(mdci)%iStab,dc(mdci)%nStab,dc(mdcj)%iStab,dc(mdcj)%nStab,iStabM,nStabM)

    call DCR(LambdT,iStabM,nStabM,iStabO,nStabO,iDCRT,nDCRT)

    if (iPrint >= 19) then
      write(u6,*)
      write(u6,*) ' g      =',nIrrep
      write(u6,*) ' u      =',dc(mdci)%nStab
      write(u6,'(9A)') '(U)=',(ChOper(dc(mdci)%iStab(ii)),ii=0,dc(mdci)%nStab-1)
      write(u6,*) ' v      =',dc(mdcj)%nStab
      write(u6,'(9A)') '(V)=',(ChOper(dc(mdcj)%iStab(ii)),ii=0,dc(mdcj)%nStab-1)
      write(u6,*) ' LambdaR=',LmbdR
      write(u6,*) ' r      =',nDCRR
      write(u6,'(9A)') '(R)=',(ChOper(iDCRR(ii)),ii=0,nDCRR-1)
      write(u6,*) ' m      =',nStabM
      write(u6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Compute normalization factor

    iuv = dc(mdci)%nStab*dc(mdcj)%nStab
    if (MolWgh == 1) then
      Fact = real(nStabO,kind=wp)/real(LambdT,kind=wp)
    else if (MolWgh == 0) then
      Fact = real(iuv*nStabO,kind=wp)/real(nIrrep**2*LambdT,kind=wp)
    else
      Fact = sqrt(real(iuv,kind=wp))*real(nStabO,kind=wp)/real(nirrep*LambdT,kind=wp)
    end if
    Fact = One/Fact
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    !       Loops over symmetry operations acting on the basis.

    nOp(1) = NrOpr(0)
    !do lDCRR=0,nDCRR-1
    do lDCRR=0,0
      call OA(iDCRR(lDCRR),B,RB)
      nOp(2) = NrOpr(iDCRR(lDCRR))
      if (iPrint >= 49) write(u6,'(A,3F6.2,2X,3F6.2)') '*',(A(i),i=1,3),(RB(i),i=1,3)
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! Pick up epsilon from memory
      call FZero(Fnl,iBas*jBas*iCmp*jCmp*nIC)
      do iB=1,iBas
        do jB=1,iBas
          ijB = (jB-1)*iBas+iB
          do iC=1,iCmp
            ijC = (iC-1)*iCmp+iC
            iTo = +(ijC-1)*iBas**2+ijB
#           ifdef _DEBUGPRINT_
            write(u6,*) 'ijB,ijC=',ijB,ijC
            write(u6,*) 'Fnl(iTo),Shells(iShll)%FockOp(iB,jB)=',Fnl(iTo),Shells(iShll)%FockOp(iB,jB)
#           endif
            Fnl(iTo) = Shells(iShll)%FockOp(iB,jB)
          end do
        end do
      end do
#     ifdef _DEBUGPRINT_
      call RecPrt('EOrb',' ',Shells(iShll)%FockOp,iBas,1)
      call RecPrt('EOrb',' ',Shells(iShll)%FockOp,iBas,iBas)
      call RecPrt('FckInt',' ',Fnl,iBas*jBas,iCmp*jCmp*nIC)
#     endif
      !                                                                *
      !*****************************************************************
      !                                                                *
      ! At this point accumulate the batch of integrals onto the
      ! final symmetry adapted integrals.

      if (iPrint >= 99) then
        call RecPrt(' Accumulated SO integrals, so far...',' ',SO,iBas*jBas,nSO)
      end if

      !------Symmetry adapt component by component

      iSOBlk = 1
      iIC = 1
      do iComp=1,nComp
        iSmLbl = lOper(iComp)
        mSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
        if (mSO == 0) then
          do iIrrep=0,nIrrep-1
            if (btest(lOper(iComp),iIrrep)) iIC = iIC+1
          end do
        else
          call SymAd1(iSmLbl,iAng,jAng,iCmp,jCmp,iShell,jShell,iShll,jShll,iAO,jAO,Fnl,iBas,jBas,nIC,iIC,SO(iSOBlk),mSO,nOp)
          iSOBlk = iSOBlk+mSO*iBas*jBas
        end if
      end do

    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Multiply with factors due to projection operators

    if (Fact /= One) call DScal_(nSO*iBas*jBas,Fact,SO,1)
    if (iPrint >= 99) then
      write(u6,*) ' Scaling SO''s',Fact
      call RecPrt(' Accumulated SO integrals',' ',SO,iBas*jBas,nSO)
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Scatter the SO's on to the non-zero blocks of the
    ! lower triangle.

    iSOBlk = 1
    do iComp=1,nComp
      iSmLbl = lOper(iComp)
      if (n2Tri(iSmLbl) /= 0) then
        mSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
      else
        mSO = 0
      end if
      if (mSO /= 0) then
        call SOSctt(SO(iSOBlk),iBas,jBas,mSO,Int1El(ip(iComp)),n2Tri(iSmLbl),iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO,nComp,Label, &
                    lOper,rHrmt)
        iSOBlk = iSOBlk+mSO*iBas*jBas
      end if
    end do
!                                                                      *
!***********************************************************************
!                                                                      *
    call mma_deallocate(Fnl)
    call mma_deallocate(SO)
!                                                                      *
!***********************************************************************
!                                                                      *
  end do
end do

call mma_deallocate(ZI)
call mma_deallocate(Zeta)

return

end subroutine Drv_Fck_Inner
