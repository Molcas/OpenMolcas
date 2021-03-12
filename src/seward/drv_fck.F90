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

use PAM2
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
#include "stdalloc.fh"
#include "print.fh"
#include "real.fh"
#include "warnings.fh"
character Label*8
real*8 CCoor(3,nComp), rNuc(nComp), PtChrg(nGrid)
integer ip(nComp), lOper(nComp), iChO(nComp), iStabO(0:7)
dimension opmol(*), opnuc(*), iopadr(*)
real*8, dimension(:), allocatable :: Int1El
integer iTwoj(0:7)
data iTwoj/1,2,4,8,16,32,64,128/

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 112
iPrint = nPrint(iRout)
if (iPrint >= 19) then
  write(6,*) ' In OneEl: Label',Label
  write(6,*) ' In OneEl: nComp'
  write(6,'(1X,8I5)') nComp
  write(6,*) ' In OneEl: lOper'
  write(6,'(1X,8I5)') lOper
  write(6,*) ' In OneEl: n2Tri'
  do iComp=1,nComp
    ip(iComp) = n2Tri(lOper(iComp))
  end do
  write(6,'(1X,8I5)')(ip(iComp),iComp=1,nComp)
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
    if (iand(lOper(iComp),iTwoj(iIrrep)) /= 0) nIC = nIC+1
  end do
end do
if (iPrint >= 20) write(6,*) ' nIC =',nIC
if (nIC == 0) Go To 999
call SOS(iStabO,nStabO,llOper)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate memory for symmetry adapted one electron integrals.
! Will just store the unique elements, i.e. low triangular blocks
! and lower triangular elements in the diagonal blocks.

call ICopy(nComp,[-1],0,ip,1)
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

call Drv_Fck_Inner(Label,ip,Int1El,LenTot,lOper,nComp,CCoor,nOrdOp,rHrmt,iChO,opmol,opnuc,ipad,iopadr,idirect,isyop,iStabO,nStabO, &
                   nIC,PtChrg,nGrid,iAddPot)
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
  !write(6,*) ' oneel *',Label,'*'

  call WrOne(iRC,iOpt,Label,iComp,Int1El(ip(iComp)),iSmLbl)

  if (Label(1:3) == 'PAM') call WrOne(iRC,iOpt,Label,1,Int1El(ip(iComp)),iSmLbl)
  iPAMcount = iPAMcount+1

  if (iRC /= 0) then
    write(6,*) ' *** Error in subroutine ONEEL ***'
    write(6,*) '     Abend in subroutine WrOne'
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
999 continue

return

end subroutine Drv_Fck

subroutine Drv_Fck_Inner(Label,ip,Int1El,LenTot,lOper,nComp,CCoor,nOrdOp,rHrmt,iChO,opmol,opnuc,ipad,iopadr,idirect,isyop,iStabO, &
                         nStabO,nIC,PtChrg,nGrid,iAddPot)
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

use Real_Spherical
use iSD_data
use Basis_Info
use Center_Info
use Sizes_of_Seward, only: S
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
#include "angtp.fh"
#include "real.fh"
#include "rmat_option.fh"
#include "stdalloc.fh"
#include "print.fh"
#include "nsd.fh"
#include "setup.fh"
real*8 A(3), B(3), RB(3), CCoor(3,nComp), PtChrg(nGrid)
character ChOper(0:7)*3, Label*8
integer nOp(2), ip(nComp), lOper(nComp), iChO(nComp), iDCRR(0:7), iDCRT(0:7), iStabM(0:7), iStabO(0:7)
integer iTwoj(0:7)
dimension opmol(*), opnuc(*), iopadr(*)
real*8, dimension(:), allocatable :: Zeta, ZI, SO, Fnl
real*8 Int1El(LenTot)
data iTwoj/1,2,4,8,16,32,64,128/
data ChOper/'E  ','x  ','y  ','xy ','z  ','xz ','yz ','xyz'/

!     Statement functions
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

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
    if (iPrint >= 29) write(6,*) ' nSO=',nSO
    if (nSO == 0) Go To 131
    call mma_allocate(SO,nSO*iBas*jBas)
    call DCopy_(nSO*iBas*jBas,[Zero],0,SO,1)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    if (iPrint >= 19) write(6,'(A,A,A,A,A)') ' ***** (',AngTp(iAng),',',AngTp(jAng),') *****'
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Allocate memory for the final integrals all in the
    ! primitive basis.
    lFinal = nIC*S%MaxPrm(iAng)*S%MaxPrm(jAng)*nElem(iAng)*nElem(jAng)
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
      write(6,*)
      write(6,*) ' g      =',nIrrep
      write(6,*) ' u      =',dc(mdci)%nStab
      write(6,'(9A)') '(U)=',(ChOper(dc(mdci)%iStab(ii)),ii=0,dc(mdci)%nStab-1)
      write(6,*) ' v      =',dc(mdcj)%nStab
      write(6,'(9A)') '(V)=',(ChOper(dc(mdcj)%iStab(ii)),ii=0,dc(mdcj)%nStab-1)
      write(6,*) ' LambdaR=',LmbdR
      write(6,*) ' r      =',nDCRR
      write(6,'(9A)') '(R)=',(ChOper(iDCRR(ii)),ii=0,nDCRR-1)
      write(6,*) ' m      =',nStabM
      write(6,'(9A)') '(M)=',(ChOper(iStabM(ii)),ii=0,nStabM-1)
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Compute normalization factor

    iuv = dc(mdci)%nStab*dc(mdcj)%nStab
    if (MolWgh == 1) then
      Fact = dble(nStabO)/dble(LambdT)
    else if (MolWgh == 0) then
      Fact = dble(iuv*nStabO)/dble(nIrrep**2*LambdT)
    else
      Fact = sqrt(dble(iuv))*dble(nStabO)/dble(nirrep*LambdT)
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
      if (iPrint >= 49) write(6,'(A,3F6.2,2X,3F6.2)') '*',(A(i),i=1,3),(RB(i),i=1,3)
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
            write(6,*) 'ijB,ijC=',ijB,ijC
            write(6,*) 'Fnl(iTo),Shells(iShll)%FockOp(iB,jB)=',Fnl(iTo),Shells(iShll)%FockOp(iB,jB)
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
            if (iand(lOper(iComp),iTwoj(iIrrep)) /= 0) iIC = iIC+1
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
      write(6,*) ' Scaling SO''s',Fact
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
    131 continue
  end do
end do

call mma_deallocate(ZI)
call mma_deallocate(Zeta)

return

! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(CCoor)
  call Unused_integer(nOrdOp)
  call Unused_integer_array(iChO)
  call Unused_real_array(opmol)
  call Unused_real_array(opnuc)
  call Unused_integer(ipad)
  call Unused_integer_array(iopadr)
  call Unused_integer(idirect)
  call Unused_integer(isyop)
  call Unused_real_array(PtChrg)
  call Unused_integer(iAddPot)
end if

end subroutine Drv_Fck_Inner
