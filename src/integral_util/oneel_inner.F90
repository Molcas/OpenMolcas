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

#include "compiler_features.h"
#ifdef _IN_MODULE_

!#define _DEBUGPRINT_
subroutine OneEl_Inner(Kernel,KrnlMm,Label,ip,lOper,nComp,CoorO,nOrdOp,rHrmt,iChO,opmol,opnuc,ipad,iopadr,idirect,isyop,iStabO, &
                       nStabO,nIC,PtChrg,nGrid,iAddPot,Array,LenTot)
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
!             Roland Lindh, Dept. of Theoretical Chemistry, University *
!             of Lund, SWEDEN.                                         *
!             Modified for general kernel routines January  91         *
!             Modified for nonsymmetrical operators February  91       *
!             Modified for better symmetry treatement October  93      *
!             Modified loop structure April 99                         *
!***********************************************************************

use iSD_data, only: iSD
use Basis_Info, only: dbsc
use Sizes_of_Seward, only: S
use stdalloc, only: mma_allocate, mma_deallocate
use rmat, only: RMat_Type_Integrals
use property_label, only: PLabel
use Constants, only: Zero, One
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
procedure(int_kernel) :: Kernel
procedure(int_mem) :: KrnlMm
character(len=8) Label
integer nComp, nOrdOp, ipad, idirect, isyop, nIC, iAddPot, LenTot
integer ip(nComp), lOper(nComp), iChO(nComp), iStabO(0:7), nGrid
real*8 CoorO(3,nComp), PtChrg(nGrid)
real*8 opmol(*), opnuc(*)
real*8 rHrmt
integer iopadr(nComp,*)
real*8 Array(LenTot)
logical, external :: Rsv_Tsk
real*8, allocatable, target :: Kern(:)
integer, allocatable :: Ind_ij(:,:)
logical Do_PGamma
real*8, dimension(:), allocatable :: Zeta, ZI, Kappa, PCoor, SOInt, Scrtch, ScrSph
real*8, allocatable, target :: FArray(:)
integer, external :: n2Tri, MemSO1
integer ixyz, nElem, iPrint, nSkal, nIJS, iS, jS, i, lFinal, lScrt1, lScrt2, MemKrn, ijS, iPrim, jPrim, iBas, jBas, iAng, jAng, &
        mFinal, mScrt1, mScrt2, lA0, lB0, MemBux, MemCux, MemKer, ijSh, iCmp, iAO, iShell, iCnttp, jCmp, jAO, jShell, jCnttp, nSO, &
        iComp, iSmLbl, ipSO, nStabO, iSOBlk, mSO, MemAux, lA1, lB1, l_SOInt, id_Tsk, nOrder
real*8 rHrmt_Save
! Statement function
nElem(ixyz) = (ixyz+1)*(ixyz+2)/2

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
iPrint = 99
#else
iPrint = 5
#endif
RMat_type_integrals = .false.
Do_PGamma = .true.

! Auxiliary memory allocation.

call mma_allocate(Zeta,S%m2Max,label='Zeta')
call mma_allocate(ZI,S%m2Max,label='ZI')
call mma_allocate(Kappa,S%m2Max,label='Kappa')
call mma_allocate(PCoor,S%m2Max*3,label='PCoor')
!                                                                      *
!***********************************************************************
!                                                                      *
call Nr_Shells(nSkal)
!                                                                      *
!***********************************************************************
!                                                                      *
! Double loop over shells. These loops decide the integral type
!                                                                      *
!***********************************************************************
!                                                                      *
! Create list of non-vanishing pairs

call mma_allocate(Ind_ij,2,nskal*(nSkal+1)/2,label='Ind_ij')
nijS = 0
is = 0
js = 0
do I=1,nSkal*(nSkal+1)/2
  nijS = nijS+1
  js = js+1
  if (jS > iS) then
    iS = jS
    jS = 1
  end if
  Ind_ij(1,nijS) = iS
  Ind_ij(2,nijS) = jS
end do
call Init_Tsk(id_Tsk,nijS)
!                                                                      *
!***********************************************************************
!                                                                      *
! Allocate scratch for the integral evaluation.

lFinal = 1
lScrt1 = 1
lScrt2 = 1
MemKrn = 1
do ijS=1,nijS
  iS = Ind_ij(1,ijS)
  jS = Ind_ij(2,ijS)
  iPrim = iSD(5,iS)
  jPrim = iSD(5,jS)
  iBas = iSD(3,iS)
  jBas = iSD(3,jS)
  iAng = iSD(1,iS)
  jAng = iSD(1,jS)

  mFinal = nIC*iPrim*jPrim*nElem(iAng)*nElem(jAng)
  lFinal = max(lFinal,mFinal)

  if (Label(1:3) == 'MAG') cycle
  mScrt1 = nIC*max(iPrim,jBas)*max(iBas,jPrim)*nElem(iAng)*nElem(jAng)
  lScrt1 = max(mScrt1,lScrt1)

  mScrt2 = nIC*iBas*jBas*nElem(iAng)*nElem(jAng)
  lScrt2 = max(mScrt2,lScrt2)

  call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)

  if (PLabel /= ' ') then
    la0 = iAng
    lb0 = jAng
    MemAux = 1+3*nElem(la0)*nElem(lb0+1)*nIC
    la1 = la0
    lb1 = lb0+1
    MemBux = 1+3*nElem(la1+1)*nElem(lb1)*nIC
    if (la1 /= 0) MemBux = MemBux+3*nElem(la1-1)*nElem(lb1)*nIC
    if (lb0 /= 0) then
      lb1 = lb0-1
      MemAux = MemAux+3*nElem(la0)*nElem(lb0-1)*nIC
      MemCux = 1+3*nElem(la1+1)*nElem(lb1)*nIC
      if (la1 /= 0) MemCux = MemCux+3*nElem(la1-1)*nElem(lb1)*nIC
    else
      MemCux = 0
    end if
    MemAux = MemAux+max(MemBux,MemCux)
    MemKer = MemKer+MemAux
  end if

  MemKrn = max(MemKer*iPrim*jPrim,MemKrn)
end do

call mma_Allocate(FArray,lFinal,label='Final')
call mma_allocate(Scrtch,lScrt1,label='Scrtch')
call mma_allocate(ScrSph,lScrt2,label='ScrSph')
call mma_allocate(Kern,MemKrn,label='Kern')
!                                                                      *
!***********************************************************************
!                                                                      *
! big loop over individual tasks, distributed over individual nodes
ijSh = 0
10 continue
! make reservation of a task on global task list and get task range
! in return. Function will be false if no more tasks to execute.
if (.not. Rsv_Tsk(id_Tsk,ijSh)) Go To 11
iS = Ind_ij(1,ijSh)
jS = Ind_ij(2,ijSh)

iCmp = iSD(2,iS)
iBas = iSD(3,iS)
iAO = iSD(7,iS)
iShell = iSD(11,iS)
iCnttp = iSD(13,iS)

jCmp = iSD(2,jS)
jBas = iSD(3,jS)
jAO = iSD(7,jS)
jShell = iSD(11,jS)
jCnttp = iSD(13,jS)

nSO = 0
do iComp=1,nComp
  iSmLbl = lOper(iComp)
  nSO = nSO+MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
end do
#ifdef _DEBUGPRINT_
write(u6,*) ' nSO=',nSO
#endif

! Do not compute matrix elements in which electronic and
! muonic basis sets are mixed.

if ((nSO > 0) .and. (dbsc(iCnttp)%fMass == dbsc(jCnttp)%fMass)) then
  l_SOInt = iBas*jBas*nSO
  call mma_allocate(SOInt,l_SOInt,label='SOInt')
  SOInt(:) = Zero
  ipSO = 1
  call OneEl_IJ(iS,jS,iPrint,Do_PGamma,Zeta,ZI,Kappa,PCoor,Kernel,KrnlMm,Label,lOper,nComp,CoorO,nOrdOp,iChO,iStabO,nStabO,nIC, &
                PtChrg,nGrid,iAddPot,SOInt,l_SOInt,FArray,lFinal,Scrtch,lScrt1,ScrSph,lScrt2,Kern,MemKrn)
  iSOBlk = ipSO
  do iComp=1,nComp
    iSmLbl = lOper(iComp)
    if (n2Tri(iSmLbl) /= 0) then
      mSO = MemSO1(iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO)
    else
      mSO = 0
    end if

    ! Special trick for integrals over electromagnetic field
    ! radiation integrals.

    rHrmt_Save = rHrmt

    if ((Label(1:5) == 'EMFR ') .or. (Label(1:5) == 'TMOM ')) then
      if (mod((iComp+5),6) < 3) then
        rHrmt = One
      else
        rHrmt = -One
      end if
    end if
    !write(u6,*) 'Label,iComp,rHrmt=',Label,iComp,rHrmt
    if (mSO /= 0) then
      call SOSctt(SOInt(iSOBlk),iBas,jBas,mSO,Array(ip(iComp)),n2Tri(iSmLbl),iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO,nComp,Label, &
                  lOper,rHrmt)
      iSOBlk = iSOBlk+mSO*iBas*jBas
    end if
    rHrmt = rHrmt_Save
  end do
  call mma_deallocate(SOInt)
end if
goto 10
11 continue
call Free_Tsk(id_Tsk)
do iComp=1,nComp
  iSmLbl = lOper(iComp)
  call GADSum(Array(ip(iComp)),n2Tri(iSmLbl))
end do

call mma_deallocate(Kern)
call mma_deallocate(ScrSph)
call mma_deallocate(Scrtch)
call mma_deallocate(FArray)
call mma_deallocate(Ind_ij)
call mma_deallocate(PCoor)
call mma_deallocate(Kappa)
call mma_deallocate(ZI)
call mma_deallocate(Zeta)

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real_array(opmol)
  call Unused_real_array(opnuc)
  call Unused_integer(ipad)
  call Unused_integer_array(iopadr)
  call Unused_integer(idirect)
  call Unused_integer(isyop)
end if

end subroutine OneEl_Inner

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(OneEl_inner)

#endif
