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

!#define _DEBUGPRINT_
subroutine OneEl_Inner(Kernel,KrnlMm,Label,ip,lOper,nComp,CoorO,nOrdOp,rHrmt,iChO,iStabO,nStabO,nIC,PtChrg,nGrid,iAddPot,Array, &
                       LenTot)
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

use Index_Functions, only: nTri_Elem, nTri_Elem1
use iSD_data, only: iSD
use Basis_Info, only: dbsc
use Sizes_of_Seward, only: S
use rmat, only: RMat_Type_Integrals
use property_label, only: PLabel
use Integral_interfaces, only: int_kernel, int_mem, OneEl_IJ
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
procedure(int_kernel) :: Kernel
procedure(int_mem) :: KrnlMm
character(len=8), intent(in) :: Label
integer(kind=iwp), intent(in) :: nComp, ip(nComp), lOper(nComp), nOrdOp, iChO(nComp), iStabO(0:7), nStabO, nIC, nGrid, iAddPot, &
                                 LenTot
real(kind=wp), intent(in) :: CoorO(3,nComp), rHrmt, PtChrg(nGrid)
real(kind=wp), intent(inout) :: Array(LenTot)
integer(kind=iwp) :: i, iAng, iAO, iBas, iCmp, iCnttp, iComp, id_Tsk, ijS, ijSh, iPrim, iPrint, ipSO, iS, iShell, iSmLbl, iSOBlk, &
                     jAng, jAO, jBas, jCmp, jCnttp, jPrim, jS, jShell, l_SOInt, lA0, lA1, lB0, lB1, lFinal, lScrt1, lScrt2, &
                     MemAux, MemBux, MemCux, MemKer, MemKrn, mFinal, mScrt1, mScrt2, mSO, nIJS, nOrder, nSkal, nSO
real(kind=wp) :: rHrmt_in
logical(kind=iwp) :: Do_PGamma
integer(kind=iwp), allocatable :: Ind_ij(:,:)
real(kind=wp), allocatable :: Kappa(:), PCoor(:), ScrSph(:), Scrtch(:), SOInt(:), Zeta(:), ZI(:)
real(kind=wp), allocatable, target :: FArray(:), Kern(:)
integer(kind=iwp), external :: MemSO1, n2Tri
logical(kind=iwp), external :: Rsv_Tsk

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

call mma_allocate(Ind_ij,2,nTri_Elem(nSkal),label='Ind_ij')
nijS = 0
is = 0
js = 0
do I=1,nTri_Elem(nSkal)
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

  mFinal = nIC*iPrim*jPrim*nTri_Elem1(iAng)*nTri_Elem1(jAng)
  lFinal = max(lFinal,mFinal)

  if (Label(1:3) == 'MAG') cycle
  mScrt1 = nIC*max(iPrim,jBas)*max(iBas,jPrim)*nTri_Elem1(iAng)*nTri_Elem1(jAng)
  lScrt1 = max(mScrt1,lScrt1)

  mScrt2 = nIC*iBas*jBas*nTri_Elem1(iAng)*nTri_Elem1(jAng)
  lScrt2 = max(mScrt2,lScrt2)

  call KrnlMm(nOrder,MemKer,iAng,jAng,nOrdOp)

  if (PLabel /= ' ') then
    la0 = iAng
    lb0 = jAng
    MemAux = 1+3*nTri_Elem1(la0)*nTri_Elem1(lb0+1)*nIC
    la1 = la0
    lb1 = lb0+1
    MemBux = 1+3*nTri_Elem1(la1+1)*nTri_Elem1(lb1)*nIC
    if (la1 /= 0) MemBux = MemBux+3*nTri_Elem1(la1-1)*nTri_Elem1(lb1)*nIC
    if (lb0 /= 0) then
      lb1 = lb0-1
      MemAux = MemAux+3*nTri_Elem1(la0)*nTri_Elem1(lb0-1)*nIC
      MemCux = 1+3*nTri_Elem1(la1+1)*nTri_Elem1(lb1)*nIC
      if (la1 /= 0) MemCux = MemCux+3*nTri_Elem1(la1-1)*nTri_Elem1(lb1)*nIC
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
do
  ! make reservation of a task on global task list and get task range
  ! in return. Function will be false if no more tasks to execute.
  if (.not. Rsv_Tsk(id_Tsk,ijSh)) exit
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
# ifdef _DEBUGPRINT_
  write(u6,*) ' nSO=',nSO
# endif

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

      rHrmt_in = rHrmt

      if ((Label(1:5) == 'EMFR ') .or. (Label(1:5) == 'TMOM ')) then
        if (mod((iComp+5),6) < 3) then
          rHrmt_in = One
        else
          rHrmt_in = -One
        end if
      end if
      !write(u6,*) 'Label,iComp,rHrmt_in=',Label,iComp,rHrmt_in
      if (mSO /= 0) then
        call SOSctt(SOInt(iSOBlk),iBas,jBas,mSO,Array(ip(iComp)),n2Tri(iSmLbl),iSmLbl,iCmp,jCmp,iShell,jShell,iAO,jAO,rHrmt_in)
        iSOBlk = iSOBlk+mSO*iBas*jBas
      end if
    end do
    call mma_deallocate(SOInt)
  end if
end do
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

end subroutine OneEl_Inner
