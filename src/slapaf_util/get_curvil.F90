!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine Get_CurviL(nq,nqRF,nqB,nqA,nqT,nqO,nsAtom,iIter,nIter,Cx,Process,Val,nB,qLbl,iRef,fconst,rMult,LuIC,Indq,iPrv, &
                      Grad_all,iGlow,iGhi,Proc_dB,iTabBonds,iTabAtoms,nBonds,nMax,iTabAI,mAtoms,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM, &
                      nB_Tot,ndB_Tot,mqB,Thr_small)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(out) :: nq, nqRF, nqB, nqA, nqT, nqO, mB_Tot, mdB_Tot
integer(kind=iwp), intent(in) :: nsAtom, iIter, nIter, nB, iRef, LuIC, iPrv, iGlow, iGhi, nBonds, iTabBonds(3,nBonds), nMax, &
                                 mAtoms, iTabAtoms(2,0:nMax,mAtoms), iTabAI(2,mAtoms), nB_Tot, ndB_Tot
real(kind=wp), intent(in) :: Cx(3,nsAtom,nIter), Thr_small
logical(kind=iwp), intent(in) :: Process, Proc_dB
real(kind=wp), intent(inout) :: Val(nB,nIter), fconst(nB), rMult(nb), Grad_all(9,iGlow:iGhi,nIter), BM(nB_Tot), dBM(ndB_Tot)
integer(kind=iwp), intent(inout) :: Indq(3,nB), iBM(nB_Tot), idBM(2,ndB_Tot), mqB(nB)
character(len=14), intent(inout) :: qLbl(nB)
integer(kind=iwp) :: nq_

!                                                                      *
!***********************************************************************
!                                                                      *
nq = 0
mB_Tot = 0
mdB_Tot = 0
if (Process) BM(:) = Zero
if (Proc_dB) dBM(:) = Zero
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
if (Process) call RecPrt('BM initial',' ',BM,1,size(BM))
#endif
nq_ = nq
call RF_Coord(nq,nsAtom,iIter,nIter,Cx,Process,Val,nB,qLbl,iRef,fconst,rMult,LuIC,Indq,Proc_dB,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM, &
              nB_Tot,ndB_Tot,mqB)
nqRF = nq-nq_
#ifdef _DEBUGPRINT_
if (Process) call RecPrt('BM after RF_Coord',' ',BM,1,size(BM))
#endif

nq_ = nq
call Bond_List(nq,nsAtom,iIter,nIter,Cx,Process,Val,nB,qLbl,fconst,rMult,LuIC,Indq,Proc_dB,iTabBonds,nBonds,iTabAI,mAtoms,mB_Tot, &
               mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,mqB)
nqB = nq-nq_
#ifdef _DEBUGPRINT_
if (Process) call RecPrt('BM after Bond_List',' ',BM,1,size(BM))
#endif

nq_ = nq
call Angle_List(nq,nsAtom,iIter,nIter,Cx,Process,Val,nB,qLbl,iRef,fconst,rMult,LuIC,Indq,Grad_all,iGlow,iGhi,iPrv,Proc_dB, &
                iTabBonds,nBonds,iTabAI,mAtoms,iTabAtoms,nMax,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,mqB,Thr_small)
nqA = nq-nq_
#ifdef _DEBUGPRINT_
if (Process) call RecPrt('BM after Angle_List',' ',BM,1,size(BM))
#endif

nq_ = nq
call Torsion_List(nq,nsAtom,iIter,nIter,Cx,Process,Val,nB,qLbl,iRef,fconst,rMult,LuIC,Indq,iPrv,Proc_dB,iTabBonds,nBonds,iTabAI, &
                  mAtoms,iTabAtoms,nMax,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,mqB)
nqT = nq-nq_
#ifdef _DEBUGPRINT_
if (Process) call RecPrt('BM after Torsion_List',' ',BM,1,size(BM))
#endif

nq_ = nq
call OutOfPlane_List(nq,nsAtom,iIter,nIter,Cx,Process,Val,nB,qLbl,iRef,fconst,rMult,LuIC,Indq,iPrv,Proc_dB,iTabBonds,nBonds, &
                     iTabAI,mAtoms,iTabAtoms,nMax,mB_Tot,mdB_Tot,BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,mqB)
nqO = nq-nq_
#ifdef _DEBUGPRINT_
if (Process) call RecPrt('BM after OutOfPlane_List',' ',BM,1,size(BM))
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Get_CurviL
