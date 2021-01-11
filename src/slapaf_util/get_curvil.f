************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Get_CurviL(
     &              nq,nqRF,nqB,nqA,nqT,nqO,
     &              nsAtom,iIter,nIter,Cx,
     &              Process,Value,nB,
     &              qLbl,iRef,
     &              fconst,rMult,LuIC,Indq,iPrv,
     &              Grad_all,iGlow,iGhi,Proc_dB,
     &              iTabBonds,iTabAtoms,nBonds,nMax,iTabAI,mAtoms,
     &              mB_Tot,mdB_Tot,
     &              BM,dBM,iBM,idBM,
     &              nB_Tot,ndB_Tot,mqB,Thr_small)
      Implicit Real*8 (a-h,o-z)
      Real*8 Cx(3,nsAtom,nIter), fconst(nB), Value(nB,nIter), rMult(nB),
     &       Grad_all(9,iGlow:iGhi,nIter), BM(nB_Tot),dBM(ndB_Tot)
      Integer Indq(3,nB),
     &        iTabBonds(3,nBonds), iTabAtoms(2,0:nMax,mAtoms),
     &        iTabAI(2,mAtoms), iBM(nB_Tot), idBM(2,ndB_Tot), mqB(nB)
      Logical Process, Proc_dB
      Character(LEN=14) qLbl(nB)
#include "Molcas.fh"
*                                                                      *
************************************************************************
*                                                                      *
      nq = 0
      mB_Tot = 0
      mdB_Tot = 0
      If (Process) Call FZero(BM,nB_Tot)
      If (Proc_dB) Call FZero(dBM,ndB_Tot)
*                                                                      *
************************************************************************
*                                                                      *
      nq_=nq
      Call RF_Coord(
     &              nq,
     &              nsAtom,iIter,nIter,Cx,
     &              Process,Value,nB,
     &              qLbl,iRef,
     &              fconst,rMult,LuIC,Indq,
     &              Proc_dB,mB_Tot,mdB_Tot,
     &              BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,mqB)
      nqRF=nq-nq_
*
      nq_=nq
      Call Bond_List(
     &              nq,
     &              nsAtom,iIter,nIter,Cx,
     &              Process,Value,nB,
     &              qLbl,fconst,
     &              rMult,LuIC,Indq,
     &              Proc_dB,iTabBonds,nBonds,iTabAI,mAtoms,
     &              mB_Tot,mdB_Tot,
     &              BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,mqB)
      nqB=nq-nq_
*
      nq_=nq
      Call Angle_List(
     &              nq,
     &              nsAtom,iIter,nIter,Cx,
     &              Process,Value,nB,
     &              qLbl,iRef,
     &              fconst,rMult,LuIC,Indq,
     &              Grad_all,iGlow,iGhi,iPrv,Proc_dB,
     &              iTabBonds,nBonds,iTabAI,mAtoms,iTabAtoms,nMax,
     &              mB_Tot,mdB_Tot,
     &              BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,mqB,Thr_small)
      nqA=nq-nq_
*
      nq_=nq
      Call Torsion_List(
     &              nq,
     &              nsAtom,iIter,nIter,Cx,
     &              Process,Value,nB,
     &              qLbl,iRef,
     &              fconst,rMult,LuIC,Indq,iPrv,Proc_dB,
     &              iTabBonds,nBonds,iTabAI,mAtoms,iTabAtoms,nMax,
     &              mB_Tot,mdB_Tot,
     &              BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,mqB)
      nqT=nq-nq_
*
      nq_=nq
      Call OutOfPlane_List(
     &              nq,
     &              nsAtom,iIter,nIter,Cx,
     &              Process,Value,nB,
     &              qLbl,iRef,
     &              fconst,rMult,LuIC,Indq,iPrv,Proc_dB,
     &              iTabBonds,nBonds,iTabAI,mAtoms,iTabAtoms,nMax,
     &              mB_Tot,mdB_Tot,
     &              BM,dBM,iBM,idBM,nB_Tot,ndB_Tot,mqB)
      nqO=nq-nq_
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
