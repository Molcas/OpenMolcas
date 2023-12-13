************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2021, Paul B Calio                                     *
************************************************************************
* ****************************************************************
* history:                                                       *
* Based on cmsbxbp.f from Jie J. Bao                             *
* ****************************************************************
      subroutine CalcbXbP_CMSNAC(bX,bP,FMO1t,FMO2t,R,H,E_Final,nTri)
      use stdalloc, only : mma_allocate, mma_deallocate
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"
#include "crun_mclr.fh"


****** Output
       Real*8,DIMENSION((nRoots-1)*nRoots/2)::bX
       Real*8,DIMENSION(nConf1*nRoots)::bP
****** Input
       INTEGER nTri
       Real*8,DIMENSION(nRoots*nTri)::FMO1t
       Real*8,DIMENSION(nRoots*nacpr2)::FMO2t
       Real*8,DIMENSION(nRoots**2)::R,H
       Real*8,DIMENSION(nRoots)::E_Final
****** Auxiliaries
       Real*8,DIMENSION(:),Allocatable::LOK,CSFOK

       CALL mma_allocate(CSFOK,nRoots*nConf1)
       CALL mma_allocate(LOK,nRoots**2)
****** Using CalcOMat in original CalcbXbP
       CALL CalcOMat(CSFOK,LOK,FMO1t,FMO2t,nTri)
       CALL CalcbP_CMSNAC(bP,CSFOK,LOK,R)
       CALL CalcbX_CMSNAC(bX,LOK,R,H,E_Final)
       CALL mma_deallocate(CSFOK)
       CALL mma_deallocate(LOK)
       return
       end subroutine
******************************************************

      Subroutine CalcbX_CMSNAC(bX,LOK,R,H,E_Final)
#include "stdalloc.fh"
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"
****** Output
      Real*8,DIMENSION((nRoots-1)*nRoots/2)::bX
****** Input
      Real*8,DIMENSION(nRoots**2)::R,H
      Real*8,DIMENSION(nRoots)::E_Final
      Real*8,DIMENSION(nRoots**2)::LOK
***** Auxiliaries
      INTEGER I,K,L,M,N,IKL,IIN,IJM,IKOL,IIK,IJK,IIL,IJL
      Real*8 TempD, dE_IJ

      CALL FZero(bX,(nRoots-1)*nRoots/2)
      I=NACstates(1)
      J=NACstates(2)
      dE_IJ = E_Final(I)-E_Final(J)

****** R_JK*HKL*RIL
      DO K=2,nRoots
       IIK=(I-1)*nRoots+K
       IJK=(J-1)*nRoots+K
      DO L=1,K-1
       IIL=IIK-K+L
       IJL=IJK-K+L
       IKL=(K-2)*(K-1)/2+L
       IKOL=(L-1)*nRoots+K
       ILOK=(K-1)*nRoots+L
*******Diagonal elements of R_JK * H_KL * R_IK
       bX(IKL)=2.0d0*(R(IJK)*R(IIK)*LOK(ILOK)-R(IJL)*R(IIL)*LOK(IKOL))
******* Additional NAC term (Requires only one line)
******* R_JK * <K|L> * R_IL =
******* R_JL * R_IK - R_JK * R_IL
       if(.not.isMECIMSPD) bX(IKL) = bX(IKL) +
     & dE_IJ * ( R(IJL)*R(IIK)-R(IJK)*R(IIL) )

*******Off-Diagonal elements of R_JK * H_KL * R_IK
       Do M=1,nRoots
        IJM=IJK-K+M
       Do N=1,nRoots
        if(M.eq.N) cycle
        TempD=0.0d0
        IIN=IIK-K+N
        IF(M.eq.K) TempD=TempD+H((L-1)*nRoots+N)
        IF(N.eq.K) TempD=TempD+H((M-1)*nRoots+L)
        IF(M.eq.L) TempD=TempD-H((K-1)*nRoots+N)
        IF(N.eq.L) TempD=TempD-H((M-1)*nRoots+K)
        bX(IKL)=bX(IKL)+TempD*R(IJM)*R(IIN)
       End Do
       End Do
      END DO
      END DO
      RETURN
      END SUBROUTINE
******************************************************


******************************************************
      subroutine CalcbP_CMSNAC(bP,CSFOK,LOK,R)
      use ipPage, only: W
#include "stdalloc.fh"
#include "Input.fh"
#include "disp_mclr.fh"
#include "Pointers.fh"
#include "Files_mclr.fh"
#include "detdim.fh"
#include "cicisp_mclr.fh"
#include "incdia.fh"
#include "spinfo_mclr.fh"
#include "real.fh"
#include "sa.fh"
***** Output
      Real*8,DIMENSION(nConf1*nRoots)::bP
***** Input
      Real*8,DIMENSION(nRoots*nConf1)::CSFOK
      Real*8,DIMENSION(nRoots**2)::LOK
      Real*8,DIMENSION(nRoots**2)::R
***** Kind quantities that help
      INTEGER I,L,K,iLoc1,iLoc2
      Real*8 tempd

      I=NACstates(1)
      J=NACstates(2)
      DO K=1,nRoots
       iLoc1=(K-1)*nConf1+1
       CALL DCopy_(nConf1,CSFOK(iLoc1),1,bP(iLoc1),1)
       Do L=1, nRoots
        tempd=-LOK((K-1)*nRoots+L)
        iLoc2=(L-1)*nConf1+1
        CALL dAXpY_(nConf1,tempd,W(ipci)%Vec(iLoc2),1,bP(iLoc1),1)
       End Do

       CALL DScal_(nConf1,2*R((J-1)*nRoots+K)*R((I-1)*nRoots+K),
     & bP(iLoc1),1)

      END DO
      RETURN
      End Subroutine
******************************************************

