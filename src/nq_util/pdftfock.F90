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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 22, 2021, created this file.               *
! ****************************************************************

      Subroutine PDFTFock(FI,FA,D1,mGrid,ActMO)
      use nq_pdft
      use nq_Info

!*****Input
      INTEGER mGrid
      Real*8,DIMENSION(mGrid*NASHT)::ActMO
      Real*8,DIMENSION(NASHT**2)::D1
!*****Output
      Real*8,DIMENSION(nPot1)::FI,FA
!*****Intermediate
      Real*8,DIMENSION(mGrid)::Fact2
      Real*8,DIMENSION(mGrid*NASHT)::SumDX,dEdPiAct
      Real*8 TempD1
      INTEGER iGrid,iIrrep,ik,k,iOff1,nGOrb
      Real*8 ddot_
      External DDot_

      nGOrb=mGrid*nOrbt

!*****calculate FI. FI=2*dEdPi*pq*sum_k{kk*Fact1}
!*****TempD1: sum_k{kk}
!*****Fact2 : 2sum_k{kk}
      DO iGrid=1,mGrid
       TempD1=0.0d0
       Do iIrrep=0,mIrrep-1
        do ik=1,nIsh(iIrrep)
         k=ik+OffOrb(iIrrep)
         IOff1=(k-1)*mGrid+iGrid
         TempD1=TempD1+MOas(IOff1)**2
        end do
       End Do
       Fact2(iGrid)=TempD1
      END DO

      CALL DScal_(mGrid,2.0d0,Fact2,1)
      IF(lft.and.lGGA) THEN
!***** In the end (drvnq_inner.f) there will be a process
!***** in which FI_pq=0.5(FI_pq+FI_qp)
!***** However, we do not want the factor of 0.5 for p'qrs part
!***** So we add another copy of GdEdPiMO to pick up the factor
!***** of 0.5
       CALL DAXpY_(nGOrb,1.0d0,GdEdPiMO,1,dEdPiMO,1)
!***** Also the pqr's part is needed with the help of the following
!****  array
       CALL TransActMO2(dEdPiAct,GdEdPiMO,mGrid)
      END IF
      CALL PDFTFock_Inner(FI,Fact2,dEdPiMO,MOas,mGrid)

      IF(lft.and.lGGA) THEN
       DO iGrid=1,mGrid
        TempD1=0.0d0
        Do iIrrep=0,mIrrep-1
         do ik=1,nIsh(iIrrep)
          k=ik+OffOrb(iIrrep)
          IOff1=(k-1)*mGrid+iGrid
          TempD1=TempD1+MOas(IOff1)*GdEdPiMO(iOff1)
         end do
        End Do
        Fact2(iGrid)=TempD1
       END DO
       CALL DScal_(mGrid,4.0d0,Fact2,1)
       CALL PDFTFock_Inner(FI,Fact2,MOas,MOas,mGrid)
      END IF
!*****calculate FA. FA=pq*sum_vx{vx*Dvx*Fact1}
!*****First calcualte sum_x{Dvx*x}
      DO iGrid=1,mGrid
       IOff1=(iGrid-1)*NASHT+1
       CALL DGEMM_('T','N',NASHT,1,NASHT,1.0d0,                         &
     & D1,NASHT,ActMO(IOff1),NASHT,0.0d0,SumDX(iOff1),NASHT)
       Fact2(iGrid)=ddot_(NASHT,ActMO(iOff1),1,SumDX(iOff1),1)
      END DO
      CALL PDFTFock_Inner(FA,Fact2,dEdPiMO,MOas,mGrid)
      IF(lft.and.lGGA) THEN
       DO iGrid=1,mGrid
        IOff1=(iGrid-1)*NASHT+1
        CALL DGEMM_('T','N',NASHT,1,NASHT,1.0d0,                        &
     &  D1,NASHT,ActMO(IOff1),NASHT,0.0d0,SumDX(iOff1),NASHT)
        Fact2(iGrid)=ddot_(NASHT,dEdPiAct(iOff1),1,SumDX(iOff1),1)
       END DO
       CALL DScal_(mGrid,2.0d0,Fact2,1)
       CALL PDFTFock_Inner(FA,Fact2,MOas,MOas,mGrid)
      END IF
      RETURN
      End Subroutine
