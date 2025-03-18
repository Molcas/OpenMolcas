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
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************
      SUBROUTINE CalcWop(Wop,D,PUVX,NPUVX,IndTUVX,Coeff,Off_Ash)
      use MCLR_Data, only: nNA, nDens2, ipMat
      use input_mclr, only: nSym,nAsh,nBas,nIsh
      Implicit None
!*****Input
      INTEGER NPUVX
      Real*8 Coeff
      REAL*8,DIMENSION(nnA**2)::D
      REAL*8,DIMENSION(NPUVX) ::PUVX
      INTEGER,DIMENSION(nnA,nnA,nnA,nnA)::IndTUVX
      INTEGER,DIMENSION(nSym)::  Off_Ash
!*****Output
      REAL*8,DIMENSION(nDens2)::Wop
!*****Auxiliaries
      INTEGER jSym,it,iu,t,u,v,x,pt,qu,iLoc1,iLoc2,jAsh
      REAL*8 tempd1

      DO jSym=1,nSym
       jAsh=nAsh(jSym)
       IF(jAsh.eq.0) Cycle
       Do iu=1,jAsh
        u=iu+off_Ash(jSym)
        qu=iu+nIsh(jSym)
        iLoc1=(qu-1)*nBas(jSym)+ipMat(jSym,jSym)-1
       Do it=1,jAsh
        t=it+off_Ash(jSym)
        pt=it+nIsh(jSym)
        tempd1=0.0d0
        do v=1,nnA
         iLoc2=(v-1)*nnA
        do x=1,nnA
         IF(IndTUVX(t,u,v,x).ne.0)                                      &
     &    tempd1=tempd1+D(iLoc2+x)*PUVX(IndTUVX(t,u,v,x))
        end do
        end do
        Wop(iLoc1+pt)=tempd1
       End Do
       End Do
      END DO

      CALL DScal_(nDens2,Coeff,Wop,1)

      END SUBROUTINE CalcWop
