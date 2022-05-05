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
! Jie J. Bao, on Dec. 08, 2021, created this file.               *
! ****************************************************************
      Subroutine CalcP2MOCube(P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,   &
     &                        nPMO3p,MOs,MOx,MOy,MOz,TabMO,P2Unzip,     &
     &                        mAO,mGrid,nMOs,do_grad)
      use nq_pdft, only: lft, lGGA
      use nq_Info
      Implicit Real*8 (A-H,O-Z)
#include "stdalloc.fh"

!*****Input
      INTEGER mAO,mGrid,nMOs,nPMO3p
      REAL*8,DIMENSION(mAO,mGrid,nMOs)::TabMO
      Real*8,DIMENSION(NASHT4)::P2Unzip
      Logical do_grad
!*****Output
      REAL*8,DIMENSION(mGrid*NASHT)::P2MOCube,MOs,MOx,MOy,MOz
      REAL*8,DIMENSION(nPMO3p)::P2MOCubex,P2MOCubey,P2MOCubez

!*****Auxiliary
      INTEGER iOff1,IOff2,IOff3,IIrrep,nGridPi,NASHT2,NASHT3,icount
      Real*8,DIMENSION(NASHT**3)::P2MO1
      Real*8,DIMENSION(NASHT**2)::P2MOSquare
      Logical lftGGA

      lftGGA=.false.
      IF(lft.and.lGGA) lftGGA=.true.
      nGridPi=mAO*mGrid
      DO iGrid=1,mGrid
       IOff1=(iGrid-1)*NASHT
       Do iIrrep=0,mIrrep-1
        IOff2=IOff_Ash(iIrrep)+1
        IOff3=IOff_BasAct(iIrrep)+1
        CALL DCopy_(nAsh(iIrrep),TabMO(1,iGrid,IOff3),nGridPi,          &
     &                             MOs(IOff1+IOff2)  ,1)
        do icount=1,nAsh(iIrrep)
        end do
       End Do
      END DO


      IF (lGGA) THEN
       DO iGrid=1,mGrid
        IOff1=(iGrid-1)*NASHT
        Do iIrrep=0,mIrrep-1
         IOff2=IOff_Ash(iIrrep)+1
         IOff3=IOff_BasAct(iIrrep)+1
         CALL DCopy_(nAsh(iIrrep),TabMO(2,iGrid,IOff3),nGridPi,         &
     &                              MOx(IOff1+IOff2)  ,1)
         CALL DCopy_(nAsh(iIrrep),TabMO(3,iGrid,IOff3),nGridPi,         &
     &                              MOy(IOff1+IOff2)  ,1)
         CALL DCopy_(nAsh(iIrrep),TabMO(4,iGrid,IOff3),nGridPi,         &
     &                              MOz(IOff1+IOff2)  ,1)
        End Do
       END DO
      END IF

      NASHT2=NASHT**2
      NASHT3=NASHT2*NASHT
      DO iGrid=1,mGrid
       IOff1=(iGrid-1)*NASHT+1

!       write(6,*) 'MOs array'
!       CALL RecPrt(' ','(10(F9.5,1X))',MOs(IOff1),1,NASHT)
!
!       write(6,*) '2RDM array'
!       CALL RecPrt(' ','(10(F9.5,1X))',P2Unzip,NASHT3,NASHT)

       CALL DGEMM_('T','N',NASHT3,1,NASHT,1.0d0,                        &
     & P2UnZip,NASHT,MOs(IOff1),NASHT,                                  &
     & 0.0d0,P2MO1,NASHT3)

!       write(6,*) 'P2MO1 array'
!       CALL RecPrt(' ','(10(F9.5,1X))',P2MO1,NASHT2,NASHT)

       CALL DGEMM_('T','N',NASHT2,1,NASHT,1.0d0,                        &
     & P2MO1,NASHT,MOs(IOff1),NASHT,                                    &
     & 0.0d0,P2MOSquare,NASHT2)

!       write(6,*) 'P2MOSquare array'
!       CALL RecPrt(' ','(10(F9.5,1X))',P2MOSquare,NASHT,NASHT)

       CALL DGEMM_('T','N',NASHT,1,NASHT,1.0d0,                         &
     & P2MOSquare,NASHT,MOs(IOff1),NASHT,                               &
     & 0.0d0,P2MOCube(iOff1),NASHT)

       IF(lftGGA.and.Do_Grad) THEN
        CALL DGEMM_('T','N',NASHT,1,NASHT,1.0d0,                        &
     &  P2MOSquare,NASHT,MOx(IOff1),NASHT,                              &
     &  0.0d0,P2MOCubex(iOff1),NASHT)
        CALL DGEMM_('T','N',NASHT,1,NASHT,1.0d0,                        &
     &  P2MOSquare,NASHT,MOy(IOff1),NASHT,                              &
     &  0.0d0,P2MOCubey(iOff1),NASHT)
        CALL DGEMM_('T','N',NASHT,1,NASHT,1.0d0,                        &
     &  P2MOSquare,NASHT,MOz(IOff1),NASHT,                              &
     &  0.0d0,P2MOCubez(iOff1),NASHT)

        CALL DGEMM_('T','N',NASHT2,1,NASHT,1.0d0,                       &
     &  P2MO1,NASHT,MOx(IOff1),NASHT,                                   &
     &  0.0d0,P2MOSquare,NASHT2)
        CALL DGEMM_('T','N',NASHT,1,NASHT,2.0d0,                        &
     &  P2MOSquare,NASHT,MOs(IOff1),NASHT,                              &
     &  1.0d0,P2MOCubex(iOff1),NASHT)

        CALL DGEMM_('T','N',NASHT2,1,NASHT,1.0d0,                       &
     &  P2MO1,NASHT,MOy(IOff1),NASHT,                                   &
     &  0.0d0,P2MOSquare,NASHT2)
        CALL DGEMM_('T','N',NASHT,1,NASHT,2.0d0,                        &
     &  P2MOSquare,NASHT,MOs(IOff1),NASHT,                              &
     &  1.0d0,P2MOCubey(iOff1),NASHT)

        CALL DGEMM_('T','N',NASHT2,1,NASHT,1.0d0,                       &
     &  P2MO1,NASHT,MOz(IOff1),NASHT,                                   &
     &  0.0d0,P2MOSquare,NASHT2)
        CALL DGEMM_('T','N',NASHT,1,NASHT,2.0d0,                        &
     &  P2MOSquare,NASHT,MOs(IOff1),NASHT,                              &
     &  1.0d0,P2MOCubez(iOff1),NASHT)
       END IF

!       write(6,*) 'P2MOCube array'
!       CALL RecPrt(' ','(10(F9.5,1X))',P2MOCube(IOff1),1,NASHT)
      END DO


      RETURN
      END SUBROUTINE
