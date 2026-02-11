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
! Copyright (C) 2022, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Apr. 12, 2022, created this file.               *
! ****************************************************************
      Subroutine CMSScaleX(X,R,DeltaR,Qnew,Qold,                        &
     &                     RCopy,GDCopy,DgCopy,                         &
     &                     GDstate,GDOrbit,Dgstate,DgOrbit,DDg,         &
     &                     nSPair,lRoots2,nGD,NAC2,nDDg,Saved)
      use CMS, only: NCMSScale
      use rasscf_global, only: CMSThreshold, lRoots
      Implicit None
#include "warnings.h"
      INTEGER nSPair,lRoots2,nGD,NAC2,nDDg
      Real*8 X(nSPair),R(lRoots2),DeltaR(lRoots2),RCopy(lRoots2),       &
     &       GDCopy(nGD),DgCopy(nGD),GDState(nGD),Dgstate(nGD),         &
     &       GDOrbit(nGD),DgOrbit(nGD),DDg(nDDg)
      Real*8 Qnew,Qold
      Logical Saved

      INTEGER nScaleMax

      Saved=.true.

      NScaleMax=5
      DO WHILE ((Qold-Qnew).gt.CMSThreshold)
       NCMSScale=NCMSScale+1
       IF(NCMSScale.eq.nScaleMax) THEN
        write(6,'(6X,A)')                                               &
     &  'Scaling does not save Qaa from decreasing.'
        write(6,'(6X,A)')                                               &
     &  'Q_a-a decreases for this step.'
        Saved=.false.
        Exit
       END IF
       CALL DCopy_(lRoots2,RCopy,1,R,1)
       CALL DCopy_(nGD,GDCopy,1,GDState,1)
       CALL DCopy_(nGD,DgCopy,1,DgState,1)
       CALL DScal_(nSPair,0.1d0,X,1)

       CALL UpDateRotMat(R,DeltaR,X,lRoots,nSPair)
       CALL RotGD(GDstate,DeltaR,nGD,lRoots,NAC2)
       CALL RotGD(Dgstate,DeltaR,nGD,lRoots,NAC2)
       CALL TransposeMat(Dgorbit,Dgstate,nGD,lRoots2,NAC2)
       CALL TransposeMat(GDorbit,GDstate,nGD,lRoots2,NAC2)
       CALL CalcDDg(DDg,GDorbit,Dgorbit,nDDg,nGD,lRoots2,NAC2)
       CALL CalcQaa(Qnew,DDg,lRoots,nDDg)

      END DO

      RETURN
      End Subroutine
