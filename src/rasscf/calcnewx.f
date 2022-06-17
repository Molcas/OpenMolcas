***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2022, Jie J. Bao                                       *
************************************************************************
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Apr. 12, 2022, created this file.               *
* ****************************************************************
      Subroutine CalcNewX(X,H,G,nSPair,XScr,GScr,EigVal,
     &                    ScrDiag,nScr)

C     Commented lines are options under development and may be used in
C     future
      use CMS, only:CMSThres,PosHess
      INTEGER nSPair,INFO,iPair,nScr
      Real*8 X(nSPair),G(nSPair),XScr(nSPair)
      Real*8 H(nSPair**2),ScrDiag(nScr),GScr(nSPair)
      Real*8 EigVal(nSPair**2)
      Real*8 MinGrad,ThreG,ThreH
******Thanks to Matthew R. Hermes for this algorithm

******Solve for x in hx=-g.
******1. diagonalize h to h' with rotation matrix C
******2. (C^T)h(C)(C^T)x=(C^T)g
******3. define x'=(C^T)x; g'=(C^T)g, we have (h')(x')=(g')
******4. (x')_p = (g')_p / (h')_pp
******5. x=(C)(C^T)x=(C)(x')

*      write(6,*) 'gradient before diag'
*      CALL RecPrt(' ',' ',G,1,nSPair)
*      write(6,*) 'hessian before diag'
*      CALL RecPrt(' ',' ',EigVal,1,nSPair)
******Step 1
      CALL DSYEV_('V','U',nSPair,H,nSPair,EigVal,ScrDiag,nScr,INFO)

******Step 3, g'=(C^T)g
      CALL DGEMM_('n','n',1,nSPair,nSPair,1.0d0,G,1,H,nSPair,
     &                                    0.0d0,GScr,1)

******Step 4
      ThreG=CMSThres*1.0d-2
      ThreH=CMSThres*1.0d-4
      MinGrad=CMSThres*1.0d5
      PosHess=.false.
      write(6,*) 'gradient'
      CALL RecPrt(' ','(1X,15(F9.6,1X))',GScr,1,nSPair)
      write(6,*) 'hessian'
      CALL RecPrt(' ','(1X,15(F9.6,1X))',EigVal,1,nSPair)

      DO iPair=1,nSPair
C       IF(abs(GScr(iPair)).lt.ThreG)  THEN
C        If(EigVal(iPair).gt.ThreH) Then
C        write(6,*) 'encounter local minimum',iPair
C        write(6,*) 'Gradient=',GScr(iPair)
C        write(6,*) 'Hessian =',EigVal(iPair)
C        XScr(iPair)=MinStep
C        Else If(abs(EigVal(iPair)).lt.ThreH) Then
C         write(6,*) 'encounter independent variable',iPair
C         XScr(iPair)=0.0d0
C        Else
C         XScr(iPair)=GScr(iPair)/Abs(EigVal(iPair))
C        End If
C       ELSE


        IF(EigVal(iPair).gt.ThreH) THEN
         PosHess=.true.
        END IF

        IF(      (Abs(GScr(iPair))  .lt.ThreG)
     &      .and.(Abs(EigVal(iPair)).lt.ThreH)) THEN
C         write(6,*) 'constant Qaa for pair',ipair
         XScr(iPair)=0.0d0
        ELSE
         XScr(iPair)=GScr(iPair)/Abs(EigVal(iPair))
        END IF

        IF(      (Abs(GScr(iPair)) .lt. ThreG)
     &      .and.(   EigVal(iPair) .gt. ThreH)) THEN
C         write(6,*) 'local minimum for pair',ipair
         XScr(iPair)=MinGrad/Abs(EigVal(iPair))
        END IF

C       END IF
      END DO

      write(6,*) 'steps taken'
      CALL RecPrt(' ','(1X,15(F9.6,1X))',XScr,1,nSPair)
******Step 5
      CALL DGEMM_('n','t',1,nSPair,nSPair,
     &            1.0d0,XScr,1,H,nSPair,
     &            0.0d0,X,1)
      RETURN
      End Subroutine




************************************************************************

      Subroutine CMSScaleX(X,R,DeltaR,Qnew,Qold,
     &                     RCopy,GDCopy,DgCopy,
     &                     GDstate,GDOrbit,Dgstate,DgOrbit,DDg,
     &                     nSPair,lRoots2,nGD,NAC2,nDDg,Saved)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      INTEGER nSPair,lRoots2,nGD,NAC2,nDDg
      Real*8 X(nSPair),R(lRoots2),DeltaR(lRoots2),RCopy(lRoots2),
     &       GDCopy(nGD),DgCopy(nGD),GDState(nGD),Dgstate(nGD),
     &       GDOrbit(nGD),DgOrbit(nGD),DDg(nDDg)
      Real*8 Qnew,Qold
      Logical Saved

      INTEGER iterscale,nScaleMax

      Saved=.true.
      ITERscale=0

*      NScaleMax=ICMSIterMax
      NScaleMax=max(20,ICMSIterMax)
      NScaleMax=min(80,NScaleMax)
      DO WHILE ((Qold-Qnew).gt.CMSThreshold)
       ITERscale=ITERscale+1
*       write(6,*) 'rescaling',ITERscale
       IF(ITERscale.eq.nScaleMax) THEN
        CALL FZero(X,nSPair)
        write(6,*) 'Scaling does not save Qaa from decreasing'
        Saved=.false.
        Exit
       END IF
       CALL DCopy_(lRoots2,RCopy,1,R,1)
       CALL DCopy_(nGD,GDCopy,1,GDState,1)
       CALL DCopy_(nGD,DgCopy,1,DgState,1)
       CALL DScal_(nSPair,0.5d0,X,1)

       CALL UpDateRotMat(R,DeltaR,X,lRoots,nSPair)
       CALL RotGD(GDstate,DeltaR,nGD,lRoots,NAC2)
       CALL RotGD(Dgstate,DeltaR,nGD,lRoots,NAC2)
       CALL TransposeMat(Dgorbit,Dgstate,nGD,lRoots2,NAC2)
       CALL TransposeMat(GDorbit,GDstate,nGD,lRoots2,NAC2)
       CALL CalcDDg(DDg,GDorbit,Dgorbit,nDDg,nGD,lRoots2,NAC2)
       CALL CalcQaa(Qnew,DDg,lRoots,nDDg)

C       write(6,*) 'rescaling',ITERscale,Qold,Qnew,CMSThreshold,
C     &            ICMSIterMax

*       IF(ITERscale.eq.iCMSIterMax) THEN
*        Exit
*       END IF
      END DO

      RETURN
      End Subroutine
