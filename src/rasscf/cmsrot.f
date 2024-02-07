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
* Copyright (C) 2020, Jie J. Bao                                       *
************************************************************************
      Subroutine CMSRot(TUVX)
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Aug. 06, 2020, created this file.               *
* ****************************************************************
      use stdalloc, only : mma_allocate, mma_deallocate
      use CMS, only: CMSNotConverged
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"

      Real*8,DIMENSION(NACPR2)::TUVX
      CHARACTER(len=16)::VecName
      Real*8,DIMENSION(:,:,:,:),Allocatable::Gtuvx
      Real*8,DIMENSION(:,:,:,:),Allocatable::DDG
      Real*8,DIMENSION(:,:,:),Allocatable::GDMat
      Real*8,DIMENSION(:,:),Allocatable::RotMat
C     Allocating Memory
      CALL mma_allocate(GDMat,LRoots*(LRoots+1)/2,NAC,NAC)
      CALL mma_allocate(RotMat,lRoots,lRoots)
      CALL mma_allocate(Gtuvx,NAC,NAC,NAC,NAC)
      CALL mma_allocate(DDG,lRoots,lRoots,lRoots,lRoots)

*     printing header
      write(6,*)
      write(6,*)
      write(6,*) '    CMS INTERMEDIATE-STATE OPTIMIZATION'
      IF(trim(CMSStartMat).eq.'XMS') THEN
       CALL ReadMat('ROT_VEC',VecName,RotMat,lroots,lroots,7,16,'N')
      ELSE
       CALL ReadMat(trim(CMSStartMat),VecName,RotMat,lroots,lroots,
     &              len_trim(CMSStartMat),16,'N')
      END IF
      CALL CMSHeader(trim(CMSStartMat),len_trim(CMSStartMat))


      CALL LoadGtuvx(TUVX,Gtuvx)

      CMSNotConverged=.false.
      CALL GetGDMat(GDMat)
      IF(lRoots.lt.NAC) THEN
*       write(6,*)"Optimization Approach 1"
       CALL GetDDgMat(DDg,GDMat,Gtuvx)
       CALL NStateOpt(RotMat,DDg)
      ELSE
*       write(6,*)"Optimization Approach 2"
       CALL NStateOpt2(RotMat,GDMat,Gtuvx)
      END IF
      VecName='CMS-PDFT'
      CALL PrintMat('ROT_VEC',VecName,RotMat,lroots,lroots,7,16,'N')

C     Deallocating Memory
      CALL mma_deallocate(GDMat)
      CALL mma_deallocate(RotMat)
      CALL mma_deallocate(Gtuvx)
      CALL mma_deallocate(DDg)
      IF(CMSNotConverged) THEN
       Call WarningMessage(2,'CMS Intermediate States Not Converged')
       Call Quit(_RC_NOT_CONVERGED_)
      END IF
      RETURN
      End Subroutine

************************************************************************

************************************************************************
      Subroutine NStateOpt(RotMat,DDg)
      use stdalloc, only : mma_allocate, mma_deallocate
      use CMS, only: CMSNotConverged
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      Real*8,DIMENSION(lRoots,lRoots,lRoots,lRoots)::DDG
      Real*8,DIMENSION(lroots,lroots)::RotMat

      INTEGER IState,JState,NPairs,IPair,ICMSIter
      Real*8 VeeSumOld,VeeSumNew,Threshold
      INTEGER,DIMENSION(:,:),Allocatable::StatePair
      Real*8,DIMENSION(:),Allocatable::theta
      Real*8,DIMENSION(:,:),Allocatable::FRot
      Logical Converged
      Real*8 CalcNSumVee
      External CalcNSumVee

      CALL mma_allocate(StatePair,LRoots*(LRoots-1)/2,2)
      CALL mma_allocate(theta,LRoots*(LRoots-1)/2)
      CALL mma_allocate(FRot,lRoots,lRoots)
      Threshold=CMSThreshold
      NPairs=lRoots*(lRoots-1)/2
      IPair=0
      DO IState=1,lRoots
       Do JState=1,IState-1
        IPair=IPair+1
        StatePair(IPair,1)=IState
        StatePair(IPair,2)=JState
       End Do
      END DO
      Converged=.false.
      CALL Copy2DMat(FRot,RotMat,lRoots,lRoots)
      VeeSumOld=CalcNSumVee(RotMat,DDg)
      ICMSIter=0
      DO WHILE(.not.Converged)
       Do IPair=1,NPairs
        theta(IPair)=0.0d0
       End Do
       ICMSIter=ICMSIter+1
       CALL ThetaOpt(FRot,theta,VeeSumNew,StatePair,NPairs,DDg)
       IF(lRoots.gt.2) THEN
       write(6,'(6X,I4,8X,F16.8,8X,ES16.4E3)')
     & ICMSIter,VeeSumNew,VeeSumNew-VeeSumOld
       ELSE
       write(6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)')
     & ICMSIter,asin(FRot(2,1))/atan(1.0d0)*45.0d0,VeeSumNew
     & ,VeeSumNew-VeeSumOld
       END IF
       IF(ABS(VeeSumNew-VeeSumOld).lt.Threshold) THEN
        If(ICMSIter.ge.ICMSIterMin) Then
         Converged=.true.
         write(6,'(4X,A)')'CONVERGENCE REACHED'
        End If
       ELSE
        if(ICMSIter.ge.ICMSIterMax) then
         Converged=.true.
         CMSNotConverged=.true.
         write(6,'(4X,A)')'NOT CONVERGED AFTER MAX NUMBER OF CYCLES'
         write(6,'(4X,A)')'TEMPORARY ROTATION MATRIX SAVED'
        end if
       END IF
       VeeSumOld=VeeSumNew
      END DO
      write(6,*)('=',i=1,71)
      CALL Copy2DMat(RotMat,FRot,lRoots,lRoots)
      CALL mma_deallocate(StatePair)
      CALL mma_deallocate(theta)
      CALL mma_deallocate(FRot)
      RETURN
      END SUBROUTINE
************************************************************************

************************************************************************
      Subroutine ThetaOpt(FRot,theta,SumVee,StatePair,NPairs,DDg)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      INTEGER NPairs
      Real*8 SumVee
      Real*8,DIMENSION(lRoots,lRoots,lRoots,lRoots)::DDG
      Real*8,DIMENSION(lroots,lroots)::FRot
      INTEGER,DIMENSION(NPairs,2)::StatePair
      Real*8,DIMENSION(NPairs)::theta

      INTEGER IPair,IState,JState
C      Real*8,DIMENSION(NPairs)::thetanew

      DO IPair=1,NPairs
       IState=StatePair(IPair,1)
       JState=StatePair(IPair,2)
       CALL
     & OptOneAngle(theta(iPair),SumVee,FRot,DDg,IState,JState,lRoots)
      END DO
      DO IPair=NPairs-1,1,-1
       IState=StatePair(IPair,1)
       JState=StatePair(IPair,2)
       CALL
     & OptOneAngle(theta(iPair),SumVee,FRot,DDg,IState,JState,lRoots)
      END DO
      RETURN
      END SUBROUTINE
************************************************************************
************************************************************************
      SUBROUTINE OptOneAngle(Angle,SumVee,RotMat,DDg,I1,I2,lRoots)
      use stdalloc, only : mma_allocate, mma_deallocate
      real*8 Angle,SumVee
      INTEGER I1,I2,lRoots
      Real*8,DIMENSION(lRoots,lRoots)::RotMat
      Real*8,DIMENSION(lRoots,lRoots,lRoots,lRoots)::DDG

      Logical Converged
      INTEGER Iter,IterMax,IA,IMax
      Real*8 Threshold,StepSize,SumOld
      Real*8,DIMENSION(:),Allocatable::Angles,Sums
      Real*8,DIMENSION(:),Allocatable::ScanA,ScanS
      Real*8,DIMENSION(:,:),Allocatable::RTmp

      Real*8 CalcNSumVee
      External CalcNSumVee
      INTEGER RMax
      External RMax

      CALL mma_allocate(Angles,4)
      CALL mma_allocate(Sums,4)
      CALL mma_allocate(ScanA,31)
      CALL mma_allocate(ScanS,31)
      CALL mma_allocate(RTmp,lRoots,lRoots)

      Converged=.false.
      stepsize=dble(atan(1.0d0))/15
      Threshold=1.0d-8

C       write(6,'(A,2(I2,2X))')
C     &'scanning rotation angles for ',I1,I2
      Angles(2)=0.0d0
      DO Iter=1,31
       ScanA(Iter)=(Iter-16)*stepsize*2
       CALL Copy2DMat(RTmp,RotMat,lRoots,lRoots)
       CALL CMSMatRot(RTmp,ScanA(Iter),I1,I2,lRoots)
       ScanS(Iter)=CalcNSumVee(RTmp,DDg)
C       IF(I2.eq.1) write(6,*) Iter,ScanA(Iter),ScanS(Iter)
      END DO

      IMax=RMax(ScanS,31)

      Iter=0
      IterMax=100
      SumOld=ScanS(IMax)
      Angles(2)=ScanA(IMax)
      DO WHILE(.not.Converged)
       Iter=Iter+1
       Angles(1)=Angles(2)-stepsize
       Angles(3)=Angles(2)+stepsize
       Do iA=1,3
        CALL Copy2DMat(RTmp,RotMat,lRoots,lRoots)
        CALL CMSMatRot(RTmp,Angles(iA),I1,I2,lRoots)
        Sums(iA)=CalcNSumVee(RTmp,DDg)
       End Do
       CALL CMSFitTrigonometric(Angles,Sums)
       CALL Copy2DMat(RTmp,RotMat,lRoots,lRoots)
       CALL CMSMatRot(RTmp,Angles(4),I1,I2,lRoots)
       Sums(4)=CalcNSumVee(RTmp,DDg)
       IF(ABS(Sums(4)-SumOld).lt.Threshold) THEN
        Converged=.true.
        Angle=Angles(4)
        CALL CMSMatRot(RotMat,Angle,I1,I2,lRoots)
        SumVee=CalcNSumVee(RotMat,DDg)
C        write(6,'(A,I3,A)')
C     &'Convergence reached after ',Iter,' micro cycles'
       ELSE
        If(Iter.eq.IterMax) Then
         Converged=.true.
        write(6,'(A,I3,A)')
     &'No convergence reached after ',Iter,' micro cycles'
        Else
         Angles(2)=Angles(4)
         SumOld=Sums(4)
        End If
       END IF
      END DO
      CALL mma_deallocate(Angles)
      CALL mma_deallocate(Sums)
      CALL mma_deallocate(ScanA)
      CALL mma_deallocate(ScanS)
      CALL mma_deallocate(RTmp)

      RETURN
      END SUBROUTINE
************************************************************************
************************************************************************
      Function RMax(A,N)
       INTEGER N,RMax
       Real*8,DIMENSION(N)::A
       INTEGER I
       RMax=1
       DO I=2,N
        IF (A(I).gt.A(RMax)) RMax=I
       END DO
       RETURN
      End Function

************************************************************************
      Function CalcNSumVee(RotMat,DDg)
      use stdalloc, only : mma_allocate, mma_deallocate
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      Real*8,DIMENSION(lRoots,lRoots,lRoots,lRoots)::DDG
      Real*8,DIMENSION(lroots,lroots)::RotMat
      Real*8,DIMENSION(:),Allocatable::Vee
      Real*8 CalcNSumVee
      INTEGER IState

      CALL mma_allocate(Vee,lRoots)
      CalcNSumVee=0.0d0
      CALL CalcVee(Vee,RotMat,DDg)
      DO IState=1,lRoots
       CalcNSumVee=CalcNSumVee+Vee(IState)
      END DO
      CALL mma_deallocate(Vee)
      RETURN
      END Function
************************************************************************
************************************************************************
      Subroutine Copy2DMat(A,B,NRow,NCol)
      INTEGER NRow,NCol,IRow,ICol
      Real*8,DIMENSION(NRow,NCol)::A,B
      DO ICol=1,NCol
       Do IRow=1,NRow
        A(IRow,ICol)=B(IRow,ICol)
       End Do
      END DO
      RETURN
      END SUBROUTINE
************************************************************************
************************************************************************
      Subroutine CMSMatRot(Mat,A,I,J,N)
      INTEGER I,J,N
      Real*8 A
      Real*8,DIMENSION(N,N)::Mat,TM
      DO K=1,N
       TM(I,K)=Mat(I,K)
       TM(J,K)=Mat(J,K)
      END DO
      DO K=1,N
       Mat(J,K)= cos(A)*TM(J,K)+sin(A)*TM(I,K)
       Mat(I,K)=-sin(A)*TM(J,K)+cos(A)*TM(I,K)
      END DO
      RETURN
      END SUBROUTINE
************************************************************************
************************************************************************
      Subroutine CMSFitTrigonometric(x,y)
      real*8,DIMENSION(4)::x,y
      real*8 s12,s23,c12,c23,d12,d23,k,a,b,c,phi,psi1,psi2,val1,val2
      real*8 atan1
      s12=sin(4.0d0*x(1))-sin(4.0d0*x(2))
      s23=sin(4.0d0*x(2))-sin(4.0d0*x(3))
      c12=cos(4.0d0*x(1))-cos(4.0d0*x(2))
      c23=cos(4.0d0*x(2))-cos(4.0d0*x(3))
      d12=y(1)-y(2)
      d23=y(2)-y(3)
      k=s12/s23
      c=(d12-k*d23)/(c12-k*c23)
      b=(d12-c*c12)/s12
      a=y(1)-b*sin(4.0d0*x(1))-c*cos(4.0d0*x(1))
      phi=atan(b/c)
      atan1=atan(1.0d0)
      psi1=phi/4.0d0
      if(psi1.gt.atan1)then
        psi2=psi1-atan1
       else
        psi2=psi1+atan1
      end if
      val1=b*sin(4.0d0*psi1)+c*cos(4.0d0*psi1)
      val2=b*sin(4.0d0*psi2)+c*cos(4.0d0*psi2)
      if (val1.gt.val2) then
       x(4)=psi1
C       y(4)=val1
      else
       x(4)=psi2
C       y(4)=val2
      end if
      y(4)=a+sqrt(b**2+c**2)
C      write(6,*)a,b,c,x(4),y(4)
      return
      END Subroutine
************************************************************************

************************************************************************
      Subroutine CalcVee(Vee,RMat,DDg)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      Real*8,DIMENSION(lRoots,lRoots,lRoots,lRoots)::DDG
      Real*8,DIMENSION(lroots,lroots)::RMat
      Real*8,DIMENSION(lroots)::Vee
      INTEGER IState,iJ,iK,iL,iM
      DO IState=1,lRoots
       Vee(IState)=0.0d0
       Do iJ=1,lRoots
        Do iK=1,lRoots
         Do iL=1,lRoots
          Do iM=1,lRoots
          Vee(Istate)=Vee(IState)+RMat(IState,iJ)*RMat(IState,iK)*
     &RMat(IState,iL)*RMat(IState,iM)*DDG(iJ,iK,iL,iM)
          End Do
         End Do
        End Do
       End Do
       Vee(IState)=Vee(IState)/2
C       write(6,'(A,I2,A,F10.6)')'The classic coulomb energy for state ',
C     & IState,' is ',Vee(IState)
      END DO
      RETURN
      END SUBROUTINE
************************************************************************
************************************************************************
      Subroutine GetDDgMat(DDg,GDMat,Gtuvx)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      Real*8,DIMENSION(lRoots,lRoots,lRoots,lRoots)::DDG
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::Gtuvx
      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GDMat

      INTEGER iI,iJ,iK,iL,it,iu,iv,ix,iII,iJJ,iKK,iLL
      DO iI=1,lRoots
       DO iJ=1,lRoots
        IF(iJ.gt.iI) THEN
         iJJ=iI
         iII=iJ
        ELSE
         iII=iI
         iJJ=iJ
        END IF
        DO iK=1,lRoots
         DO iL=1,lRoots
          IF(iL.gt.iK) THEN
           iLL=iK
           iKK=iL
          ELSE
           iLL=iL
           iKK=iK
          END IF
          DDG(iI,iJ,iK,iL)=0.0d0
          do it=1,NAC
           do iu=1,NAC
            do iv=1,NAC
             do ix=1,NAC
              DDG(iI,iJ,iK,iL)=DDG(iI,iJ,iK,iL)
     & +GDMat(iII*(iII-1)/2+iJJ,it,iu)*GDMat(iKK*(iKK-1)/2+iLL,iv,ix)
     & *Gtuvx(it,iu,iv,ix)
             end do
            end do
           end do
          end do
         END DO
        END DO
       END DO
      END DO
      RETURN
      End Subroutine
************************************************************************

      Subroutine LoadGtuvx(TUVX,Gtuvx)
* ****************************************************************
* Purpose:                                                       *
* Loading TUVX array to a 4-D tensor.                            *
* Copyied from src/molcas_ci_util/david5.f                       *
* ****************************************************************
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"

      Real*8,DIMENSION(NACPR2)::TUVX
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::Gtuvx
      INTEGER it,iu,iv,ix,ituvx,ixmax
      ituvx=0
      DO it=1,NAC
       Do iu=1,it
        dO iv=1,it
         ixmax=iv
         if (it==iv) ixmax=iu
         do ix=1,ixmax
          ituvx=ituvx+1
          Gtuvx(it,iu,iv,ix)=TUVX(ituvx)
          Gtuvx(iu,it,iv,ix)=TUVX(ituvx)
          Gtuvx(it,iu,ix,iv)=TUVX(ituvx)
          Gtuvx(iu,it,ix,iv)=TUVX(ituvx)
          Gtuvx(iv,ix,it,iu)=TUVX(ituvx)
          Gtuvx(ix,iv,it,iu)=TUVX(ituvx)
          Gtuvx(iv,ix,iu,it)=TUVX(ituvx)
          Gtuvx(ix,iv,iu,it)=TUVX(ituvx)
         end do
        eND dO
       End Do
      END DO
      RETURN
      End Subroutine
************************************************************************
      Subroutine NStateOpt2(RotMat,GDMat,Gtuvx)
      use stdalloc, only : mma_allocate, mma_deallocate
      use CMS, only: CMSNotConverged
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"

      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GDMat
      Real*8,DIMENSION(lRoots,lRoots)::RotMat
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::Gtuvx

      INTEGER IState,JState,NPairs,IPair,ICMSIter
      Real*8 VeeSumOld,VeeSumNew,Threshold,VeeSumChange
      INTEGER,DIMENSION(:,:),Allocatable::StatePair
      Real*8,DIMENSION(:),Allocatable::theta
      Real*8,DIMENSION(:),Allocatable::Vee
      Real*8,DIMENSION(:,:),Allocatable::FRot
      Logical Converged
      Real*8 SumArray
      External SumArray

      CALL mma_allocate(StatePair,LRoots*(LRoots-1)/2,2)
      CALL mma_allocate(theta,LRoots*(LRoots-1)/2)
      CALL mma_allocate(Vee,LRoots)
      CALL mma_allocate(FRot,lRoots,lRoots)
      Threshold=CMSThreshold
      NPairs=lRoots*(lRoots-1)/2
      IPair=0
      DO IState=1,lRoots
       Do JState=1,IState-1
        IPair=IPair+1
        StatePair(IPair,1)=IState
        StatePair(IPair,2)=JState
       End Do
      END DO
      Converged=.false.
      CALL Copy2DMat(FRot,RotMat,lRoots,lRoots)
      CALL RotGDMat(FRot,GDMat)
      CALL CalcVee2(Vee,GDMat,Gtuvx)
      VeeSumOld=SumArray(Vee,lRoots)
      ICMSIter=0
*        write(6,'(6X,I4,8X,F16.8,8X,ES16.4E3)')
*     &  ICMSIter,VeeSumOld,0.0d0
      DO WHILE(.not.Converged)
       Do IPair=1,NPairs
        theta(IPair)=0.0d0
       End Do
       ICMSIter=ICMSIter+1
       CALL ThetaOpt2
     & (FRot,theta,VeeSumChange,StatePair,NPairs,GDMat,Vee,Gtuvx)
       VeeSumNew=VeeSumOld+VeeSumChange
       IF(lRoots.gt.2) THEN
        write(6,'(6X,I4,8X,F16.8,8X,ES16.4E3)')
     &  ICMSIter,VeeSumNew,VeeSumChange
*        CALL RecPrt(' ',' ',Vee,lRoots,1)
*        write(6,*) SumArray(Vee,lRoots)
       ELSE
       write(6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)')
     & ICMSIter,asin(FRot(2,1))/atan(1.0d0)*45.0d0,VeeSumNew
     & ,VeeSumChange
*       CALL RecPrt(' ',' ',Vee,lRoots,1)
*       write(6,*) SumArray(Vee,lRoots)
       END IF
       IF(ABS(VeeSumChange).lt.Threshold) THEN
        If(ICMSIter.ge.ICMSIterMin) Then
         Converged=.true.
         write(6,'(4X,A)')'CONVERGENCE REACHED'
        End If
       ELSE
        if(ICMSIter.ge.ICMSIterMax) then
         Converged=.true.
         CMSNotConverged=.true.
         write(6,'(4X,A)')'NOT CONVERGED AFTER MAX NUMBER OF CYCLES'
         write(6,'(4X,A)')'TEMPORARY ROTATION MATRIX SAVED'
        end if
       END IF
*         Converged=.true.
       VeeSumOld=VeeSumNew
      END DO
      write(6,*)('=',i=1,71)

      CALL Copy2DMat(RotMat,FRot,lRoots,lRoots)
      CALL mma_deallocate(StatePair)
      CALL mma_deallocate(theta)
      CALL mma_deallocate(Vee)
      CALL mma_deallocate(FRot)
      RETURN
      END SUBROUTINE

************************************************************************

      SubRoutine OptOneAngle2(ang,change,R,GD,I1,I2,Vee,G)
      use stdalloc, only : mma_allocate, mma_deallocate
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      Real*8 ang,change
      Integer I1,I2
      Real*8,DIMENSION(lRoots,lRoots)::R
      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GD
      Real*8,DIMENSION(lRoots)::Vee
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::G

      Logical Converged
      INTEGER Itera,Itermax,IA,IMax
      Real*8 Threshold,StepSize,SumOld,Vee1,Vee2,SumOld2
      Real*8,DIMENSION(:),Allocatable::Angles,Sums
      Real*8,DIMENSION(:),Allocatable::ScanA,ScanS

      INTEGER RMax
      External RMax

      CALL mma_allocate(Angles,4)
      CALL mma_allocate(Sums,4)
      CALL mma_allocate(ScanA,31)
      CALL mma_allocate(ScanS,31)

      Converged=.false.
      stepsize=dble(atan(1.0d0))/15
      Threshold=1.0d-8

      Angles(2)=0.0d0
      DO Itera=1,31
       ScanA(Itera)=(Itera-16)*stepsize*2
       CALL
     &SumVeeNew(ScanS(Itera),ScanA(Itera),GD,I1,I2,G,Vee1,Vee2,.false.)
C       IF(I2.eq.1) write(6,*) Iter,ScanA(Iter),ScanS(Iter)
      END DO

      IMax=RMax(ScanS,21)

      Itera=0
      IterMax=100
      SumOld=Vee(I1)+Vee(I2)
      SumOld2=SumOld
      Angles(2)=ScanA(IMax)
      DO WHILE(.not.Converged)
       Itera=Itera+1
       Angles(1)=Angles(2)-stepsize
       Angles(3)=Angles(2)+stepsize
       Do iA=1,3
       CALL SumVeeNew(Sums(iA),Angles(iA),GD,I1,I2,G,Vee1,Vee2,.false.)
       End Do
       CALL CMSFitTrigonometric(Angles,Sums)
       CALL SumVeeNew(Sums(4),Angles(4),GD,I1,I2,G,Vee1,Vee2,.false.)
       change=Sums(4)-SumOld
       IF(ABS(change).lt.Threshold) THEN
        Converged=.true.
        Ang=Angles(4)
        Vee(I1)=Vee1
        Vee(I2)=Vee2
        CALL SumVeeNew(Sums(4),Ang,GD,I1,I2,G,Vee1,Vee2,.true.)
        CALL CMSMatRot(R,Ang,I1,I2,lRoots)
       ELSE
        If(Itera.eq.IterMax) Then
         Converged=.true.
        write(6,'(A,I3,A)')
     &'No convergence reached after ',Itera,' micro cycles'
        Else
         Angles(2)=Angles(4)
         SumOld=Sums(4)
        End If
       END IF
      END DO
      change=Vee(I1)+Vee(I2)-SumOld2
      CALL mma_deallocate(Angles)
      CALL mma_deallocate(Sums)
      CALL mma_deallocate(ScanA)
      CALL mma_deallocate(ScanS)
      RETURN
      End Subroutine

************************************************************************
      Subroutine SumVeeNew(SV,A,GD,I1,I2,G,V1,V2,Update)
      use stdalloc, only : mma_allocate, mma_deallocate
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"

      Real*8 SV,A,V1,V2
      INTEGER I1,I2
      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GD
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::G
      Real*8,DIMENSION(:,:),Allocatable::D11,D22
      Real*8,DIMENSION(:,:,:),Allocatable::D1J,D2J
      Logical Update

      INTEGER t,u,v,x,i11,i22,i12
      INTEGER J,I1J,I2J
      IF(Update) THEN
       CALL mma_allocate(D1J,lRoots,NAC,NAC)
       CALL mma_allocate(D2J,lRoots,NAC,NAC)
*      calculating
       DO J=1,I2-1                           !(J<I2<I1)
        I1J=(I1-1)*I1/2+J
        I2J=(I2-1)*I2/2+J
        Do t=1,NAC
        Do u=1,NAC
         D2J(J,t,u)= cos(A)*GD(I2J,t,u)+sin(A)*GD(I1J,t,u)
         D1J(J,t,u)=-sin(A)*GD(I2J,t,u)+cos(A)*GD(I1J,t,u)
        End Do
        End Do
       END DO
       J=I2                                  !(J=I2<I1)
      i11=(I1+1)*I1/2
      i22=(I2+1)*I2/2
      i12=(I1-1)*I1/2+I2
        Do t=1,NAC
        Do u=1,NAC
         D2J(i2,t,u)=GD(i11,t,u)*sin(A)**2+GD(i22,t,u)*cos(A)**2
     &+cos(A)*sin(A)*(GD(i12,u,t)+GD(i12,t,u))
         D1J(i2,t,u)=cos(A)*sin(A)*(GD(i11,t,u)-GD(i22,t,u))
     &+GD(i12,t,u)*cos(A)**2-GD(i12,u,t)*sin(A)**2
        End Do
        End Do
       DO J=I2+1,I1-1                        !(I2<J<I1)
        I1J=(I1-1)*I1/2+J
        I2J=(J-1)*J/2+I2
        Do t=1,NAC
        Do u=1,NAC
         D2J(J,t,u)= cos(A)*GD(I2J,u,t)+sin(A)*GD(I1J,t,u)
         D1J(J,t,u)=-sin(A)*GD(I2J,u,t)+cos(A)*GD(I1J,t,u)
        End Do
        End Do
       END DO
       J=I1                                  !(I2<J=I1)
      i11=(I1+1)*I1/2
      i22=(I2+1)*I2/2
      i12=(I1-1)*I1/2+I2
        Do t=1,NAC
        Do u=1,NAC
         D1J(i1,t,u)=GD(i11,t,u)*cos(A)**2+GD(i22,t,u)*sin(A)**2
     &   -cos(A)*sin(A)*(GD(i12,t,u)+GD(i12,u,t))
        End Do
        End Do

       DO J=I1+1,lRoots                      !(I2<I1<J)
        I1J=(J-1)*J/2+I1
        I2J=(J-1)*J/2+I2
        Do t=1,NAC
        Do u=1,NAC
         D2J(J,t,u)= cos(A)*GD(I2J,u,t)+sin(A)*GD(I1J,u,t)
         D1J(J,t,u)=-sin(A)*GD(I2J,u,t)+cos(A)*GD(I1J,u,t)
        End Do
        End Do
       END DO
*      updating
       DO J=1,I2-1                           !(J<I2<I1)
        I1J=(I1-1)*I1/2+J
        I2J=(I2-1)*I2/2+J
        Do t=1,NAC
        Do u=1,NAC
         GD(I2J,t,u)=D2J(J,t,u)
         GD(I1J,t,u)=D1J(J,t,u)
        End Do
        End Do
       END DO
       J=I2                                  !(J=I2<I1)
*      i11=(I1+1)*I1/2
      i22=(I2+1)*I2/2
      i12=(I1-1)*I1/2+I2
        Do t=1,NAC
        Do u=1,NAC
         GD(I22,t,u)=D2J(I2,t,u)
         GD(I12,t,u)=D1J(I2,t,u)
        End Do
        End Do
       DO J=I2+1,I1-1                        !(I2<J<I1)
        I1J=(I1-1)*I1/2+J
        I2J=(J-1)*J/2+I2
        Do t=1,NAC
        Do u=1,NAC
         GD(I2J,t,u)=D2J(J,u,t)
         GD(I1J,t,u)=D1J(J,t,u)
        End Do
        End Do
       END DO
       J=I1                                  !(I2<J=I1)
       i11=(I1+1)*I1/2
        Do t=1,NAC
        Do u=1,NAC
         GD(i11,t,u)=D1J(I1,t,u)
        End Do
        End Do

       DO J=I1+1,lRoots                      !(I2<I1<J)
        I1J=(J-1)*J/2+I1
        I2J=(J-1)*J/2+I2
        Do t=1,NAC
        Do u=1,NAC
         GD(I2J,t,u)=D2J(J,u,t)
         GD(I1J,t,u)=D1J(J,u,t)
        End Do
        End Do
       END DO
       Call mma_deallocate(D1J)
       Call mma_deallocate(D2J)
      ELSE
       CALL mma_allocate(D11,NAC,NAC)
       CALL mma_allocate(D22,NAC,NAC)
       i11=(I1+1)*I1/2
       i22=(I2+1)*I2/2
       i12=(I1-1)*I1/2+I2
       V1=0.0d0
       V2=V1
       DO t=1,NAC
        Do u=1,NAC
         D11(t,u)=GD(i11,t,u)*cos(A)**2+GD(i22,t,u)*sin(A)**2
     &   -cos(A)*sin(A)*(GD(i12,t,u)+GD(i12,u,t))
         D22(t,u)=GD(i11,t,u)*sin(A)**2+GD(i22,t,u)*cos(A)**2
     &   +cos(A)*sin(A)*(GD(i12,u,t)+GD(i12,t,u))
        End Do
       END DO
       DO t=1,NAC
       DO u=1,NAC
        Do v=1,NAC
        Do x=1,NAC
         V1=V1+D11(t,u)*D11(v,x)*G(t,u,v,x)
         V2=V2+D22(t,u)*D22(v,x)*G(t,u,v,x)
        End Do
        End Do
       END DO
       END DO
       V1=V1/2.0d0
       V2=V2/2.0d0
       SV=V1+V2
       Call mma_deallocate(D11)
       Call mma_deallocate(D22)
      END IF
      RETURN
      End Subroutine
************************************************************************

************************************************************************


      Subroutine ThetaOpt2(R,theta,deltaQ,SPair,NP,GD,Vee,G)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      INTEGER NP
      Real*8,DIMENSION(NP)::theta
      Real*8 Change,deltaQ
      INTEGER,DIMENSION(NP,2)::SPair
      Real*8,DIMENSION(lroots,lroots)::R
      Real*8,DIMENSION(lroots)::Vee
      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GD
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::G

      INTEGER IP,I,J
      deltaQ=0.0d0
      DO IP=1,NP
       I=SPair(IP,1)
       J=SPair(IP,2)
       CALL OptOneAngle2(theta(iP),change,R,GD,I,J,Vee,G)
       deltaQ=deltaQ+change
      END DO

      DO IP=NP-1,1,-1
       I=SPair(IP,1)
       J=SPair(IP,2)
       CALL OptOneAngle2(theta(iP),change,R,GD,I,J,Vee,G)
       deltaQ=deltaQ+change
      END DO

      RETURN
      END SUBROUTINE

************************************************************************
      Function SumArray(A,N)
      INTEGER N,I
      Real*8,DIMENSION(N)::A
      Real*8 SumArray
      SumArray=0.0d0
      DO I=1, N
       SumArray=SumArray+A(I)
      END DO
      RETURN
      End Function


************************************************************************
      Subroutine CalcVee2(Vee,GD,Gtuvx)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GD
      Real*8,DIMENSION(lRoots)::Vee
      Real*8,DIMENSION(NAC,NAC,NAC,NAC)::Gtuvx

      INTEGER I,t,u,v,x,III
      DO I=1,lRoots
       Vee(I)=0.0d0
       III=I*(I+1)/2
       Do t=1,nac
       Do u=1,nac
       Do v=1,nac
       Do x=1,nac
        Vee(I)=Vee(I)+GD(III,t,u)*GD(III,v,x)*Gtuvx(t,u,v,x)
       End Do
       End Do
       End Do
       End Do
       Vee(I)=Vee(I)/2.0d0
      END DO
      RETURN
      END SUBROUTINE

************************************************************************
      Subroutine RotGDMat(R,GD)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"

      Real*8,DIMENSION(LRoots*(LRoots+1)/2,NAC,NAC)::GD,GD2
      Real*8,DIMENSION(lRoots,lRoots)::R

      INTEGER I,J,K,L,p,q,iIJ,iKL,ip,iq

      DO p=1,nac
      DO q=1,nac
       Do I=1,lRoots
       Do J=1,I
        iIJ=(I-1)*I/2+J
        GD2(iIJ,p,q)=0.0d0
        do K=1,lRoots
        do L=1,lRoots
         IF(K.gt.L) THEN
          iKL=(K-1)*K/2+L
          ip=p
          iq=q
         ELSE
          iKL=(L-1)*L/2+K
          ip=q
          iq=p
         END IF
         GD2(iIJ,p,q)=GD2(iIJ,p,q)+GD(iKL,ip,iq)*R(I,K)*R(J,L)
        end do
        end do
       End Do
       End Do
      END DO
      END DO

      DO p=1,nac
      DO q=1,nac
       Do I=1,lRoots
       Do J=1,I
        iIJ=(I-1)*I/2+J
        GD(iIJ,p,q)=GD2(iIJ,p,q)
       End Do
       End Do
      END DO
      END DO
      RETURN
      END SUBROUTINE

