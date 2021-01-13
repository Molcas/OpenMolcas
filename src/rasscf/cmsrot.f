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
* Copyright (C) 2020, Jie J. Bao                                       *
************************************************************************
      Subroutine CMSRot(TUVX)
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Aug. 06, 2020, created this file.               *
* ****************************************************************
      use stdalloc, only : mma_allocate, mma_deallocate
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.fh"
#include "rasscf_lucia.fh"

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

      CALL ReadMat('ROT_VEC',VecName,RotMat,lroots,lroots,7,16,'T')

      CALL LoadGtuvx(TUVX,Gtuvx)

      CALL GetGDMat(GDMat)

      CALL GetDDgMat(DDg,GDMat,Gtuvx)

      CALL NStateOpt(RotMat,DDg)

      VecName='CMS-PDFT'
      CALL PrintMat('ROT_VEC',VecName,RotMat,lroots,lroots,7,16,'T')

C     Deallocating Memory
      CALL mma_deallocate(GDMat)
      CALL mma_deallocate(RotMat)
      CALL mma_deallocate(Gtuvx)
      CALL mma_deallocate(DDg)
      RETURN
      End Subroutine

***********************************************************************

***********************************************************************
      Subroutine NStateOpt(RotMat,DDg)
      use stdalloc, only : mma_allocate, mma_deallocate
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.fh"
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
      write(6,*)
      write(6,*)
      write(6,*) '    CMS INTERMEDIATE STATES OPTIMIZATION'
      write(6,'(4X,A12,2X,ES8.2E2)')
     &'THRESHOLD',Threshold
      write(6,'(4X,A12,2X,I8)')
     &'MAX CYCLES',ICMSIterMax
      write(6,'(4X,A12,2X,I8)')
     &'MIN CYCLES',ICMSIterMin
      write(6,*)('=',i=1,71)
      IF(lRoots.gt.2) THEN
      write(6,'(4X,A8,2X,2(A16,11X))')
     &'Cycle','Q_a-a','Difference'
      ELSE
      write(6,'(4X,A8,2X,A18,6X,A8,12X,A12)')
     &'Cycle','Rot. Angle (deg.)','Q_a-a','Q_a-a Diff.'
      END IF
      write(6,*)('-',i=1,71)
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
         write(6,'(4X,A)')'NOT CONVERGED AFTER MAX NUMBER OF CYCLES'
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
***********************************************************************

***********************************************************************
      Subroutine ThetaOpt(FRot,theta,SumVee,StatePair,NPairs,DDg)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.fh"
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
***********************************************************************
***********************************************************************
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
***********************************************************************
***********************************************************************
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

***********************************************************************
      Function CalcNSumVee(RotMat,DDg)
      use stdalloc, only : mma_allocate, mma_deallocate
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.fh"
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
***********************************************************************
***********************************************************************
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
***********************************************************************
***********************************************************************
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
***********************************************************************
***********************************************************************
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
***********************************************************************

***********************************************************************
      Subroutine CalcVee(Vee,RMat,DDg)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.fh"
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
***********************************************************************
***********************************************************************
      Subroutine GetDDgMat(DDg,GDMat,Gtuvx)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.fh"
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
***********************************************************************

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
#include "warnings.fh"

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
***********************************************************************
