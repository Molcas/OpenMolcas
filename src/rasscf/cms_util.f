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
* Copyright (C) 2022, Jie J. Bao                                       *
************************************************************************
******************************************************************
* history:                                                       *
* Jie J. Bao, on Apr. 07, 2022, created this file.               *
******************************************************************

* This file contains simple codes called in CMSNewton. Complicated ones
* are written in files with the name as the subroutine name.

      Subroutine PrintCMSIter(iStep,Qnew,Qold,RMat,lRoots)
      use CMS, only: iCMSOpt,NPosHess,LargestQaaGrad,NCMSScale
      INTEGER iStep,lRoots
      Real*8 Qnew,Qold,Diff
      Real*8 RMat(lRoots**2)

*      write(6,*) 'iteration information'
      Diff=Qnew-Qold
      IF(iCMSOpt.eq.2) THEN


       If(lRoots.eq.2) Then
        write(6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)')
     &  iStep,asin(RMat(3))/atan(1.0d0)*45.0d0,Qnew,Diff
       Else
         write(6,'(6X,I4,2X,F14.8,2X,ES14.4E3)')
     &   iStep, Qnew,Diff
       End If


      ELSE


C       If(lRoots.eq.2) Then
C        write(6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)')
C     &  iStep,asin(RMat(3))/atan(1.0d0)*45.0d0,Qnew,Diff
C       Else
        if (NCMSScale.gt.0) then
      write(6,'(6X,I4,2X,F14.8,2X,ES12.2E3,2X,I5,2X,ES14.4E3,3X,A3,'//
     &        'I1)')
     &   iStep, Qnew,Diff,nPosHess,LargestQaaGrad,'1E-',NCMSScale
        else
       write(6,'(6X,I4,2X,F14.8,2X,ES12.2E3,2X,I5,2X,ES14.4E3,3X,A3)')
     &   iStep, Qnew,Diff,nPosHess,LargestQaaGrad,'1.0'
        end if
C       End If


      END IF
      RETURN
      End Subroutine
************************************************************************

      Subroutine UnzipTUVX(TUVX,gtuvx,nTUVX)
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      INTEGER nTUVX
      Real*8 gtuvx(nTUVX),TUVX(NACPR2)
      INTEGER it,iu,iv,ix,ituvx,ixmax,
     &        jtuvx,jtuxv,jutvx,jutxv,
     &        jvxtu,jvxut,jxvtu,jxvut,
     &        NAC3,NAC2

*      CALL FZero(gtuvx,nTUVX)

      NAC2=NAC**2
      NAC3=NAC2*NAC
      ituvx=0
      DO it=1,NAC
       Do iu=1,it
        dO iv=1,it
         ixmax=iv
         if (it==iv) ixmax=iu
         do ix=1,ixmax
          ituvx=ituvx+1
          jtuvx=(it-1)*NAC3+(iu-1)*NAC2+(iv-1)*NAC+ix
          jtuxv=(it-1)*NAC3+(iu-1)*NAC2+(ix-1)*NAC+iv
          jutvx=(iu-1)*NAC3+(it-1)*NAC2+(iv-1)*NAC+ix
          jutxv=(iu-1)*NAC3+(it-1)*NAC2+(ix-1)*NAC+iv
          jvxtu=(iv-1)*NAC3+(ix-1)*NAC2+(it-1)*NAC+iu
          jvxut=(iv-1)*NAC3+(ix-1)*NAC2+(iu-1)*NAC+it
          jxvtu=(ix-1)*NAC3+(iv-1)*NAC2+(it-1)*NAC+iu
          jxvut=(ix-1)*NAC3+(iv-1)*NAC2+(iu-1)*NAC+it
          Gtuvx(jtuvx)=TUVX(ituvx)
          Gtuvx(jtuxv)=TUVX(ituvx)
          Gtuvx(jutvx)=TUVX(ituvx)
          Gtuvx(jutxv)=TUVX(ituvx)
          Gtuvx(jvxtu)=TUVX(ituvx)
          Gtuvx(jvxut)=TUVX(ituvx)
          Gtuvx(jxvtu)=TUVX(ituvx)
          Gtuvx(jxvut)=TUVX(ituvx)
         end do
        eND dO
       End Do
      END DO
      RETURN
      End Subroutine
************************************************************************


      Subroutine CMSTail()
      write(6,*)('=',i=1,71)
      RETURN
      End Subroutine
************************************************************************



      Subroutine CMSHeader(CMSSFile,LenCMSS)
      use CMS, only: iCMSOpt, CMSGuessFile
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "WrkSpc.fh"
#include "SysDef.fh"
#include "input_ras.fh"
#include "warnings.h"
      INTEGER LenCMSS
      CHARACTER(len=LenCMSS)::CMSSFile
      write(6,*)
      write(6,*)
      write(6,'(4X,A35)')
     & 'CMS INTERMEDIATE-STATE OPTIMIZATION'
      IF(CMSSFile.eq.'XMS') THEN
       write(6,'(5X,A11,9X,A25)')
     &'START MATRX','XMS INTERMEDIATE STATES'
      ELSE
       write(6,'(5X,A11,9X,A25)')
     &'START MATRX',CMSGuessFile
      END IF
      IF(iCMSOpt.eq.1) THEN
       write(6,'(5X,A8,12X,A25)')
     & 'OPT ALGO','NEWTON'
      ELSE IF(iCMSOpt.eq.2) THEN
       write(6,'(5X,A8,12X,A25)')
     & 'OPT ALGO','JACOBI'
      END IF
      write(6,'(5X,A15,5X,16X,ES9.2E2)')
     &'Q_a-a THRESHOLD',CMSThreshold
      IF(iCMSOpt.eq.1) THEN
        write(6,'(5X,A15,5X,16X,ES9.2E2)')
     &  'GRAD  THRESHOLD',CMSThreshold*1.0d-2
      END IF
      write(6,'(5X,A10,10X,I25)')
     &'MAX CYCLES',ICMSIterMax
      write(6,'(5X,A10,10X,I25)')
     &'MIN CYCLES',ICMSIterMin
      write(6,*)('=',i=1,71)
      IF(iCMSOpt.eq.2) THEN
       If(lRoots.gt.2) Then
       write(6,'(4X,A8,2X,2(A16,11X))')
     & 'Cycle','Q_a-a','Difference'
       Else
       write(6,'(4X,A8,2X,A18,6X,A8,12X,A12)')
     & 'Cycle','Rot. Angle (deg.)','Q_a-a','Q_a-a Diff.'
       End If
      ELSE
       write(6,'(6X,A5,7X,A5,8X,A10,2X,A6,5X,A7,4X,A4)')
     & 'Cycle','Q_a-a','Difference','# Pos.','Largest','Step'
       write(6,'(43X,A7,4X,A8,3X,A6)')
     & 'Hessian','Gradient','Scaled'
      END IF
      write(6,*)('-',i=1,71)

      RETURN
      End Subroutine
************************************************************************



      Subroutine AntiOneDFoil(TwoD,OneD,m,n)
      INTEGER M,N,I,J,iLoc
      Real*8,DIMENSION(m,n)::TwoD
      Real*8,DIMENSION(m*n)::OneD
      iLoc=1
      DO J=1,N
       Do I=1,M
        TwoD(I,J)=OneD(iLoc)
        iLoc=iLoc+1
       End Do
      END DO
      RETURN
      End Subroutine

      Subroutine OneDFoil(OneD,TwoD,m,n)
      INTEGER M,N,I,J,iLoc
      Real*8,DIMENSION(m,n)::TwoD
      Real*8,DIMENSION(m*n)::OneD
      iLoc=1
      DO J=1,N
       Do I=1,M
        OneD(iLoc)=TwoD(I,J)
        iLoc=iLoc+1
       End Do
      END DO
      RETURN
      End Subroutine
************************************************************************



      Subroutine InitRotMat(RotMat,lRoots,CMSSFile,LenCMSS)
      INTEGER LenCMSS,lRoots
      CHARACTER(Len=LenCMSS)::CMSSFile
      Real*8,DIMENSION(lRoots,lRoots)::RotMat
      CHARACTER(Len=16)::ScrChar

      IF(CMSSFile.eq.'XMS') THEN
        CALL ReadMat('ROT_VEC',ScrChar,RotMat,lroots,lroots,7,16,'T')
      ELSE
        CALL ReadMat(CMSSFile ,ScrChar,RotMat,lroots,lroots,LenCMSS,16,
     &               'T')
      END IF
      RETURN
      End Subroutine
************************************************************************

      Subroutine UpdateRotMat(RMat,ExpX,X,lRoots,nSPair)
      INTEGER lRoots,nSPair
      Real*8 X(nSPair)
      Real*8 RMat(lRoots**2),RScr(lRoots**2)
      Real*8 ExpX(lRoots**2)


      CALL ExpMat(ExpX,X,lRoots,nSPair)
      CALL DGEMM_('n','n',lRoots,lRoots,lRoots,1.0d0,RMat,lRoots,
     &                                               ExpX,lRoots,
     &                                         0.0d0,RScr,lRoots)
      CALL DCopy_(lRoots**2,RScr,1,RMat,1)
      RETURN
      End Subroutine





