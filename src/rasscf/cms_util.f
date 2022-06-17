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
******************************************************************
* history:                                                       *
* Jie J. Bao, on Apr. 07, 2022, created this file.               *
******************************************************************

* This file contains simple codes called in CMSNewton. Complicated ones
* are written in files with the name as the subroutine name.

      Subroutine PrintCMSIter(iStep,Qnew,Qold,RMat,lRoots)
      use CMS, only: iCMSOpt,CMSScaled,PosHess
      INTEGER iStep,lRoots
      Real*8 Qnew,Qold,Diff
      Real*8 RMat(lRoots**2)

*      write(6,*) 'iteration information'
      Diff=Qnew-Qold
      IF(lRoots.eq.2) THEN
       write(6,'(6X,I4,8X,F6.1,9X,F16.8,5X,ES16.4E3)')
     & iStep,asin(RMat(3))/atan(1.0d0)*45.0d0,Qnew,Diff
      ELSE
       If (iCMSOpt.eq.1) Then
        if(CMSScaled.and.PosHess) then
         write(6,'(6X,I4,8X,F16.8,8X,ES16.4E3,A1)')
     &   iStep, Qnew,Diff,'*^'
        else if(CMSScaled) then
         write(6,'(6X,I4,8X,F16.8,8X,ES16.4E3,A1)')
     &   iStep, Qnew,Diff,'*'
        else if (PosHess) then
         write(6,'(6X,I4,8X,F16.8,8X,ES16.4E3,A1)')
     &   iStep, Qnew,Diff,'^'
        else
         write(6,'(6X,I4,8X,F16.8,8X,ES16.4E3)')
     &   iStep, Qnew,Diff
        end if
       Else
        write(6,'(6X,I4,8X,F16.8,8X,ES16.4E3)')
     &  iStep, Qnew,Diff
       End If
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
      write(6,*) '    CMS INTERMEDIATE-STATE OPTIMIZATION'
      IF(CMSSFile.eq.'XMS') THEN
       write(6,'(5X,A12,8X,A23)')
     &'START MATRX','XMS INTERMEDIATE STATES'
      ELSE
       write(6,'(5X,A12,8X,A23)')
     &'START MATRX',CMSGuessFile
      END IF
      IF(iCMSOpt.eq.1) THEN
      write(6,'(4X,A12,8X,A8)')
     &'OPT ALGO  ','NEWTON'
      ELSE IF(iCMSOpt.eq.2) THEN
      write(6,'(4X,A12,8X,A8)')
     &'OPT ALGO  ','JACOBI'
      END IF
      write(6,'(4X,A12,8X,ES9.2E2)')
     &'THRESHOLD ',CMSThreshold
      write(6,'(4X,A12,8X,I8)')
     &'MAX CYCLES',ICMSIterMax
      write(6,'(4X,A12,8X,I8)')
     &'MIN CYCLES',ICMSIterMin
      write(6,'(6X,A)')
     &'A ^ sign means Q_a-a is at a saddle point.'
      write(6,'(6X,A)')
     &"A * sign means a scaled step is taken in the Newton's method."
      write(6,*)('=',i=1,71)
      IF(lRoots.gt.2) THEN
      write(6,'(4X,A8,2X,2(A16,11X))')
     &'Cycle','Q_a-a','Difference'
      ELSE
      write(6,'(4X,A8,2X,A18,6X,A8,12X,A12)')
     &'Cycle','Rot. Angle (deg.)','Q_a-a','Q_a-a Diff.'
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


      Subroutine GetDiagScr(nScr,Mat,EigVal,nDim)
      INTEGER nScr,nDim,INFO
      Real*8 Mat(nDim**2)
      Real*8 EigVal(nDim)
      Real*8 Scr(2)

      CALL DSYEV_('V','U',nDim,Mat,nDim,EigVal,Scr,-1,INFO)
      NScr=INT(Scr(1))
      RETURN
      End Subroutine



