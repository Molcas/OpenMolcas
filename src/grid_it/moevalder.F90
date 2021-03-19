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
! Copyright (C) 1995, Roland Lindh                                     *
!               2000, Valera Veryazov                                  *
!               2014, Thomas Dresselhaus                               *
!***********************************************************************
      Subroutine MOEvalDer(MOValue,iDir,nMOs,nCoor,CCoor,CMOs,nCMO,DoIt)

      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"
      Real*8 Ccoor(3,nCoor),MOValue(nCoor*nMOs),CMOs(nCMO)
      Integer DoIt(nMOs),mAO
      integer nDrv

      mAO  = 4
      nDrv = 1

      Call GetMem('MOTMP','Allo','Real',iMoTmp,4*nCoor*nMOs)

      Call MOEval(work(iMoTmp),nMOs,nCoor,CCoor,CMOs,nCMO,DoIt,nDrv,mAO)

! iDir = 1 then do dX
! iDir = 2 then do dY
! iDir = 3 then do dZ
      write(6,*) "iDir:",iDir
      if(iDir.gt.0.and.iDir.lt.4) then
        DO I=1,nCoor*nMOs
          IJ=iDir+1+(I-1)*4
          MOValue(I)=Work(iMoTmp-1+IJ)
        END DO
      else ! do gradient
        DO I=1,nCoor*nMOs
          IJX=2+(I-1)*4
          IJY=3+(I-1)*4
          IJZ=4+(I-1)*4
          MOValue(I)=Work(iMoTmp-1+IJX)+                                &
     &               Work(iMoTmp-1+IJY)+                                &
     &               Work(iMoTmp-1+IJZ)
        END DO
      end if
      Call GetMem('MOTMP','Free','Real',iMoTmp,4*nCoor*nMOs)

      Return
      End
