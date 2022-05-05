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
      Subroutine CalcPUVXOff()
      use nq_Info

      INTEGER IOff1,iIrrep,jIrrep,kIrrep,lIrrep,iOrb,jAct,kAct,lAct,    &
     &        ijIrrep,klIrrep,nklAct

      IOff1=0
      DO kIrrep=0,mIrrep-1
       kAct=nAsh(kIrrep)
       Do lIrrep=0,kIrrep
        lAct=nAsh(lIrrep)
        nklAct=kAct*lAct
        If(kIrrep.eq.lIrrep) nklAct=kAct*(kAct+1)/2
        OffVX(lIrrep,kIrrep)=IOff1
        nVX(lIrrep,kIrrep)=nklAct
        IOff1=IOff1+nklAct
       End Do
      END DO
      nVXt=iOff1

      IOff1=0
      DO jIrrep=0,mIrrep-1
       jAct=nAsh(jIrrep)
       Do kIrrep=0,mIrrep-1
        kAct=nAsh(kIrrep)
        do lIrrep=0,kIrrep
         lAct=nAsh(lIrrep)
         nklAct=kAct*lAct
         If(kIrrep.eq.lIrrep) nklAct=kAct*(kAct+1)/2
          OffUVX(lIrrep,kIrrep,jIrrep)=IOff1
          nUVX(lIrrep,kIrrep,jIrrep)=jAct*nklAct
          IOff1=iOff1+jAct*nklAct
        end do
       End Do
      END DO
      nUVXt=IOff1

      IOff1=0
      DO iIrrep=0,mIrrep-1
       OffPUVX(iIrrep)=IOff1
       iOrb=mOrb(iIrrep)
       Do jIrrep=0,mIrrep-1
        jAct=nAsh(jIrrep)
        ijIrrep=1+IEOR(iIrrep,jIrrep)
        Do kIrrep=0,mIrrep-1
         kAct=nAsh(kIrrep)
         do lIrrep=0,kIrrep
          lAct=nAsh(lIrrep)
          klIrrep=1+IEOR(kIrrep,lIrrep)
          IF(ijIrrep.eq.klIrrep) THEN
           iOff1=iOff1+iOrb*nUVX(lIrrep,kIrrep,jIrrep)
          END IF
         end do
        End Do
       End Do
      END DO
      nPot2=IOff1

!      write(6,*)'OffPUVX new method',nPot2,MaxUVX
!      write(6,'(8(I5,2X))')(OffPUVX(iIrrep),iIrrep=0,mIrrep-1)
      RETURN
      End Subroutine
