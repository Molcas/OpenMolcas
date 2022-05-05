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
      Subroutine TransferMO(MOas,TabMO,mAO,mGrid,nMOs,iAO)
      use nq_Info

!*****Purpose:
!*****Transferring MO information to MOas to be used in dgemm.
!*****It records the MO values on each grid point, too.
!*****But the difference from TransActMO is that the first and
!*****the second elements are the values of the first MO at grid
!*****point 1 and grid point 2.

!*****Input
      INTEGER mAO,mGrid,nMOs,iAO
      Real*8,DIMENSION(mAO,mGrid,nMOs)::TabMO
!*****Output
      Real*8,DIMENSION(mGrid*nOrbt)::MOas

!*****Auxiliary
      INTEGER iIrrep,IOff1,iOff2,iOff3,nCP
      IOff3=0
      DO iIrrep=0,mIrrep-1
       IOff1=OffBasFro(iIrrep)+1
       IOff2=IOff3*mGrid+1
       nCP=mOrb(iIrrep)*mGrid
       CALL DCopy_(nCP,TabMO(iAO,1,IOff1),mAO,MOas(IOff2),1)
       IOff3=IOff3+mOrb(iIrrep)
      END DO
      RETURN
      End Subroutine
