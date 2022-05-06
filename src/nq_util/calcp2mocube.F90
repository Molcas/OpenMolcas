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
subroutine CalcP2MOCube(P2MOCube,P2MOCubex,P2MOCubey,P2MOCubez,nPMO3p,MOs,MOx,MOy,MOz,TabMO,P2Unzip,mAO,mGrid,nMOs,do_grad)

use nq_pdft, only: lft, lGGA
use nq_Info

implicit real*8(A-H,O-Z)
#include "stdalloc.fh"
! Input
integer mAO, mGrid, nMOs, nPMO3p
real*8, dimension(mAO,mGrid,nMOs) :: TabMO
real*8, dimension(NASHT4) :: P2Unzip
logical do_grad
! Output
real*8, dimension(mGrid*NASHT) :: P2MOCube, MOs, MOx, MOy, MOz
real*8, dimension(nPMO3p) :: P2MOCubex, P2MOCubey, P2MOCubez
! Auxiliary
integer iOff1, IOff2, IOff3, IIrrep, nGridPi, NASHT2, NASHT3, icount
real*8, dimension(NASHT**3) :: P2MO1
real*8, dimension(NASHT**2) :: P2MOSquare
logical lftGGA

lftGGA = .false.
if (lft .and. lGGA) lftGGA = .true.
nGridPi = mAO*mGrid
do iGrid=1,mGrid
  IOff1 = (iGrid-1)*NASHT
  do iIrrep=0,mIrrep-1
    IOff2 = IOff_Ash(iIrrep)+1
    IOff3 = IOff_BasAct(iIrrep)+1
    call DCopy_(nAsh(iIrrep),TabMO(1,iGrid,IOff3),nGridPi,MOs(IOff1+IOff2),1)
    do icount=1,nAsh(iIrrep)
    end do
  end do
end do

if (lGGA) then
  do iGrid=1,mGrid
    IOff1 = (iGrid-1)*NASHT
    do iIrrep=0,mIrrep-1
      IOff2 = IOff_Ash(iIrrep)+1
      IOff3 = IOff_BasAct(iIrrep)+1
      call DCopy_(nAsh(iIrrep),TabMO(2,iGrid,IOff3),nGridPi,MOx(IOff1+IOff2),1)
      call DCopy_(nAsh(iIrrep),TabMO(3,iGrid,IOff3),nGridPi,MOy(IOff1+IOff2),1)
      call DCopy_(nAsh(iIrrep),TabMO(4,iGrid,IOff3),nGridPi,MOz(IOff1+IOff2),1)
    end do
  end do
end if

NASHT2 = NASHT**2
NASHT3 = NASHT2*NASHT
do iGrid=1,mGrid
  IOff1 = (iGrid-1)*NASHT+1

  !call RecPrt('MOs array','(10(F9.5,1X))',MOs(IOff1),1,NASHT)

  !call RecPrt('2RDM array','(10(F9.5,1X))',P2Unzip,NASHT3,NASHT)

  call DGEMM_('T','N',NASHT3,1,NASHT,1.0d0,P2UnZip,NASHT,MOs(IOff1),NASHT,0.0d0,P2MO1,NASHT3)

  !call RecPrt('P2MO1 array','(10(F9.5,1X))',P2MO1,NASHT2,NASHT)

  call DGEMM_('T','N',NASHT2,1,NASHT,1.0d0,P2MO1,NASHT,MOs(IOff1),NASHT,0.0d0,P2MOSquare,NASHT2)

  !call RecPrt('P2MOSquare array','(10(F9.5,1X))',P2MOSquare,NASHT,NASHT)

  call DGEMM_('T','N',NASHT,1,NASHT,1.0d0,P2MOSquare,NASHT,MOs(IOff1),NASHT,0.0d0,P2MOCube(iOff1),NASHT)

  if (lftGGA .and. Do_Grad) then
    call DGEMM_('T','N',NASHT,1,NASHT,1.0d0,P2MOSquare,NASHT,MOx(IOff1),NASHT,0.0d0,P2MOCubex(iOff1),NASHT)
    call DGEMM_('T','N',NASHT,1,NASHT,1.0d0,P2MOSquare,NASHT,MOy(IOff1),NASHT,0.0d0,P2MOCubey(iOff1),NASHT)
    call DGEMM_('T','N',NASHT,1,NASHT,1.0d0,P2MOSquare,NASHT,MOz(IOff1),NASHT,0.0d0,P2MOCubez(iOff1),NASHT)

    call DGEMM_('T','N',NASHT2,1,NASHT,1.0d0,P2MO1,NASHT,MOx(IOff1),NASHT,0.0d0,P2MOSquare,NASHT2)
    call DGEMM_('T','N',NASHT,1,NASHT,2.0d0,P2MOSquare,NASHT,MOs(IOff1),NASHT,1.0d0,P2MOCubex(iOff1),NASHT)

    call DGEMM_('T','N',NASHT2,1,NASHT,1.0d0,P2MO1,NASHT,MOy(IOff1),NASHT,0.0d0,P2MOSquare,NASHT2)
    call DGEMM_('T','N',NASHT,1,NASHT,2.0d0,P2MOSquare,NASHT,MOs(IOff1),NASHT,1.0d0,P2MOCubey(iOff1),NASHT)

    call DGEMM_('T','N',NASHT2,1,NASHT,1.0d0,P2MO1,NASHT,MOz(IOff1),NASHT,0.0d0,P2MOSquare,NASHT2)
    call DGEMM_('T','N',NASHT,1,NASHT,2.0d0,P2MOSquare,NASHT,MOs(IOff1),NASHT,1.0d0,P2MOCubez(iOff1),NASHT)
  end if

  !call RecPrt('P2MOCube','(10(F9.5,1X))',P2MOCube(IOff1),1,NASHT)
end do

return

end subroutine CalcP2MOCube
