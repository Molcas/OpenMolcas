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
use nq_Info, only: IOff_Ash, IOff_BasAct, mIrrep, nAsh, NASHT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nPMO3p, mAO, mGrid, nMOs
real(kind=wp), intent(out) :: P2MOCube(NASHT,mGrid), P2MOCubex(NASHT,nPMO3p), P2MOCubey(NASHT,nPMO3p), P2MOCubez(NASHT,nPMO3p), &
                              MOs(NASHT,mGrid), MOx(NASHT,mGrid), MOy(NASHT,mGrid), MOz(NASHT,mGrid)
real(kind=wp), intent(in) :: TabMO(mAO,mGrid,nMOs), P2Unzip(NASHT,NASHT,NASHT,NASHT)
logical(kind=iwp), intent(in) :: do_grad
integer(kind=iwp) :: IIrrep, iGrid, IOff1, IOff2, NASHT2, NASHT3
logical(kind=iwp) :: lftGGA
real(kind=wp), allocatable :: P2MO1(:), P2MOSquare(:)

lftGGA = .false.
if (lft .and. lGGA) lftGGA = .true.
do iGrid=1,mGrid
  do iIrrep=0,mIrrep-1
    IOff1 = IOff_Ash(iIrrep)
    IOff2 = IOff_BasAct(iIrrep)
    MOs(IOff1+1:IOff1+nAsh(iIrrep),iGrid) = TabMO(1,iGrid,IOff2+1:IOff2+nAsh(iIrrep))
  end do
end do

if (lGGA) then
  do iGrid=1,mGrid
    do iIrrep=0,mIrrep-1
      IOff1 = IOff_Ash(iIrrep)
      IOff2 = IOff_BasAct(iIrrep)
      MOx(IOff1+1:IOff1+nAsh(iIrrep),iGrid) = TabMO(2,iGrid,IOff2+1:IOff2+nAsh(iIrrep))
      MOy(IOff1+1:IOff1+nAsh(iIrrep),iGrid) = TabMO(3,iGrid,IOff2+1:IOff2+nAsh(iIrrep))
      MOz(IOff1+1:IOff1+nAsh(iIrrep),iGrid) = TabMO(4,iGrid,IOff2+1:IOff2+nAsh(iIrrep))
    end do
  end do
end if

NASHT2 = NASHT**2
NASHT3 = NASHT2*NASHT
call mma_allocate(P2MO1,NASHT3,label='P2MO1')
call mma_allocate(P2MOSquare,NASHT2,label='P2MOSquare')
do iGrid=1,mGrid

  !call RecPrt('MOs array','(10(F9.5,1X))',MOs(:,iGrid),1,NASHT)

  !call RecPrt('2RDM array','(10(F9.5,1X))',P2Unzip,NASHT3,NASHT)

  call DGEMM_('T','N',NASHT3,1,NASHT,One,P2UnZip,NASHT,MOs(:,iGrid),NASHT,Zero,P2MO1,NASHT3)

  !call RecPrt('P2MO1 array','(10(F9.5,1X))',P2MO1,NASHT2,NASHT)

  call DGEMM_('T','N',NASHT2,1,NASHT,One,P2MO1,NASHT,MOs(:,iGrid),NASHT,Zero,P2MOSquare,NASHT2)

  !call RecPrt('P2MOSquare array','(10(F9.5,1X))',P2MOSquare,NASHT,NASHT)

  call DGEMM_('T','N',NASHT,1,NASHT,One,P2MOSquare,NASHT,MOs(:,iGrid),NASHT,Zero,P2MOCube(:,iGrid),NASHT)

  if (lftGGA .and. Do_Grad) then
    call DGEMM_('T','N',NASHT,1,NASHT,One,P2MOSquare,NASHT,MOx(:,iGrid),NASHT,Zero,P2MOCubex(:,iGrid),NASHT)
    call DGEMM_('T','N',NASHT,1,NASHT,One,P2MOSquare,NASHT,MOy(:,iGrid),NASHT,Zero,P2MOCubey(:,iGrid),NASHT)
    call DGEMM_('T','N',NASHT,1,NASHT,One,P2MOSquare,NASHT,MOz(:,iGrid),NASHT,Zero,P2MOCubez(:,iGrid),NASHT)

    call DGEMM_('T','N',NASHT2,1,NASHT,One,P2MO1,NASHT,MOx(:,iGrid),NASHT,Zero,P2MOSquare,NASHT2)
    call DGEMM_('T','N',NASHT,1,NASHT,Two,P2MOSquare,NASHT,MOs(:,iGrid),NASHT,One,P2MOCubex(:,iGrid),NASHT)

    call DGEMM_('T','N',NASHT2,1,NASHT,One,P2MO1,NASHT,MOy(:,iGrid),NASHT,Zero,P2MOSquare,NASHT2)
    call DGEMM_('T','N',NASHT,1,NASHT,Two,P2MOSquare,NASHT,MOs(:,iGrid),NASHT,One,P2MOCubey(:,iGrid),NASHT)

    call DGEMM_('T','N',NASHT2,1,NASHT,One,P2MO1,NASHT,MOz(:,iGrid),NASHT,Zero,P2MOSquare,NASHT2)
    call DGEMM_('T','N',NASHT,1,NASHT,Two,P2MOSquare,NASHT,MOs(:,iGrid),NASHT,One,P2MOCubez(:,iGrid),NASHT)
  end if

  !call RecPrt('P2MOCube','(10(F9.5,1X))',P2MOCube(:,iGrid),1,NASHT)
end do
call mma_deallocate(P2MO1)
call mma_deallocate(P2MOSquare)

return

end subroutine CalcP2MOCube
