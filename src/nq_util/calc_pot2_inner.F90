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
! Jie J. Bao, on Dec. 25, 2021, created this file.               *
! ****************************************************************
subroutine Calc_Pot2_Inner(Pot2,mGrid,MOP,MOU,MOV,MOX,lSum)

use nq_Info, only: mIrrep, mOrb, nAsh, nIsh, nOrbt, nPot2, nUVX, nUVXt, nVX, nVXt, OffOrb, OffPUVX, OffUVX, OffVX
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Pot2(nPot2)
integer(kind=iwp), intent(in) :: mGrid
real(kind=wp), intent(in) :: MOP(mGrid,nOrbt), MOU(mGrid,nOrbt), MOV(mGrid,nOrbt), MOX(mGrid,nOrbt)
logical(kind=iwp) :: lSum
! Note: when lSum is .true., calculate P(U'VX+UV'X+UVX'), otherwise calculate PUVX
integer(kind=iwp) :: iOff0, iOff1, iOff2, iOff3, ioffu, iStack, iVX, nnUVX, nporb, pIrrep, puIrrep, u, uIrrep, v, vIrrep, vorb, x, &
                     xIrrep, xMax, xorb
real(kind=wp), allocatable :: MOUVX(:,:), MOVX1(:,:), MOVX2(:,:)

call mma_allocate(MOVX1,mGrid,nVXt)
if (lSum) call mma_allocate(MOVX2,mGrid,nVXt)
call mma_allocate(MOUVX,mGrid,nUVXt)

do vIrrep=0,mIrrep-1
  do xIrrep=0,vIrrep
    do v=1,nAsh(vIrrep)
      vorb = v+OffOrb(vIrrep)+nIsh(vIrrep)
      if (xIrrep == vIrrep) then
        xMax = v
        iStack = v*(v-1)/2
      else
        xMax = nAsh(xIrrep)
        iStack = (v-1)*xMax
      end if
      do x=1,xMax
        xorb = x+OffOrb(xIrrep)+nIsh(xIrrep)
        MOVX1(:,OffVX(xIrrep,vIrrep)+iStack+x) = MOV(:,vorb)*MOX(:,xorb)
      end do
    end do
  end do
end do

if (lSum) then
  do vIrrep=0,mIrrep-1
    do xIrrep=0,vIrrep
      do v=1,nAsh(vIrrep)
        vorb = v+OffOrb(vIrrep)+nIsh(vIrrep)
        if (xIrrep == vIrrep) then
          xMax = v
          iStack = v*(v-1)/2
        else
          xMax = nAsh(xIrrep)
          iStack = (v-1)*xMax
        end if
        do x=1,xMax
          xorb = x+OffOrb(xIrrep)+nIsh(xIrrep)
          IOff3 = (OffVX(xIrrep,vIrrep)+iStack+x-1)*mGrid
          MOVX1(:,OffVX(xIrrep,vIrrep)+iStack+x) = MOVX1(:,OffVX(xIrrep,vIrrep)+iStack+x)+MOX(:,vorb)*MOV(:,xorb)
          MOVX2(:,OffVX(xIrrep,vIrrep)+iStack+x) = MOU(:,vorb)*MOV(:,xorb)
        end do
      end do
    end do
  end do
end if

do uIrrep=0,mIrrep-1
  IOffU = OffOrb(uIrrep)+nIsh(uIrrep)
  do vIrrep=0,mIrrep-1
    do xIrrep=0,vIrrep
      do iVX=1,nVX(xIrrep,vIrrep)
        IOff0 = OffUVX(xIrrep,vIrrep,uIrrep)+(iVX-1)*nAsh(uIrrep)
        do u=1,nAsh(uIrrep)
          MOUVX(:,IOff0+u) = MOU(:,iOffU+u)*MOVX1(:,Offvx(xIrrep,vIrrep)+iVX)
        end do
      end do
    end do
  end do
end do

if (lSum) then
  do uIrrep=0,mIrrep-1
    IOffU = OffOrb(uIrrep)+nIsh(uIrrep)
    do vIrrep=0,mIrrep-1
      do xIrrep=0,vIrrep
        do iVX=1,nVX(xIrrep,vIrrep)
          IOff0 = OffUVX(xIrrep,vIrrep,uIrrep)+(iVX-1)*nAsh(uIrrep)
          do u=1,nAsh(uIrrep)
            MOUVX(:,IOff0+u) = MOUVX(:,IOff0+u)+MOX(:,iOffU+u)*MOVX2(:,Offvx(xIrrep,vIrrep)+iVX)
          end do
        end do
      end do
    end do
  end do
end if

! Use dgemm to calculate PUVX at this grid point
do pIrrep=0,mIrrep-1
  nporb = mOrb(pIrrep)
  if (nporb == 0) cycle
  if (nAsh(pIrrep) == 0) cycle
  IOff1 = OffOrb(pIrrep)+1
  IOff2 = OffPUVX(pIrrep)+1
  do uIrrep=0,mIrrep-1
    puIrrep = ieor(pIrrep,uIrrep)
    do vIrrep=0,mIrrep-1
      xIrrep = ieor(puIrrep,vIrrep)
      nnUVX = nUVX(xIrrep,vIrrep,uIrrep)
      if ((xIrrep > vIrrep) .or. (nnUVX == 0)) cycle
      IOff3 = OffUVX(xIrrep,vIrrep,uIrrep)+1
      call DGEMM_('T','N',npOrb,nnUVX,mGrid,One,MOP(:,iOff1:),mGrid,MOUVX(:,IOff3:),mGrid,One,Pot2(iOff2:),npOrb)
      IOff2 = IOff2+nnUVX*npOrb
    end do
  end do
end do

call mma_deallocate(MOVX1)
if (lSum) call mma_deallocate(MOVX2)
call mma_deallocate(MOUVX)

return

end subroutine Calc_Pot2_Inner
