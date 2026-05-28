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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine CnstInt(Mode,INT1,INT2)

use CHOVEC_IO, only: NVLOC_CHOBATCH
use caspt2_global, only: FIMO_all
use caspt2_module, only: IfChol, NASH, NASHT, NBAS, NBAST, NBSQT, NBTCH, NBTCHES, NFRO, NISH, NSYM
#ifdef _MOLCAS_MPP_
use Para_Info, only: Is_Real_Par
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Mode
real(kind=wp), intent(out) :: INT1(nAshT,nAshT), INT2(nAshT,nAshT,nAshT,nAshT)
integer(kind=iwp) :: iAshI, IB, IB1, IB2, IBEND, IBGRP, IBSTA, iOrb, iSym, iSymA, iSymB, iSymI, iSymJ, iT, iTU, iU, iV, iVX, iX, &
                     jAshI, jOrb, JSYM, kAshI, lAshI, MXBGRP, MXPIQK, NADDBUF, nBasI, NBGRP, NCHOBUF, nCorI, nFroI, nIshI, nKET, NV
real(kind=wp) :: SCAL
integer(kind=iwp), allocatable :: BGRP(:,:)
real(kind=wp), allocatable :: KET(:), WRK1(:), WRK2(:)
integer(kind=iwp), parameter :: Inactive = 1, Active = 2, Virtual = 3

INT1(:,:) = Zero
Int2(:,:,:,:) = Zero

iSym = 1
nFroI = nFro(iSym)
nIshI = nIsh(iSym)
nCorI = nFroI+nIshI
nBasI = nBas(iSym)
!nOrbI = nOrb(iSym)

call mma_allocate(WRK1,NBSQT,Label='WRK1')
call mma_allocate(WRK2,NBSQT,Label='WRK2')

! --- One-Electron Integral

!! Read H_{\mu \nu}
do jAshI=1,nAsh(iSym)
  INT1(1:nAsh(iSym),jAshI) = INT1(1:nAsh(iSym),jAshI)+FIMO_all(nCorI+nBasI*(nCorI+jAshI-1)+1:nCorI+nBasI*(nCorI+jAshI-1)+nAsh(iSym))
end do

! --- Two-Electron Integral

iSymA = 1
iSymI = 1
iSymB = 1
iSymJ = 1

if (IfChol) then
  !nSh(1:nSym,Inactive) = NISH(1:nSym)
  !nSh(1:nSym,Active) = NASH(1:nSym)
  !nSh(1:nSym,Virtual) = NSSH(1:nSym)
  do JSYM=1,NSYM
    IB1 = NBTCHES(JSYM)+1
    IB2 = NBTCHES(JSYM)+NBTCH(JSYM)

    MXBGRP = IB2-IB1+1
    if (MXBGRP <= 0) cycle
    call mma_allocate(BGRP,2,MXBGRP,Label='BGRP')
    IBGRP = 1
    do IB=IB1,IB2
      BGRP(:,IBGRP) = IB
      IBGRP = IBGRP+1
    end do
    NBGRP = MXBGRP

    call MEMORY_ESTIMATE(JSYM,BGRP,NBGRP,NCHOBUF,MXPIQK,NADDBUF)
    call mma_allocate(KET,NCHOBUF,Label='KETBUF')
    do IBGRP=1,NBGRP

      IBSTA = BGRP(1,IBGRP)
      IBEND = BGRP(2,IBGRP)

      NV = sum(NVLOC_CHOBATCH(IBSTA:IBEND))

      !! int2(tuvx) = (tu|vx)/2
      !! This can be computed without frozen orbitals
      call Get_Cholesky_Vectors(Active,Active,JSYM,KET,size(KET),nKet,IBSTA,IBEND)

      if (IBGRP == 1) SCAL = Zero
      if (IBGRP /= 1) SCAL = One
      call DGEMM_('N','T',NASH(JSYM)**2,NASH(JSYM)**2,NV,Half,KET,NASH(JSYM)**2,KET,NASH(JSYM)**2,SCAL,INT2,NASH(JSYM)**2)
    end do
    call mma_deallocate(KET)
    call mma_deallocate(BGRP)
  end do
else
  do iAshI=1,nAsh(iSym)
    iOrb = nCorI+iAshI
    do jAshI=1,nAsh(iSym)
      jOrb = nCorI+jAshI

      call Coul(iSymA,iSymI,iSymB,iSymJ,iOrb,jOrb,WRK1,WRK2)
      !! Put in INT2
      do lAshI=1,nAsh(iSym)
        do kAshI=1,nAsh(iSym)
          INT2(iAshI,jAshI,1:nAsh(iSym),lAshI) = INT2(iAshI,jAshI,1:nAsh(iSym),lAshI)+ &
                                                 WRK1(nCorI+nBasT*(nCorI+lAshI-1)+1:nCorI+nBasT*(nCorI+lAshI-1)+nAsh(iSym))*Half
        end do
      end do
    end do
  end do
end if
call mma_deallocate(WRK1)
call mma_deallocate(WRK2)
if (Mode == 0) then
  do IT=1,nAshT
    do iU=1,nAshT
      iTU = iT+nAshT*(iU-1)
      do iV=1,nAshT
        do iX=1,nAshT
          iVX = iV+nAshT*(iX-1)
          if (iVX > iTU) then
            INT2(iT,iU,iV,IX) = INT2(iT,iU,iV,iX)+INT2(iV,iX,iT,iU)
            INT2(iV,iX,iT,iU) = Zero
          end if
        end do
      end do
    end do
  end do
end if

#ifdef _MOLCAS_MPP_
if (is_real_par()) call GADGOP(INT2,nAshT**4,'+')
#endif

return

end subroutine CnstInt
