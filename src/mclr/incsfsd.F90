!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine InCSFSD(iState,State_sym)

use Symmetry_Info, only: Mul
use Str_Info, only: CNSM
use MCLR_Data, only: i1, iALLO, iAnders, lConf, llDet, LuCSF2SD
use stdalloc, only: mma_allocate
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: iState, State_sym
integer(kind=iwp) :: i, iad, iAdr, idum(1), iSym

! Place pointer

iSym = Mul(State_Sym,iState)

if ((isym == 1) .and. (i1 == 1)) return
if (isym == iAnders) return

iAdr = 2
if (iSym == 1) iAdr = 1
iad = 0
do i=1,iState-1
  call iDafile(LUCSF2SD,0,idum,lldet,iad)
  call iDafile(LUCSF2SD,0,idum,lconf,iad)
end do

if (iSym /= 1) then
  if (iAnders == -9) then
    call mma_allocate(CNSM(2)%icts,lldet,Label='ICTS')
    call mma_allocate(CNSM(2)%iconf,lConf,Label='ICONF')
    iAllo = 1
  end if
  iAnders = isym
end if
if (iSym == 1) then
  if (i1 == -9) then
    call mma_allocate(CNSM(1)%icts,lldet,Label='ICTS')
    call mma_allocate(CNSM(1)%iconf,lConf,Label='ICONF')
    i1 = 1
  end if
end if

!open(unit=1422,file='det.index') ! yma
!do i=1,lldet
!  write(1422,*)CNSM(iAdr)%icts(i)
!end do
!close(1422)

! calculated from zoo.f, the GUGA number for determinent
call iDafile(LUCSF2SD,2,CNSM(iAdr)%icts,lldet,iad)
call iDafile(LUCSF2SD,2,CNSM(iAdr)%iconf,lconf,iad)

end subroutine InCSFSD
