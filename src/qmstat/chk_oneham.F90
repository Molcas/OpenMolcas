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

subroutine Chk_OneHam(nBas)

use qmstat_global, only: MxSymQ
use Index_Functions, only: nTri_Elem
use OneDat, only: sNoNuc, sNoOri
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nBas(MxSymQ)
integer(kind=iwp) :: iComp, iopt, irc, iSmLbl, Lu_One, nBT
real(kind=wp) :: dNorm
character(len=8) :: Label_Pure, Label_Read
real(kind=wp), allocatable :: OneP(:), OneR(:)
integer(kind=iwp), external :: IsFreeUnit
real(kind=wp), external :: dnrm2_

Lu_One = IsFreeUnit(49)
Label_Read = 'OneHam  '
Label_Pure = 'OneHam 0'
nBT = nTri_Elem(nBas(1))
iopt = 0
call OpnOne(irc,iopt,'ONEINT',Lu_One)
call mma_allocate(OneR,nBT,label='Read')
call mma_allocate(OneP,nBT,label='Pure')

irc = -1
iopt = ibset(ibset(0,sNoOri),sNoNuc)
iSmLbl = 0
iComp = 1
call RdOne(irc,iopt,Label_Read,iComp,OneR,iSmLbl)
irc = -1
call RdOne(irc,iopt,Label_Pure,iComp,OneP,iSmLbl)
call ClsOne(irc,Lu_One)

OneP(:) = OneP-OneR

dNorm = dnrm2_(nBT,OneP,1)

if (dNorm > 1.0e-8_wp) then
  write(u6,*)
  write(u6,*)
  write(u6,*) ' WARNING!'
  write(u6,*)
  write(u6,*) '   Your one-electron hamiltonian is not purely vacuum. This means that the Hamiltonian'
  write(u6,*) '   in QmStat can be contaminated. Is this intentional? If not, then make sure that the ONEINT'
  write(u6,*) '   file comes directly from a Seward calculation without any calls from'
  write(u6,*) '   FFPT (or similar) in between.'
  write(u6,*)
  write(u6,*)
end if

call mma_deallocate(OneP)
call mma_deallocate(OneR)

return

end subroutine Chk_OneHam
