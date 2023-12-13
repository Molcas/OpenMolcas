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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine mkbiks_cvb()

use casvb_global, only: aikcof, bikcof, ikcoff, ipr, kbasiscvb, nel
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: i2s1, ibeg, iend, ifns, nalf1, ndet, nel1, ntmp
logical(kind=iwp) :: share
character(len=*), parameter :: basis(7) = ['Kotani    ','Serber    ','Rumer     ','Rumer (LT)','projected ','Determ    ', &
                                           'Determ    ']
integer(kind=iwp), external :: ifns_cvb

aikcof(0) = real(kbasiscvb,kind=wp)
bikcof(0) = real(kbasiscvb,kind=wp)
if (kbasiscvb == 6) return

if (ipr(1) >= 1) write(6,6100) trim(basis(kbasiscvb))

share = associated(bikcof,aikcof)
do nel1=0,nel
  do nalf1=0,nel
    do i2s1=0,nel
      ibeg = ikcoff(nel1,nalf1,i2s1)+1
      if (ibeg /= 0) then
        ntmp = (nel1+i2s1)/2
        ifns = ifns_cvb(nel1,ntmp,kbasiscvb)
        call icomb_cvb(nel1,nalf1,ndet)
        iend = ibeg+ndet*ifns-1
        call bikset_cvb(aikcof(ibeg:iend),bikcof(ibeg:iend),nel1,nalf1,i2s1,ndet,ifns,kbasiscvb,share,ipr(1))
      end if
    end do
  end do
end do

return
6100 format(/,' Generate ',a,' spin functions.')

end subroutine mkbiks_cvb
