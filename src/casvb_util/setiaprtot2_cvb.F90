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

subroutine setiaprtot2_cvb(civec,iapr,ixapr,ibpr,ixbpr,npvb,nda,ndb)

implicit real*8(a-h,o-z)
dimension civec(nda,ndb)
dimension iapr(npvb), ixapr(nda+1), ibpr(npvb), ixbpr(ndb+1)

idetvb = 0
ixapr(1) = 1
do ia=1,nda
  do ib=1,ndb
    if (civec(ia,ib) == 1d0) then
      idetvb = idetvb+1
      if (idetvb > npvb) then
        write(6,*) ' Error in setiaprtot!',npvb
        call abend_cvb()
      end if
      iapr(idetvb) = ib
    end if
  end do
  ixapr(ia+1) = idetvb+1
end do

idetvb = 0
ixbpr(1) = 1
do ib=1,ndb
  do ia=1,nda
    if (civec(ia,ib) == 1d0) then
      idetvb = idetvb+1
      if (idetvb > npvb) then
        write(6,*) ' Error in setiaprtot!',npvb
        call abend_cvb()
      end if
      ibpr(idetvb) = ia
    end if
  end do
  ixbpr(ib+1) = idetvb+1
end do

return

end subroutine setiaprtot2_cvb
