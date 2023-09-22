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

subroutine vbgendet2_cvb(iapr,ixapr,ibpr,ixbpr,iconfs,idetvb,nconf,nconfion,nda,ndb,ndetvb,nel,noe,nalf,nbet,norb,idetavb,idetbvb, &
                         iwrk1,iwrk2)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: ndetvb, iapr(ndetvb), nda, ixapr(nda+1), ibpr(ndetvb), ndb, ixbpr(ndb+1), nconf, noe, iconfs(noe,nconf), &
                     idetvb(ndetvb), nel, nconfion(0:nel), nalf, nbet, norb, idetavb(ndetvb), idetbvb(ndetvb), iwrk1(ndetvb), &
                     iwrk2(ndetvb)
integer(kind=iwp) :: i, j
logical(kind=iwp), parameter :: debug = .false.

if (debug) then
  write(u6,*) ' Generate determinant information :'
  write(u6,*) ' ----------------------------------'
end if

! vbgenabdet gives all VB alpha and beta strings, to use in CASSCF
! space we sort in A/B strings to get increasing order:
call vbgenabdet_cvb(idetavb,idetbvb,iconfs,nconf,nconfion,ndetvb,nel,noe,nalf,nbet,norb)

call sortindxi_cvb(ndetvb,idetbvb,idetvb)
do i=1,ndetvb
  iwrk1(i) = idetbvb(idetvb(i))
  ibpr(i) = idetavb(idetvb(i))
end do
ixbpr(1) = 1
do i=1,ndb
  do j=ixbpr(i),ndetvb
    if (iwrk1(j) /= i) exit
  end do
  ixbpr(i+1) = j
end do
do i=1,ndb
  call sortindxi_cvb(ixbpr(i+1)-ixbpr(i),ibpr(ixbpr(i)),iwrk2)
  do j=1,ixbpr(i+1)-ixbpr(i)
    iwrk1(j) = ibpr(iwrk2(j)+ixbpr(i)-1)
  end do
  call imove_cvb(iwrk1,ibpr(ixbpr(i)),ixbpr(i+1)-ixbpr(i))
end do
if (debug) then
  write(u6,*) ' ixbpr='
  write(u6,'(10i6)') ixbpr
  write(u6,*) ' ibpr='
  write(u6,'(10i6)') ibpr
end if

call sortindxi_cvb(ndetvb,idetavb,idetvb)
do i=1,ndetvb
  iwrk1(i) = idetavb(idetvb(i))
  iapr(i) = idetbvb(idetvb(i))
end do
ixapr(1) = 1
do i=1,nda
  do j=ixapr(i),ndetvb
    if (iwrk1(j) /= i) exit
  end do
  ixapr(i+1) = j
end do
do i=1,nda
  call sortindxi_cvb(ixapr(i+1)-ixapr(i),iapr(ixapr(i)),iwrk2)
  do j=1,ixapr(i+1)-ixapr(i)
    iwrk1(j) = iapr(iwrk2(j)+ixapr(i)-1)
  end do
  call imove_cvb(iwrk1,iapr(ixapr(i)),ixapr(i+1)-ixapr(i))
  do j=1,ixapr(i+1)-ixapr(i)
    iwrk1(j) = idetvb(iwrk2(j)+ixapr(i)-1)
  end do
  call imove_cvb(iwrk1,idetvb(ixapr(i)),ixapr(i+1)-ixapr(i))
end do
if (debug) then
  write(u6,*) ' ixapr='
  write(u6,'(10i6)') ixapr
  write(u6,*) ' iapr='
  write(u6,'(10i6)') iapr
end if
do i=1,ndetvb
  idetavb(idetvb(i)) = i
end do
call imove_cvb(idetavb,idetvb,ndetvb)
if (debug) then
  write(u6,*) ' idetvb='
  write(u6,'(10i6)') idetvb
end if

return

end subroutine vbgendet2_cvb
