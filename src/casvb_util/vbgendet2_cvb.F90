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

implicit real*8(a-h,o-w,y-z),integer(x)
dimension iapr(ndetvb), ixapr(nda+1), ibpr(ndetvb), ixbpr(ndb+1)
dimension iconfs(noe,nconf), idetvb(ndetvb)
dimension nconfion(0:nel)
dimension idetavb(ndetvb), idetbvb(ndetvb)
dimension iwrk1(ndetvb), iwrk2(ndetvb)
logical debug
data debug/.false./

if (debug) then
  write(6,*) ' Generate determinant information :'
  write(6,*) ' ----------------------------------'
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
    if (iwrk1(j) /= i) goto 400
  end do
  j = ndetvb+1
400 ixbpr(i+1) = j
end do
do i=1,ndb
  call sortindxi_cvb(ixbpr(i+1)-ixbpr(i),ibpr(ixbpr(i)),iwrk2)
  do j=1,ixbpr(i+1)-ixbpr(i)
    iwrk1(j) = ibpr(iwrk2(j)+ixbpr(i)-1)
  end do
  call imove_cvb(iwrk1,ibpr(ixbpr(i)),ixbpr(i+1)-ixbpr(i))
end do
if (debug) then
  write(6,*) ' ixbpr='
  write(6,'(10i6)') ixbpr
  write(6,*) ' ibpr='
  write(6,'(10i6)') ibpr
end if

call sortindxi_cvb(ndetvb,idetavb,idetvb)
do i=1,ndetvb
  iwrk1(i) = idetavb(idetvb(i))
  iapr(i) = idetbvb(idetvb(i))
end do
ixapr(1) = 1
do i=1,nda
  do j=ixapr(i),ndetvb
    if (iwrk1(j) /= i) goto 1000
  end do
  j = ndetvb+1
1000 ixapr(i+1) = j
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
  write(6,*) ' ixapr='
  write(6,'(10i6)') ixapr
  write(6,*) ' iapr='
  write(6,'(10i6)') iapr
end if
do i=1,ndetvb
  idetavb(idetvb(i)) = i
end do
call imove_cvb(idetavb,idetvb,ndetvb)
if (debug) then
  write(6,*) ' idetvb='
  write(6,'(10i6)') idetvb
end if

return

end subroutine vbgendet2_cvb
