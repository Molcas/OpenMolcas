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

subroutine diis(wrk,wrksize,diispointt,diispointr,key)
! 1) increment key
! 2) if key >=firstext do:
! Tn = DIIS (previous cycext)
! Tn=(T21,T22,T23,T13,T14)
! if key < firstext
! Tn=Tn(prev)
!
! diispointt - pointer of T stack (I)
! diispointr - pointer of R stack (I)
! key        - manipulation key (I/O)

use ccsd_global, only: cycext, firstext, fullprint, t13, t14, t21, t22, t23, v1, v2, v3, v4
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: wrksize, diispointt(4), diispointr(4), key
real(kind=wp) :: wrk(wrksize)
integer(kind=iwp) :: i, j, lun1, nhelp, rc
real(kind=wp) :: cdiis(4), rdiis1(4,4)

!ulf
do i=1,4
  cdiis(i) = Zero
  do j=1,4
    rdiis1(i,j) = Zero
  end do
end do
!1 increment key
key = key+1

if (key < firstext) then

  ! get Tn from last position

  ! get lun number
  lun1 = diispointt(1)
  ! rewind lun1
  call filemanager(2,lun1,rc)
  ! T2aaaa
  call getmediate(wrk,wrksize,lun1,t21%pos0,t21%d,t21%i,rc)
  ! T2bbbb
  call getmediate(wrk,wrksize,lun1,t22%pos0,t22%d,t22%i,rc)
  ! T2abab
  call getmediate(wrk,wrksize,lun1,t23%pos0,t23%d,t23%i,rc)
  ! T1aa
  call getmediate(wrk,wrksize,lun1,t13%pos0,t13%d,t13%i,rc)
  ! T1bb
  call getmediate(wrk,wrksize,lun1,t14%pos0,t14%d,t14%i,rc)
  ! rewind lun1
  call filemanager(2,lun1,rc)

  return
end if

!2.1 make overlap matrix

!2.1.1 rewind R-files
call diisrf(diispointr,cycext)

!2.1.2.1 read R-T21
call diisra(wrk,wrksize,diispointr,cycext,v1%d,v1%i,v1%pos0,v2%d,v2%i,v2%pos0,v3%d,v3%i,v3%pos0,v4%d,v4%i,v4%pos0)

!2.1.2.2 add overlap mtx
call diish1(wrk,wrksize,4,rdiis1,v1%d,v2%d,v3%d,v4%d,v1%i,v2%i,v3%i,v4%i,cycext,1)

!2.1.3.1 read R-T22
call diisra(wrk,wrksize,diispointr,cycext,v1%d,v1%i,v1%pos0,v2%d,v2%i,v2%pos0,v3%d,v3%i,v3%pos0,v4%d,v4%i,v4%pos0)

!2.1.3.2 add overlap mtx
call diish1(wrk,wrksize,4,rdiis1,v1%d,v2%d,v3%d,v4%d,v1%i,v2%i,v3%i,v4%i,cycext,0)

!2.1.4.1 read R-T23
call diisra(wrk,wrksize,diispointr,cycext,v1%d,v1%i,v1%pos0,v2%d,v2%i,v2%pos0,v3%d,v3%i,v3%pos0,v4%d,v4%i,v4%pos0)

!2.1.4.2 add overlap mtx
call diish1(wrk,wrksize,4,rdiis1,v1%d,v2%d,v3%d,v4%d,v1%i,v2%i,v3%i,v4%i,cycext,0)

!2.1.5.1 read R-T13
call diisra(wrk,wrksize,diispointr,cycext,v1%d,v1%i,v1%pos0,v2%d,v2%i,v2%pos0,v3%d,v3%i,v3%pos0,v4%d,v4%i,v4%pos0)

!2.1.5.2 add overlap mtx
call diish1(wrk,wrksize,2,rdiis1,v1%d,v2%d,v3%d,v4%d,v1%i,v2%i,v3%i,v4%i,cycext,0)

!2.1.6.1 read R-T14
call diisra(wrk,wrksize,diispointr,cycext,v1%d,v1%i,v1%pos0,v2%d,v2%i,v2%pos0,v3%d,v3%i,v3%pos0,v4%d,v4%i,v4%pos0)

!2.1.6.2 add overlap mtx
call diish1(wrk,wrksize,2,rdiis1,v1%d,v2%d,v3%d,v4%d,v1%i,v2%i,v3%i,v4%i,cycext,0)

!2.2.1 calc DIIS coefficients
call diish2(rdiis1,cycext,cdiis)

!2.2.2 write DIIS coefficients
if (fullprint > 1) write(u6,'(6X,A,4(F9.5,2X))') 'DIIS coefficients   :',(cdiis(nhelp),nhelp=1,cycext)

!2.3 make new vector

!2.3.1 rewind T-files
call diisrf(diispointt,cycext)

!2.3.2.1 read T21
call diisra(wrk,wrksize,diispointt,cycext,v1%d,v1%i,v1%pos0,v2%d,v2%i,v2%pos0,v3%d,v3%i,v3%pos0,v4%d,v4%i,v4%pos0)

!2.3.2.2 make new T21
call diish3(wrk,wrksize,t21%d,v1%d,v2%d,v3%d,v4%d,cdiis,cycext)

!2.3.3.1 read T22
call diisra(wrk,wrksize,diispointt,cycext,v1%d,v1%i,v1%pos0,v2%d,v2%i,v2%pos0,v3%d,v3%i,v3%pos0,v4%d,v4%i,v4%pos0)

!2.3.3.2 make new T22
call diish3(wrk,wrksize,t22%d,v1%d,v2%d,v3%d,v4%d,cdiis,cycext)

!2.3.4.1 read T23
call diisra(wrk,wrksize,diispointt,cycext,v1%d,v1%i,v1%pos0,v2%d,v2%i,v2%pos0,v3%d,v3%i,v3%pos0,v4%d,v4%i,v4%pos0)

!2.3.4.2 make new T23
call diish3(wrk,wrksize,t23%d,v1%d,v2%d,v3%d,v4%d,cdiis,cycext)

!2.3.5.1 read T13
call diisra(wrk,wrksize,diispointt,cycext,v1%d,v1%i,v1%pos0,v2%d,v2%i,v2%pos0,v3%d,v3%i,v3%pos0,v4%d,v4%i,v4%pos0)

!2.3.5.2 make new T13
call diish3(wrk,wrksize,t13%d,v1%d,v2%d,v3%d,v4%d,cdiis,cycext)

!2.3.6.1 read T14
call diisra(wrk,wrksize,diispointt,cycext,v1%d,v1%i,v1%pos0,v2%d,v2%i,v2%pos0,v3%d,v3%i,v3%pos0,v4%d,v4%i,v4%pos0)

!2.3.6.2 make new T14
call diish3(wrk,wrksize,t14%d,v1%d,v2%d,v3%d,v4%d,cdiis,cycext)

return

end subroutine diis
