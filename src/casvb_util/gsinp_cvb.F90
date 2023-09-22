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

subroutine gsinp_cvb(orbs,irdorbs,ip_cvb,nvbinp,kbasiscvb_inp,mxaobf,mxorb,kbasis,strtvb)

use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: mxorb, irdorbs(mxorb), ip_cvb, nvbinp, kbasiscvb_inp, mxaobf, kbasis
real(kind=wp) :: orbs(mxaobf,mxorb), strtvb
#include "WrkSpc.fh"
integer(kind=iwp) :: idum(1), iorb, istr, mouse, mxread, nread
integer(kind=iwp), parameter :: ncmp = 4, ngs = 7
character(len=*), parameter :: guess(ngs) = ['ORB     ','STRUC   ','READ    ','AOBASIS ','MOBASIS ','END     ','ENDGUESS']
integer(kind=iwp), external :: mavailr_cvb, mheapr_cvb
logical(kind=iwp), external :: firsttime_cvb

if (firsttime_cvb()) call touch_cvb('INPGS')
mouse = 1
do
  call fstring_cvb(guess,ngs,istr,ncmp,2)
  if (istr == 1) then
    ! 'ORB'
    call int_cvb(idum,1,nread,0)
    iorb = idum(1)
    if ((iorb <= 0) .or. (iorb > mxorb)) then
      write(u6,*) ' Illegal orbital number read :',iorb
      call abend_cvb()
    end if
    if (nread == 0) then
      write(u6,*) ' Orbital label in ORB keyword not found!'
      call abend_cvb()
    end if
    irdorbs(iorb) = mouse
    call fzero(orbs(1,iorb),mxaobf)
    call real_cvb(orbs(1,iorb),mxaobf,nread,0)
  else if (istr == 2) then
    ! 'STRUC'
    !  If previous orb. permutation disable:
    call mhpfreer_cvb(ip_cvb)
    mxread = mavailr_cvb()/2
    ip_cvb = mheapr_cvb(mxread)
    call realz_cvb(work(ip_cvb),mxread,nvbinp,0)
    call mreallocr_cvb(ip_cvb,nvbinp)
    kbasiscvb_inp = kbasis
  !else if (istr == 3) then
  !  ! 'READ'
  !  call fstring_cvb(readgs,nrdgs,istr2,ncmp,1)
  !  if (istr2 == 1) then
  !    ! 'ORB'
  !    iorb1 = 0
  !    iorb2 = 0
  !    call int_cvb(idum,1,nread,0)
  !    iorb1 = idum(1)
  !    if (nread == 0) then
  !      write(u6,*) ' No orbital number in READ,ORB keyword!'
  !      call abend_cvb()
  !    end if
  !    call fstring_cvb('TO      ',1,jstr,ncmp,1)
  !    if (jstr /= 0) then
  !      call int_cvb(idum,1,nread,0)
  !      iorb2 = idum(1)
  !      if (nread == 0) then
  !        write(u6,*) ' No orbital number after READ,...,TO, !'
  !        call abend_cvb()
  !      end if
  !    else
  !      iorb2 = iorb1
  !    end if
  !    jorb1 = iorb1
  !    jorb2 = iorb2
  !    call setstrtvb_cvb(strtvb)
  !    recordnm = strtvb
  !    do
  !      call fstring_cvb(readgs2,nrdgs2,istr3,ncmp,1)
  !      if (istr3 == 1) then
  !        call real_cvb(recordnm,1,nread,0)
  !        if (nread == 0) then
  !          write(u6,*) ' No identifier after READ,...,FROM, !'
  !          call abend_cvb()
  !        end if
  !      else if (istr3 == 2) then
  !        call int_cvb(idum,1,nread,0)
  !        jorb1 = idum(1)
  !        if (nread == 0) then
  !          write(u6,*) ' No orbital number after READ,...,AS, !'
  !          call abend_cvb()
  !        end if
  !        call fstring_cvb('TO      ',1,jstr,ncmp,1)
  !        if (jstr /= 0) then
  !          call int_cvb(idum,1,nread,0)
  !          jorb2 = idum(1)
  !          if (nread == 0) then
  !            write(u6,*) ' No orbital number after READ,...,TO, !'
  !            call abend_cvb()
  !          end if
  !        else
  !          jorb2 = jorb1
  !        end if
  !      end if
  !      if (istr3 == 0) exit
  !    end do
  !    call othergs_cvb(orbs,work(ip_cvb),recordnm,1,iorb1,iorb2,jorb1,jorb2)
  !  else if (istr2 == 2) then
  !    ! 'STRUC'
  !    istruc1 = 0
  !    istruc2 = 0
  !    call int_cvb(idum,1,nread,0)
  !    istruc1 = idum(1)
  !    if (nread == 0) then
  !      write(u6,*) ' No structure number in READ,STRUC keyword!'
  !      call abend_cvb()
  !    end if
  !    call fstring_cvb('TO      ',1,jstr,ncmp,1)
  !    if (jstr /= 0) then
  !      call int_cvb(idum,1,nread,0)
  !      istruc2 = idum(1)
  !      if (nread == 0) then
  !        write(u6,*) ' No structure number after READ,...,TO, !'
  !        call abend_cvb()
  !      end if
  !    else
  !      istruc2 = istruc1
  !    end if
  !    jstruc1 = istruc1
  !    jstruc2 = istruc2
  !    call setstrtvb_cvb(strtvb)
  !    recordnm = strtvb
  !    do
  !      call fstring_cvb(readgs2,nrdgs2,istr3,ncmp,1)
  !      if (istr3 == 1) then
  !        call real_cvb(recordnm,1,nread,0)
  !        if (nread == 0) then
  !          write(u6,*) ' No identifier after READ,...,FROM, !'
  !          call abend_cvb()
  !        end if
  !      else if (istr3 == 2) then
  !        call int_cvb(idum,1,nread,0)
  !        jstruc1 = idum(1)
  !        if (nread == 0) then
  !          write(u6,*) ' No structure number after READ,...,AS, !'
  !          call abend_cvb()
  !        end if
  !        call fstring_cvb('TO      ',1,jstr,ncmp,1)
  !        if (jstr /= 0) then
  !          call int_cvb(idum,1,nread,0)
  !          jstruc2 = idum(1)
  !          if (nread == 0) then
  !            write(u6,*) ' No structure number after READ,...,TO, !'
  !            call abend_cvb()
  !          end if
  !        else
  !          jstruc2 = jstruc1
  !        end if
  !      end if
  !      if (istr3 == 0) exit
  !    end do
  !    call othergs_cvb(orbs,work(ip_cvb),recordnm,2,istruc1,istruc2,jstruc1,jstruc2)
  !  else if (istr2 == 3) then
  !    ! 'ALL'
  !    call setstrtvb_cvb(strtvb)
  !    recordnm = strtvb
  !    call fstring_cvb('FROM    ',1,jstr,ncmp,1)
  !    if (jstr /= 0) then
  !      call real_cvb(recordnm,1,nread,0)
  !      if (nread == 0) then
  !        write(u6,*) ' No identifier after READ,...,FROM, !'
  !        call abend_cvb()
  !      end if
  !    end if
  !    call getguess_cvb(orbs,work(ip_cvb),recordnm,kbasiscvb_inp)
  !  end if
  else if (istr == 4) then
    ! 'AOBASIS'
    mouse = 2
  else if (istr == 5) then
    ! 'MOBASIS'
    mouse = 1
  end if
  ! 'END' , 'ENDGUESS' or unrecognized keyword -- end GUESS input:
  if ((istr == 0) .or. (istr == 6) .or. (istr == 7)) exit
end do

return
! Avoid unused argument warnings
if (.false.) call Unused_real(strtvb)

end subroutine gsinp_cvb
