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

subroutine orthcon_cvb(ipairs,ipair,igroups,ngroup,iorthlst,mxortl,mxpair)

implicit real*8(a-h,o-z)
#include "main_cvb.fh"
#include "optze_cvb.fh"
#include "files_cvb.fh"
#include "print_cvb.fh"
parameter(nstrin=7,ncmp=4,mxgroup=40)
character*8 string(nstrin)
character*3 glabel(mxgroup)
dimension ipairs(2,mxpair), ipair(mxorb_cvb,mxorb_cvb)
dimension igroups(mxorb_cvb,mxgroup), ngroup(mxgroup)
dimension iorthlst(mxortl)
save string
data string/'GROUP   ','ORTH    ','PAIRS   ','STRONG  ','FULL    ','END     ','ENDORTHC'/

call izero(ipair,mxorb_cvb*mxorb_cvb)
ngrp = 0
2000 call fstring_cvb(string,nstrin,istr,ncmp,2)
if (istr == 1) then
  ! 'GROUP'
  ngrp = ngrp+1
  if (ngrp > mxgroup) then
    write(6,*) ' Too many GROUP keywords in input!',mxgroup
    call abend_cvb()
  end if
  glabel(ngrp) = ' '
  call string_cvb(glabel(ngrp),1,nread,1)
  if ((glabel(ngrp)(1:1) < 'A') .or. (glabel(ngrp)(1:1) > 'Z')) then
    write(6,*) ' Group label must begin with a character A-Z: ',glabel(ngrp)
    call abend_cvb()
  end if
  call int_cvb(igroups(1,ngrp),mxorb_cvb,ngroup(ngrp),0)
  if (ngroup(ngrp) == -1) then
    write(6,*) ' Too many elements for group ',glabel(ngrp)
    call abend_cvb()
  end if
  do i=1,ngroup(ngrp)
    if ((igroups(i,ngrp) < 1) .or. (igroups(i,ngrp) > mxorb_cvb)) then
      write(6,*) ' Illegal orbital number in group ',glabel(ngrp),' :',igroups(i,ngrp)
      call abend_cvb()
    end if
  end do
  do i=1,ngrp-1
    if (glabel(ngrp) == glabel(i)) then
      write(6,*) ' Repeated label in GROUP keywords : ',glabel(ngrp)
      call abend_cvb()
    end if
  end do
else if (istr == 2) then
  ! 'ORTH'
  nsp = 0
175 continue
  call int_cvb(iorthlst(1+nsp),mxortl-nsp,nread,0)
  nsp = nsp+nread
  if (mxortl-nsp > 0) then
    call fstring_cvb(glabel,ngrp,igrp,3,0)
    if (igrp > 0) then
      nsp = nsp+1
      iorthlst(nsp) = -igrp
      if (mxortl-nsp > 0) goto 175
    end if
  end if
  do isp=1,nsp
    do jsp=isp+1,nsp
      ior = iorthlst(isp)
      jor = iorthlst(jsp)
      if ((ior > 0) .and. (jor > 0)) then
        ipair(ior,jor) = 1
      else if ((ior > 0) .and. (jor < 0)) then
        jor2 = -jor
        do jo=1,ngroup(jor2)
          ipair(ior,igroups(jo,jor2)) = 1
        end do
      else if ((ior < 0) .and. (jor > 0)) then
        ior2 = -ior
        do io=1,ngroup(ior2)
          ipair(jor,igroups(io,ior2)) = 1
        end do
      else if ((ior < 0) .and. (jor < 0)) then
        ior2 = -ior
        jor2 = -jor
        do io=1,ngroup(ior2)
          do jo=1,ngroup(jor2)
            ipair(igroups(io,ior2),igroups(jo,jor2)) = 1
          end do
        end do
      end if
    end do
  end do
else if (istr == 3) then
  ! 'PAIRS'
  nsp = 0
975 continue
  call int_cvb(ipairs(1+nsp,1),2*mxpair-nsp,nread,0)
  nsp = nsp+nread
  if (2*mxpair-nsp > 0) then
    call fstring_cvb(glabel,ngrp,igrp,3,0)
    if (igrp > 0) then
      nsp = nsp+1
      ipairs(nsp,1) = -igrp
      if (mxortl-nsp > 0) goto 975
    end if
  end if
  if (mod(nsp,2) == 1) then
    write(6,*) ' Odd number of orthogonalization numbers in PAIRS!'
    call abend_cvb()
  end if
  npairs = nsp/2
  do ipar=1,npairs
    ior = ipairs(1,ipar)
    jor = ipairs(2,ipar)
    if ((ior > 0) .and. (jor > 0)) then
      ipair(ior,jor) = 1
    else if ((ior > 0) .and. (jor < 0)) then
      jor2 = -jor
      do jo=1,ngroup(jor2)
        ipair(ior,igroups(jo,jor2)) = 1
      end do
    else if ((ior < 0) .and. (jor > 0)) then
      ior2 = -ior
      do io=1,ngroup(ior2)
        ipair(jor,igroups(io,ior2)) = 1
      end do
    else if ((ior < 0) .and. (jor < 0)) then
      ior2 = -ior
      jor2 = -jor
      do io=1,ngroup(ior2)
        do jo=1,ngroup(jor2)
          ipair(igroups(io,ior2),igroups(jo,jor2)) = 1
        end do
      end do
    end if
  end do
else if (istr == 4) then
  ! 'STRONG'
  do i=1,mxorb_cvb
    do j=i+1,mxorb_cvb
      if (.not. ((mod(i,2) == 1) .and. (j == i+1))) ipair(i,j) = 1
    end do
  end do
else if (istr == 5) then
  ! 'FULL'
  do i=1,mxorb_cvb
    do j=i+1,mxorb_cvb
      ipair(i,j) = 1
    end do
  end do
end if
! 'END' , 'ENDORTHC' or unrecognized keyword -- end of ORTHCON input:
if (.not. ((istr == 6) .or. (istr == 7) .or. (istr == 0))) goto 2000
call izero(ipairs,2*mxpair)
nort = 0
do i=1,mxorb_cvb
  do j=i+1,mxorb_cvb
    if ((ipair(i,j) == 1) .or. (ipair(j,i) == 1)) then
      nort = nort+1
      ipairs(1,nort) = i
      ipairs(2,nort) = j
    end if
  end do
end do

return

end subroutine orthcon_cvb
