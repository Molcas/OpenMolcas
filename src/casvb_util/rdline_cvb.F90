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

subroutine rdline_cvb(nfield)

implicit real*8(a-h,o-z)
! BLANKDELIM signifies whether blanks are used to delimit fields:
logical blankdelim
#include "luinp_cvb.fh"
#include "rdline.fh"
parameter(nblank=2,ncomeol=3,neol=4,neofield=1,neof=2,nalias=2)
character*1 blanks(nblank)
! COMEOL are comments that comment out the rest of the line
! (might add 'bracketing' comments (e.g. /* ... */ ) later):
character*3 comeol(ncomeol)
character*1 eol(neol)
character*1 eofield(neofield)
character*10 eof(neof)
character*5 alias(nalias,2)
save blankdelim
save blanks, comeol, eof, eol, eofield
data blankdelim/.true./
data blanks/' ',','/,comeol/'***','!  ','*  '/
data eof/'---       ','ENDOFINPUT'/
data eol/';','=','{','}'/
data eofield/' '/
data alias/'={   ',';    ','}    ',';END;'/

50 continue
if (nline == -1) then
  nfield = -1
  return
end if
if (iline < nline) then
  iline = iline+1
  goto 899
else
  iline = 1
end if
! Read full input line from file and make preparations:
100 read(inp,'(a)',end=200) line
lenline = len_trim_cvb(line)
call strip_blanks_cvb(line,lenline,blanks,nblank,blankdelim)
call upper_case_cvb(line,lenline)
! Check for "end-of-file" character sequence:
do ieof=1,neof
  ilength = len_trim_cvb(eof(ieof))
  if (line(1:ilength) == eof(ieof)(1:ilength)) goto 200
end do
! Comment strings (skip rest of line):
indmin = lenline+1
do icom=1,ncomeol
  ind = index(line(1:lenline),comeol(icom)(1:len_trim_cvb(comeol(icom))))
  if (ind /= 0) indmin = min(indmin,ind)
end do
lenline = len_trim_cvb(line(1:indmin-1))
if (lenline == 0) goto 100
! Aliases:
do ialias=1,nalias
560 ind = index(line(1:lenline),alias(ialias,1)(1:len_trim_cvb(alias(ialias,1))))
  if (ind /= 0) then
    call charinsert_cvb(alias(ialias,2),len_trim_cvb(alias(ialias,2)),line,lenline,ind,len_trim_cvb(alias(ialias,1)))
    goto 560
  end if
end do
lenline = len_trim_cvb(line(1:indmin-1))
if (lenline == 0) goto 100
!  Split into lines:
call izero(ilv,lenline)
do ieol=1,neol
  if = 0
  ilength = len_trim_cvb(eol(ieol))
650 ind = index(line(if+1:lenline),eol(ieol)(1:ilength))
  if (ind /= 0) then
    ilv(ind+if) = 1
    if = ind+if
    goto 650
  end if
end do
nlold = nline
nline = 1
do ich=1,lenline
  if (ilv(ich) == 1) nline = nline+1
end do
! Split into fields:
do ieofield=1,neofield
  if = 0
  ilength = max(1,len_trim_cvb(eofield(ieofield)))
850 ind = index(line(if+1:lenline),eofield(ieofield)(1:ilength))
  if (ind /= 0) then
    ilv(ind+if) = 2
    if = ind+if
    goto 850
  end if
end do
! Eliminate field separators at end of lines:
do ieofield=1,neofield
  ihadchar = 0
  do ich=lenline,1,-1
    if (ilv(ich) == 1) then
      ihadchar = 0
    else if (ilv(ich) == 2) then
      if (ihadchar == 0) ilv(ich) = 0
      ilength = len_trim_cvb(eofield(ieofield))
      do i=0,ilength-1
        line(i+ich:i+ich) = ' '
      end do
      ihadchar = 1
    else
      ihadchar = 1
    end if
  end do
end do
! Count the number of fields in present card (ILINE).
! Also make sure line is not all empty:
899 continue
nfield = 1
jline = 1
ilinebeg = 0
ilineend = -1
do ich=1,lenline
  if (jline == iline-1) ilinebeg = ich+1
  if (ilv(ich) == 1) jline = jline+1
  if ((jline == iline+1) .and. (ilineend == -1)) ilineend = ich-1
  if ((jline == iline) .and. (ilv(ich) == 2)) nfield = nfield+1
end do
if (iline == 1) ilinebeg = 1
if (ilineend == -1) ilineend = lenline
! Go back and read a new line if this one is empty:
if (ilinebeg > ilineend) goto 50
if (len_trim_cvb(line(ilinebeg:ilineend)) == 0) goto 50
return
200 nline = -1
nfield = -1

return

end subroutine rdline_cvb

block data rdline_bd
#include "rdline.fh"
data iline/0/,nline/0/,nlold/0/
end block data rdline_bd
