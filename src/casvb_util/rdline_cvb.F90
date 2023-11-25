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

use casvb_global, only: iline, ilv, inp, lenline, line, nline, nlold
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(out) :: nfield
integer(kind=iwp) :: ialias, ich, icom, ieof, ieofield, ieol, iff, ihadchar, ilength, ilinebeg, ilineend, ind, indmin, istatus, &
                     jline
integer(kind=iwp), parameter :: nalias = 2, nblank = 2, ncomeol = 3, neof = 2, neofield = 1, neol = 4
character(len=*), parameter :: alias(nalias,2) = reshape(['={   ',';    ','}    ',';END;'],[nalias,2]), &
                               blanks(nblank) = [' ',','], &
                               comeol(ncomeol) = ['***','!  ','*  '], & ! COMEOL are comments that comment out the rest of the line
                                                                        ! (might add 'bracketing' comments (e.g. /* ... */ ) later)
                               eof(neof) = ['---       ','ENDOFINPUT'], &
                               eofield(neofield) = [' '], &
                               eol(neol) = [';','=','{','}']
logical(kind=iwp), parameter :: blankdelim = .true. ! BLANKDELIM signifies whether blanks are used to delimit fields

do
  if (nline == -1) then
    nfield = -1
    return
  end if
  if (iline < nline) then
    iline = iline+1
  else
    iline = 1
    ! Read full input line from file and make preparations:
    do
      read(inp,'(a)',iostat=istatus) line
      if (istatus < 0) then
        nline = -1
        nfield = -1
        return
      end if
      call strip_blanks_cvb(line,blanks,nblank,blankdelim)
      call upcase(line)
      ! Check for "end-of-file" character sequence:
      do ieof=1,neof
        ilength = len_trim(eof(ieof))
        if (line(1:ilength) == trim(eof(ieof))) then
          nline = -1
          nfield = -1
          return
        end if
      end do
      ! Comment strings (skip rest of line):
      lenline = len_trim(line)
      indmin = lenline+1
      do icom=1,ncomeol
        ind = index(line(1:lenline),trim(comeol(icom)))
        if (ind /= 0) indmin = min(indmin,ind)
      end do
      lenline = len_trim(line(1:indmin-1))
      if (lenline /= 0) then
        ! Aliases:
        do ialias=1,nalias
          do
            ind = index(line(1:lenline),trim(alias(ialias,1)))
            if (ind == 0) exit
            call charinsert_cvb(alias(ialias,2),len_trim(alias(ialias,2)),line,lenline,ind,len_trim(alias(ialias,1)))
          end do
        end do
        lenline = len_trim(line(1:indmin-1))
        if (lenline /= 0) exit
      end if
    end do
    ! Split into lines:
    ilv(1:lenline) = 0
    do ieol=1,neol
      iff = 0
      do
        ind = index(line(iff+1:lenline),trim(eol(ieol)))
        if (ind == 0) exit
        ilv(ind+iff) = 1
        iff = ind+iff
      end do
    end do
    nlold = nline
    nline = 1
    do ich=1,lenline
      if (ilv(ich) == 1) nline = nline+1
    end do
    ! Split into fields:
    do ieofield=1,neofield
      iff = 0
      ilength = max(1,len_trim(eofield(ieofield)))
      do
        ind = index(line(iff+1:lenline),eofield(ieofield)(1:ilength))
        if (ind == 0) exit
        ilv(ind+iff) = 2
        iff = ind+iff
      end do
    end do
    ! Eliminate field separators at end of lines:
    do ieofield=1,neofield
      ihadchar = 0
      do ich=lenline,1,-1
        if (ilv(ich) == 1) then
          ihadchar = 0
        else if (ilv(ich) == 2) then
          if (ihadchar == 0) ilv(ich) = 0
          ilength = len_trim(eofield(ieofield))
          line(ich:ich+ilength-1) = ''
          ihadchar = 1
        else
          ihadchar = 1
        end if
      end do
    end do
  end if
  ! Count the number of fields in present card (ILINE).
  ! Also make sure line is not all empty:
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
  if (ilinebeg <= ilineend) then
    if (len_trim(line(ilinebeg:ilineend)) /= 0) exit
  end if
end do

return

end subroutine rdline_cvb
