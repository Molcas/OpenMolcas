************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      logical function next_non_comment(lu,line)
      implicit none
      integer, intent(in) :: lu
      character(*), intent(out) :: line
      next_non_comment = .false.
100   continue
      read(lu,'(A)',End=900) line
      line = adjustl(line)
      if ( line(1:1).eq.'*' ) goto 100
      if ( line.eq.' ' ) goto 100
      next_non_comment = .true.
900   continue
      return
      end function
