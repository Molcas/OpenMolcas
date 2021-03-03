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

      subroutine Untested(label)
      implicit none
      character(len=*), intent(in) :: label

      write(6,*)
      write(6,*) label
      write(6,*) 'This code is untested and should be carefully '//
     &           'verified.'
#ifndef _DEVEL_
      call abend()
#endif

      end subroutine Untested
