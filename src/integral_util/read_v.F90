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
      subroutine Read_v(lunit,work,istrt,iend,inc,ierr)
      implicit None
      Integer lUnit, iStrt, iEnd, Inc, iErr
      real*8 work(iend)

      Integer i
      ierr=0
      read(lunit,*,err=100) (work(istrt+i),i=0,iend-istrt,inc)
      goto 110
100   ierr=1
110   continue
      return
      end subroutine Read_v
