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
      subroutine cvbmn_rvb(icode)
      use definitions, only: iwp
      use casvb_global, only: esym, n_iter
      implicit none
      integer(kind=iwp), intent(in):: icode

!  ICODE=0 standard casvb calculation
!  ICODE=1 variational calculation
!  ICODE=2 end of variational calculation (print summary)

      call cvbstart_rvb_lt9(icode)
      call main_cvb()
      call setretvals_cvb(esym,n_iter)

      end subroutine cvbmn_rvb
