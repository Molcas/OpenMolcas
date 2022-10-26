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
      subroutine small_svd(m,n,amat,umat,vmat,svals)
      implicit real*8 (a-h,o-z)
      dimension amat(m,*)
      dimension umat(m,*)
      dimension vmat(n,*)
      dimension svals(*)

      if (m.ge.n) then
       call svd_mgen(m,n,amat,umat,vmat,svals)
      else
       call svd_mlen(m,n,amat,umat,vmat,svals)
      end if
      return
      end
