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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************
      subroutine getritrfinfo(nnbstr_,maxvec_,n2_)
      use Cholesky, only: infvec_N2, MaxVec, nnBstR
      use definitions, only: iwp
      implicit none
      integer(kind=iwp), intent(out) :: nnbstr_(8,3), maxvec_, n2_
      maxvec_ = maxvec
      n2_     = infvec_n2
      nnbstr_ = nnbstr
      end subroutine getritrfinfo
