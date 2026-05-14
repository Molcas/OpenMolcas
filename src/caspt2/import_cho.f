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
      subroutine import_cho(numcho_pt2,infvec_n2_pt2,maxvec_pt2)
      use Cholesky, only: infvec_N2, MaxVec, nSym, NumCho
      use definitions, only: iwp
      implicit none
      integer(kind=iwp), intent(out):: infvec_n2_pt2, maxvec_pt2
      integer(kind=iwp), intent(out):: numcho_pt2(8)

      numcho_pt2(:)=0
      numcho_pt2(1:nSym)=numcho(1:nSym)
      maxvec_pt2=maxvec
      infvec_N2_pt2=infvec_N2
      end subroutine import_cho
