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

* Legacy code that is not used by current LAPACK but could be "needed"
* by other libraries. Define here just as dummy procedures
*
* If you modify this, the files f[1-5].fh should be updated too, use delayed.py

      module LEGACY_MOD

      contains

      subroutine dgetf2( m, n, a, lda, ipiv, info )
        implicit none
        integer :: info, lda, m, n
        integer :: ipiv( * )
        real*8 :: a( lda, * )
      end subroutine dgetf2

      subroutine dpotf2( uplo, n, a, lda, info )
        implicit none
        character :: uplo
        integer :: info, lda, n
        real*8 :: a( lda, * )
      end subroutine dpotf2

      end module LEGACY_MOD
