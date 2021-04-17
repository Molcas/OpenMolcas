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
        subroutine MkV_K22 (W1,W2,dim)
c
c       this routine do
c       W1(p) = -W2(p)
c
c       N.B. toto je iba plytke copy s opacnym znamienkom,
c            da sa aj s blasmi spravit
c       N.B. Kvajt odflaknute
c
c
        implicit none
        integer dim
        real*8 W1(1:dim)
        real*8 W2(1:dim)
c
c       help variables
        integer p
c
c
        do p=1,dim
          W1(p)=-W2(p)
        end do
c
c
        return
        end
