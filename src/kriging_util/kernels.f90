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
! Copyright (C) 2019, Gerardo Raggi                                    *
!***********************************************************************
        SUBROUTINE kernels(iter,nInter)
            use globvar
            integer i,z,j,iter,nInter,lm
            real*8 value
!
            call miden(iter)
!
        END

        SUBROUTINE setlkriging(lv,iter,nInter)
            use globvar
            integer iter,nInter,i
            real*8 lv
            do i = 1,nInter
                l(i)=lv
            enddo
            call covarmatrix(iter,nInter)
            call k(iter)
            write (6,*) 'set l value, lh:',l(1)
        END

        subroutine miden(iter)
            use globvar
            integer j,iter
            iden=0
            forall(j=1:iter) iden(j,j)=1
        end subroutine
