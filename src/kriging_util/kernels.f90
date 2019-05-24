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
        SUBROUTINE setlkriging(lv)
            use globvar
            integer i
            real*8 lv
            call miden()
            do i = 1,nInter_save
                l(i)=lv
            enddo
            call covarmatrix(nPoints_Save,nInter_save)
            call k(nPoints_save)
!           write (6,*) 'set l value, lh:',l(1)
        END SUBROUTINE setlkriging

        subroutine miden()
            use globvar
            integer j
            iden=0
            forall(j=1:nPoints_save) iden(j,j)=1
        end subroutine
