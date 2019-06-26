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
        SUBROUTINE setlkriging(lv,nInter)
            use globvar
            integer i, nInter
            real*8 lv(nInter)
            If (nInter.eq.nInter_Save) Then
               do i = 1,nInter_save
                l(i)=lv(i)
               enddo
            Else If (nInter.eq.1) Then
               do i = 1,nInter_save
                l(i)=lv(1)
               enddo
            Else
               Write (6,*) "setlkriging: illegal nInter value."
               Call Abend()
            End If
            call covarMatrix(nPoints_Save,nInter_save)
            call k(nPoints_save)
!           write (6,*) 'set l value, lh:',l(1)
        END SUBROUTINE setlkriging

        subroutine miden(iden,msize)
            !use globvar
            integer j,msize
            real*8 iden(msize,msize)
            iden=0
            forall(j=1:msize) iden(j,j)=1
        end subroutine