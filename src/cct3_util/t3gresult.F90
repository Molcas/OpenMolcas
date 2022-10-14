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
        subroutine t3gresult (symi,symj,i,j,eaaa,eaab,eabb,ebbb)
!
!       this routine get:
!       1) value of SymI,SymJ,I,J
!       2) present stage of energies
!       from T3tEne file
!
        integer symi,symj,i,j
        real*8 eaaa,eaab,eabb,ebbb
!
!       help variable
!
        integer lun
!
        lun=1
        Call Molcas_Open(Lun,'T3tEne')
!       open (unit=lun,file='T3tEne')
!
!       read blank, since there is Symimin,imin,Symjmin,jmin
!       on first line
!
        read (lun,*) i
!
        read (lun,*) symi,symj
        read (lun,*) i,j
        read (lun,*) eaaa
        read (lun,*) eaab
        read (lun,*) eabb
        read (lun,*) ebbb
!
!98      format (2x,2(i4,2x))
!99      format (2x,f22.16)
!
        close (lun)
!
        return
        end
