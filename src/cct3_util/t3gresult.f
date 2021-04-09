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
        subroutine t3gresult (symi,symj,i,j,eaaa,eaab,eabb,ebbb)
c
c       this routine get:
c       1) value of SymI,SymJ,I,J
c       2) present stage of energies
c       from T3tEne file
c
        integer symi,symj,i,j
        real*8 eaaa,eaab,eabb,ebbb
c
c       help variable
c
        integer lun
c
        lun=1
        Call Molcas_Open(Lun,'T3tEne')
*       open (unit=lun,file='T3tEne')
c
c       read blank, since there is Symimin,imin,Symjmin,jmin
c       on first line
c
        read (lun,*) i
c
        read (lun,*) symi,symj
        read (lun,*) i,j
        read (lun,*) eaaa
        read (lun,*) eaab
        read (lun,*) eabb
        read (lun,*) ebbb
c
c98      format (2x,2(i4,2x))
c99      format (2x,f22.16)
c
        close (lun)
c
        return
        end
