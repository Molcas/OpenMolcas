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
        subroutine UrobInt (W,maxa,maxbe,n,lun)
c
c       vyraba subor LunInt so simulovanymi (a"be"|b"ga") integralmi
c       kde je za sebou N blokov
c
        implicit none
#include "chcc1.fh"
#include "o2v4.fh"
c
        integer maxa,maxbe,n,lun
        real*8 W(1)
c
c       help variables
        integer length,i
c
c       open (unit=Lun,file='IntFil',form='unformatted')
c
c1      def legth
c
        length=maxa*maxa*maxbe*maxbe
c
c2      cycle over N
c
        do i=1,n
c
c2.1    full W with random numbers
c
        call RNFill (length,W(1),1.0d-2)
c
c2.2    write block
c
        write (6,*) 'Vint',i,length
        call wri_chcc (lun,length,W(1))
c
        end do
c
c
        rewind(lun)
c
        return
        end
