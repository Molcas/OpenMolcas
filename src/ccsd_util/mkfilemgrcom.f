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
       subroutine mkfilemgrcom
c
c     this routine make names of temp files and define initial values
c     of filestatus
c
#include "filemgr.fh"
#include "ccsd1.fh"
c
c     help variable
c
       integer nhelp
c
c1    def filestatus
c
       do 100 nhelp=10,maxfiles
       filestatus(nhelp)=0
 100    continue
c
c2    def disk addresses
c
       do 200 nhelp=10,maxfiles
       daddr(nhelp)=0
 200    continue
c
c3    def filenames
c
       filename(10)='INTAB'
       filename(11)='INTA1'
       filename(12)='INTA2'
       filename(13)='INTA3'
       filename(14)='INTA4'
       filename(15)='INTSTA'
       filename(16)=filerst
       filename(17)='Temp17'
       filename(18)='Temp18'
       filename(19)='Temp19'
       filename(20)='Temp20'
       filename(21)='Temp21'
       filename(22)='Temp22'
       filename(23)='Temp23'
       filename(24)='Temp24'
       filename(25)='Temp25'
       filename(26)='Temp26'
       filename(27)='Temp27'
       filename(28)='Temp28'
       filename(29)='Temp29'
       filename(30)='Temp30'
       filename(31)='Temp31'
       filename(32)='Temp32'
       filename(33)='Temp33'
       filename(34)='Temp34'
       filename(35)='Temp35'
       filename(36)='Temp36'
       filename(37)='Temp37'
       filename(38)='Temp38'
       filename(39)='Temp39'
       filename(40)='Temp40'
       filename(41)='Temp41'
       filename(42)='Temp42'
       filename(43)='Temp43'
       filename(44)='Temp44'
       filename(45)='Temp45'
       filename(46)='Temp46'
       filename(47)='Temp47'
       filename(48)='Temp48'
       filename(49)='Temp49'
       filename(50)='Temp50'
c

       return
       end
