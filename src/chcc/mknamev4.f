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
        subroutine MkNameV4 (i,j,k,l,Schem,Nomen)
c
c       help routine to ReaW4, producing name of V4file
c       ex: Schem='XY', i=1, j=3, k=5, l=07 ->  Nomen='XY01030507'
c
        implicit none
        integer i,j,k,l
        character*2 Schem
        character*10 Nomen
c
c       help variables
        character*1 Chr(1:10)
        character*2 digit(1:64)
        character*2 ichr,jchr,kchr,lchr
        character*2 baza
        character*10 meno
c
        equivalence (Chr(1),meno)
        equivalence (Chr(1),baza)
        equivalence (Chr(3),ichr)
        equivalence (Chr(5),jchr)
        equivalence (Chr(7),kchr)
        equivalence (Chr(9),lchr)
c
c
c        quite a porno this piece
        digit(1)='01'
        digit(2)='02'
        digit(3)='03'
        digit(4)='04'
        digit(5)='05'
        digit(6)='06'
        digit(7)='07'
        digit(8)='08'
        digit(9)='09'
        digit(10)='10'
        digit(11)='11'
        digit(12)='12'
        digit(13)='13'
        digit(14)='14'
        digit(15)='15'
        digit(16)='16'
        digit(17)='17'
        digit(18)='18'
        digit(19)='19'
        digit(20)='20'
        digit(21)='21'
        digit(22)='22'
        digit(23)='23'
        digit(24)='24'
        digit(25)='25'
        digit(26)='26'
        digit(27)='27'
        digit(28)='28'
        digit(29)='29'
        digit(30)='30'
        digit(31)='31'
        digit(32)='32'
        digit(33)='33'
        digit(34)='34'
        digit(35)='35'
        digit(36)='36'
        digit(37)='37'
        digit(38)='38'
        digit(39)='39'
        digit(40)='40'
        digit(41)='41'
        digit(42)='42'
        digit(43)='43'
        digit(44)='44'
        digit(45)='45'
        digit(46)='46'
        digit(47)='47'
        digit(48)='48'
        digit(49)='49'
        digit(50)='50'
        digit(51)='51'
        digit(52)='52'
        digit(53)='53'
        digit(54)='54'
        digit(55)='55'
        digit(56)='56'
        digit(57)='57'
        digit(58)='58'
        digit(59)='59'
        digit(60)='60'
        digit(61)='61'
        digit(62)='62'
        digit(63)='63'
        digit(64)='64'
c
c
        baza=Schem
        ichr=digit(i)
        jchr=digit(j)
        kchr=digit(k)
        lchr=digit(l)
        Nomen=meno
c
        return
        end
