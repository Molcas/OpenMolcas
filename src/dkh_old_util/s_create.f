************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2004, Alexander Wolf                                   *
*               2004,2007, Markus Reiher                               *
************************************************************************
      subroutine s_create (dkhscfflg,sused,snumber,scounter,stimes,
     *                     wstimes,sorder,smult,s,scrleng,scrchar,
     *                     stimestot,wopscounter,no_s)
c
c******************************************************************************
c
c   This SR belongs to dkhparser_symbolic.
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.1.0
c
c   last modified: 15.01.2007 (M. Reiher, ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      logical dkhscfflg,no_s
      integer sused,snumber,scounter(maxsnumber),stimes(maxsnumber),
     *        wstimes(maxuops,maxsnumber),sorder(maxsnumber,3),
     *        smult(maxsnumber),scrleng(maxsnumber),sinfo(maxsnumber,2),
     *        stimestot,wopscounter
c
      character*(4) s(maxsnumber)
      character*(9) scrchar(maxsnumber)
c
C     integer j,k,dkh_char2int,idum1,idum2
      integer j,k,dkh_char2int
      character*(3) dkh_int2char
c
*** use flag, which determines the memory available
      if (no_s) then
        sused=0
      else
        sused=216
      end if
c
      do 10 j=1,sused
        scounter(j)=1
        stimes(j)=0
        do 15 k=1,wopscounter
          wstimes(k,j)=0
  15    continue
        do 18 k=1,3
          sorder(j,k)=0
  18    continue
        sinfo(j,1)=0
        sinfo(j,2)=0
        smult(j)=0
        s(j)='    '
        scrchar(j)='         '
  10  continue
      snumber=0
      stimestot=0
c
      do 20 j=1,sused
        s(j)(1:4)='S'//dkh_int2char(j)
  20  continue
      scrchar(1)  ='W01O01   '
      scrchar(2)  ='O01W01   '
      scrchar(3)  ='W01W01   '
      scrchar(4)  ='W01E01W01'
      scrchar(5)  ='W01CO0   '
      scrchar(6)  ='CO0W01   '
      scrchar(7)  ='W01CE0W01'
      scrchar(8)  ='W02O01   '
      scrchar(9)  ='O01W02   '
      scrchar(10) ='W02W01   '
      scrchar(11) ='W01W02   '
      scrchar(12) ='W02E01W01'
      scrchar(13) ='W01E01W02'
      scrchar(14) ='W02W02   '
      scrchar(15) ='W02E01W02'
      scrchar(16) ='W02CO0   '
      scrchar(17) ='CO0W02   '
      scrchar(18) ='W02CE0W01'
      scrchar(19) ='W01CE0W02'
      scrchar(20) ='W02CE0W02'
      scrchar(21) ='W03O01   '
      scrchar(22) ='O01W03   '
      scrchar(23) ='W03W01   '
      scrchar(24) ='W01W03   '
      scrchar(25) ='W03E01W01'
      scrchar(26) ='W01E01W03'
      scrchar(27) ='W03W02   '
      scrchar(28) ='W02W03   '
      scrchar(29) ='W03E01W02'
      scrchar(30) ='W02E01W03'
      scrchar(31) ='W03W03   '
      scrchar(32) ='W03E01W03'
      scrchar(33) ='W03CO0   '
      scrchar(34) ='CO0W03   '
      scrchar(35) ='W03CE0W01'
      scrchar(36) ='W01CE0W03'
      scrchar(37) ='W03CE0W02'
      scrchar(38) ='W02CE0W03'
      scrchar(39) ='W03CE0W03'
      scrchar(40) ='W04O01   '
      scrchar(41) ='O01W04   '
      scrchar(42) ='W04W01   '
      scrchar(43) ='W01W04   '
      scrchar(44) ='W04E01W01'
      scrchar(45) ='W01E01W04'
      scrchar(46) ='W04W02   '
      scrchar(47) ='W02W04   '
      scrchar(48) ='W04E01W02'
      scrchar(49) ='W02E01W04'
      scrchar(50) ='W04W03   '
      scrchar(51) ='W03W04   '
      scrchar(52) ='W04E01W03'
      scrchar(53) ='W03E01W04'
      scrchar(54) ='W04W04   '
      scrchar(55) ='W04E01W04'
      scrchar(56) ='W04CO0   '
      scrchar(57) ='CO0W04   '
      scrchar(58) ='W04CE0W01'
      scrchar(59) ='W01CE0W04'
      scrchar(60) ='W04CE0W02'
      scrchar(61) ='W02CE0W04'
      scrchar(62) ='W04CE0W03'
      scrchar(63) ='W03CE0W04'
      scrchar(64) ='W04CE0W04'
      scrchar(65) ='W05O01   '
      scrchar(66) ='O01W05   '
      scrchar(67) ='W05W01   '
      scrchar(68) ='W01W05   '
      scrchar(69) ='W05E01W01'
      scrchar(70) ='W01E01W05'
      scrchar(71) ='W05W02   '
      scrchar(72) ='W02W05   '
      scrchar(73) ='W05E01W02'
      scrchar(74) ='W02E01W05'
      scrchar(75) ='W05W03   '
      scrchar(76) ='W03W05   '
      scrchar(77) ='W05E01W03'
      scrchar(78) ='W03E01W05'
      scrchar(79) ='W05W04   '
      scrchar(80) ='W04W05   '
      scrchar(81) ='W05E01W04'
      scrchar(82) ='W04E01W05'
      scrchar(83) ='W05W05   '
      scrchar(84) ='W05E01W05'
      scrchar(85) ='W05CO0   '
      scrchar(86) ='CO0W05   '
      scrchar(87) ='W05CE0W01'
      scrchar(88) ='W01CE0W05'
      scrchar(89) ='W05CE0W02'
      scrchar(90) ='W02CE0W05'
      scrchar(91) ='W05CE0W03'
      scrchar(92) ='W03CE0W05'
      scrchar(93) ='W05CE0W04'
      scrchar(94) ='W04CE0W05'
      scrchar(95) ='W05CE0W05'
      scrchar(96) ='W06O01   '
      scrchar(97) ='O01W06   '
      scrchar(98) ='W06W01   '
      scrchar(99) ='W01W06   '
      scrchar(100)='W06E01W01'
      scrchar(101)='W01E01W06'
      scrchar(102)='W06W02   '
      scrchar(103)='W02W06   '
      scrchar(104)='W06E01W02'
      scrchar(105)='W02E01W06'
      scrchar(106)='W06W03   '
      scrchar(107)='W03W06   '
      scrchar(108)='W06E01W03'
      scrchar(109)='W03E01W06'
      scrchar(110)='W06W04   '
      scrchar(111)='W04W06   '
      scrchar(112)='W06E01W04'
      scrchar(113)='W04E01W06'
      scrchar(114)='W06W05   '
      scrchar(115)='W05W06   '
      scrchar(116)='W06E01W05'
      scrchar(117)='W05E01W06'
      scrchar(118)='W06W06   '
      scrchar(119)='W06E01W06'
      scrchar(120)='W06CO0   '
      scrchar(121)='CO0W06   '
      scrchar(122)='W06CE0W01'
      scrchar(123)='W01CE0W06'
      scrchar(124)='W06CE0W02'
      scrchar(125)='W02CE0W06'
      scrchar(126)='W06CE0W03'
      scrchar(127)='W03CE0W06'
      scrchar(128)='W06CE0W04'
      scrchar(129)='W04CE0W06'
      scrchar(130)='W06CE0W05'
      scrchar(131)='W05CE0W06'
      scrchar(132)='W06CE0W06'
      scrchar(133)='W07O01   '
      scrchar(134)='O01W07   '
      scrchar(135)='W07W01   '
      scrchar(136)='W01W07   '
      scrchar(137)='W07E01W01'
      scrchar(138)='W01E01W07'
      scrchar(139)='W07W02   '
      scrchar(140)='W02W07   '
      scrchar(141)='W07E01W02'
      scrchar(142)='W02E01W07'
      scrchar(143)='W07W03   '
      scrchar(144)='W03W07   '
      scrchar(145)='W07E01W03'
      scrchar(146)='W03E01W07'
      scrchar(147)='W07W04   '
      scrchar(148)='W04W07   '
      scrchar(149)='W07E01W04'
      scrchar(150)='W04E01W07'
      scrchar(151)='W07W05   '
      scrchar(152)='W05W07   '
      scrchar(153)='W07E01W05'
      scrchar(154)='W05E01W07'
      scrchar(155)='W07W06   '
      scrchar(156)='W06W07   '
      scrchar(157)='W07E01W06'
      scrchar(158)='W06E01W07'
      scrchar(159)='W07W07   '
      scrchar(160)='W07E01W07'
      scrchar(161)='W07CO0   '
      scrchar(162)='CO0W07   '
      scrchar(163)='W07CE0W01'
      scrchar(164)='W01CE0W07'
      scrchar(165)='W07CE0W02'
      scrchar(166)='W02CE0W07'
      scrchar(167)='W07CE0W03'
      scrchar(168)='W03CE0W07'
      scrchar(169)='W08O01   '
      scrchar(170)='O01W08   '
      scrchar(171)='W08W01   '
      scrchar(172)='W01W08   '
      scrchar(173)='W08E01W01'
      scrchar(174)='W01E01W08'
      scrchar(175)='W08W02   '
      scrchar(176)='W02W08   '
      scrchar(177)='W08E01W02'
      scrchar(178)='W02E01W08'
      scrchar(179)='W08W03   '
      scrchar(180)='W03W08   '
      scrchar(181)='W08E01W03'
      scrchar(182)='W03E01W08'
      scrchar(183)='W08W04   '
      scrchar(184)='W04W08   '
      scrchar(185)='W08E01W04'
      scrchar(186)='W04E01W08'
      scrchar(187)='W08W05   '
      scrchar(188)='W05W08   '
      scrchar(189)='W08E01W05'
      scrchar(190)='W05E01W08'
      scrchar(191)='W08W06   '
      scrchar(192)='W06W08   '
      scrchar(193)='W08E01W06'
      scrchar(194)='W06E01W08'
      scrchar(195)='W08W07   '
      scrchar(196)='W07W08   '
      scrchar(197)='W08E01W07'
      scrchar(198)='W07E01W08'
      scrchar(199)='W08W08   '
      scrchar(200)='W08E01W08'
      scrchar(201)='W08CO0   '
      scrchar(202)='CO0W08   '
      scrchar(203)='W08CE0W01'
      scrchar(204)='W01CE0W08'
      scrchar(205)='W08CE0W02'
      scrchar(206)='W02CE0W08'
      scrchar(207)='W09O01   '
      scrchar(208)='O01W09   '
      scrchar(209)='W09W01   '
      scrchar(210)='W01W09   '
      scrchar(211)='W09CO0   '
      scrchar(212)='CO0W09   '
      scrchar(213)='W09CE0W01'
      scrchar(214)='W01CE0W09'
      scrchar(215)='W10CO0   '
      scrchar(216)='CO0W10   '
c
      do 25 j=1,sused
        if (scrchar(j)(7:9).eq.'   ') then
          scrleng(j)=6
        else
          scrleng(j)=9
        endif
  25  continue
c
      do 30 j=1,2
        do 31 k=1,sused
          if (scrleng(k).eq.6) then
            if (index(scrchar(k)(3*j-2:3*j),'C').eq.0) then
              sinfo(k,j)=dkh_char2int(2,scrchar(k)(3*j-1:3*j))
            else
              sinfo(k,j)=0
            endif
          endif
          if (scrleng(k).eq.9) then
            sinfo(k,j)=dkh_char2int(2,scrchar(k)(6*j-4:6*j-3))
          endif
  31    continue
  30  continue
c
      if (dkhscfflg) then
        do 130 j=1,sused
          sorder(j,3)=sinfo(j,1)+sinfo(j,2)
          if (scrleng(j).eq.9) sorder(j,3)=sorder(j,3)+1
          if (scrleng(j).eq.6 .and. index(scrchar(j)(1:6),
     *                            'CO0').gt.0) sorder(j,3)=sorder(j,3)+1
 130    continue
      else
        do 131 j=1,sused
          sorder(j,1)=sinfo(j,1)+sinfo(j,2)
          if (index(scrchar(j)(1:scrleng(j)),'E01').gt.0)
     *                                         sorder(j,1)=sorder(j,1)+1
          if (index(scrchar(j)(1:scrleng(j)),'C').gt.0) sorder(j,2)=1
          sorder(j,3)=sorder(j,1)+sorder(j,2)
 131    continue
      endif
c
      if (dbgflg.ge.1) then
        write (dbgunit,1001) sused
1001    format (/2X,'Creation and initialization of scratch arrays "s"',
     *          ' successful.',5X,'(sused=',I3,')')
      endif
c
      if (dbgflg.ge.3) then
         write (dbgunit,1156)sused
 1156    format (/2X,'The following ',I3,' Sxxx-matrices have been set',
     *           ' up:',//2X,'Sxxx',2X,'scrchar',2X,'sinfo',1X,
     *           'sinfo',2X,'sorder',1X,'sorder',1X,'sorder',2X,
     *           'smult',2X,'scrleng',/18X,'(,1)',2X,'(,2)',3X,'(V)',
     *           4X,'(X)',3X,'(tot)',/)
         do 6298 k=1,sused
           write (dbgunit,8877) s(k),scrchar(k),sinfo(k,1),sinfo(k,2),
     *                          sorder(k,1),sorder(k,2),sorder(k,3),
     *                          smult(k),scrleng(k)
8877       format (2X,A4,2X,A9,2X,I2,4X,I2,4X,I2,5X,I2,4X,I3,5X,I3,5X,
     *             I3)
6298     continue
      endif
c
      return
      end
