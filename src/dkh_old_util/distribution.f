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
* Copyright (C) 2004,2005, Alexander Wolf                              *
*               2004,2005, Markus Reiher                               *
************************************************************************
      subroutine distribution (unit,dkhorder,xorder,dkhscfflg,
     *               ordercounter,opcounter,operleng,oporder,evenodd,
     *               doperators,operators,xordercounter,xopcounter,
     *               xoperleng,xoporder,xevenodd,xdoperators,xoperators)
c
c******************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 18.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer unit,dkhorder,xorder,ordercounter(0:maxorder),opcounter,
     *        operleng(maxoperators),oporder(maxoperators,3),
     *        evenodd(maxoperators),
     *        xordercounter(0:maxorder),xopcounter,
     *        xoperleng(maxoperators),xoporder(maxoperators),
     *        xevenodd(maxoperators)
      real*8 doperators(maxoperators),
     *                 xdoperators(maxoperators)
      character*(maxlength) operators(maxoperators),
     *                      xoperators(maxoperators)
c
      logical dkhscfflg
      integer i,j,idum,ddistribution(0:100),xddistribution(0:100)
C     real*8 DABS,DBLE,dscr
      real*8 dscr
c
      if (dkhscfflg) then
        write (unit,1023)
1023    format (2X,88('-'),//9X,'order',6X,'ordercounter',/9X,'(tot)')
        do 20 i=0,dkhorder
          write (unit,1025) i,ordercounter(i)
1025      format (10X,I2,8X,I7)
  20    continue
        write(unit,1027) opcounter
1027    format(20X,7('-'),/20X,I7)
      else
        write (unit,1033)
1033    format (2X,88('-'),//9X,'order',6X,'ordercounter',5X,
     *          'xordercounter',/10X,'(V)')
        do 30 i=0,max(dkhorder,xorder)
          write (unit,1035) i,ordercounter(i),xordercounter(i)
1035      format (10X,I2,8X,I7,11X,I7)
  30    continue
        write(unit,1045) opcounter,xopcounter
1045    format(20X,7('-'),10X,8('-'),/20X,I7,11X,I7)
      endif
c
      do 40 i=0,100
        ddistribution(i)=0
        xddistribution(i)=0
  40  continue
      idum=0
      dscr=1.0d0
      do 45 i=0,20
        if (abs(dscr-dkhzero).le.dkhzero) then
          idum=i
          goto 47
        endif
        dscr=dscr/10.0d0
  45  continue
  47  continue
      do 50 i=1,opcounter
        if (abs(doperators(i)).lt.dkhzero)
     *        ddistribution(idum+1)=ddistribution(idum+1)+1
        dscr=dkhzero
        do 70 j=idum,1,-1
          if (abs(doperators(i)).ge.dscr .and.
     *          abs(doperators(i)).lt.dscr*10.0d0)
     *             ddistribution(j)=ddistribution(j)+1
          dscr=dscr*10.0d0
  70    continue
        if (abs(doperators(i)).eq.1.0d0)
     *                    ddistribution(0)=ddistribution(0)+1
  50  continue
c
c
      if (.not.dkhscfflg) then
        do 51 i=1,xopcounter
          if (abs(xdoperators(i)).lt.dkhzero)
     *          xddistribution(idum+1)=xddistribution(idum+1)+1
          dscr=dkhzero
          do 71 j=idum,1,-1
            if (abs(xdoperators(i)).ge.dscr .and.
     *            abs(xdoperators(i)).lt.dscr*10.0d0)
     *               xddistribution(j)=xddistribution(j)+1
            dscr=dscr*10.0d0
  71      continue
          if (abs(xdoperators(i)).eq.1.0d0)
     *                      xddistribution(0)=xddistribution(0)+1
  51    continue
      endif
c
      if (dkhscfflg) then
        write (unit,1043)
1043    format (/2X,'Distribution of coefficients (doperators):',
     *          3X,'(Absolute values!)',/2X,63('-'),/2X)
        write (unit,1046) dkhzero,ddistribution(idum+1)
1046    format (5X,'[ 0.0     ,',D8.1,' [',5X,I6)
        do 77 i=idum,1,-1
          write (unit,1048) 10.0d0**(-i),10.0d0**(-i+1),ddistribution(i)
1048      format (5X,'[',D8.1,' ,',D8.1,' [',5X,I6)
  77    continue
        write(unit,1049) ddistribution(0),opcounter
1049    format(14X,' = 1.00',10X,I6,/30X,7('-'),/30X,I7)
      else
        write (unit,1053)
1053    format (/2X,'Distribution of coefficients (doperators) and ',
     *          '(xdoperators):',3X,'(Absolute values!)',/2X,60('-'),
     *          /2X)
        write (unit,1056) dkhzero,ddistribution(idum+1),
     *                    xddistribution(idum+1)
1056    format (5X,'[ 0.0     ,',D8.1,' [',5X,I6,12X,I6)
        do 80 i=idum,1,-1
          write (unit,1085) 10.0d0**(-i),10.0d0**(-i+1),
     *                         ddistribution(i),x ddistribution(i)
1085      format (5X,'[',D8.1,' ,',D8.1,' [',5X,I6,12X,I6)
  80    continue
        write(unit,1095) ddistribution(0),xddistribution(0),opcounter,
     *                   xopcounter
1095    format(14X,' = 1.00',10X,I6,12X,I6/30X,7('-'),11X,7('-')/30X,I7,
     *         11X,I7)
      endif
c
      do 140 i=0,100
        ddistribution(i)=0
        xddistribution(i)=0
 140  continue
c
      idum=(maxlength/10)+1
      do 150 i=1,opcounter
        do 170 j=1,idum
          if (operleng(i).ge.10*(j-1) .and. operleng(i).lt.10*j)
     *             ddistribution(j)=ddistribution(j)+1
 170    continue
 150  continue
c
      do 152 i=1,xopcounter
        do 172 j=1,idum
          if (xoperleng(i).ge.10*(j-1) .and. xoperleng(i).lt.10*j)
     *             xddistribution(j)=xddistribution(j)+1
 172    continue
 152  continue
c
      if (dkhscfflg) then
        write (unit,2052)
2052    format (/2X,'Distribution of length of operators (operleng):',
     *          2X,/2X,47('-'),/2X)
        do 179 j=1,idum
          write (unit,2086) 10*(j-1),10*j,ddistribution(j)
2086      format (10X,'[ ',I3,' , ',I3,' [',8X,I6,12X,I6)
 179    continue
        write(unit,2096) opcounter
2096    format(30X,7('-'),/30X,I7)
      else
        write (unit,2053)
2053    format (/2X,'Distribution of length of operators (operleng) ',
     *          'and xoperators (xoperleng):',2X,/2X,74('-'),/2X)
        do 180 j=1,idum
          write (unit,2085) 10*(j-1),10*j,ddistribution(j),
     *                      xddistribution(j)
2085      format (10X,'[ ',I3,' , ',I3,' [',8X,I6,12X,I6)
 180    continue
        write(unit,2095) opcounter,xopcounter
2095    format(30X,7('-'),11X,7('-')/30X,I7,11X,I7)
      endif
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer_array(oporder)
        call Unused_integer_array(evenodd)
        call Unused_character(operators)
        call Unused_integer_array(xoporder)
        call Unused_integer_array(xevenodd)
        call Unused_character(xoperators)
      end if
      end
