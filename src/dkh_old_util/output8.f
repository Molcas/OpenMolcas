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
      subroutine output8 (unit,dkhscfflg,uused,unumber,utimes,uorder,
     *                    umult,u,uscrchar,uscrleng,utimestot,umulttot)
c
c***********************************************************************
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
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      logical dkhscfflg
      integer unit,uused,unumber,utimes(maxunumber),
     *        uorder(maxunumber,3),umult(maxunumber),
     *        uscrleng(maxunumber),utimestot,umulttot
      character*(4) u(maxunumber)
      character*(3) uscrchar(maxunumber)
      integer j,idum
c
      if (dkhscfflg) then
        write (unit,1070)
1070    format (/2X,'Statistic about Uxxx substitution procedure in ',
     *        'steps 3--5:',/2X,58('-'),//5X,'Note:  * All values ',
     *        'refer to low-level situation of steps 3--5.',/12X,
     *        '* All substitutions for operators, Sxxx, and Txxx ',
     *        'have been counted.')
      else
        write (unit,1071)
1071    format (/2X,'Statistic about Uxxx substitution procedure in ',
     *        'steps 3--5:',/2X,58('-'),//5X,'Note:  * All values ',
     *        'refer to low-level situation of steps 3--5.',/12X,
     *        '* All substitutions for (x)operators, Sxxx, and Txxx ',
     *        'have been counted.')
      endif
c
        write (unit,1072)
1072    format (//5X,'#',4X,'steps',3X,'-->',3X,'steps',4X,
     *         'order',1X,'order',1X,'order',3X,'no. of subs.',3X,
     *         'no. of mat.mult.',/10X,'1--2',10X,'3--5',6X,
     *         '(V)',3X,'(X)',2X,'(tot)',
     *         5X,'(utimes)',9X,'(umult)',/)
c
      idum=0
      do 25 j=1,uused
        if (utimes(j).gt.0) then
          idum=idum+1
          write (dbgunit,1081) idum,uscrchar(j),u(j),uorder(j,1),
     *                       uorder(j,2),uorder(j,3),utimes(j),umult(j)
1081      format (3X,I3,4X,A3,11X,A4,6X,I2,4X,I2,4X,I2,5X,I8,9X,I5)
        endif
  25  continue
c
      write (unit,1091) utimestot,umulttot
1091  format (53X,8('-'),9X,5('-'),/53X,I8,9X,I5,/2X)
c
      return
c Avoid unused argument warnings
      if (.false.) then
        call Unused_integer(unumber)
        call Unused_integer_array(uscrleng)
      end if
      end
