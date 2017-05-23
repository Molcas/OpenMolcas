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
      subroutine output7 (unit,dkhscfflg,sused,tcounter,ttimes,sorder,
     *                    tmult,t,tscrchar,ttimestot,tmulttot)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 19.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      logical dkhscfflg
      integer unit,sused,tcounter(maxsnumber),ttimes(maxsnumber),
     *        sorder(maxsnumber,3),tmult(maxsnumber),ttimestot,tmulttot
      character*(4) t(maxsnumber)
      character*(11) tscrchar(maxsnumber)
      integer j,idum
c
      write (unit,1071)
1071  format (//2X,'Statistic about evaluation of Txxx-expressions in ',
     *        'step 5:',/2X,57('-'),//5X,'Note: All values refer to',
     *        ' low-level situation of step5.')
c
      if (dkhscfflg) then
        write (unit,1072)
1072    format (//5X,'#',7X,'step1',4X,'-->',3X,'step5',3X,'order',3X,
     *          'no. of subs.',3X,'no. of terms',3X,'no. of mat.mult.',
     *          /36X,'(tot)',5X,'(ttimes)',5X,'(tcounter)',6X,'(tmult)',
     *          10X,/)
      else
       write (unit,1073)
1073    format (//5X,'#',6X,'step1',3X,'-->',3X,'step5',4X,'order',1X,
     *        'order',1X,'order',3X,'no. of subs.',3X,'no. of terms',3X,
     *        'no. of mat.mult.',/36X,'(V)',3X,'(X)',2X,'(tot)',5X,
     *        '(ttimes)',7X,'(tcounter)',7X,'(tmult)',/)
      endif
c
      idum=0
      do 25 j=1,sused
        if (ttimes(j).gt.0) then
          idum=idum+1
          if (dkhscfflg) then
            write (dbgunit,1081) idum,tscrchar(j),t(j),sorder(j,3),
     *                         ttimes(j),tcounter(j),tmult(j)
1081        format (3X,I3,5X,A11,6X,A4,5X,I2,4X,I8,8X,I6,8X,I8)
          else
            write (dbgunit,1082) idum,tscrchar(j),t(j),sorder(j,1),
     *                           sorder(j,2),sorder(j,3),ttimes(j),
     *                           tcounter(j),tmult(j)
1082        format (3X,I3,5X,A11,4X,A4,6X,I2,4X,I2,4X,I2,6X,I8,6X,I8,7X,
     *              I9)
          endif
        endif
  25  continue
c
      if (dkhscfflg) then
        write (unit,1091) ttimestot,tmulttot
1091    format (43X,8('-'),21X,9('-'),/43X,I8,21X,I9,/2X)
      else
        write (unit,1092) ttimestot,tmulttot
1092    format (56X,8('-'),21X,9('-'),/56X,I8,21X,I9,/2X)
      endif
c
      return
      end
