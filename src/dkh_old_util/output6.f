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
      subroutine output6 (unit,dkhscfflg,sused,scounter,stimes,sorder,
     *                    smult,s,scrchar,stimestot,smulttot)
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
      integer unit,sused,scounter(maxsnumber),stimes(maxsnumber),
     *        sorder(maxsnumber,3),smult(maxsnumber),stimestot,smulttot
      character*(4) s(maxsnumber)
      character*(9) scrchar(maxsnumber)
      integer j,idum
c
      write (unit,1071)
1071  format (//2X,'Statistic about evaluation of Sxxx-expressions in ',
     *        'step 4:',/2X,58('-'),//5X,'Note: All values refer to',
     *        ' low-level situation of step4.')
c
      if (dkhscfflg) then
        write (unit,1072)
1072    format (//5X,'#',6X,'step1',3X,'-->',3X,'step4',4X,'order',3X,
     *          'no. of subs.',3X,'no. of terms',3X,'no. of mat.mult.',
     *          /35X,'(tot)',5X,'(stimes)',7X,'(scounter)',7X,'(smult)',
     *          10X,/)
      else
        write (unit,1073)
1073    format (//5X,'#',6X,'step1',3X,'-->',3X,'step4',4X,'order',1X,
     *        'order',1X,'order',3X,'no. of subs.',3X,'no. of terms',3X,
     *        'no. of mat.mult.',/36X,'(V)',3X,'(X)',2X,'(tot)',5X,
     *        '(stimes)',7X,'(scounter)',7X,'(smult)',/)
      endif
c
      idum=0
      do 25 j=1,sused
        if (stimes(j).gt.0) then
          idum=idum+1
          if (dkhscfflg) then
            write (dbgunit,1081) idum,scrchar(j),s(j),sorder(j,3),
     *                           stimes(j),scounter(j),smult(j)
1081        format (3X,I3,5X,A9,6X,A4,6X,I2,6X,I6,7X,I8,8X,I9)
          else
            write (dbgunit,1082) idum,scrchar(j),s(j),sorder(j,1),
     *                           sorder(j,2),sorder(j,3),stimes(j),
     *                           scounter(j),smult(j)
1082        format (3X,I3,5X,A9,6X,A4,6X,I2,4X,I2,4X,I2,6X,I8,6X,I8,7X,
     *              I9)
          endif
        endif
  25  continue
c
      if (dkhscfflg) then
        write (unit,1090) stimestot,smulttot
1090    format (42X,8('-'),23X,9('-'),/42X,I8,23X,I9,/2X)
      else
        write (unit,1091) stimestot,smulttot
1091    format (56X,8('-'),21X,9('-'),/56X,I8,21X,I9,/2X)
      endif
c
      return
      end
