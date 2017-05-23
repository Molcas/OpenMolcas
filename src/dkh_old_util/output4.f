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
      subroutine output4 (dkhscfflg,unit,sused,stimes,sorder,s,scrchar,
     *                    stimestot)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 16.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      logical dkhscfflg
      integer unit,sused,stimes(maxsnumber),sorder(maxsnumber,3),
     *        stimestot
      character*(4) s(maxsnumber)
      character*(9) scrchar(maxsnumber)
      integer j,idum
c
      idum=0
      do 25 j=1,sused
        if (stimes(j).gt.0) then
          idum=idum+1
          if (dkhscfflg) then
            write (unit,1081) idum,scrchar(j),s(j),sorder(j,3),stimes(j)
1081        format (3X,I3,5X,A9,6X,A4,5X,I2,7X,I5)
          else
            write (unit,1083) idum,scrchar(j),s(j),sorder(j,1),
     *                        sorder(j,2),sorder(j,3),stimes(j)
1083        format (3X,I3,5X,A9,6X,A4,5X,I2,4X,I2,4X,I2,4X,I5)
          endif
        endif
  25  continue
c
      if (dkhscfflg) then
        write (unit,1091) stimestot
1091    format (42X,7('-'),/42X,I7,/2X)
      else
        write (unit,1092) stimestot
1092    format (51X,7('-'),/51X,I7,/2X)
      endif
c
      return
      end
