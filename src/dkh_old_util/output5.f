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
      subroutine output5 (dkhscfflg,unit,sused,ttimes,sorder,t,tscrchar,
     *                    ttimestot)
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
      integer unit,sused,ttimes(maxsnumber),sorder(maxsnumber,3),
     *        ttimestot
      character*(4) t(maxsnumber)
      character*(11) tscrchar(maxsnumber)
      integer j,idum
c
      idum=0
      do 25 j=1,sused
        if (ttimes(j).gt.0) then
          idum=idum+1
          if (dkhscfflg) then
            write (unit,1081) idum,tscrchar(j),t(j),sorder(j,3),
     *                        ttimes(j)
1081        format (3X,I3,5X,A11,6X,A4,5X,I2,7X,I5)
          else
            write (unit,1083) idum,tscrchar(j),t(j),sorder(j,1),
     *                        sorder(j,2),sorder(j,3),ttimes(j)
1083        format (3X,I3,5X,A11,6X,A4,5X,I2,4X,I2,4X,I2,4X,I7)
          endif
        endif
  25  continue
c
      if (dkhscfflg) then
        write (unit,1091) ttimestot
1091    format (44X,7('-'),/44X,I7,/2X)
      else
        write (unit,1092) ttimestot
1092    format (54X,8('-'),/54X,I8,/2X)
      endif
c
      return
      end
