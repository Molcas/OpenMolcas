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
      subroutine output2b (unit,dkhscfflg,ucounter,uopsleng,uoporder,
     *                     eouops,duops,uops)
c
c*******************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 15.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      logical dkhscfflg
      integer unit,ucounter,uopsleng(maxuops),uoporder(maxuops),
     *        eouops(maxuops)
      REAL*8 duops(maxuops)
      character*(maxlength) uops(maxuops)
      integer j
c
      if (dkhscfflg) then
        write (unit,1012)
1012    format (/6X,'#',2X,'length',92X,'order',2X,'symm.',2X,'coeff.',
     *          /107X,'(tot)',/)
      else
        write (unit,1013)
1013    format (/6X,'#',2X,'length',92X,'order',2X,'symm.',2X,'coeff.',
     *          /108X,'(V)',/)
      endif
      do 10 j=1,ucounter
        write (unit,1015) j,uopsleng(j),uops(j)(1:uopsleng(j)),
     *                       uoporder(j),eouops(j),duops(j)
1015    format (I7,2X,I3,2X,A90,4X,I2,5X,I2,3X,F17.14)
  10  continue
c
      return
      end
