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
      subroutine output3 (unit,opcounter,operleng,oporder,evenodd,
     *                    doperators,operators)
c
c*******************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 12.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer unit,opcounter,operleng(maxoperators),
     *        oporder(maxoperators,3),evenodd(maxoperators)
      REAL*8 doperators(maxoperators)
      character*(maxlength) operators(maxoperators)
      integer j
c
      write (unit,1012)
1012  format (/6X,'#',2X,'length',90X,'order',1X,'order',1X,
     *        'order',2X,'symm.',2X,'coeff.'/,106X,'(V)',3X,'(X)',
     *        2X,'(tot)',/)
      do 10 j=1,opcounter
        write (unit,1015) j,operleng(j),operators(j)(1:operleng(j)),
     *                    oporder(j,1),oporder(j,2),oporder(j,3),
     *                    evenodd(j),doperators(j)
1015    format (I7,2X,I3,2X,A90,2X,I2,4X,I2,4X,I2,5X,I2,3X,F17.14)
  10  continue
c
      return
      end
