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
      subroutine output12 (unit,ttcounter,ttorder1,ttorder2,ttorder3,
     *                     tt,ttchar,termcounter,termleng,termorder,
     *                     dtcoeff,term)
c
c****************************************************************
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
c*****************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer unit,ttcounter,ttorder1,ttorder2,ttorder3,termcounter,
     *        termleng(maxoperators),termorder(maxoperators,3)
      REAL*8 dtcoeff(maxoperators)
#if defined(_MOLCAS_) || defined(MOLPRO)
      integer term(*)
      character*(maxlength) termstr
#else
      character*(maxlength) term(maxoperators)
#endif
      character*(4) tt
      character*(11) ttchar
      integer j
c
      write (unit,1001) tt,ttchar,ttorder1,ttorder2,ttorder3,ttcounter
1001  format ('***',/A4,1X,A11,1X,I2,1X,I2,1X,I2,2X,I7)
      do 10 j=1,ttcounter
#if defined(_MOLCAS_) || defined(MOLPRO)
        call get_dkoperators_i(j,termstr,term)
        write (unit,2001) j,termleng(j),termorder(j,1),termorder(j,2),
     *                    termorder(j,3),termstr(1:termleng(j)),
     *                    dtcoeff(j)
#else
        write (unit,2001) j,termleng(j),termorder(j,1),termorder(j,2),
     *                    termorder(j,3),term(j)(1:termleng(j)),
     *                    dtcoeff(j)
#endif
2001    format (I7,2X,I3,5X,I2,1X,I2,1X,I2,1X,A90,4X,F17.14)
  10  continue
      write (unit,3001)
3001  format (/2X)
c
      return
c Avoid unused argument warnings
      if (.false.) call Unused_integer(termcounter)
      end
