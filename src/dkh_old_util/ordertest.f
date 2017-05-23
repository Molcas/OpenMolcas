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
      subroutine ordertest (dkhorder,paramtype,scrtxt,ctrflg)
c
c*************************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2)
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 25.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c*************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer dkhorder,ctrflg,idum,dkh_char2int
      character*(3) paramtype
      character*(maxlength) scrtxt
c
      idum=0
      if (scrtxt(14:14).eq.' ') then
        idum=dkh_char2int(1,scrtxt(15:15))
      else
        idum=dkh_char2int(2,scrtxt(14:15))
      endif
      if (ctrflg.eq.1) then
        if (dkhorder.gt.idum) then
          write (stdout,1010) dkhorder,idum,idum
1010      format ('SR ordertest (1): The desired dkhorder = ',I2,
     *            ' is larger than dkhorder = ',I2,/,'stored in ',
     *            'dkhops.13.',/,'--> Reduce dkhorder:',
     *            ' dkhorder = ',I2,'.',/2X)
          dkhorder=idum
        endif
      else
        if (dkhorder.gt.idum) then
          write (stdout,1014) dkhorder,idum,idum
1014      format ('SR ordertest (2): The desired xorder = ',I2,
     *            ' is larger than xorder = ',I2,/,'stored in ',
     *            'dkhops.13.',/,'--> Reduce xorder:',
     *            ' xorder = ',I2,'.',/2X)
          dkhorder=idum
        endif
      endif
c
      if (ctrflg.eq.1) paramtype(1:3)=scrtxt(26:28)
c
      return
      end
