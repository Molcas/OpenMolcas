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
* Copyright (C) 2004, Alexander Wolf                                   *
*               2004,2006, Markus Reiher                               *
************************************************************************
      subroutine finalize2 (length,operator,uused,utimes,u,uscrleng,
     *                      uscrchar)
c
c******************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 11.10.2006 (MR, ETH Zurich)
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c******************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer length,uused,utimes(maxunumber),uscrleng(maxunumber),
     *        i,k,pos,idum,dummyleng
      character*(maxlength) operator
      character*(3) uscrchar(maxunumber)
      character*(4) u(maxunumber)
c
c--------------------------------------------------------------------------
c
1010  continue
      do 10 i=1,uused
        idum=0
        dummyleng=3
        if (uscrleng(i).eq.2) dummyleng=2
        pos=0
        pos=index(operator(1:length),uscrchar(i)(1:dummyleng))
        if (pos.gt.0) then
          idum=4-uscrleng(i)
          if (length+idum.gt.maxlength) then
            write (stdout,1050)
1050        format (2X,'SR finalize2: Parameter maxlength not large',
     *              ' enough!',//2X,'STOP.',/2X)
            CALL Abend
          endif
          do 20 k=length,pos+uscrleng(i),-1
            operator(k+idum:k+idum)=operator(k:k)
  20      continue
          operator(pos:pos+3)=u(i)
          length=length+idum
          utimes(i)=utimes(i)+1
          goto 1010
        endif
  10  continue
c
      return
      end
