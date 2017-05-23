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
      subroutine replace_Sxxx (idum,poss,length,term,scrleng,
     *                         scrchar)
c
c*************************************************************************
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
c*************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer idum,poss,length,scrleng(maxsnumber),dummyleng,reslengl,
     *        reslengr,j
      character*(maxlength) term,rescharl,rescharr,dummychar
      character*(9) scrchar(maxsnumber)
c
      do 10 j=1,maxlength
        dummychar(j:j)=' '
  10  continue
c
      dummyleng=4
      call store_reschar (length,term,poss,dummyleng,reslengl,reslengr,
     *                    rescharl,rescharr)
c
      length=length-4+scrleng(idum)
      if (scrleng(idum).eq.6) then
        dummychar(1:6)=scrchar(idum)(1:6)
        dummyleng=6
      endif
      if (scrleng(idum).eq.9) then
        dummychar(1:9)=scrchar(idum)(1:9)
        dummyleng=9
      endif
c
      call concatenate (length,term,reslengl,rescharl,
     *                  dummyleng,dummychar,reslengr,rescharr)
c
      return
      end
