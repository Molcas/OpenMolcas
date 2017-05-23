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
      subroutine store_reschar (termleng,term,hit,scrleng,
     *                          reslengl,reslengr,rescharl,rescharr)
c
c***************************************************************************
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
c***************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer termleng,hit,scrleng,reslengl,reslengr
      character*(maxlength) term,rescharl,rescharr
      integer k
      do 10 k=1,maxlength
        rescharl(k:k)=' '
        rescharr(k:k)=' '
  10  continue
      reslengl=hit-1
      reslengr=termleng-(reslengl+scrleng)
      if (reslengl.ne.0) rescharl(1:reslengl)=term(1:reslengl)
      if (reslengr.ne.0) rescharr(1:reslengr)=term(hit+scrleng:termleng)
c
      return
      end
