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
      logical function equalstring (reslengl1,reslengr1,rescharl1,
     *                rescharr1,reslengl2,reslengr2,rescharl2,rescharr2)
c
c***************************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 26.06.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c***************************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      integer reslengl1,reslengr1,reslengl2,reslengr2
      character*(maxlength) rescharl1,rescharr1,rescharl2,rescharr2
      logical equaldummy
c
      equaldummy=.true.
c
      if (reslengl1.ne.reslengl2) equaldummy=.false.
      if (reslengr1.ne.reslengr2) equaldummy=.false.
      if (rescharl1.ne.rescharl2) equaldummy=.false.
      if (rescharr1.ne.rescharr2) equaldummy=.false.
c
      equalstring=equaldummy
c
      return
      end
