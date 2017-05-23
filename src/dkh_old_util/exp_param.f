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
      subroutine exp_param (ducoeffs)
c
c***********************************************************************
c
c   This SR belongs to dkhparser_symbolic (dkhparser1).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 26.06.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Jena)
c
c***********************************************************************
c
      implicit none
#include "dkhparameters.fh"
c
      REAL*8 ducoeffs(maxorder)
      integer i
      intrinsic DBLE
      ducoeffs(1)=1.0d0
      do 10 i=2,maxorder
        ducoeffs(i)=ducoeffs(i-1)/DBLE(i)
  10  continue
c
      return
      end
