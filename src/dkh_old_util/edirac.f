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
      real*8 function edirac (n,j,Z,c)
c
c*************************************************************************
c
c   This SR belongs to the test suite 'parser2'
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   last modified: 21.07.2005
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c*************************************************************************
c
      implicit none
c
      real*8 n,j,Z,c
c
      edirac = 0.0d0
      edirac = (Z/c) * 1.d0/(n-(j+0.5d0) +
     *                   sqrt((j+0.5d0)*(j+0.5d0)-(Z*Z)/(c*c)))
      edirac = c*c / sqrt(1.d0+edirac*edirac)
      edirac = edirac - c*c
c
      return
      end
