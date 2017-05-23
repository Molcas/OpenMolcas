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
*               2004, Markus Reiher                                    *
************************************************************************
      character*(2) function int2char2 (number)
c
c*************************************************************************
c
c   This SR belongs to dkhparser_numeric (dkhparser2).
c
c   written by:  Alexander Wolf and Markus Reiher  (Univ. Jena)
c
c   version:  2.0.1
c
c   THIS ORIGINAL SUB IS KEPT BECAUSE THE ROUTINE IS CALLED BY sew_molpro.f
c   IT IS IDENTICAL WITH dkh_int2char2 AND MAY BE REMOVED FROM THE PACKAGE
c   IN PROGRAM ENVIRONMENTS OTHER THAN MOLPRO OR MOLCAS
c
c   first version: 10.04.2004  (Theoretical Chemistry, Univ. Bonn)
c
c*************************************************************************
c
      implicit none
c
      integer number,idum(2),j
      character*(2) string
c
c   Translate integer number --> character string of length 2
c
      idum(1)=INT(number/10)
      idum(2)=INT(number-idum(1)*10)
c
      do 100 j=1,2
        if (idum(j).eq.0)  string(j:j)='0'
        if (idum(j).eq.1)  string(j:j)='1'
        if (idum(j).eq.2)  string(j:j)='2'
        if (idum(j).eq.3)  string(j:j)='3'
        if (idum(j).eq.4)  string(j:j)='4'
        if (idum(j).eq.5)  string(j:j)='5'
        if (idum(j).eq.6)  string(j:j)='6'
        if (idum(j).eq.7)  string(j:j)='7'
        if (idum(j).eq.8)  string(j:j)='8'
        if (idum(j).eq.9)  string(j:j)='9'
 100  continue
c
      int2char2(1:2)=string(1:2)
c
      return
      end
