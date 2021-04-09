************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
       subroutine setb (wrk,wrksize,
     & mapda,mapdb,factor)
c
c     this routine do
c     B = factor . A
c
c     mapda  - direct map of A (I)
c     mapdb  - direct map of B (I)
c     factor - numerical factor (I)
c
c     mediate B must have defined maps, and their must be
c     of identicat type as those for A. If they are not defined,
c     use grc0 before setb
c
c     N.B. this routine should be done using matrix operations

#include "wrk.fh"
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       real*8 factor
c
c     help variables
c
       integer possa0,possb0,length,nhelp
c
c
c1    def the length of the mediate
       nhelp=mapda(0,5)
       length=mapda(nhelp,1)+mapda(nhelp,2)-mapda(1,1)
       if (length.eq.0) return
c
c2    def initial possitions
       possa0=mapda(1,1)
       possb0=mapdb(1,1)
c
c3    set B=f.A
       do 100 nhelp=0,length-1
       wrk(possb0+nhelp)=factor*wrk(possa0+nhelp)
 100    continue
c
       return
       end
