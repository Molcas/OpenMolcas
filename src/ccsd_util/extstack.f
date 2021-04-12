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
* Copyright (C) 2006, Pavel Neogrady                                   *
************************************************************************
        subroutine extstack (wrk,wrksize,
     &                 mapda,mapdb,b,dimb)
c
c       This routine do:
c       A(ij) <- B(ij,_b) for given b
c
c       A special routine used only in sumoverab for stacking
c       case.
c
c       Yet it is assumed,that blocks in A and B are in the same
c       order. To je pomerne odflaknuty predpoklad, a moze to
c       byt bugous
c
#include "wrk.fh"
c
       integer mapda(0:512,1:6)
       integer mapdb(0:512,1:6)
       integer b,dimb
c
c     help variables
c
       integer ii,dimij,possa,possb
c
c
        do ii=1,mapda(0,5)
          dimij=mapda(ii,2)
          possa=mapda(ii,1)
          possb=mapdb(ii,1)
          call extstackhlp1 (wrk(possa),wrk(possb),dimij,dimb,b)
        end do
c
        return
        end
