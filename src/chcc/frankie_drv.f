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
c
c this file contains :
c
c frankie_drv
c frankie_drv_fake
c
c ----------------------------------------------------------------
c
        subroutine frankie_drv (NChHere)
c
        implicit none
c
#include "chcc1.fh"
#include "cholesky.fh"
#include "WrkSpc.fh"
c
        integer NChHere
c
c ----------------------------------------------------------------
c
c - transform AO Cholesky vectors to MO basis localy on each node
c   with fragmented cholesky index. _AI1 -> _CDtmp1
c
        call frankie(nfr,no,nv,printkey)
c
c        take local # of Cholesky Vectors on this node
        NChHere=NumCho(1)
c
        return
        end
c
c -------------------------------------------------------------------
c
        subroutine frankie_drv_fake (NChHere)
c
        implicit none
c
#include "chcc1.fh"
#include "cholesky.fh"
#include "choptr.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
c
        integer NChHere,rc
        real*8 FracMem
c
c ----------------------------------------------------------------
        FracMem=0.0d0
        Call Cho_X_init(rc,FracMem) ! initialize cholesky info
c
c        take local # of Cholesky Vectors on this node
        NChHere=NumCho(1)

        Call Cho_X_final(rc)
c
        return
        end
