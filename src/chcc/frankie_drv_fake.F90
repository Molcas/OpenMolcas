!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
        subroutine frankie_drv_fake (NChHere)
!
        implicit none
!
#include "chcc1.fh"
#include "cholesky.fh"
#include "choorb.fh"
#include "WrkSpc.fh"
!
        integer NChHere,rc
        real*8 FracMem
!
! ----------------------------------------------------------------
        FracMem=0.0d0
        Call Cho_X_init(rc,FracMem) ! initialize cholesky info
!
!        take local # of Cholesky Vectors on this node
        NChHere=NumCho(1)

        Call Cho_X_final(rc)
!
        return
        end
