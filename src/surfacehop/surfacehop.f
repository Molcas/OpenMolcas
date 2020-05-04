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
* Copyright (C) 2015, Luis Manuel Frutos                               *
*               2015, Ignacio Fdez. Galvan                             *
*               2015, Alessio Valentini                                *
************************************************************************
      subroutine surfacehop(return)
      use Tully_variables
      implicit none
#include "warnings.fh"
#include "stdalloc.fh"
      integer return
      integer NSTATE, LUIPH, IAD, ITOC15(15)
      integer NCI, IDISK, I
      real*8, allocatable :: CIBigArray(:)

      CALL qEnter('SurfaceHop')

      call initial_surfacehop()
      call rdinp_surfacehop()

      LUIPH=20
      call DANAME(LUIPH,"JOBIPH")
      IAD=0
      call IDAFILE(LUIPH,2,ITOC15,15,IAD)
      call getIphInfo(LUIPH,NCI,NSTATE,ITOC15)
      call mma_allocate(CIBigArray,NCI*NSTATE)
      CIBigArray(:)=0.0D0

      IDISK=ITOC15(4)
      do I=1, NSTATE
         call DDAFILE(LUIPH,2,CIBigArray(NCI*(I-1)+1),NCI,IDISK)
      end do
      call DACLOS(LUIPH)

*      call recprt('CI coefficients',"",CIBigArray,NCI,NSTATE)

      call tully(CIBigArray,NSTATE,NCI)

      call mma_deallocate(CIBigArray)
      CALL qExit('SurfaceHop')
      return=_RC_ALL_IS_WELL_
      end
