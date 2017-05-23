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
      SUBROUTINE initial_surfacehop()
      use Tully_variables
      implicit none
#include "surfacehop.fh"
      logical Found
      real*8  DT
      tullyL=.false.
      decoherence=.false.
      DECO=0.0
      Ethreshold=huge(Ethreshold)
      RandThreshold=0.0
      tullySubVerb=.false.
      fixedrandL=.false.
      FixedRand=-1.0
      lH5Restart= .False.

      call qpg_dscalar('Timestep',Found)
      if (Found) then
         call Get_dScalar('Timestep',DT)
         NSUBSTEPS=Int(200*DT/41.3)
      else
         NSUBSTEPS=0
      end if
      end subroutine

