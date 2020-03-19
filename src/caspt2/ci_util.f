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
* Copyright (C) 2019, Stefano Battaglia                                *
************************************************************************

* Load the CI vector of state Istate from LUCIEX into memory
      subroutine loadCI(CI, Istate)
      implicit real(8) (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
      real(8) CI(Nconf)
      integer ID, Istate

* Skip over states
      ID=IDCIEX
      do I=1,Istate-1
        call ddafile(LUCIEX,0,CI,Nconf,ID)
      end do

* Load the CI array
      call ddafile(LUCIEX,2,CI,Nconf,ID)

      return
      end


************************************************************************


* write the CI vector of state Istate from memory into LUCIEX
      subroutine writeCI(CI, Istate)
      implicit real(8) (A-H,O-Z)
#include "rasdim.fh"
#include "caspt2.fh"
      real(8) CI(Nconf)
      integer ID, Istate

* Skip over states
      ID=IDCIEX
      do I=1,Istate-1
        call ddafile(LUCIEX,0,CI,Nconf,ID)
      end do

* Write the CI array
      call ddafile(LUCIEX,1,CI,Nconf,ID)

      return
      end


