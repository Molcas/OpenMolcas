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
      subroutine rdminit

      implicit real(8) (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "pt2_guga.fh"
#include "WrkSpc.fh"

      integer offset

      call QENTER('rdminit')

      if (IPRGLB.GE.DEBUG) then
        write(6,*)' Entered rdminit.'
      end if

* Get CASSCF MO coefficients
      call getmem('LCMO','ALLO','REAL',LCMO,NCMO)
      IDISK=IAD1M(1)
      call ddafile(LUONEM,2,WORK(LCMO),NCMO,IDISK)

* Allocate memory for CI vector
      call getmem('LCI','ALLO','REAL',LCI,Nconf)

* Initialize array of 1-RDMs with zeros
      call dcopy_(Nstate*NDREF,[0.0D0],0,WORK(LDMIX),1)

* Start long loop over all states and compute the weighted density
* of each state using the weights in WORK(LDWGT)
      do I=1,Nstate

        if (ISCF.NE.0) then
* Then we still need the "CI array": It is used in subroutine calls
          WORK(LCI)=1.0D0
        else
* Get the CI array
          call loadCI(WORK(LCI),I)
        end if

* Compute 1-particle active density matrix GAMMA1
        call POLY1(WORK(LCI))
* Restructure GAMMA1 as DREF array, but keep it in DMIX
        call GETDREF(WORK(LDREF))

* Loop over states to compute the contribution of state I to states J
        do J=1,Nstate
* Retrieve the weight of the contribution of state I to the density
* of state J
          wij = WORK(LDWGT+(I-1) + NSTATE*(J-1))
* Compute offset to access the density of state J in WORK, needed to
* know where to store the weighted density of state I
          offset = NDREF*(J-1)
* Multiply density of state I with weight wij and add it to whatever
* is already in LDMIX (contributions of other states already computed)
* and store it in LDMIX
          call daxpy_(NDREF,wij,WORK(LDREF),1,WORK(LDMIX+offset),1)
        end do

* End of long loop over states
      end do

* Deallocate everything
      call getmem('LCMO','FREE','REAL',LCMO,NCMO)
      call getmem('LCI','FREE','REAL',LCI,NCONF)

      call QEXIT('rdminit')
      return
      end