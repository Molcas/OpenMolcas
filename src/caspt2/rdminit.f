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
      subroutine rdminit()

      use caspt2_global, only:iPrGlb
      use caspt2_global, only: CMO, CMO_Internal, DREF, DMIX, DWGT, NCMO
      use caspt2_global, only: LUONEM
      use PrintLevel, only: debug
      use stdalloc, only: mma_allocate, mma_deallocate
      implicit real(8) (A-H,O-Z)

#include "caspt2.fh"
#include "pt2_guga.fh"

      REAL*8, ALLOCATABLE:: CI(:)

      if (IPRGLB.GE.DEBUG) then
        write(6,*)' Entered rdminit.'
      end if

* Get CASSCF MO coefficients
      call mma_allocate(CMO_Internal,NCMO,Label='CMO_Internal')
      CMO=>CMO_Internal
      IDISK=IAD1M(1)
      call ddafile(LUONEM,2,CMO,NCMO,IDISK)

* Allocate memory for CI vector
      call mma_allocate(CI,Nconf,Label='CI')

* Initialize array of 1-RDMs with zeros
      DMIX(:,:)=0.0D0

* Start long loop over all states and compute the weighted density
* of each state using the weights in DWGT
      do I=1,Nstate

        if (ISCF.NE.0) then
* Then we still need the "CI array": It is used in subroutine calls
          CI(1)=1.0D0
        else
* Get the CI array
          call loadCI(CI,I)
        end if

* Compute 1-particle active density matrix GAMMA1
        call POLY1(CI,nConf)
* Restructure GAMMA1 as DREF array, but keep it in DMIX
        call GETDREF(DREF,SIZE(DREF))

* Loop over states to compute the contribution of state I to states J
        do J=1,Nstate
* Retrieve the weight of the contribution of state I to the density
* of state J
          wij = DWGT(I,J)
* Multiply density of state I with weight wij and add it to whatever
* is already in DMIX (contributions of other states already computed)
* and store it in DMIX
          call daxpy_(SIZE(DREF),wij,DREF,1,DMIX(:,J),1)
        end do

* End of long loop over states
      end do

* Deallocate everything
      call mma_deallocate(CMO_Internal)
      nullify(CMO)
      call mma_deallocate(CI)

      return
      end
