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
      SUBROUTINE Backtransform(Heff,Ueff,U0)
      USE REFWFN
      IMPLICIT REAL*8 (A-H,O-Z)
C Back-transform Heff and Ueff to the basis of the original
C CASSCF states.
#include "rasdim.fh"
#include "caspt2.fh"
#include "stdalloc.fh"
      real(8) Heff(Nstate,Nstate),Ueff(Nstate,Nstate),U0(Nstate,Nstate)
      real(8),allocatable :: U0transpose(:,:),Utmp(:,:)


      if (IFXMS) then

* First we need to back-transform the effective Hamiltonian in the
* basis of original CASSCF states by U0 * Heff * U0^T
* Note that in the case of a normal MS-CASPT2 this and the next step
* do not have any effect on Heff and Ueff
        call mma_allocate(U0transpose,Nstate,Nstate)
        call trnsps(Nstate,Nstate,U0,U0transpose)
        call transmat(Heff,U0transpose,Nstate)
        call mma_deallocate(U0transpose)

* Compute transformation matrix that diagonalizes the effective
* Hamiltonian expressed in the basis of original CASSCF states,
* i.e. simply combine the two transf matrices: Ueff = U0 * Ueff
        call mma_allocate(Utmp,Nstate,Nstate)
        call dgemm_('N','N',Nstate,Nstate,Nstate,
     &               1.0d0,U0,Nstate,Ueff,Nstate,
     &               0.0d0,Utmp,Nstate)
        Ueff=Utmp
        call mma_deallocate(Utmp)

      end if


      RETURN
      END
