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
* Copyright (C) 2021, Ignacio Fdez. Galvan                             *
************************************************************************
!  ExpKap
!
!> @brief Compute an orbital rotation matrix from the rotation parameters
!> @author Ignacio Fdez. Galv&aacute;n
!>
!> @details
!> Computes the orbital rotation matrix corresponding to the parametrized
!> form as an antisymmetric matrix, i.e. its exponential.
!> The input parameters \p kapOV are the unique elements of the occupied-virtual
!> block of each symmetry.
!>
!> @param[in]  kapOV  Parameters of the antisymmetric matrix
!> @param[out] U      Unitary matrix to transform old CMOs
!> @param[in]  mynOcc Number of occupied orbitals (including frozen) in each symmetry
!***********************************************************************

      SubRoutine ExpKap(kapOV,U,mynOcc)
*
      Implicit None
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
*
*     declaration subroutine parameters
      Real*8 kapOV(nOV),U(nOFS)
      Integer mynOcc(8)
*
      Integer iKap,iSym,iU,j,jU,mOcc,mOrb,mVir
      Real*8 Cpu1,Cpu2,Tim1,Tim2,Tim3,theta
*
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
*
      iU = 1
      iKap = 1
      U(:) = Zero
      Do iSym=1,nSym
        mOrb = nOrb(iSym)-nFro(iSym)
        mOcc = mynOcc(iSym)-nFro(iSym)
        mVir = mOrb-mOcc
        If (mVir*mOcc == 0) Cycle
        ! Put the non-zero values in the occ-vir nondiagonal block
        jU = iU+mOcc
        Do j=1,mOcc
          U(jU:jU+mVir-1) = kapOV(iKap:iKap+mVir-1)
          iKap = iKap+mVir
          jU = jU+mOrb
        End Do
        ! Compute the exponential
        Call Exp_Schur(mOrb,U(iU),theta)
        iU = iU+mOrb**2
      End Do
*
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld(10) = TimFld(10) + (Cpu2 - Cpu1)
      Return
      End
