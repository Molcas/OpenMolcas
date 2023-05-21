!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
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
!> @param[in]  nkapOV number of elements in kapOV
!> @param[out] U      Unitary matrix to transform old CMOs
!> @param[in]  mynOcc Number of occupied orbitals (including frozen) in each symmetry
!***********************************************************************

      SubRoutine ExpKap(kapOV,nKapOV,U,mynOcc)
      use InfSCF, only: nOFs, nSym, nFro, TimFld, nOrb
      Use Constants, only: Pi, Zero
!
      Implicit None

#define  qnext
!
!     declaration subroutine parameters
      Integer nKapOV
      Real*8 kapOV(nkapOV),U(nOFS)
      Integer mynOcc(8)
!
      Integer iKap,iSym,iU,j,jU,mOcc,mOrb,mVir
      Real*8 Cpu1,Cpu2,Tim1,Tim2,Tim3

#ifndef qnext
      Real*8 theta
#endif

      Real*8, Parameter :: Thrs = 1.0D-14

      Do j = 1, nKapOV
         If (Abs(KapOV(j))>Pi) Then
            Write (6,*) 'ExpKap: KapOV too large:',KapOV(j)
            Call Abend()
         End If
      End Do
      Call Timing(Cpu1,Tim1,Tim2,Tim3)
!
      iU = 1
      iKap = 1
      U(:) = Zero

      Do iSym=1,nSym
        mOrb = nOrb(iSym)-nFro(iSym)
        mOcc = mynOcc(iSym)-nFro(iSym)
        mVir = mOrb-mOcc

        If (mVir*mOcc == 0) Cycle

        jU = iU+mOcc

        Do j=1,mOcc
          U(jU:jU+mVir-1) = kapOV(iKap:iKap+mVir-1)
          iKap = iKap+mVir
          jU = jU+mOrb
        End Do

#ifdef  qnext
        Call matexp(mOrb,mOcc,U(iU:iU+mOrb**2))
#else
        Call Exp_Schur(mOrb,U(iU:iU+mOrb**2),theta)
#endif

        iU = iU+mOrb**2
      End Do
!
      Do j=1,nOFS
        If (abs(U(j)).lt.Thrs) U(j) = Zero
      End do
!
      Call Timing(Cpu2,Tim1,Tim2,Tim3)
      TimFld(10) = TimFld(10) + (Cpu2 - Cpu1)

      Return
      End SubRoutine ExpKap
