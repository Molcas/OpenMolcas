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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_IntegralErrorStat(Mode,AB,CD,Stat)
C
C=======================================================================
C          WARNING: THIS ROUTINE MAY REQUIRE A LOT OF MEMORY !!!
C=======================================================================
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: compute integral error statistics for LDF.
C              On exit, array Stat contains the following info:
C
C              Stat(1): min error
C              Stat(2): max error
C              Stat(3): min abs error
C              Stat(4): max abs error
C              Stat(5): average error
C              Stat(6): avrage abs error
C              Stat(7): std dev with respect to average error
C              Stat(8): std dev with respect to average abs error
C              Stat(9): RMS error
C
C     Input:
C       Mode   --- 1: robust, 2: nonrobust, 3: half-and-half
C       AB, CD --- atom pairs for which the statistics is to be computed
C
C     It is assumed that all integrals (AB|CD) can be stored in core (!)
C     LDF info must be properly set up before calling this routine.
C
      Implicit None
      Integer Mode
      Integer AB
      Integer CD
      Real*8  Stat(9)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Logical Add

      Integer A, B, C, D
      Integer nAB, nCD
      Integer ip_Int, l_Int

      Real*8 xmin, xmax, xamin, xamax, average, abs_average
      Real*8 std, abs_std, rms

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Init Stat
      Call Cho_dZero(Stat,9)

      ! Get atoms
      A=AP_Atoms(1,AB)
      B=AP_Atoms(2,AB)
      C=AP_Atoms(1,CD)
      D=AP_Atoms(2,CD)

      ! Get pair dimensions
      nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
      nCD=LDF_nBas_Atom(C)*LDF_nBas_Atom(D)

      ! Return if nothing to do
      l_Int=nAB*nCD
      If (l_Int.lt.1) Return

      ! Allocate integrals (u_A v_B | k_C l_D)
      Call GetMem('IESInt','Allo','Real',ip_Int,l_Int)

      ! Compute integrals difference, exact - approximate
      Call LDF_ComputeValenceIntegrals(AB,CD,l_Int,Work(ip_Int))
      Call dScal_(l_Int,-1.0d0,Work(ip_Int),1)
      Add=.True.
      Call LDF_ComputeApproximateIntegrals(Mode,Add,AB,CD,l_Int,
     &                                     Work(ip_Int))

      ! Compute statistics
      xmin=9.9d99
      xmax=-9.9d99
      xamin=9.9d99
      xamax=-9.9d99
      average=0.0d0
      abs_average=0.0d0
      rms=0.0d0
      Do i=0,l_Int-1
         xmin=min(xmin,Work(ip_int+i))
         xmax=max(xmax,Work(ip_Int+i))
         xamin=min(xamin,abs(Work(ip_int+i)))
         xamax=max(xamax,abs(Work(ip_Int+i)))
         average=average+Work(ip_Int+i)
         abs_average=abs_average+abs(Work(ip_Int+i))
         rms=rms+Work(ip_Int+i)**2
      End Do
      average=average/dble(l_Int)
      abs_average=abs_average/dble(l_Int)
      rms=sqrt(rms/dble(l_Int))
      std=0.0d0
      abs_std=0.0d0
      Do i=0,l_int-1
         std=std+(Work(ip_Int+i)-average)**2
         abs_std=abs_std+(abs(Work(ip_Int+i))-abs_average)**2
      End Do
      std=sqrt(std/dble(l_int))
      abs_std=sqrt(abs_std/dble(l_int))

      Stat(1)=xmin
      Stat(2)=xmax
      Stat(3)=xamin
      Stat(4)=xamax
      Stat(5)=average
      Stat(6)=abs_average
      Stat(7)=std
      Stat(8)=abs_std
      Stat(9)=rms

      ! Deallocate
      Call GetMem('IESInt','Free','Real',ip_Int,l_Int)

      End
