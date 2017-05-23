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
* Copyright (C) 1997, Jeppe Olsen                                      *
************************************************************************
*
* Testing new routines with active/passive divisions
*                          Jeppe, sept. 20 97
*
* Splitting the strings into active and passive parts allow us to
* write the a b loop  as
*
* Sigma(I_pa_a,I_pa_b,I_ac_a,I_ac_b) =
* sum(ijkl) <I_ac_a!Ea(ij)!J_ac_b><I_ac_b!Eb(kl)!I_ac_b>(ij!kl)
*           C(J_pa_a,J_pa_b,I_ac_a,I_pa_b)
*. One possibility is to use the above directly, and vectorize
*  over the passive strings.
*
* Another possibility is to use N-1 resolution in the alpha space
* as
*   1) Construct C(Ka_a,I_pa_a,I_pa_b,j,J_ac_b)
*   2) Obtain    S(Ka_a,I_pa_a,I_pa_b,I,I_ac_b) =
*      Sum(kl)   <I_ac_b!Eb(kl)!J_ac_b> Sum(ij)(ij!kl)
*      C(Ka_a,I_Pa_a,I_Pa_b,j,J_ac_b)
*
* The latter is the straightforward extension of LUCILLE's normal route.
*
* About reorganization of strings
*
* I will construct an array IREO(IA,IP)=I, that for given active
* and passive parts gives original number.
* Invoking symmetry :
* We are reordering strings of given supergroup and symmetry.
* The readressing goes thus as
* Loop over symmetry of active part => symmetry of passive part
*   Obtain offset for reordered strings
*          Number of active and passive strings of appropriate sym
*          Reorder array.
*
      SUBROUTINE MAT_P_MATT(A,B,NR,NC,COEF)
*
* A(I,J) = A(I,J) + Coef*B(J,I)
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input and output
      DIMENSION A(NR,NC)
*. Input
      DIMENSION B(NC,NR)
*
      DO J = 1, NC
      DO I = 1, NR
        A(I,J) = A(I,J) + COEF*B(J,I)
      END DO
      END DO
*
      RETURN
      END
