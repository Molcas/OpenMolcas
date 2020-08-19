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
* Copyright (C) 2020, Ignacio Fdez. Galvan                             *
************************************************************************
*  Random_Vector
*
*> @brief
*>   Generate a random vector in an \f$ N \f$ dimensional hypersphere
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Generate a random \f$ N \f$ dimensional vector of unit length or less,
*> by using random normal-distributed numbers. \cite Mul1959-CACM-2-19
*> The normal-distributed random numbers are generated with the Box--Muller
*> transform. \cite Box1958-AMS-29-610
*>
*> @param[in]  N    Dimension of the generated vector
*> @param[out] Vec  Generated random vector
*> @param[in]  UVec Specifies whether the vector should be of unit length
************************************************************************
      SUBROUTINE Random_Vector(N,Vec,UVec)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: N
      REAL*8, INTENT(OUT) :: Vec(N)
      LOGICAL, INTENT(IN) :: UVec
#include "real.fh"
      INTEGER :: i
      INTEGER, SAVE :: iSeed=0
      REAL*8 :: U, V, X, Y, m, sm, tot_m
      REAL*8, EXTERNAL :: Random_Molcas
      REAL*8, PARAMETER :: Thr=1.0D-8
*
*     Initialize random seed
      IF (iSeed==0) CALL GetSeed(iSeed)
*
*     To reduce numerical errors, repeat until the size is reasonable
      tot_m = Zero
      DO WHILE ((tot_m < Thr) .OR. (tot_m > One/Thr))
        tot_m = Zero
        DO i=1,N,2
*         Get two independent normal distributed-variales, X and Y
*         See doi:10.1214/aoms/1177706645
          U = Random_Molcas(iSeed)
          V = Two*Pi*Random_Molcas(iSeed)
          m = -Two*LOG(U)
          sm = SQRT(m)
          X = sm*COS(V)
          Y = sm*SIN(V)
*         Add them to the vector,
*         being careful with the last one if N is odd
*         See doi:10.1145/377939.377946
          IF (i==N) THEN
            Vec(i) = X
            tot_m = tot_m + X*X
          ELSE
            Vec(i:i+1) = [X, Y]
            tot_m = tot_m + m
          END IF
        END DO
      END DO
*     Normalize the vector
*     and scale by a random (0,1) number if no unit vector desired
      IF (UVec) THEN
        sm = One
      ELSE
        sm = Random_Molcas(iSeed)
      END IF
      Vec(:) = sm/SQRT(tot_m)*Vec(:)
*
      RETURN
      END
