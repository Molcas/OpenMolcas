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
* Copyright (C) 2013, Ignacio Fdez. Galvan                             *
************************************************************************
*     Find the "alpha" value that produces a given step length in the
*     RS-RFO family of algorithms. This routine assumes that a larger
*     alpha produces a smaller step.
*
*     Use a "regula falsi"-like method, instead of diagonalization
*     as in Eqs. (16) and (20)
*
      SUBROUTINE Find_RFO_Root(x1,y1,x2,y2,x3,y3,Val)
      IMPLICIT NONE
#include "real.fh"
      REAL*8 x1,x2,x3,y1,y2,y3,Val,xx1,yy1,xx2,yy2
      REAL*8 new,nnew,delta,denom,coefA,coefB,coefC,discr,Thr
      PARAMETER ( Thr = 1.0D-16 )
*
*---- First, we must find a value of alpha that gives a small enough step
*
      IF (y2 .GT. Val) THEN
        y2=y3
*
*       In the first iteration, simply try with alpha+1
        IF (x2 .EQ. Zero) THEN
          new=x1+One
          x2=new
*
*       If the step is below the threshold, a bracketing pair has been found.
*       Suggest an intermediate point (with secant method, or midpoint)
        ELSE IF (y3 .LT. Val) THEN
          new=x1+(Val-y1)/(y2-y1)*(x2-x1)
          IF ((new.LE.x1).OR.(new.GE.x2)) new=Half*(x1+x2)
*
*       If it is not there yet, extrapolate using a straight line (times an
*       arbitrary factor of 1.5, to overcome systematic undershooting), and
*       update the other end point
        ELSE
*
*         If y2 is very close to y1, or actually lower, just double the step
*
          IF (y1-y2 .GT. Thr) THEN
            delta=(Val-y2)*(x1-x2)/(y1-y2)
            new=x2+1.5D0*delta
          ELSE
            new=x2+Two*(x2-x1)
          END IF
          x1=x2
          y1=y2
          x2=new
        END IF
*
*---- Once a bracketing pair has been found, we arrive here with two end points
*     (long and short), and a middle point, for which the step length has just
*     been computed
*
      ELSE
*
*       Use the midpoint (bisection) or secant (regula falsi) as safety net
*       Do this in the subinterval where the root should be
        xx1=x1
        yy1=y1
        xx2=x2
        yy2=y2
        IF (y3 .LT. Val) THEN
          xx2=x3
          yy2=y3
        ELSE
          xx1=x3
          yy1=y3
        END IF
        new=xx1+(Val-yy1)/(yy2-yy1)*(xx2-xx1)
        IF ((new.LE.xx1).OR.(new.GE.xx2)) new=Half*(xx1+xx2)
        nnew=new
*
*       Define a quadratic fitting for the 3 points
        denom=(x1-x2)*(x1-x3)*(x2-x3)
        IF (ABS(denom) .GT. Thr) THEN
          coefA=(x3*(y2-y1)+x2*(y1-y3)+x1*(y3-y2))/denom
          coefB=(x3**2*(y1-y2)+x2**2*(y3-y1)+x1**2*(y2-y3))/denom
          coefC=(x2*x3*(x2-x3)*y1+
     &           x3*x1*(x3-x1)*y2+x1*x2*(x1-x2)*y3)/denom
          discr=coefB**2-Four*coefA*(coefC-Val)
        ELSE
          discr=Zero
        END IF
*
*       Find the point between x1 and x2 where the step would be exactly Val
        IF (discr .GT. Zero) THEN
*         Choose root depending on sign of y1-y2
          IF (y1-y2 .GT. Zero) THEN
            nnew=(-coefB-SQRT(discr))/(Two*coefA)
          ELSE IF (y1-y2 .LT. Zero) THEN
            nnew=(-coefB+SQRT(discr))/(Two*coefA)
          END IF
        END IF
*       Last check: make sure the new point is in the correct interval
        IF ((nnew.GT.xx1).AND.(nnew.LT.xx2)) new=nnew
*
*       Update the end points
        x1=xx1
        y1=yy1
        x2=xx2
        y2=yy2
      END IF
*
*     Update the suggested value
      x3=new
*
      RETURN
      END
