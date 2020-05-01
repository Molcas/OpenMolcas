************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Real*8 function dev(N, Fcal, Fexp)
c this function Returns the standard deviation between experimental
c data points and the computed ones
c     N --- number of data points ( Integer, input)
c  Fcal --- calculated array of size (N), Real(kind=8) ::, input
c  Fexp --- experimental array of size (N), Real(kind=8) ::, input
c   dev --- standard deviation, Real(kind=8) ::, output;
      Implicit None
      Integer, parameter        :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)       :: N
      Real(kind=8), intent(in) :: Fcal(N), Fexp(N)
      Integer                   :: i
      Real(kind=8)             :: diff, X
      dev=0.0_wp
      X=0.0_wp
      Do i=1,N
        diff=0.0_wp
        diff=Fcal(i)-Fexp(i)
           X=X+diff*diff/dble(N)
      End Do
      dev=sqrt(X)
      Return
      End
