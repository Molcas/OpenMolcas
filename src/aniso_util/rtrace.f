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
      Subroutine rtrace(N,A,B)
      ! removes the trace of a Real array
      Implicit None
      Integer, parameter         :: wp=SELECTED_REAL_KIND(p=15,r=307)
      Integer, intent(in)        :: N ! size of the array
      Real(kind=wp), intent(in)  :: A(N) ! input
      Real(kind=wp), intent(out) :: B(N) ! output
      ! local variables
      Integer       :: i
      Real(kind=wp) :: AS
      AS=0.0_wp
      Call dcopy_(N,0.0_wp,0, B,1)
      !compute the equal-weighted average
      Do i=1,N
        AS=AS+A(i)/dble(N)
      End Do
      !translate each of the elements
      Do i=1,N
        B(i)=A(i)-AS
      End Do
      Return
      End

