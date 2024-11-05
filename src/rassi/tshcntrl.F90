!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
Module tshcntrl
! This file preserves the values of the variables for the
! trajectory surface hopping algorithm.
INTEGER  ISTATE1,ISTATE2,NCI1,NCI2,nHop
LOGICAL  ChkHop,lHop
! ISTATE1 - Number of the current relaxed state
! ISTATE2 - Number of the state that interacts with the current state
! NCI1    - Configuration's number of state 1
! ChkHop  - If .TRUE. switches on the TSH algorithm
! lHop    - Is .TRUE. if nHop is set in the RunFile
! nHop    - Number of transitions (Hops) already occured
End Module tshcntrl
