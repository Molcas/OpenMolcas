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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

function Cho_LK_MaxVecPerBatch()
! Thomas Bondo Pedersen, May 2013.
!
! Return the max number of vectors to be handled in core in the LK
! algorithm.
!
! Implemented as a function to make all LK flavors (SCF, RASSCF,
! RASSI) behave the same way. The value 25 has been chosen to mimic
! the approximate number of vectors in each reduced set of full CD
! with the one-step algorithm. This should make the LK algorithm
! (nearly) independent of available memory.

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp) :: Cho_LK_MaxVecPerBatch

Cho_LK_MaxVecPerBatch = 25
#ifdef _DEBUGPRINT_
write(u6,'(A,I10)') 'Cho_LK_MaxVecPerBatch=',Cho_LK_MaxVecPerBatch
#endif

end function Cho_LK_MaxVecPerBatch
