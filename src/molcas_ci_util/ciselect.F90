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
! Copyright (C) 1996, Markus P. Fuelscher                              *
!***********************************************************************

subroutine CiSelect(S1,S2)
!***********************************************************************
!                                                                      *
!     purpose:                                                         *
!     Select CI_vector which matches best the test vectors             *
!                                                                      *
!     calling arguments:                                               *
!     S1      : array of real*8, input/output                          *
!               overlap matrix with test vectors                       *
!     S2      : array of real*8, input/output                          *
!               norm of the test configurations in the CI vector       *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     written by:                                                      *
!     M.P. Fuelscher                                                   *
!     University of Lund, Sweden, 1996                                 *
!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
!     history: none                                                    *
!                                                                      *
!***********************************************************************

implicit real*8(A-H,O-Z)
#include "rasdim.fh"
#include "general.fh"
#include "rasscf.fh"
dimension S1(lRoots,lRoots)
dimension S2(lRoots,lRoots)
dimension iTemp(mxRoot)

if (ITER == 1) return

! make a local copy of the present selection vector of weights
call iCopy(mxRoot,[0],1,iTemp,1)

! make a new choice using the overlap ov the etst vector with
! the CI vector in the subspace of the test vector
do kRoot=1,nRoots
  ! search for the element of largest overlap
  maxS1 = 1
  S1max = S1(1,kRoot)
  do jRoot=1,lRoots
    if (S1(jRoot,kRoot) > S1max) then
      maxS1 = jRoot
      S1max = S1(jRoot,kRoot)
    end if
  end do
  iTemp(kRoot) = maxS1
  ! cleanup,to ensure that we don't pick the same root twice
  do jRoot=1,nRoots
    S1(maxS1,jRoot) = S1(maxS1,jRoot)-999999.0d0
  end do
end do

! Restore S1
do kRoot=1,nRoots
  do jRoot=1,nRoots
    S1(iTemp(kRoot),jRoot) = S1(iTemp(kRoot),jRoot)+999999.0d0
  end do
end do

! print
if (nRoots == 1) then
  write(6,'(6X,A,T45,10I6)') 'new root selected:',(iTemp(kRoot),kRoot=1,nRoots)
  write(6,'(6X,A,T45,10F6.3)') 'overlap           ',(S1(iTemp(kRoot),kRoot),kRoot=1,nRoots)
else
  write(6,'(6X,A,T45,10I6)') 'new roots selected:',(iTemp(kRoot),kRoot=1,nRoots)
  write(6,'(6X,A,T45,10F6.3)') 'overlap           ',(S1(iTemp(kRoot),kRoot),kRoot=1,nRoots)
end if

! Compare the overlap elemets <1|2> and <1|1> where
! |1> is the CI vector in the subspace of the test vector and
! |2> is the test vector.
! The program breaks execution if the following roules apply:
! If <1|2> is less than 0.5*<1|1>
! If the weight of <1|2> or <1|1> is less than 0.3
istop = 0
do kRoot=1,nRoots
  S1jk = S1(iTemp(kRoot),kRoot)
  S2jk = S2(iTemp(kRoot),kRoot)
  if (S1jk < 0.5d0*S2jk) istop = ibset(istop,0)
  if (sqrt(S1jk) < 0.316d0) istop = ibset(istop,1)
  if (sqrt(S2jk) < 0.3d0) istop = ibset(istop,2)
end do

! If the stop flag has been set write an approriate message
! and to stop execution change the iteration counter.
if (istop >= 1) then
  write(6,*)
  write(6,'(6X,120A1)') ('=',i=1,120)
  if (btest(istop,0)) then
    write(6,'(6X,A)') 'The projection of the CI vector(s) onto the model vector(s)'
    write(6,'(6X,A)') 'is smaller than half the norm of the subspace.'
  end if
  if (btest(istop,1)) then
    write(6,'(6X,A)') 'The overlap of the projected CI vector(s) and the model vector(s) is smaller than 0.1'
  end if
  if (btest(istop,2)) then
    write(6,'(6X,A)') 'The weight(s) of the subspace is(are) smaller than 30% of the total wave function(s)'
  end if
  write(6,'(6X,A)') 'Please, check your model space'
  write(6,'(6X,A)') 'The program stops after the next iteration'
  write(6,'(6X,120A1)') ('=',i=1,120)
  write(6,*)
  write(6,*)
  MAXIT = ITER
else
  call iCopy(nRoots,iTemp,1,iRoot,1)
end if

return

end subroutine CiSelect
