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
! Copyright (C) Jeppe Olsen                                            *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine GRAPW(W,Y,MINEL,MAXEL,NORB,NEL)
! A graph of strings has been defined from
!
!  MINEL(I) is the smallest allowed number of electrons in
!  orbitals 1 through I
!
!  MAXEL(I) is the largest allowed number of electrons in
!  orbitals 1 through I
!
! Set up vertex weights W
! Set up arc weights    Y
!
! Reverse lexical ordering is used with
! weights of unoccupied orbitals set to 0
!
! Jeppe Olsen

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: NORB, MINEL(NORB), MAXEL(NORB), NEL
integer(kind=iwp), intent(out) :: W(NORB+1,NEL+1), Y(NORB,NEL)
integer(kind=iwp) :: IEL, IORB

W(:,:) = 0
Y(:,:) = 0

!================
!  Vertex weights
!================

! (Weight for vertex(IEL,IORB) is stored in W(IORB+1,IEL+1))
W(1,1) = 1
do IEL=0,NEL
  do IORB=1,NORB
    if ((MINEL(IORB) <= IEL) .and. (IEL <= MAXEL(IORB))) then
      if (IEL > 0) then
        W(IORB+1,IEL+1) = W(IORB,IEL+1)+W(IORB,IEL)
      else
        W(IORB+1,1) = W(IORB,1)
      end if
    end if
  end do
end do

! Weight for arc connecting vertices (IORB-1,IEL-1) and(IORB,IEL)
! is stored in Y(IORB,IEL)
! Y(IORB,IEL) = W(IORB-1,IEL)

do IEL=1,NEL
  do IORB=1,NORB
    if ((MINEL(IORB) <= IEL) .and. (IEL <= MAXEL(IORB))) Y(IORB,IEL) = W(IORB-1+1,IEL+1)
  end do
end do

#ifdef _DEBUGPRINT_
write(u6,*) ' vertex weights'
call IWRTMA(W,NORB+1,NEL+1,NORB+1,NEL+1)
write(u6,*) ' arc weights'
call IWRTMA(Y,NORB,NEL,NORB,NEL)
#endif

end subroutine GRAPW
