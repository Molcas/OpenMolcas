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

subroutine FTWOI(DLT,DSQ,FLT,nFLT,FSQ,LBUF,X1,X2)

#include "intent.fh"

use motra_global, only: Debug, FnTwoAO, iPrint, LuTwoAO, nBas, nFro, nSym
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nFLT, LBUF
real(kind=wp), intent(in) :: DLT(*), DSQ(*)
real(kind=wp), intent(inout) :: FLT(nFLT), FSQ(*)
real(kind=wp), intent(_OUT_) :: X1(*), X2(*)
integer(kind=iwp) :: IOPT, IRC, ISTLTT, ISYM, KEEP(8), NB, NB1, NB2, NBSX(8), NSYM2
real(kind=wp) :: ExFac
logical(kind=iwp) :: FoundTwoEls, ISQUAR

call f_Inquire(FnTwoAo,FoundTwoEls)
if (.not. FoundTwoEls) then
  write(u6,*) 'FTwoi: OrdInt not found!'
  call Abend()
end if

IOPT = 0
call OPNORD(IRC,IOPT,FNTWOAO,LUTWOAO)
call GETORD(IRC,ISQUAR,NSYM2,NBSX,KEEP)

! Compare content of 1el and 2el integral file

if (NSYM2 /= NSYM) then
  write(u6,*) 'FTwoi: NSYM2.NE.NSYM'
  write(u6,*) 'NSYM2=',NSYM2
  write(u6,*) 'NSYM=',NSYM
  call Abend()
end if
do ISYM=1,NSYM
  NB1 = NBAS(ISYM)
  NB2 = NBSX(ISYM)
  if (NB1 /= NB2) then
    write(u6,*) 'FTwoi: NB1.NE.NB2'
    write(u6,*) 'NB1=',NB1
    write(u6,*) 'NB2=',NB2
    call Abend()
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
ExFac = One
call FockTwo(nSym,nBas,nFro,Keep,DLT,DSQ,FLT,nFLT,FSQ,LBUF,X1,X2,ExFac)
!                                                                      *
!***********************************************************************
!                                                                      *
! Close electron repulsion integral file

call CLSORD(IRC)
!                                                                      *
!***********************************************************************
!                                                                      *
! Print the Fock-matrix

if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(u6,'(6X,A)') 'Fock matrix in AO basis'
  ISTLTT = 1
  do ISYM=1,NSYM
    NB = NBAS(ISYM)
    if (NB > 0) then
      write(u6,'(6X,A,I2)') 'symmetry species:',ISYM
      call TRIPRT(' ',' ',FLT(ISTLTT),NB)
      ISTLTT = ISTLTT+NB*(NB+1)/2
    end if
  end do
end if

return

end subroutine FTWOI
