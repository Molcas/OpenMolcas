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

implicit real*8(A-H,O-Z)
#include "motra_global.fh"
#include "files_motra.fh"
dimension FSQ(*), FLT(nFLT), DSQ(*), DLT(*), X1(*), X2(*)
dimension NBSX(8), KEEP(8)
logical FoundTwoEls, ISQUAR

call f_Inquire(FnTwoAo,FoundTwoEls)
if (.not. FoundTwoEls) then
  write(6,*) 'FTwoi: OrdInt not found!'
  call Abend()
end if

call OPNORD(IRC,0,FNTWOAO,LUTWOAO)
call GETORD(IRC,ISQUAR,NSYM2,NBSX,KEEP)

! Compare content of 1el and 2el integral file

if (NSYM2 /= NSYM) then
  write(6,*) 'FTwoi: NSYM2.NE.NSYM'
  write(6,*) 'NSYM2=',NSYM2
  write(6,*) 'NSYM=',NSYM
  call Abend()
end if
do ISYM=1,NSYM
  NB1 = NBAS(ISYM)
  NB2 = NBSX(ISYM)
  if (NB1 /= NB2) then
    write(6,*) 'FTwoi: NB1.NE.NB2'
    write(6,*) 'NB1=',NB1
    write(6,*) 'NB2=',NB2
    call Abend()
  end if
end do
!                                                                      *
!***********************************************************************
!                                                                      *
ExFac = 1.0d0
call FockTwo(nSym,nBas,nFro,Keep,DLT,DSQ,FLT,nFLT,FSQ,LBUF,X1,X2,ExFac)
!                                                                      *
!***********************************************************************
!                                                                      *
! Close electron repulsion integral file

call CLSORD(IRC,0)
!                                                                      *
!***********************************************************************
!                                                                      *
! Print the Fock-matrix

if ((IPRINT >= 5) .or. (DEBUG /= 0)) then
  write(6,'(6X,A)') 'Fock matrix in AO basis'
  ISTLTT = 1
  do ISYM=1,NSYM
    NB = NBAS(ISYM)
    if (NB > 0) then
      write(6,'(6X,A,I2)') 'symmetry species:',ISYM
      call TRIPRT(' ',' ',FLT(ISTLTT),NB)
      ISTLTT = ISTLTT+NB*(NB+1)/2
    end if
  end do
end if

return

end subroutine FTWOI
