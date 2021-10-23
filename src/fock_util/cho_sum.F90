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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine CHO_SUM(rc,nSym,nBas,iUHF,DoExchange,FLT,FSQ)
!****************************************************************
!  Author : F. Aquilante
!
!  Purpose:
!           Accumulates the Coulomb and Exchange contribution
!           to the frozen AO-Fock matrices for alpha and beta
!           spin as defined in the calling routine
!*****************************************************************

use Data_Structures, only: DSBA_Type
implicit real*8(a-h,o-z)
integer rc, nBas(8), nSym, iUHF
type(DSBA_Type) FLT(*), FSQ(*)
logical DoExchange(*)
!*************************************************
!Statement function
iTri(i,j) = max(i,j)*(max(i,j)-3)/2+i+j
!*************************************************

if (iUHF == 1) then
  nDen = 3
else
  nDen = 1
end if

! Accumulate the contributions and Square the final matrix
! FLT is in lower triangular storage
! FSQ is in squared storage
!
! the lower triangular part of FSQ is added to FLT
!
if (nDen == 1) then

  do ISYM=1,NSYM
    NB = NBAS(ISYM)
    if (NB > 0) then
      if (DoExchange(1)) then
        do IB=1,NB
          do JB=IB,NB
            IJB = iTri(JB,IB)
            FLT(1)%SB(ISYM)%A1(IJB) = FLT(1)%SB(ISYM)%A1(IJB)+FSQ(1)%SB(ISYM)%A2(JB,IB)
          end do
        end do
      end if
      call SQUARE(FLT(1)%SB(ISYM)%A1,FSQ(1)%SB(ISYM)%A2,1,NB,NB)
    end if
  end do

else ! nDen=3

  do ISYM=1,NSYM
    NB = NBAS(ISYM)
    if (NB > 0) then
      if (DoExchange(2)) then
        do IB=1,NB
          do JB=IB,NB
            IJB = iTri(JB,IB)
            FLT(1)%SB(ISYM)%A1(IJB) = FLT(1)%SB(ISYM)%A1(IJB)+FSQ(2)%SB(ISYM)%A2(JB,IB)
            FLT(2)%SB(ISYM)%A1(IJB) = FLT(2)%SB(ISYM)%A1(IJB)+FSQ(3)%SB(ISYM)%A2(JB,IB)
          end do
        end do
      end if

      call SQUARE(FLT(1)%SB(ISYM)%A1,FSQ(2)%SB(ISYM)%A2,1,NB,NB)
      call SQUARE(FLT(2)%SB(ISYM)%A1,FSQ(3)%SB(ISYM)%A2,1,NB,NB)

    end if
  end do

end if ! nDen=3

! Print the Fock-matrix
#ifdef _DEBUGPRINT_
write(6,'(6X,A)') 'TEST PRINT FROM CHO_SUM.'
write(6,'(6X,A)') 'FROZEN FOCK MATRIX IN AO BASIS.'

if (nDen > 1) then

  do jDen=1,2
    if (jDen == 1) write(6,'(6X,A)') 'SPIN ALPHA'
    if (jDen == 2) write(6,'(6X,A)') 'SPIN BETA'
    do ISYM=1,NSYM
      NB = NBAS(ISYM)
      if (NB > 0) then
        write(6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
        call TRIPRT(' ',' ',FLT(jDen)%SB(ISYM)%A1,NB)
      end if
    end do
  end do

else ! nDen=1

  do ISYM=1,NSYM
    NB = NBAS(ISYM)
    if (NB > 0) then
      write(6,'(6X,A,I2)') 'SYMMETRY SPECIES:',ISYM
      call TRIPRT(' ',' ',FLT(1)%SB(ISYM)%A1,NB)
    end if
  end do

end if
#endif

rc = 0

return

end subroutine CHO_SUM
