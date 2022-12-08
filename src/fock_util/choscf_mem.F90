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

subroutine CHOSCF_MEM(nSym,nBas,nD,DoExchange,pNocc,ALGO,REORD,MinMem,lOff1)
!****************************************************************
!  Author : F. Aquilante
!
!  Purpose:
!            Returns an array MinMem that contains the minimum
!            amount of memory required for each Cholesky
!            vector for the building of the frozen AO-Fock matrix
!
!            The first vector is stored in a larger memory block
!            in some cases in order to re-use the allocated memory
!            for computing the exchange intermediate.
!            The value of LOFF1 specifies the amount of memory
!            reserved for the 1st vector
!
!            Note that this is a specialized code and can be used
!                 only in connection with the truth table defined
!                 in the calling routine
!*****************************************************************

use ChoArr, only: nDimRS
use Symmetry_Info, only: Mul
use Data_Structures, only: Integer_Pointer
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: nSym, nBas(nSym), nD, ALGO
logical(kind=iwp), intent(in) :: DoExchange(*), REORD
type(Integer_Pointer), intent(in) :: pNocc(*)
integer(kind=iwp), intent(out) :: MinMem(nSym), lOff1
integer(kind=iwp) :: i, iSym, iSymp, j, jSym, ksym, Nab, nDen, Nmax, NSab, Mabmx(nSym), Moccmx(nSym), MxBas(nSym)
logical(kind=iwp) :: xToDo

!*************************************************
nDen = nD*(nD+1)/2
xToDo = DoExChange(nD)

!============================
lOff1 = 0
do i=1,nDen
  do j=1,nSym
    lOff1 = max(lOff1,pNocc(i)%I1(j))
  end do
end do

do j=1,nSym
  Moccmx(j) = 0
  do i=1,nDen
    Moccmx(j) = max(Moccmx(j),pNocc(i)%I1(j))
  end do
end do

! Max dimension of a read symmetry block
Nmax = 0
do iSym=1,nSym
  if ((NBAS(iSym) > Nmax) .and. (Moccmx(isym) /= 0)) then
    Nmax = NBAS(iSym)
  end if
end do

lOff1 = Nmax*lOff1

!============================
do jSym=1,nSym

  ! Total length of a vector
  Mabmx(jSym) = 0
  MxBas(jSym) = 0
  Nab = 0
  NSab = 0
  do ksym=1,nSym
    iSymp = Mul(ksym,jSym)
    if ((iSymp > ksym) .and. ((Moccmx(iSymp) /= 0) .or. (Moccmx(ksym) /= 0))) then
      Nab = Nab+nBas(ksym)*nBas(iSymp)
      Mabmx(jSym) = max(Mabmx(jSym),max(nBas(ksym)*Moccmx(iSymp),nBas(iSymp)*Moccmx(ksym)))
      MxBas(jSym) = max(MxBas(jSym),nBas(ksym)*nBas(iSymp))
    else
      if (iSymp == ksym) then ! all are needed for the Coulomb
        Nab = Nab+nBas(ksym)*(nBas(ksym)+1)/2
        NSab = NSab+nBas(ksym)**2
      end if
    end if
  end do

  !============================================================
  ! The following memory management is bound to the truth table
  ! assigned in the calling routine
  !============================================================

  if (.not. xToDo) then

    MinMem(jSym) = Nab+1 ! to store 1 vector of L(rJ),s + V(J)
    if (.not. REORD) then
      MinMem(jSym) = Nab+nDimRS(jSym,1)
      ! L + read 1 vect reduced set1
    end if

  else ! Memory for "off-diagonal" exchange only (jSym /= 1)

    if (REORD) then
      MinMem(jSym) = 2*Nab ! 1 vector of L(rJ),s + W(rJ),s
    else
      if (ALGO == 2) then
        MinMem(jSym) = Nab+max(nDimRS(jSym,1),Mabmx(jSym))
      else
        MinMem(jSym) = Nab+max(nDimRS(jSym,1),MxBas(jSym))
      end if
    end if

  end if

  !-----
  if (xToDo .and. (jSym == 1)) then ! jSym=1 is a special case

    if (nSym == 1) then

      if (ALGO == 2) then
        if (lOff1 > Nab) then
          MinMem(jSym) = lOff1+NSab
        else
          MinMem(jSym) = Nab+NSab
          lOff1 = Nab
        end if
      else

        MinMem(jSym) = 2*NSab
        lOff1 = NSab

      end if

    else  ! nSym > 1   jSym=1

      if (ALGO == 2) then
        if (lOff1 > nBas(1)*(nBas(1)+1)/2) then
          MinMem(jSym) = Nab-nBas(1)*(nBas(1)+1)/2+lOff1+Nmax**2
        else
          MinMem(jSym) = Nab+Nmax**2
          lOff1 = nBas(1)*(nBas(1)+1)/2
        end if
      else

        MinMem(jSym) = 2*(Nmax**2)+(Nab-nBas(1)*(nBas(1)+1)/2)
        lOff1 = Nmax**2

      end if

    end if

  end if ! JSYM=1  special case

! ====================================================================

end do ! loop over jSym

return

end subroutine CHOSCF_MEM
