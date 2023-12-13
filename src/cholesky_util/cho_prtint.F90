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
! Copyright (C) 2010, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_PrtInt(iSCD,iSAB,xInt,lInt)
!
! Thomas Bondo Pedersen, April 2010.
!
! Purpose: Print integral shell quadruple (IfcSew=2 or 3).

use Index_Functions, only: nTri_Elem
use Cholesky, only: IFCSEW, iOff_col, iShP2Q, iShP2RS, iSP2F, LuPri, nBstSh, nDim_Batch, nnBstR, nSym
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iSCD, iSAB, lInt
real(kind=wp), intent(in) :: xInt(lInt)
integer(kind=iwp) :: AB, CD, iAB, iCD, iSA, iSB, iSC, iSD, iSym, kOffI, nAB, nCD, nRow(8)
real(kind=wp) :: xNorm
character(len=*), parameter :: SecNam = 'Cho_PrtInt'

! Set row dimension
if (IfcSew == 2) then
  nRow(1:nSym) = nnBstR(1:nSym,2)
else if (IfcSew == 3) then
  nRow(1:nSym) = nDim_Batch(1:nSym)
else
  call Cho_Quit(SecNam//': Illegal IfcSew',103)
  nRow(1:nSym) = 0 ! avoid compiler warnings
end if

! Get full shell pair dimensions
call Cho_InvPck(iSP2F(iSCD),iSC,iSD,.true.)
if (iSC == iSD) then
  nCD = nTri_Elem(nBstSh(iSC))
else
  nCD = nBstSh(iSC)*nBstSh(iSD)
end if
call Cho_InvPck(iSP2F(iSAB),iSA,iSB,.true.)
if (iSA == iSB) then
  nAB = nTri_Elem(nBstSh(iSA))
else
  nAB = nBstSh(iSA)*nBstSh(iSB)
end if

! Loop through integral shell quadruple
write(LuPri,'(//,A,I4,A,I4,A,I4,A,I4,A)') 'Shell Quadruple (',iSC,',',iSD,'|',iSA,',',iSB,'):'
do AB=1,nAB
  iAB = iShP2Q(1,AB)
  if (iAB > 0) then
    iSym = iShP2Q(2,AB)
    kOffI = iOff_Col(iSym)+nRow(iSym)*(iAB-1)
    xNorm = Zero
    do CD=1,nCD
      iCD = iShP2RS(1,CD)
      if (iCD > 0) then
        if (iShP2RS(2,CD) == iSym) then
          write(Lupri,'(2X,A,I4,A,I4,A,ES15.6)') '(',CD,'|',AB,') =',xInt(kOffI+iCD)
          xNorm = xNorm+xInt(kOffI+iCD)**2
        end if
      end if
    end do
    write(Lupri,'(A,I4,A,ES15.6)') '**Norm of column',AB,':',sqrt(xNorm)
  end if
end do

end subroutine Cho_PrtInt
