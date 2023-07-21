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

use ChoArr, only: iSP2F, nBstSh, iShP2RS, iShP2Q, nDim_Batch

implicit none
integer iSCD, iSAB
integer lInt
real*8 xInt(lInt)
#include "cholesky.fh"
character*10 SecNam
parameter(SecNam='Cho_PrtInt')
integer nRow(8)
integer iSC, iSD, iSA, iSB
integer nCD, nAB
integer AB, CD
integer iAB, iCD
integer iSym
integer kOffI
real*8 xNorm

! Set row dimension
if (IfcSew == 2) then
  do iSym=1,nSym
    nRow(iSym) = nnBstR(iSym,2)
  end do
else if (IfcSew == 3) then
  do iSym=1,nSym
    nRow(iSym) = nDim_Batch(iSym)
  end do
else
  call Cho_Quit(SecNam//': Illegal IfcSew',103)
  do iSym=1,nSym ! avoid compiler warnings
    nRow(iSym) = 0
  end do
end if

! Get full shell pair dimensions
call Cho_InvPck(iSP2F(iSCD),iSC,iSD,.true.)
if (iSC == iSD) then
  nCD = nBstSh(iSC)*(nBstSh(iSC)+1)/2
else
  nCD = nBstSh(iSC)*nBstSh(iSD)
end if
call Cho_InvPck(iSP2F(iSAB),iSA,iSB,.true.)
if (iSA == iSB) then
  nAB = nBstSh(iSA)*(nBstSh(iSA)+1)/2
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
    xNorm = 0.0d0
    do CD=1,nCD
      iCD = iShP2RS(1,CD)
      if (iCD > 0) then
        if (iShP2RS(2,CD) == iSym) then
          write(Lupri,'(2X,A,I4,A,I4,A,1P,D15.6)') '(',CD,'|',AB,') =',xInt(kOffI+iCD)
          xNorm = xNorm+xInt(kOffI+iCD)**2
        end if
      end if
    end do
    write(Lupri,'(A,I4,A,1P,D15.6)') '**Norm of column',AB,':',sqrt(xNorm)
  end if
end do

end subroutine Cho_PrtInt
