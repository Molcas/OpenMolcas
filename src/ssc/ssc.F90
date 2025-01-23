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
! Copyright (C) 1989-1992, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************

subroutine ssc(iRC)
use definitions, only: wp, iwp, u6
use spool, only: SpoolInp, Close_LuSpool
use Breit, only: D_tensor
use Constants, only: Zero, Two, Three, gelectron, c_in_au, Four
Implicit None
integer(kind=iwp), Intent(out) :: iRC

integer(kind=iwp) :: LuSpool, nDiff, i, j
integer(kind=iwp), parameter :: nH=3
logical(kind=iwp) :: DoRys
integer(kind=iwp), external :: IsFreeUnit
character(len=8) :: Method
real(kind=wp) R2, EVec(nH,nH), EVal(nH*(nH+1)/2), D, E


LuSpool = 37
LuSpool=IsFreeUnit(LuSpool)
call SpoolInp(LuSpool)

nDiff = 2
DoRys = .True.
call IniSew(DoRys,nDiff)

call Close_LuSpool(LuSpool)

D_tensor(:,:) = Zero
Call Drv2El_BP()

Call ClsSew()
D_tensor(:,:)=-(gelectron**2 * Three)/(Four * c_in_au**2) * D_tensor(:,:)

call Get_cArray('Relax Method',Method,8)
Write (u6,*)
Write (u6,'(2A)') 'Method:', Method
Write (u6,*)

call RecPrt('The D tensor',' ',D_tensor,3,3)
Write (u6,*)
Write (u6,*)

R2=D_tensor(1,1)+D_tensor(2,2)+D_tensor(3,3)
D_tensor(1,1)=D_tensor(1,1)-(R2/Three)
D_tensor(2,2)=D_tensor(2,2)-(R2/Three)
D_tensor(3,3)=D_tensor(3,3)-(R2/Three)

call RecPrt('The D tensor, in the traceless form ',' ',D_tensor,3,3)
Write (u6,*)
Write (u6,*)

Do i = 1, nH
   Do j = 1, i
      EVal(i*(i-1)/2+j)=D_tensor(i,j)
   End Do
End Do
Call TriPrt('The D tensor in triangular form',' ',EVal,nH)

call unitmat(EVec,nH)

! Compute eigenvalues and eigenvectors

call Jacob(EVal,EVec,nH,nH)
call Jacord(EVal,EVec,nH,nH)

Call TriPrt('The diagonal D tensor',' ',EVal,nH)
Write (u6,*)
Write (u6,*)
Call RecPrt('The eigenvectors of the D tensor',' ',EVec,nH,nH)
Write (u6,*)
Write (u6,*)

D=(Three/Two)*EVal(1)
E=(EVal(6)-EVal(3))/Two
Write (u6,'(2X,A,E15.5,A)') 'D=(3/2)*D_zz      =',D,'cm-1'
Write (u6,'(2X,A,E15.5,A)') 'E=(1/2)*(D_xx-Dyy)=',E,'cm-1'
Write (u6,'(2X,A,E15.5  )') 'E/D               =',E/D
Write (u6,*)
Write (u6,*)

! add info to the check file.
Call Add_Info('D Tensor',D_tensor,6,8)

iRC=0

end subroutine ssc
