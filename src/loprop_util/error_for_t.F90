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

real*8 function Error_for_t(t,rMP,xrMP,xxrMP,xnrMP,EC,A,R_ij,C_o_C,ij,l,nij,lMax,nElem,nAtoms,nPert,Scratch_New,Scratch_Org, &
                            iPrint_Errors)

implicit real*8(A-H,O-Z)
#include "real.fh"
dimension EC(3,nij), A(3,nij), C_o_C(3), R_ij(3)
dimension rMP(nij,0:nElem-1,0:nPert-1), xnrMP(nij,nElem)
dimension xrMP(nij,nElem), xxrMP(nij,nElem)
dimension Scratch_Org(nij*(2+lMax+1)), Scratch_New(nij*(2+lMax+1))

iDim = nij*nElem
A(1,ij) = EC(1,ij)+t*R_ij(1)
A(2,ij) = EC(2,ij)+t*R_ij(2)
A(3,ij) = EC(3,ij)+t*R_ij(3)
call dCopy_(iDim,rMP,1,xrMP,1)
do ij_temp=1,nij
  call ReExpand(xrMP,nij,nElem,EC(1,ij_temp),C_o_C,ij_temp,lMax)
end do
call daxpy_(iDim,One,xnrMP,1,xrMP,1)
call dCopy_(iDim,xrMP,1,xxrMP,1)
call CutOff_Error(l,lMax,xrMP,xxrMP,nij,A,C_o_C,nElem,Scratch_New,Scratch_Org,nAtoms,iPrint_Errors,Error)

Error_for_t = Error

return

end function Error_for_t
