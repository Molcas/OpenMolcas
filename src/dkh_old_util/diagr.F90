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
! Copyright (C) 1986,1995, Bernd Artur Hess                            *
!               2005, Jesper Wisborg Krogh                             *
!***********************************************************************

subroutine Diagr(A,N,EIG,EW,SINV,AUX,AUXI)

implicit real*8(A-H,O-Z)
dimension A(*), AUX(N,N), SINV(N,N), EIG(N,N), EW(*), AUXI(*)
#include "WrkSpc.fh"

if (n == 0) return
#ifdef MOLPRO
call Square(Aux,A,n,n)
#else
call Square(A,Aux,n,1,n)
#endif
call DGEMM_('N','N',n,n,n,1.0d0,Aux,n,Sinv,n,0.0d0,Eig,n)
call dGemm_tri('T','N',n,n,n,1.0d0,SINV,n,EIG,n,0.0d0,AUXI,n)

call dCopy_(n*n,[0.0d0],0,Eig,1)
call dCopy_(n,[1.0d0],0,Eig,n+1)
call dCopy_(N*(N+1)/2,AUXI,1,AUX,1)
!call NIDiag(AUXI,EIG,N,N)
call NIDiag_New(AUXI,EIG,N,N)
call vEig(N,AUXI,EW)
call JacOrd2(EW,Eig,n,n)

return

end subroutine Diagr
