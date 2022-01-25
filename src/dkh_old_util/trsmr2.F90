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

subroutine TrSmr2(A,B,C,N,H,G,W)
! Performs the equivalent operations as two calls to TrSmr
! RESULT IS IN C = G^T * (B^T * A * B) * G
! and where C and A are triangular packed.

implicit real*8(A-H,O-Z)
dimension A(N*(N+2)/2), B(N,N), C(N*(N+1)/2), H(N,N), W(N,N), G(N,N)

#ifdef MOLPRO
call Square(W,A,n,n)
#else
call Square(A,W,n,1,n)
#endif
call DGEMM_('T','N',n,n,n,1.0d0,B,n,W,n,0.0d0,H,n)
call DGEMM_('N','N',n,n,n,1.0d0,H,n,B,n,0.0d0,W,n)
call DGEMM_('T','N',n,n,n,1.0d0,G,n,W,n,0.0d0,H,n)
call dGemm_tri('N','N',n,n,n,1.0d0,H,n,G,n,0.0d0,C,n)

return

end subroutine TrSmr2
