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
       subroutine unpckhelp3 (a,b,dimp,dimq,dime,dimf,eadd,noe,fadd,nof)
!
!     this routine do:
!     b(e,f) = a(pe,qf)
!
       integer dimp,dimq,dime,dimf,eadd,noe,fadd,nof
       real*8 a(1:dimp,1:dimq)
       real*8 b(1:dime,1:dimf)
!
!     help variables
       integer pe,qf,f
!
       do 100 qf=fadd+1,fadd+nof
       f=qf-fadd
       do 101 pe=eadd+1,eadd+noe
       b(pe-eadd,f)=a(pe,qf)
 101    continue
 100    continue
!
       return
       end
