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

subroutine dkh_prop(n,s,t,v,w,X,pXp,clight,dkord,xord,dkparam)
! Apply the arbitrary order DKH transformation to property integral

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
! Input :
!   X     matrix of property operator
!   pXp   matrix representation of <pxXpx>+<pyXpy>+<pzXpz>
!   w     aka pVp
! Output :
!   X     store the transformed property integral
integer(kind=iwp), intent(in) :: n, dkord, xord, dkparam
real(kind=wp), intent(in) :: s(n,n), t(n,n), v(n,n), w(n,n), clight
real(kind=wp), intent(inout) :: X(n,n), pXp(n,n)
integer(kind=iwp) :: nn, vord, nz, m, iTr, iBk, iEL, iES, iOL, iOS, iEp, iE0, iKC, iCo, iSco, iM, iZ, iW
#include "WrkSpc.fh"

! Transform to free-particle FW picture

nn = n*n+4
call getmem('Tr  ','ALLOC','REAL',iTr,nn)
call getmem('Back','ALLOC','REAL',iBk,nn)
call getmem('mEL ','ALLOC','REAL',iEL,nn)
call getmem('mES ','ALLOC','REAL',iES,nn)
call getmem('mOL ','ALLOC','REAL',iOL,nn)
call getmem('mOS ','ALLOC','REAL',iOS,nn)
call getmem('Ep  ','ALLOC','REAL',iEp,n+4)
call getmem('E0  ','ALLOC','REAL',iE0,n+4)
call getmem('KC  ','ALLOC','REAL',iKC,n*3+4)
call XDR_fpFW(n,s,t,v,w,Work(iTr),Work(iBk),Work(iEL),Work(iES),Work(iOL),Work(iOS),Work(iEp),work(iE0),Work(iKC),Work(iKC+n), &
              Work(iKC+2*n),clight)

! Calculate the DKH unitary transformation ( in terms of a series of W matrices 1..xord )

vord = xord*2
m = n*n
nz = m*vord
call getmem('Wsav','ALLOC','REAL',iW,m*xord*2+4)
call getmem('Cof ','ALLOC','REAL',iCo,vord+8)
call dkh_cofu(vord,dkparam,Work(iCo))

call getmem('Cof2','ALLOC','REAL',iSCo,vord+8)
call getmem('Mat ','ALLOC','REAL',iM,m*6+4)
call getmem('Mat2','ALLOC','REAL',iZ,nz*10+4)
call dkh_ham(n,vord,xord,vord,Work(iEL),Work(iES),Work(iOL),Work(iOS),Work(iEp),Work(iE0),Work(iCo),Work(iSco),Work(iM), &
             Work(iM+m),Work(iM+m*2),Work(iM+m*3),Work(iM+m*4),Work(iM+m*5),Work(iZ),Work(iZ+nz),Work(iZ+nz*2),Work(iZ+nz*3), &
             Work(iZ+nz*4),Work(iZ+nz*5),Work(iZ+nz*6),Work(iZ+nz*7),Work(iZ+nz*8),Work(iZ+nz*9),Work(iW))

! Apply W[1..xord] determined transformation to property operator X

! convert X to fpFW picture
call XDR_fpFWprop(n,Work(iTr),X,pXp,Work(iKC),Work(iKC+n),Work(iKC+2*n),Work(iEL),Work(iES),Work(iOL),Work(iOS),Work(iM))
call dkh_xpx(n,vord,xord,vord,Work(iEL),Work(iES),Work(iOL),Work(iOS),Work(iEp),Work(iE0),Work(iCo),Work(iSco),Work(iM), &
             Work(iM+m),Work(iM+m*2),Work(iM+m*3),Work(iM+m*4),Work(iM+m*5),Work(iZ),Work(iZ+nz),Work(iZ+nz*2),Work(iZ+nz*3), &
             Work(iZ+nz*4),Work(iZ+nz*5),Work(iZ+nz*6),Work(iZ+nz*7),Work(iZ+nz*8),Work(iZ+nz*9),Work(iW))
call getmem('Cof2','FREE','REAL',iSCo,vord+8)
call getmem('Mat ','FREE','REAL',iM,m*6+4)
call getmem('Mat2','FREE','REAL',iZ,nz*10+4)

! Back transform to original non-orthogonal basis picture

call dmxma(n,'C','N',Work(iBk),Work(iEL),Work(iES),One)
call dmxma(n,'N','N',Work(iES),Work(iBk),X,One)

! Free temp memories

call getmem('Cof ','FREE','REAL',iCo,vord+8)
call getmem('Wsav','FREE','REAL',iW,m*xord*2+4)
call getmem('Tr  ','FREE','REAL',iTr,nn)
call getmem('Back','FREE','REAL',iBk,nn)
call getmem('mEL ','FREE','REAL',iEL,nn)
call getmem('mES ','FREE','REAL',iES,nn)
call getmem('mOL ','FREE','REAL',iOL,nn)
call getmem('mOS ','FREE','REAL',iOS,nn)
call getmem('Ep  ','FREE','REAL',iEP,n+4)
call getmem('E0  ','FREE','REAL',iE0,n+4)
call getmem('KC  ','FREE','REAL',iKC,n*3+4)

return
! Avoid unused argument warnings
if (.false.) call Unused_integer(dkord)

end subroutine dkh_prop
