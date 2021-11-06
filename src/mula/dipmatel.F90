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

subroutine DipMatEl(Dij,W,L,U,FC00,nMat,nInc,nDec,D0,D1,D2,D3,D4,max_term,base,m_ord,nosc,max_mOrd,max_nOrd2)
!  Purpose:
!    Calculate matrix elements of the transition dipole.
!
!  Input:
!    D0       : Real*8 variable - the zero order term of the
!               transition dipole.
!    D1       : Real*8 array - the first order term of the
!               transition dipole.
!    D2       : Real*8 two dimensional array - the second order
!               term of the transition dipole.
!    D3       : Real*8 three dimensional array - the third order
!               term of the transition dipole.
!    D4       : Real*8 four dimensional array - the fourth order
!               term of the transition dipole.
!    W        : Real*8 two dimensional array
!    L,U      : Real*8 two dimensional array
!    FC00     : Real*8 variable
!    nMat     : Two dimensional integer array.
!    max_term : Integer - maximum order of the transition dipole terms.
!
!  Output:
!    Dij      : Real*8 two dimensional array - contains the
!               matrix elements of the transition dipole.
!
!  Uses:
!    MatElMod

!use PotKin
implicit real*8(a-h,o-z)
#include "dims.fh"
real*8 Base(nOsc,nOsc)
real*8 Dij(0:max_mOrd,0:max_mOrd)
real*8 W(nOsc,nOsc)
real*8 L(0:max_mOrd,0:max_mOrd)
real*8 U(0:max_nOrd2,0:max_nord2)
integer nMat(0:ndim1,ndim2), ndec(0:ndim1,ndim2), ninc(0:ndim1,ndim2)
real*8 D1(nosc)
real*8 D2(nosc,nosc)
real*8 D3(nosc,nosc,nosc)
real*8 D4(nosc,nosc,nosc,nosc)
#include "WrkSpc.fh"

! Initialize.
max_nOrd = max_mOrd
mPlus = max_mOrd+1
nPlus = max_nOrd+1
nOscOld = nOsc
l_A = (max_mOrd+1)*(max_nOrd+1)
call GetMem('A','Allo','Real',ipA,l_A)
call GetMem('Wtemp','Allo','Real',ipWtemp,nOscOld*nOsc)
call DGEMM_('N','N',nOscOld,nOsc,nOsc,1.0d0,Base,nOscOld,W,nOsc,0.0d0,Work(ipWtemp),nOscOld)
call dcopy_(l_A,[0.0d0],0,Work(ipA),1)
call PotEnergy(Work(ipA),nMat,nInc,nDec,D0,D1,D2,D3,D4,max_term,Work(ipWTemp),m_ord,nosc,nOscOld)

call GetMem('Wtemp','Free','Real',ipWtemp,nOscOld*nOsc)
call GetMem('Temp','Allo','Real',ipTemp,l_A)
call DGEMM_('N','T',mplus,mplus,nplus,1.0d0,Work(ipA),mplus,U,mplus,0.0d0,Work(ipTemp),mplus)
call DGEMM_('N','N',mPlus,nPlus,mPlus,1.0d0,L,mPlus,Work(ipTemp+max_mOrd+1+(max_mOrd+1)*max_nOrd2-1),mPlus,0.0d0,Dij,mPlus)
call GetMem('Temp','Free','Real',ipTemp,l_A)
call GetMem('A','Free','Real',ipA,l_A)

!Dij = FC00*Dij
call dscal_((max_mOrd+1)*(max_mOrd+1),FC00,Dij,1)

end subroutine DipMatEl
