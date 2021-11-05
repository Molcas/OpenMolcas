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
!!-----------------------------------------------------------------------!
!!
      Subroutine DipMatEl(Dij,W,L,U,FC00,nMat,nInc,nDec,D0,D1,          &
     &       D2,D3,D4,max_term,base,m_ord,nosc,max_mOrd,max_nOrd2)
!!
!!  Purpose:
!!    Calculate matrix elements of the transition dipole.
!!
!!  Input:
!!    D0       : Real*8 variable - the zero order term of the
!!               transition dipole.
!!    D1       : Real*8 array - the first order term of the
!!               transition dipole.
!!    D2       : Real*8 two dimensional array - the second order
!!               term of the transition dipole.
!!    D3       : Real*8 three dimensional array - the third order
!!               term of the transition dipole.
!!    D4       : Real*8 four dimensional array - the fourth order
!!               term of the transition dipole.
!!    W        : Real*8 two dimensional array
!!    L,U      : Real*8 two dimensional array
!!    FC00     : Real*8 variable
!!    nMat     : Two dimensional integer array.
!!    max_term : Integer - maximum order of the transition dipole terms.
!!
!!  Output:
!!    Dij      : Real*8 two dimensional array - contains the
!!               matrix elements of the transition dipole.
!!
!!  Uses:
!!    MatElMod
!!
!       Use PotKin
      Implicit Real*8 ( a-h,o-z )
#include "dims.fh"
      Real*8 Base(nOsc,nOsc)
      Real*8 Dij( 0:max_mOrd,0:max_mOrd )
      Real*8 W(nOsc,nOsc)
      Real*8 L( 0:max_mOrd,0:max_mOrd )
      Real*8 U( 0:max_nOrd2,0:max_nord2 )
      Integer nMat( 0:ndim1,ndim2 ),ndec( 0:ndim1,ndim2 ),              &
     &         ninc( 0:ndim1,ndim2 )
      Real*8 D1( nosc )
      Real*8 D2( nosc,nosc )
      Real*8 D3( nosc,nosc,nosc )
      Real*8 D4( nosc,nosc,nosc,nosc )
#include "WrkSpc.fh"
!!
!!---- Initialize.
      max_nOrd = max_mOrd
      mPlus = max_mOrd+1
      nPlus  = max_nOrd +1
      nOscOld = nOsc
      l_A=(max_mOrd+1)*(max_nOrd+1)
      Call GetMem('A','Allo','Real',ipA,l_A)
      Call GetMem('Wtemp','Allo','Real',ipWtemp,nOscOld*nOsc)
      Call DGEMM_('N','N',                                              &
     &            nOscOld,nOsc,nOsc,                                    &
     &            1.0d0,Base,nOscOld,                                   &
     &            W,nOsc,                                               &
     &            0.0d0,Work(ipWtemp),nOscOld)
      call dcopy_(l_A,[0.0d0],0,Work(ipA),1)
      Call PotEnergy(                                                   &
     &      Work(ipA),nMat,nInc,nDec,D0,D1,D2,D3,D4,max_term,           &
     &      Work(ipWTemp),m_ord,nosc,nOscOld)
!!
      Call GetMem('Wtemp','Free','Real',ipWtemp,nOscOld*nOsc)
      Call GetMem('Temp','Allo','Real',ipTemp,l_A)
      Call DGEMM_('N','T',                                              &
     &            mplus,mplus,nplus,                                    &
     &            1.0d0,Work(ipA),mplus,                                &
     &            U,mplus,                                              &
     &            0.0d0,Work(ipTemp),mplus)
      Call DGEMM_('N','N',                                              &
     &            mPlus,nPlus,mPlus,                                    &
     &            1.0d0,L,mPlus,                                        &
     &            Work(ipTemp+max_mOrd+1+(max_mOrd+1)*max_nOrd2-1),     &
     &            mPlus,                                                &
     &            0.0d0,Dij,mPlus)
      Call GetMem('Temp','Free','Real',ipTemp,l_A)
      Call GetMem('A','Free','Real',ipA,l_A)

!       Dij = FC00*Dij
      call dscal_((max_mOrd+1)*(max_mOrd+1),FC00,Dij,1)
!!
      End
