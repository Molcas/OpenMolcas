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
! Copyright (C) 1996, Niclas Forsberg                                  *
!***********************************************************************

subroutine AnharmonicFreq(x_anharm,harmfreq,anharmfreq,nOsc)
!  Purpose:
!    Calculate the anharmonic frequencies.
!
!  Input:
!    x_anharm   : Real*8 two dimensional array - anharmonicity
!                 constants.
!    harmfreq   : Real*8 array - harmonic frequencies.
!
!  Output:
!    anharmfreq : Real*8 array - anharmonic frequencies.
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1996.

implicit real*8(a-h,o-z)
#include "Constants_mula.fh"
real*8 G1, G2
real*8 x_anharm(nosc,nosc)
real*8 harmfreq(nosc)
real*8 anharmfreq(nosc)
#include "WrkSpc.fh"

!nDim = nOsc
call GetMem('level1','Allo','Inte',iplevel1,nOsc)
call GetMem('level2','Allo','Inte',iplevel2,nOsc)

G1 = 0.0d0
G2 = G1
do iv=0,nOsc-1
  iWork(iplevel1+iv) = 0
end do
!level1 = 0
do istate=1,nOsc
  do iv=0,nOsc-1
    iWork(iplevel2+iv) = 0
  end do
  !level2 = 0
  iWork(iplevel2+istate-1) = 1
  !l_harm = nOsc
  call TransEnergy(G1,x_anharm,harmfreq,iWork(iplevel1),G2,x_anharm,harmfreq,iWork(iplevel2),anharmfreq(istate),nOsc)
end do

call GetMem('level1','Free','Inte',iplevel1,nOsc)
call GetMem('level2','Free','Inte',iplevel2,nOsc)

end subroutine AnharmonicFreq
