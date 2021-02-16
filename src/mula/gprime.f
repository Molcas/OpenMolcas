************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1996, Niclas Forsberg                                  *
************************************************************************
C!-----------------------------------------------------------------------!
C!
      Subroutine CalcGprime(Gprime,Mass,xvec,InterVec,
     &       AtCoord,NumOfAt,h,NumInt)
C!
C!  Purpose:
C!    Calculate first derivatives of G.
C!
C!  Input:
C!    Mass     : Real*8 array - the mass of the atoms.
C!    xvec     : Real*8 array - the geometry in internal
C!               coordinates.
C!    InterVec : Integer array.
C!    NumOfAt  : Integer - the number of atoms.
C!
C!  Output:
C!    Gprime   : Real*8 two dimensional array - first
C!               derivative of the inverse mass tensor G.
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1996.
C!
c       Implicit None
#include "Constants_mula.fh"
#include "dims.fh"
      Real*8  h
      Real*8 Gprime(ngdim,ngdim,ngdim)
      Integer   NumInt,NumOfAt
      Real*8 Mass  (NumOfAt)
      Real*8 xvec  (NumInt)
      Real*8 xtmp  (NumInt)
      Integer InterVec(*)
      Real*8 AtCoord (3,NumOfAt)
      Integer   icoord,iterm
      Integer   ih,k,j
#include "WrkSpc.fh"
C!
C!---- Initialize.
      Call GetMem('Stemp','Allo','Real',ipStemp,3*NumOfAt*NumInt)
      Call GetMem('Gtemp','Allo','Real',ipGtemp,NumInt*NumInt*4)
C!
      Do icoord = 1,NumInt
      do iv=1,NumInt
      xtmp(iv) = xvec(iv)
      enddo
c              Gtemp = 0.0d0
      call dcopy_(NumInt*NumInt*4,0.0d0,0,Work(ipGtemp),1)
      iterm = 1
      Do ih = -3,3,2
      xtmp(icoord) = xvec(icoord)+dble(ih)*h
C!Call Int_To_Cart(InterVec,xtmp,AtCoord,NumOfAt,NumInt,Mass)
      Call Int_to_Cart1(InterVec,xtmp,AtCoord,
     &          NumOfAt,NumInt       )
c             Stemp = 0.0d0
      call dcopy_(3*NumOfAt*NumInt,0.0d0,0,Work(ipStemp),1)

      Call CalcS(AtCoord,InterVec,Work(ipStemp),NumInt,NumOfAt)
c Gtemp(1,1,iterm)=1+NumInt*( 0 + NumInt*(iterm-1))-1
c             Call CalcG(Gtemp(1,1,iterm),Mass,Work(ipStemp))
      Call CalcG(Work(ipGtemp+NumInt*NumInt*(iterm-1)),
     &       Mass,Work(ipStemp),NumInt,NumOfAt)

      iterm = iterm+1
      End Do
      Do k = 1,NumInt
      Do j = 1,NumInt
c Gtemp(j,k,1)=j+NumInt*(k-1+NumInt*(1-1))-1
c Gtemp(j,k,2)=j+NumInt*(k-1+NumInt*(2-1))-1
c Gtemp(j,k,3)=j+NumInt*(k-1+NumInt*(3-1))-1

c                Gprime(j,k,icoord) = (Gtemp(j,k,1)-27.0d0*Gtemp(j,k,2)+
c     &                     27.0d0*Gtemp(j,k,3)-Gtemp(j,k,4))/(48.0d0*h)
      Gprime(j,k,icoord) =
     &          (Work(ipGtemp+j+NumInt*(k-1)-1)-
     &             27.0d0*Work(ipGtemp+j+NumInt*(k-1+NumInt)-1)+
     &             27.0d0*Work(ipGtemp+j+NumInt*(k-1+NumInt*2)-1)-
     &                    Work(ipGtemp+j+NumInt*(k-1+NumInt*3)-1))/
     &                    (48.0d0*h)
      End Do
      End Do
      End Do
C!Call Int_To_Cart(InterVec,xvec,AtCoord,NumOfAt,NumInt,Mass)
      Call Int_to_Cart1(InterVec,xtmp,AtCoord,NumOfAt,NumInt)
C!
      Call GetMem('Stemp','Free','Real',ipStemp,3*NumOfAt*NumInt)
      Call GetMem('Gtemp','Free','Real',ipGtemp,NumInt*NumInt*4)
C!
      End
C!
C!-----------------------------------------------------------------------!
C!-----------------------------------------------------------------------!
C!
      Subroutine CalcGdbleprime(Gdbleprime,Mass,xvec,InterVec,
     &       AtCoord,NumOfAt,h,NumInt)
C!
C!  Purpose:
C!    Calculate second derivatives of G.
C!
C!  Input:
C!    Mass       : Real*8 array - the mass of the atoms.
C!    xvec       : Real*8 array - the geometry in internal
C!                 coordinates.
C!    InterVec   : Integer array.
C!    NumOfAt    : Integer - the number of atoms.
C!
C!  Output:
C!    Gdbleprime : Real*8 two dimensional array - second
C!                 derivative of the inverse mass tensor G.
C!
C!  Written by:
C!    Niclas Forsberg,
C!    Dept. of Theoretical Chemistry, Lund University, 1996.
C!
c       Implicit None
#include "Constants_mula.fh"
#include "dims.fh"
      Real*8   h
      Integer  icoord,jcoord,k,j
      Integer  NumInt,NumOfAt
      Real*8 Gdbleprime(ngdim,ngdim,ngdim,ngdim)
      Real*8 Mass  (NumOfAt)
      Real*8 xvec  (NumInt)
      Integer InterVec(*)
      Real*8 AtCoord (3,NumOfAt)
#include "WrkSpc.fh"
C!
C!---- Initialize.
      Call GetMem('xtmp','Allo','Real',ipxtmp,NumInt)
      NumInt3=NumInt*NumInt*NumInt
      Call GetMem('Gprime1','Allo','Real',ipGprime1,NumInt3)
      Call GetMem('Gprime2','Allo','Real',ipGprime2,NumInt3)
      Call GetMem('Gprime3','Allo','Real',ipGprime3,NumInt3)
      Call GetMem('Gprime4','Allo','Real',ipGprime4,NumInt3)
C!
      Do jcoord = 1,NumInt
      call dcopy_(NumInt,xvec,1,Work(ipxtmp),1)
c          xtmp = xvec
      Work(ipxtmp+jcoord-1) = xvec(jcoord)-3.0d0*h
      Call CalcGprime(Work(ipGprime1),Mass,
     &    Work(ipxtmp),InterVec,AtCoord,NumOfAt,h,NumInt)
      Work(ipxtmp+jcoord-1) = xvec(jcoord)-1.0d0*h
      Call CalcGprime(Work(ipGprime2),Mass,
     &    Work(ipxtmp),InterVec,AtCoord,NumOfAt,h,NumInt)
      Work(ipxtmp+jcoord-1) = xvec(jcoord)+1.0d0*h
      Call CalcGprime(Work(ipGprime3),Mass,
     &    Work(ipxtmp),InterVec,AtCoord,NumOfAt,h,NumInt)
      Work(ipxtmp+jcoord-1) = xvec(jcoord)+3.0d0*h
      Call CalcGprime(Work(ipGprime4),Mass,
     &    Work(ipxtmp),InterVec,AtCoord,NumOfAt,h,NumInt)
      Do icoord = 1,NumInt
      Do k = 1,NumInt
      Do j = 1,NumInt
c Gprime(j,k,icoord)= j+NumInt*(k-1+NumInt*(icoord-1))-1
      ivv=j+NumInt*(k-1+NumInt*(icoord-1))-1
      Gdbleprime(j,k,icoord,jcoord) =
     &                        (Work(ipGprime1+ivv)-
     &                         27.0d0*Work(ipGprime2+ivv)+
     &                         27.0d0*Work(ipGprime3+ivv)-
     &                                Work(ipGprime4+ivv))/(48.0d0*h)
      End Do
      End Do
      End Do
      End Do
C!Call Int_To_Cart(InterVec,xvec,AtCoord,NumOfAt,NumInt,Mass)
      Call Int_To_Cart1(InterVec,xvec,AtCoord,NumOfAt,NumInt)
C!
      Call GetMem('xtmp','Free','Real',ipxtmp,NumInt)
      Call GetMem('Gprime1','Free','Real',ipGprime1,NumInt3)
      Call GetMem('Gprime2','Free','Real',ipGprime2,NumInt3)
      Call GetMem('Gprime3','Free','Real',ipGprime3,NumInt3)
      Call GetMem('Gprime4','Free','Real',ipGprime4,NumInt3)
C!
      End
