************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine  pgamma
************************************************************************
*                                                                      *
* Object: to compute the arrays gammath and gammaph for the angular    *
*         integration within a R-matrix run                            *
*                                                                      *
************************************************************************
      Implicit real*8 (a-h,o-z)
*     external gammat,gammaf
#include "real.fh"
#include "rmat.fh"
#include "gam.fh"
*
      iPrint=0
*
*
* initialize arrays
      do m=-2,2*lgamma+2
       do n=-2,2*lgamma+2
        gammath(m,n)=0.0d0
        gammaph(m,n)=0.0d0
       End Do
      End Do
*
*   Set intitial values for recursion of gammath
      gammath(0,0)=2.0d0
      gammath(1,0)=Pi/two
*
*  m=0
       m=0
       do n=0,2*lgamma+2,2
         gammath(0,n+2)=DBLE(n+1)/DBLE(m+n+3)*gammath(0,n)
*        hgt=gammat(x)
*        write(*,*)' m,n,gammath,hgt',m,n,gammath(m,n),hgt
       End Do
       do n=1,2*lgamma+2,2
         gammath(0,n)=0.0d0
*        hgt=gammat(x)
*        write(*,*)' m,n,gammath,hgt',m,n,gammath(m,n),hgt
       End Do
*
*  m.gt.0
      do m=1,2*lgamma+2
       do n=0,2*lgamma+2,2
         gammath(m,n+2)=DBLE(n+1)/DBLE(m+n+3)*gammath(m,n)
*        hgt=gammat(x)
*        write(*,*)' m,n,gammath,hgt',m,n,gammath(m,n),hgt
       End Do
       do n=1,2*lgamma+2,2
         gammath(m,n)=0.0d0
*        hgt=gammat(x)
*        write(*,*)' m,n,gammath,hgt',m,n,gammath(m,n),hgt
       End Do
        gammath(m+1,0)=DBLE(m+1)/DBLE(m+2)*gammath(m-1,0)
      End Do
*
*
*   Set intitial values for recursion of gammaph
      gammaph(0,0)=2.0d0*pi
      gammaph(1,0)=0.0d0
      gammaph(0,1)=0.0d0
*
*  m=0
      m=0
      do n=0,2*lgamma+2
        gammaph(0,n+2)=DBLE(n+1)/DBLE(m+n+2)*gammaph(0,n)
*       hgp=gammaf(x)
*       write(*,*)' m,n,gammaph,hgp',m,n,gammaph(m,n),hgp
      End Do
*
*  m.gt.0
      do m=1,2*lgamma+2
       do n=0,2*lgamma+2
         gammaph(m,n+2)=DBLE(n+1)/DBLE(m+n+2)*gammaph(m,n)
*        hgp=gammaf(x)
*        write(*,*)' m,n,gammaph,hgp',m,n,gammaph(m,n),hgp
       End Do
        gammaph(m+1,0)=DBLE(m)/DBLE(m+1)*gammaph(m-1,0)
      End Do
*                                                                      *
************************************************************************
*                                                                      *

      Return
      End
