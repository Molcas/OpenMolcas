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
      Subroutine CovRadT_Init()
      Implicit Real*8 (a-h,o-z)
#include "covradt_data.fh"
      Real*8 Temp(0:92)
*     From: http://wulff.mit.edu/pt/pert5.html
c     Same as below but with H,C,N,O from urocanic acid article
      Data Temp/
     & 0.00d0,
     & 3.08d0,                                              0.93d0,
     & 1.23d0,.90d0,     0.82d0,3.59d0,3.31d0,3.08d0,0.72d0,0.71d0,
     & 1.54d0,.36d0,     1.18d0,1.11d0,1.06d0,1.02d0,0.99d0,0.98d0,
     & 2.03d0,.74d0,
     &              1.44d0,1.32d0,1.22d0,1.18d0,1.17d0,
     &              1.17d0,1.16d0,1.15d0,1.17d0,1.25d0,
     &                   1.26d0,1.22d0,1.20d0,1.16d0,1.14d0,1.89d0,
     & 2.16d0,.91d0,
     &              1.62d0,1.45d0,1.34d0,1.3d0,1.27d0,
     &              1.25d0,1.25d0,1.28d0,1.34d0,1.41d0,
     &                   1.44d0,1.41d0,1.40d0,1.36d0,1.33d0,1.31d0,
     & 2.35d0,.98d0,1.25d0,
     &              1.65d0,1.65d0,1.64d0,1.63d0,1.62d0,1.85d0,1.61d0,
     &              1.59d0,1.59d0,1.58d0,1.57d0,1.56d0,1.70d0,1.56d0,
     &                     1.44d0,1.34d0,1.30d0,.28d0,
     &              1.26d0,1.27d0,1.3 d0,1.34d0,1.49d0,
     &                   1.48d0,1.47d0,1.46d0,1.53d0,1.47d0,1.3d0,
     & 2.50d0,2.05d0,1.50d0,
     &                    1.65d0,1.50d0,1.42d0/
*
       Do i = 0, 92
          CovRadT_(i)=Temp(i)
       End Do
c
c     Original radii - totally weird
c      Data CovRadT_/
c     & 0.00d0,
c     & 0.32d0,                                              0.93d0,
c     & 1.23d0,.90d0,     0.82d0,0.77d0,0.75d0,0.73d0,0.72d0,0.71d0,
c     & 1.54d0,.36d0,     1.18d0,1.11d0,1.06d0,1.02d0,0.99d0,0.98d0,
c     & 2.03d0,.74d0,
c     &              1.44d0,1.32d0,1.22d0,1.18d0,1.17d0,
c     &              1.17d0,1.16d0,1.15d0,1.17d0,1.25d0,
c     &                   1.26d0,1.22d0,1.20d0,1.16d0,1.14d0,1.89d0,
c     & 2.16d0,.91d0,
c     &              1.62d0,1.45d0,1.34d0,1.3d0,1.27d0,
c     &              1.25d0,1.25d0,1.28d0,1.34d0,1.41d0,
c     &                   1.44d0,1.41d0,1.40d0,1.36d0,1.33d0,1.31d0,
c     & 2.35d0,.98d0,1.25d0,
c     &              1.65d0,1.65d0,1.64d0,1.63d0,1.62d0,1.85d0,1.61d0,
c     &              1.59d0,1.59d0,1.58d0,1.57d0,1.56d0,1.70d0,1.56d0,
c     &                     1.44d0,1.34d0,1.30d0,.28d0,
c     &              1.26d0,1.27d0,1.3 d0,1.34d0,1.49d0,
c     &                   1.48d0,1.47d0,1.46d0,1.53d0,1.47d0,1.3d0,
c     & 2.50d0,2.05d0,1.50d0,
c     &                    1.65d0,1.50d0,1.42d0/

*     Following  elements has been faked by me: 86, 87, 88, 89, 91
c     From: H.B. Schlegel, TCA, 66, 333 (1984).
c     Data CovRadT/
c    &  0.0d0,
c    &  0.32D0,0.60D0,
c    &  1.2D0,1.05D0,0.81D0,0.77D0,2*0.74D0,2*0.72D0,
c    &  1.5D0,1.40D0,1.30D0,1.17D0,1.10D0,1.04D0,2*0.99D0,
c    &  1.8D0,1.60D0,11*1.4D0,1.3D0,2*1.2D0,2*1.1D0/
*
      Return
      End
