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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Aug. 06, 2020, created this file.               *
! ****************************************************************
      Subroutine QaaVerif(G2q,ng2,PUVX,NPUVX,IndTUVX)
      use MCLR_Data, only: nNA
      use input_mclr, only: ntAsh
      Implicit None
      INTEGER nG2,nPUVX
      Real*8,DIMENSION(nG2)::G2q
      Real*8,DIMENSION(NPUVX)::PUVX
      INTEGER,DIMENSION(ntAsh,ntAsh,ntAsh,ntAsh)::IndTUVX
      INTEGER I,J,K,L,IJKL,lMax
      Real*8 dQdX

      ijkl=0
      dQdX=0.0d0
      do i=1,nna
        do j=1,i
          do k=1,i
            if(i.eq.k) then
              lmax = j
            else
              lmax = k
            end if
            do l=1,lmax
              ijkl = ijkl + 1
              dQdX=dQdX+G2q(ijkl)*PUVX(IndTUVX(I,J,K,L))
            end do
          end do
        end do
      end do

      write(6,*) 'dQdX in QaaVerif=',dQdX

      End Subroutine QaaVerif
