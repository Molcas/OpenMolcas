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
        subroutine DefParReord (NaGrpR,maxdim)

c
c       This routine do:
c       define parameters in chcc_reord.fh using NaGrpR,maxdim
c
c       I/O parameter description:
c       NxGrpR   - # of groups in a (=b) set (I)
c       maxdim   - # maximal dimension of a (=b) Groups(O)
c
        implicit none
#include "chcc1.fh"
#include "chcc_reord.fh"
#include "chcc_files.fh"
c
        integer NaGrpR,maxdim
c
c       help variables
c
        real*8 rdim
        integer i,j
        integer Up(1:MaxGrp),Low(1:MaxGrp)
c
c
c1      define parameters of Groups of a set
c
        rdim=1.0d0*nv/(1.0d0*NaGrpR)
c
        do i=1,NaGrpR
c
           if (i.eq.1) then
             Up(i)=int(rdim*i)
             Low(i)=1
           else if (i.eq.NaGrpR) then
             Up(i)=nv
             Low(i)=Up(i-1)+1
           else
             Up(i)=int(rdim*i)
             Low(i)=Up(i-1)+1
           end if
c
           DimGrpaR(i)=(Up(i)-Low(i))+1
c
        end do
c
c
c2      find maximal dimensions of a'
c
        maxdim=DimGrpaR(1)
        do i=1,NaGrpR
          if (DimGrpaR(i).gt.maxdim) then
          maxdim=DimGrpaR(i)
          end if
        end do
c
c
c3.1    def L2Name, T2Name, I2Name,I3Name
c
        do i=1,MaxGrp
        do j=1,MaxGrp
          call DefParo3v3Hlp1(i,j,'L2',L2Name(i,j))
          call DefParo3v3Hlp1(i,j,'T2',T2Name(i,j))
          call DefParo3v3Hlp1(i,j,'I2',I2Name(i,j))
          call DefParo3v3Hlp1(i,j,'I3',I3Name(i,j))
        end do
        end do
c
c3.2    def L1Name,I1Name
c
        do i=1,MaxGrp
          call DefParo3v3Hlp2 (i,'L1vc',L1Name(i))
          call DefParo3v3Hlp2 (i,'I1in',I1Name(i))
        end do
c
c3.3    def L0Name,I0Name
c
        L0Name='L0vctr'
        I0Name='I0intg'
c
        return
        end
