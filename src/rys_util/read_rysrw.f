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
* Copyright (C) 1990,1991, Roland Lindh                                *
************************************************************************
      subroutine read_rysrw
************************************************************************
*                                                                      *
* Object: to setup the coefficients for the Rys roots and weights.     *
*                                                                      *
*     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
*             September '90                                            *
*             Modified to DaFile February '91                          *
************************************************************************
      use vRys_RW
CVV: some variables used under #ifdef are not defined.
c      implicit none
#include "SysDef.fh"
#include "itmax.fh"
#include "stdalloc.fh"
      character(*), parameter :: RYSRW_NAME = 'RYSRW'
      integer, parameter :: lu_rysrw = 22
      logical :: found_rysrw
*
      integer :: mRys, nOrder
      real*8 :: acc(maxrys)
*
      integer :: iRys
      integer :: nMap_Tot, nMem_Tot, nx0_Tot, nMem
      integer :: io
*
*     Open file for data base
*
      call f_Inquire(RYSRW_NAME,found_rysrw)
      if (.not.found_rysrw) then
        call warningmessage(2,
     &              ' the rysrw.ascii file does not exist.')
        call abend()
      end if
      call molcas_open(lu_rysrw,RYSRW_NAME)
#ifdef _DEBUG_
      Write (6,*) ' nDisk=',nDisk
#endif
*
*     Read initial data
*
      io = 1
      Do While (io .ne. 0)
        Read (lu_rysrw,*,IOStat=io) mRys,nOrder,Acc
      End Do
      If (mRys.gt.MaxRys) Then
         Call WarningMessage(2,
     &      ' Database requires new code!'//
     &      ' Database and code are at incompatible levels!')
         Call Abend()
      End If
      nMxRys=mRys
#ifdef _DEBUG_
      Write (6,*)
      Write (6,*) ' Reading tables for roots and weights of Rys poly.'
      Write (6,*) ' Highest order is:',mRys
      Write (6,*) ' Order of approximating polynomial:',nOrder
      Write (6,*) ' Relative accuracy of computed values:',(Acc(i),
     &              i=1,mRys)
      Write (6,*)
#endif
*
*     Read value of T at which asymptotic formulas will be used
*
      Call mma_allocate(TMax,mRys,label='TMax')
      Call InR(Tmax,mRys,lu_rysrw)
#ifdef _DEBUG_
      Call RecPrt(' Tmax',' ',Tmax,mRys,1)
#endif
*
*     Read increment of tables
*
      Call mma_allocate(ddx,mRys,label='ddx')
      Call InR(ddx,mRys,lu_rysrw)
#ifdef _DEBUG_
      Call RecPrt(' ddx ',' ',ddx,mRys,1)
#endif
*
*     Read size of map array
*
      Call InI(nMap,mRys,lu_rysrw)
#ifdef _DEBUG_
      Write (6,*) ' nMap=',nMap
#endif
*
*     Read number of subranges
*
      Call InI(nx0,mRys,lu_rysrw)
#ifdef _DEBUG_
      Write (6,*) ' nx0=',nx0
#endif
*
*     Read map array and x0 array for each order of Rys polynomials
*
      nMap_Tot = 0
      nx0_Tot  = 0
      Do iRys = 1, mRys
         iMap(iRys) = nMap_Tot + 1
         nMap_Tot = nMap_Tot + nMap(iRys)
         ix0(iRys) = nx0_Tot + 1
         nx0_Tot = nx0_Tot + nx0(iRys)
      End Do
      call mma_allocate(Map,nMap_Tot,label='Map')
      Call mma_allocate(x0,nx0_Tot,label='x0')
      Do iRys = 1, mRys
        Call InI(Map(iMap(iRys)),nMap(iRys),lu_rysrw)
*
        Call InR(x0(ix0(iRys)),nx0(iRys),lu_rysrw)
      End Do
*
*     Allocate memory for coefficients
*
      nMem_Tot = 0
      Do iRys = 1, mRys
         iCffR(0,iRys) = nMem_Tot + 1
         nMem=nx0(iRys)*iRys
         nMem_Tot = nMem_Tot + 14*nMem
      End Do
      Call mma_allocate(Cff,nMem_Tot,label='Cff')
      Do iRys = 1, mRys
*
*     Read coefficients from file
*
        nMem=nx0(iRys)*iRys
        iCffR(1,iRys) = iCffR(0,iRys) + nMem
        iCffR(2,iRys) = iCffR(1,iRys) + nMem
        iCffR(3,iRys) = iCffR(2,iRys) + nMem
        iCffR(4,iRys) = iCffR(3,iRys) + nMem
        iCffR(5,iRys) = iCffR(4,iRys) + nMem
        iCffR(6,iRys) = iCffR(5,iRys) + nMem
*
        ICffW(0,iRys) = iCffR(6,iRys) + nMem
        iCffW(1,iRys) = iCffW(0,iRys) + nMem
        iCffW(2,iRys) = iCffW(1,iRys) + nMem
        iCffW(3,iRys) = iCffW(2,iRys) + nMem
        iCffW(4,iRys) = iCffW(3,iRys) + nMem
        iCffW(5,iRys) = iCffW(4,iRys) + nMem
        iCffW(6,iRys) = iCffW(5,iRys) + nMem
*
        Call InR(Cff(iCffR(0,iRys)),nMem*14,lu_rysrw)
*
      End Do
*
      Close (lu_rysrw)
*
      Return
      End

      Subroutine InR(A,n,Lu)
      Implicit None
      Integer n, Lu
      Real*8 A(n)
      Integer i, iEnd, j
      Do i = 1, n, 3
        iend = Min(i+2,n)
c The numbers are actually E21.15, but some compilers
c warn about the size, and it shouldn't matter for reading
        Read (Lu,'(3E21.14)') (A(j),j=i,iend)
      End Do
      Return
      End

      Subroutine InI(A,n,Lu)
      Implicit None
      Integer n, Lu
      Integer A(n)
      Integer i, iEnd, j
      Do i = 1, n, 3
        iend = Min(i+2,n)
        Read (Lu,*) (A(j),j=i,iend)
      End Do
      Return
      End
