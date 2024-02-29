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
      Subroutine GugaCtl_dmrg()
*
      use stdalloc, only: mma_allocate, mma_deallocate
      use gugx, only: IFCAS, SGS, CIS
      Implicit Real*8 (A-H,O-Z)
*
#include "Input.fh"
#include "Pointers.fh"
#include "detdim.fh"
#include "spinfo_mclr.fh"
*
      Interface
      SUBROUTINE MKGUGA(STSYM,Skip_MKSGNUM)
      IMPLICIT None

      Integer STSYM
      Logical, Optional:: Skip_MKSGNUM
      End SUBROUTINE MKGUGA
      End Interface

      Associate ( nLev => SGS%nLev,
     &            LM1RAS=>SGS%LM1RAS, LM3RAS=>SGS%LM3RAS,
     &            LV1RAS=>SGS%LV1RAS, LV3RAS=>SGS%LV3RAS,
     &            IA0 => SGS%IA0, IB0 => SGS%IB0, IC0 => SGS%IC0)
*
      ntRas1=0
      ntRas2=0
      ntRas3=0
      Do iSym=1,nSym
         ntRas1=ntRas1+nRs1(iSym)
         ntRas2=ntRas2+nRs2(iSym)
         ntRas3=ntRas3+nRs3(iSym)
      End Do
*
      B0=iSpin-1
      A0=(nActEl-B0)/2
      C0=ntASh-A0-B0
      If ( (2*A0+B0).ne.nActEl ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine GUGACTL ***'
         Write (6,*) ' 2*A0+B0.ne.nActEl '
         Write (6,*)
      End If
      If ( A0.lt.0 ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine GUGACTL ***'
         Write (6,*) ' A0.lt.0'
         Write (6,*)
      End If
      If ( B0.lt.0 ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine GUGACTL ***'
         Write (6,*) ' B0.lt.0'
         Write (6,*)
      End If
      If ( C0.lt.0 ) then
         Write (6,*)
         Write (6,*) ' *** Error in subroutine GUGACTL ***'
         Write (6,*) ' C0.lt.0'
         Write (6,*)
      End If
*
      Call mma_allocate(SGS%ISM,ntAsh,Label='SGS%ISM')
      iOrb=0
      Do iSym=1,nSym
         Do iBas=1,nRs1(iSym)
            iOrb=iOrb+1
            SGS%ISM(iOrb)=iSym
         End Do
      End Do
      Do iSym=1,nSym
         Do iBas=1,nRs2(iSym)
            iOrb=iOrb+1
            SGS%ISM(iOrb)=iSym
         End Do
      End Do
      Do iSym=1,nSym
         Do iBas=1,nRs3(iSym)
            iOrb=iOrb+1
            SGS%ISM(iOrb)=iSym
         End Do
      End Do
*
      NLEV=ntASh
      LV1RAS=ntRas1
      LV3RAS=LV1RAS+ntRas2
      LM1RAS=2*LV1RAS-nHole1
      LM3RAS=nActEl-nElec3
!
      IFCAS=1
      SGS%nLev=nLev
      SGS%nSym=nSym
      Call mkGUGA(State_Sym)
      NCSF(1:nSym) = CIS%NCSF(1:nSym)
      NCONF=CIS%NCSF(State_Sym)

*
      Call mkGUGA_Free()

      End Associate

      End Subroutine GugaCtl_dmrg
