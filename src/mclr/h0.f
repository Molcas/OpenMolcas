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
* Copyright (C) 1997, Anders Bernhardsson                              *
************************************************************************
      Subroutine H0(rdia,MP1,MP2,MQ,isym,nprciv,TimeDep)
      Use Exp, only: H0S, H0F, SBIDT
      use negpre
      Use Iso_C_Binding
      use Arrays, only: Int2, FIMO
*
* frontend to jeppes explicit routines
*
      implicit Real*8(a-h,o-z)
#include "detdim.fh"
#include "Pointers.fh"
#include "orbinp_mclr.fh"
#include "cstate_mclr.fh"
#include "crun_mclr.fh"
#include "cicisp_mclr.fh"
#include "spinfo_mclr.fh"
#include "incdia.fh"
#include "stdalloc.fh"
      Real*8 rdia(*)
      Logical TimeDep
      Real*8, Allocatable:: H0T(:), Vec2(:)
      Real*8, Allocatable, Target:: H0Scr(:)
      Integer, Pointer:: iH0Scr(:)
      Integer, Allocatable:: SBCNF(:)
*
      MXP1=MP1
      MXP2=MP2
      MXQ=MQ
      ispc=1
      NDET= NINT(XISPSM(ISYM,ISPC))
*
      NSBDET = MXP1 + MXP2 + MXQ
      MXP = MXP1 + MXP2
      LH0T= MXP*(MXP+1)/2 + MXP1*MXQ
      MXCSFC = 0
      MXDTFC = 0
      DO  ITYP = 1, NTYP
        MXCSFC = MAX(MXCSFC,NCPCNT(ITYP) )
        MXDTFC = MAX(MXDTFC,NDPCNT(ITYP) )
      End Do

      nactel=naelci(1)+nbelci(1)
      If (TimeDep) Then
         EnA=E2_td(FIMO,Int2,0,-1)
      Else
         EnA=E2(FIMO,Int2,0,-1)
      End If
      LH0SCR = MAX(6*NSBDET,4*NSBDET+4*NOCOB,MXP1*(MXP1+1)/2+MXP1**2)
      LVEC2 = 2 * NACTEL + MXCSFC**2
     &      + 6*MXDTFC+2*MXDTFC**2
     &      + MAX(MXDTFC*NACTEL+2*NACTEL,4*NACOB+2*NACTEL)
      LVEC2=Max(lvec2,ndet)
*
      If (isym.eq.irefsm) then
         ieaw=1
      else
         ieaw=2
      end if
      CALL mma_allocate(SBIDT,NSBDET,Label='SBIDT')
      CALL mma_allocate(H0S,MXP**2,Label='H0S')
      CALL mma_allocate(H0F,MXP,Label='H0F')
      CALL mma_allocate(H0T,LH0T,Label='H0T')
      CALL mma_allocate(SBCNF,NSBDET,Label='SBCNF')
      Call mma_allocate(H0SCR,LH0SCR,Label='H0Scr')
      Call C_F_Pointer(C_Loc(H0SCR),iH0Scr,[LH0SCR])
      Call mma_allocate(VEC2,lvec2,Label='Vec2')
*
      Call H0MAT_MCLR(H0T,SBIDT,SBCNF,
     &                MXP1,MXP2,MXQ,NACOB,NPRCIV,
     &                NOCSF,ISYM,IDC,PSSIGN,ECOREP,
     &                rDIA,Vec2,H0Scr,iH0Scr,ieaw)

*
      do i=1,nprciv
         H0T(i*(i+1)/2)=H0T(i*(i+1)/2)-ENA
      End Do
      IF (NGP)  Call mkp1(nprciv,SBIDT,H0T,rdia)
*     Call Triprt('PRECI',' ',H0T,nprciv)
*     write(*,*) (SBIDT(i),i=1,nprciv)
      Call mma_deallocate(Vec2)
      Nullify(iH0Scr)
      Call mma_deallocate(H0Scr)
      Call mma_deallocate(SBCNF)
      call square(H0T,H0S,1,NPRCIV,NPRCIV)
      Call mma_deallocate(H0T)
*
      irc=0
      call dgetrf_(NPRCIV,NPRCIV,H0S,NPRCIV,H0F,irc)
      If (irc.ne.0) Then
         Write(6,*) 'Sorry but you have an singular ci matrix'
         Write(6,*) 'Set ExpDimension and restart mclr'
         Call Abend()
      End If
*
      Return
      End
