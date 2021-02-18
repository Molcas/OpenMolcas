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

      SUBROUTINE INTDIA(DIAG,NSPC,ISPC,ISM,LSPC,IAMCMP,ecore)
      Use Str_Info
*
* CI diagonal in SD basis for the NCSPC ci spaces defined by
* ISPC,ISM
*
* if IAMCMP .ne. 0 : then it is assumed that a complex
* hermitian eigenvalued problem is being solved by
* doubling the dimensions. The diagonal is then
* constructed and written out twice
*
      IMPLICIT REAL*8(A-H,O-Z)
*
* ==============
*.Specific Input
* ==============
*
      INTEGER ISPC(NSPC),LSPC(NSPC),ISM(NSPC)
*
* ==============
*.General Input
* ==============
*
*./ORBINP/ : NACOB used
*
#include "detdim.fh"
#include "orbinp_mclr.fh"
#include "cicisp_mclr.fh"
#include "cstate_mclr.fh"
#include "csm.fh"
#include "stdalloc.fh"
#include "cprnt_mclr.fh"
#include "spinfo_mclr.fh"
#include "crun_mclr.fh"
#include "dmrginfo_mclr.fh"
      Real*8 DIAG(*)
      Integer idum(1)
      Real*8, Allocatable:: JA(:), KA(:), XA(:), XB(:), SCR(:), H1D(:)
      Integer, Allocatable:: BLTP(:), IOIO(:)
*
* ======
*.Output
* ======
*
*       OBS THIS WILL JUST WORK FOR CASSCF/RASSCF RESPONSE
      LUDIA=0

      if(doDMRG)then  ! yma
         call dmrg_dim_change_mclr(RGras2(1:8),ntoob,0)
         call dmrg_dim_change_mclr(RGras2(1:8),nacob,0)
      end if

*
*
**. Local memory
*
      Call mma_allocate(JA  ,NTOOB**2, Label='JA')
      Call mma_allocate(KA  ,NTOOB**2, Label='KA')
      Call mma_allocate(XA  ,NACOB,Label='XA')
      Call mma_allocate(XB  ,NACOB,Label='XB')
      Call mma_allocate(SCR ,2*NACOB,Label='SCR')
      if(doDMRG)then !yma
* wired Call mma_allocate(H1D ,NACOB_KLH1D,Label='H1D')
        Call mma_allocate(H1D ,NACOB,Label='H1D')
      else
        Call mma_allocate(H1D ,NACOB,Label='H1D')
      end if
      Call mma_allocate(BLTP,NSMST,Label='BLTP')

*
*. Largest NOCTPA*NOCTPB block
      MXOCOC = 0
      DO IISPC = 1, NSPC
        NOCTPA = NOCTYP(IASTFI(ISPC(IISPC)))
        NOCTPB = NOCTYP(IBSTFI(ISPC(IISPC)))
        MXOCOC = MAX(MXOCOC,NOCTPA*NOCTPB)
      END DO
      Call mma_allocate(IOIO,NOCTPA*NOCTPB,Label='IOIO')
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
      CALL GT1DIA_MCLR(H1D)
      CALL GTJK_MCLR(JA,KA)
*
*. K goes to J - K
      ONE = +1.0D0
      ONEG = -1.0D0
      CALL VECSUM(KA,KA,JA,ONEG,ONE,NTOOB **2)
*
*
*. Loop over internal CI spaces
*
      IF(IAMCMP.EQ.0) THEN
        NLOOP = 1
      ELSE
        NLOOP = 2
      END IF
      DO 200 ILOOP = 1, NLOOP
      DO 100 IISPC = 1, NSPC
        IATP = IASTFI(ISPC(IISPC))
        IBTP = IBSTFI(ISPC(IISPC))
        NAEL = NELEC(IATP)
        NBEL = NELEC(IBTP)
        NOCTPA = NOCTYP(IATP)
        NOCTPB = NOCTYP(IBTP)
        MNRS1C = MNR1IC(ISPC(IISPC))
        MXRS3C = MXR3IC(ISPC(IISPC))
*
        CALL ZBLTP(ISMOST(1,ISM(IISPC)),NSMST,IDC,BLTP,idum)
        CALL IAIBCM_MCLR(MNRS1C,MXRS3C,NOCTPA,NOCTPB,
     &                   Str(IATP)%EL1,Str(IATP)%EL3,
     &                   Str(IBTP)%EL1,Str(IBTP)%EL3,
     &                   IOIO,IPRDIA)
*
        IF(ICISTR.LE.1) THEN
          LLUDIA = 0
        ELSE
          LLUDIA = LUDIA
        END IF
        CALL CIDIA4(NAEL,Str(IATP)%OCSTR,NBEL,Str(IBTP)%OCSTR,
     &              NACOB,DIAG,NSMST,H1D,
     &              ISMOST(1,ISM(IISPC)),BLTP,
     &              XA,XB,SCR,JA,KA,Str(IATP)%NSTSO,Str(IBTP)%NSTSO,
     &              IOIO,NOCTPA,NOCTPB,Str(IATP)%ISTSO,
     &              Str(IBTP)%ISTSO,LLUDIA,ECORE,
     &              PLSIGN,PSSIGN,IPRDIA,NTOOB,ICISTR)
*
      IF(ICISTR.LE.1.AND.LUDIA.GT.0) THEN
*. Each CI space is written in one record
        CALL ITODS(LSPC(IISPC),1,0,LUDIA)
        CALL TODSC_MCLR(DIAG,LSPC(IISPC),0,LUDIA)
      END IF
  100 CONTINUE
*. Write end of vector mark
      IF(LUDIA.GT.0) CALL ITODS([-1],1,0,LUDIA)
  200 CONTINUE

*
      Call mma_deallocate(IOIO)
      Call mma_deallocate(BLTP)
      Call mma_deallocate(H1D)
      Call mma_deallocate(SCR)
      Call mma_deallocate(XB)
      Call mma_deallocate(XA)
      Call mma_deallocate(KA)
      Call mma_deallocate(JA)

      if(doDMRG)then  ! yma
         call dmrg_dim_change_mclr(LRras2(1:8),ntoob,0)
         call dmrg_dim_change_mclr(LRras2(1:8),nacob,0)
      end if

      RETURN
      END
