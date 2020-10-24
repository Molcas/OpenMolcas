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
#include "strbas_mclr.fh"
#include "cstate_mclr.fh"
#include "strinp_mclr.fh"
#include "stinf_mclr.fh"
#include "csm.fh"
#include "WrkSpc.fh"
#include "cprnt_mclr.fh"
#include "spinfo_mclr.fh"
#include "crun_mclr.fh"
#include "dmrginfo_mclr.fh"
      Real*8 DIAG(*)
*
* ======
*.Output
* ======
*
*       OBS THIS WILL JUST WORK FOR CASSCF/RASSCF RESPONSE
      ISYM=ism(1)
      LUDIA=0

      if(doDMRG)then  ! yma
         call dmrg_dim_change_mclr(RGras2(1:8),ntoob,0)
         call dmrg_dim_change_mclr(RGras2(1:8),nacob,0)
      end if

*
*
**. Local memory
*
      Call GetMem('KLJ   ','ALLO','REAL',KLJ   ,NTOOB**2)
      Call GetMem('KLK   ','ALLO','REAL',KLK   ,NTOOB**2)
      Call GetMem('KLSC2 ','ALLO','REAL',KLSCR2,2*NTOOB**2)
      Call GetMem('KLXA  ','ALLO','REAL',KLXA  ,NACOB)
      Call GetMem('KLXB  ','ALLO','REAL',KLXB  ,NACOB)
      Call GetMem('KLSCR ','ALLO','REAL',KLSCR ,2*NACOB)
      if(doDMRG)then !yma
* wired Call GetMem('KLH1D ','ALLO','REAL',KLH1D ,NACOB_KLH1D)
        Call GetMem('KLH1D ','ALLO','REAL',KLH1D ,NACOB)
      else
        Call GetMem('KLH1D ','ALLO','REAL',KLH1D ,NACOB)
      end if
*     Call GetMem('KLSMOS','ALLO','REAL',KLSMOS,NSMST)
      Call GetMem('KLSMO2','ALLO','INTE',KLBLTP,NSMST)

*
        KLSVST = 1
*. Largest NOCTPA*NOCTPB block
      MXOCOC = 0
      DO IISPC = 1, NSPC
        NOCTPA = NOCTYP(IASTFI(ISPC(IISPC)))
        NOCTPB = NOCTYP(IBSTFI(ISPC(IISPC)))
        MXOCOC = MAX(MXOCOC,NOCTPA*NOCTPB)
      END DO
      Call GetMem('KLIOIO','ALLO','INTE',KLIOIO,NOCTPA*NOCTPB)
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
      CALL GT1DIA_MCLR(WORK(KLH1D))
      CALL GTJK_MCLR(WORK(KLJ),WORK(KLK))
*
*. K goes to J - K
      ONE = +1.0D0
      ONEG = -1.0D0
      CALL VECSUM(WORK(KLK),WORK(KLK),WORK(KLJ),
     &            ONEG,ONE,NTOOB **2)
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
        CALL ZBLTP(ISMOST(1,ISM(IISPC)),NSMST,IDC,
     &       iWORK(KLBLTP),iWORK(KLSVST))
        CALL IAIBCM_MCLR(MNRS1C,MXRS3C,NOCTPA,NOCTPB,
     &       iWORK(KEL1(IATP)),iWORK(KEL3(IATP)),
     &       iWORK(KEL1(IBTP)),iWORK(KEL3(IBTP)),
     &       iWORK(KLIOIO),IPRDIA)
*
        IF(ICISTR.LE.1) THEN
          LLUDIA = 0
        ELSE
          LLUDIA = LUDIA
        END IF
        CALL CIDIA4(NAEL,Str(IATP)%OCSTR,
     &              NBEL,Str(IBTP)%OCSTR,
     &       NACOB,DIAG,NSMST,WORK(KLH1D),
     &       ISMOST(1,ISM(IISPC)),iWORK(KLBLTP),
     &       WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &       WORK(KLK),Str(IATP)%NSTSO,Str(IBTP)%NSTSO,
     &       iWORK(KLIOIO),NOCTPA,NOCTPB,Str(IATP)%ISTSO,
     &       Str(IBTP)%ISTSO,LLUDIA,ECORE,
     &       PLSIGN,PSSIGN,IPRDIA,NTOOB,ICISTR)
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
      Call GetMem('KLIOIO','FREE','INTE',KLIOIO,NOCTPA*NOCTPB)
      Call GetMem('KLJ   ','FREE','REAL',KLJ   ,NTOOB**2)
      Call GetMem('KLK   ','FREE','REAL',KLK   ,NTOOB**2)
      Call GetMem('KLSC2 ','FREE','REAL',KLSCR2,2*NTOOB**2)
      Call GetMem('KLXA  ','FREE','REAL',KLXA  ,NACOB)
      Call GetMem('KLXB  ','FREE','REAL',KLXB  ,NACOB)
      Call GetMem('KLSCR ','FREE','REAL',KLSCR ,2*NACOB)
      Call GetMem('KLH1D ','FREE','REAL',KLH1D ,NACOB)
      Call GetMem('KLSMO2','FREE','INTE',KLBLTP,NSMST)

      if(doDMRG)then  ! yma
         call dmrg_dim_change_mclr(LRras2(1:8),ntoob,0)
         call dmrg_dim_change_mclr(LRras2(1:8),nacob,0)
      end if

      RETURN
      END
