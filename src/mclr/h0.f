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
      Subroutine h0(rdia,MP1,MP2,MQ,
     &              isym,KH0S,KH0F,KSBIDT,nprciv,
     &               TimeDep)
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
#include "WrkSpc.fh"
#include "csfbas_mclr.fh"
#include "negpre.fh"
      Real*8 rdia(*)
      Logical TimeDep
*
      MXP1=MP1
      MXP2=MP2
      MXQ=MQ
      ispc=1
      NDET= nint(XISPSM(ISYM,ISPC))
      IF(NOCSF.EQ.1) THEN
        NVAR=nDET
      ELSE
        NVAR = NCSASM(ISYM)
      END  IF
*
      LBLK=NVAR
      NSBDET = MXP1 + MXP2 + MXQ
      MXP = MXP1 + MXP2
      LH0 = MXP*(MXP+1)/2 + MXP1*MXQ
      MXCSFC = 0
      MXDTFC = 0
      DO  ITYP = 1, NTYP
        MXCSFC = MAX(MXCSFC,NCPCNT(ITYP) )
        MXDTFC = MAX(MXDTFC,NDPCNT(ITYP) )
      End Do

      nactel=naelci(1)+nbelci(1)
      If (TimeDep) Then
         EnA=E2_td(Work(ipFiMo),Work(k2int),0,-1)
      Else
         EnA=E2(Work(ipFiMo),Work(k2int),0,-1)
      End If
      LH0SCR =
     &  MAX(6*NSBDET,4*NSBDET+4*NOCOB,MXP1*(MXP1+1)/2+MXP1**2)
      LVEC2=2 * NACTEL + MXCSFC**2
     &               + 6*MXDTFC+2*MXDTFC**2
     &               + MAX(MXDTFC*NACTEL+2*NACTEL,4*NACOB+2*NACTEL)
      LVEC2=Max(lvec2,ndet)
*
      If (isym.eq.irefsm) then
       ieaw=1
      else
       ieaw=2
      end if
      CALL GetMem('KSBIDT','ALLO','INTE',KSBIDT,NSBDET)
      CALL GetMem('EXPH0S','ALLO','REAL',KH0S,MXP**2)
      CALL GetMem('EXPH0F','ALLO','INTE',KH0F,MXP)
      CALL GetMem('EXPH0','ALLO','REAL',KH0,LH0)
      CALL GetMem('KSBCNF','ALLO','INTE',KSBCNF,NSBDET)
      Call GETMEM('KH0SCR','ALLO','REAL',KH0SCR,LH0SCR)
      iKH0SCR=ip_of_iWork(Work(KH0SCR))
      Call GetMem('KVEC2','ALLO','REAL',KVEC2,lvec2)
*
      Call H0MAT_MCLR(Work(KH0),
     &       iWORK(KSBIDT),iWORK(KSBCNF),
     &       MXP1,MXP2,MXQ,NACOB,NPRCIV,
     &       NOCSF,ISYM,IDC,PSSIGN,ECOREP,
     &       rDIA,Work(KVEC2),work(kh0scr),iWork(ikh0scr),ieaw)

*
      do i=1,nprciv
       Work(KH0+i*(i+1)/2-1)=Work(KH0+i*(i+1)/2-1)-ENA
      End Do
      IF (NGP)
     & Call mkp1(nprciv,iWORK(KSBIDT),Work(kh0),rdia)
*     Call Triprt('PRECI',' ',Work(KH0),nprciv)
*     write(*,*) (iWork(KSBIDT+i),i=0,nprciv-1)
      Call GetMem('KVEC2','FREE','REAL',KVEC2,lvec)
      Call GETMEM('KH0SCR','FREE','REAL',KH0SCR,LH0SCR)
      CALL GetMem('KSBCNF','FREE','INTE',KSBCNF,NSBDET)
*     CALL GetMem('KSBIDT','FREE','INTE',KSBIDT,NSBDET)
      call square(Work(KH0),Work(KH0S),1,NPRCIV,NPRCIV)
      CALL GetMem('EXPH0','FREE','REAL',KH0,LH0)
*
cVV: next dummy statement is needed for Intel MKL installation.
      if(KH0S.eq.666) write(6,*) Work(KH0S)
      irc=0
      call dgetrf_(NPRCIV,NPRCIV,Work(KH0S),NPRCIV,
     &            iWork(KH0F),irc)
      If (irc.ne.0) Then
        Write(6,*) 'Sorry but you have an singular ci matrix'
        Write(6,*) 'Set ExpDimension and restart mclr'
        call Abend
      End If
*
*
      Return
      end
