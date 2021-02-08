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
* Copyright (C) 1999, Anders Bernhardsson                              *
************************************************************************
      Subroutine Dens_IF(i_root,CA,CB,OCCA,OCCB)
*
*     A small stupid interface for creating the alpha and beta
*     occupation numbers and corresponding molecular orbitals,
*     used in MOLDEN.
*
*     EAW 990118
*
      Implicit Real*8 (a-h,o-z)
#include "SysDef.fh"
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
#include "casvb.fh"
#include "WrkSpc.fh"
      Real*8 CA(*),CB(*),OCCA(*),OCCB(*)
*
      Call GetMem('DS','ALLO','REAL',ipDS,NACPAR)
      Call GetMem('DT','ALLO','REAL',ipDT,NACPAR)
      Call GetMem('DA','ALLO','REAL',ipDA,NACPAR)
      Call GetMem('DB','ALLO','REAL',ipDB,NACPAR)
      Call GetMem('CMO','ALLO','REAL',ipC,NTOT2)
      Call GetMem('UNITY','ALLO','REAL',ipUnity,NAC**2)
*
*     READ IN ORBITALS
* Averaged...
      If ( iOrbTyp.ne.2 ) iad15=IADR15(2)
* Canonical...
      If ( iOrbTyp.eq.2 ) iad15=IADR15(9)

* Read-in orbitals from Jobiph following instructions from previous lines...
      Call DDAFile(JOBIPH,2,Work(ipC),NTOT2,IAD15)

*     COEN WANTED IT AS A BLOCKED MATRIX, SO HERE THEY COME...
      ip1=ipC
      ip2=1
      call dcopy_(nTOT**2,[0.0d0],0,CA,1)
      call dcopy_(nTOT**2,[0.0d0],0,CB,1)
      Do iS=1,nSym
         Do i=1,nbas(is)
            call dcopy_(nbas(is),Work(ip1),1,CA(ip2),1)
            call dcopy_(nbas(is),Work(ip1),1,CB(ip2),1)
            ip1=ip1+nbas(is)
            ip2=ip2+NTOT
         End Do
         ip2=ip2+nbas(is)
      End DO

*     --- FOLLOWING FOR STANDARD (PSEUDO/AVERAGED/CANONICAL) ORBITALS ---
      IF (i_root.eq.0) THEN
* Averaged...
       If ( iOrbTyp.ne.2 ) then
*        SIMPLY READ OCC NOS AS ALPHAS AND ZERO BETA OCC NOS
         iad15=IADR15(2)
         Call DDAFile(JOBIPH,0,[0.0d0],NTOT2,IAD15)
         CALL DDAFILE(JOBIPH,2,OCCA,NTOT,IAD15)
         call dcopy_(nTOT,[0.0d0],0,OCCB,1)
       end if
* Canonical...
       If ( iOrbTyp.eq.2 ) then
*        SIMPLY ZERO ALPHA and BETA OCC NOS
         call dcopy_(nTOT,[0.0d0],0,OCCA,1)
         call dcopy_(nTOT,[0.0d0],0,OCCB,1)
       end if
*
         IF (IFVB.NE.0)THEN
            Call GetMem('VB','ALLO','REAL',ipVB,NAC*NAC)
            call getvb2mo_cvb(Work(ipVB))
            Call GetMem('ACTMO1','ALLO','REAL',ipAM1,NTOT*NAC)
            Call GetMem('ACTMO2','ALLO','REAL',ipAM2,NTOT*NAC)
c           Gather active MOs ...
c           Also count no of active electrons ...
            iAC=0
            IMO=1
            IOCC=1
            OCCNO=0d0
            nAct=0
            Do iS=1,nSym
               call dcopy_(nTOT*nash(is),
     &                    CA((NISH(iS)+NFRO(IS))*NTOT+IMO),1,
     &                    Work(IAC+ipAM1),1)
               Do J=0,NASH(IS)-1
                  OCCNO=OCCNO+OCCA(J+NISH(IS)+NFRO(IS)+IOCC)
               End Do
               nAct=nAct+NASH(iS)
               IAC=IAC+NASH(iS)*NTOT
               IMO=IMO+nBas(is)*ntot
               IOCC=IOCC+NBAS(iS)
            End Do

            Call DGEMM_('N','N',NTOT,NAC,NAC,
     &                 1d0,Work(ipAM1),NTOT,Work(ipVB),NAC,
     &                 0d0,Work(ipAM2),NTOT)
*
*           Scatter active MOs ...
*           Also reset active occ nos - we choose nel/nact since that
*           gives as much meaning as anything ...
*
            iAC=0
            IMO=1
            IOCC=1
            Do iS=1,nSym
               call dcopy_(nTOT*nash(is),Work(IAC+ipAM2),1,
     &                    CA((NISH(iS)+NFRO(IS))*NTOT+IMO),1)
               call dcopy_(nash(is),[OCCNO/DBLE(nAct)],0,
     &                    OCCA(NISH(iS)+NFRO(IS)+IOCC),1)
               IAC=IAC+NASH(iS)*NTOT
               IMO=IMO+nBas(is)*ntot
               IOCC=IOCC+NBAS(iS)
            End Do
            Call GetMem('ACTMO1','FREE','REAL',ipAM1,NTOT*NAC)
            Call GetMem('ACTMO2','FREE','REAL',ipAM2,NTOT*NAC)
            Call GetMem('VB','FREE','REAL',ipVB,NAC*NAC)
         End If
*
      ELSE
*                                                                      *
************************************************************************
*                                                                      *
*     --- FOLLOWING FOR NATURAL *SPIN* ORBITALS ---
*
*     READ IN DENSITIES
*
      iad15=IADR15(3)
      Do i=1,i_root
        Call DDAFile(JOBIPH,2,Work(ipDS),NACPAR,IAD15)
        Call DDAFile(JOBIPH,2,Work(ipDT),NACPAR,IAD15)
        Call DDAFile(JOBIPH,0,[0.0d0],NACPR2,IAD15)
        Call DDAFile(JOBIPH,0,[0.0d0],NACPR2,IAD15)
      End Do
*
*     CREATE SPIN DENSITIES
*
      Do i=0,NACPAR-1
       Work(ipDA+i)=0.5d0*(Work(ipDS+i)+Work(ipDT+i))
       Work(ipDB+i)=0.5d0*(Work(ipDS+i)-Work(ipDT+i))
      End Do

*
*     DIAGONALIZE THE SPIN DENSITIES
*
*      FIRST ALPHA
*
      call dcopy_(nac**2,[0.0d0],0,Work(ipUnity),1)
      call dcopy_(nac,[1.0d0],0,Work(ipUnity),nac+1)
      Call Jacob(Work(ipDA),Work(ipUNITY),NAC,NAC)
*
* TRANSFORM THE ACTIVE ORBITALS
*
      iAC=0
      iAC2=1
      ip=ipC
      Do iS=1,nSym
         If (nBas(iS).ne.0) then
            ip=ip+nBas(iS)*nIsh(iS)
            iAC2=iAC2+nish(is)*nTot
            If (NASH(IS).gt.0)
     &      Call DGEMM_('N','N',NBAS(IS),NASH(IS),NASH(IS),
     &                  1.0d0,Work(ip),NBAS(IS),
     &                  WORK(IPUNITY+IAC),NAC,
     &                  0.0d0,CA(iAC2),NTOT)
            iAC=iAC+NASH(IS)*NAC+NASH(IS)
            iAC2=iAC2+(nbas(is)-nish(is))*NTOT+nbas(is)
            ip=ip+nbas(is)*(nbas(is)-nish(is))
         End If
      End Do
*
*     COPY OCCUPATION NUMBERS
*
      ip=1
      i=0
      ii=0
      Do iS=1,nSYM
         call dcopy_(nBAS(is),[0.0d0],0,OCCA(ip),1)
         call dcopy_((nFro(is)+nish(is)),[1.0d0],0,OCCA(ip),1)
         ip=ip+nFro(is)+nish(is)
         Do iA=1,nash(is)
            ii=ii+1
            i=i+ii
            OCCA(ip)=Work(ipDA+i-1)
            ip=ip+1
         End Do
         ip=ip+nbas(is)-nFro(is)-nish(is)-nash(is)
      End Do
*
*      THEN ONCE AGAIN FOR BETA....
*
      call dcopy_(nac**2,[0.0d0],0,Work(ipUnity),1)
      call dcopy_(nac,[1.0d0],0,Work(ipUnity),nac+1)
*
      Call Jacob(Work(ipDB),Work(ipUNITY),NAC,NAC)
*
      iAC=0
      iAC2=1
      ip=ipC
      Do iS=1,nSym
      if (nbas(is).ne.0) then
       ip=ip+nBas(is)*nIsh(is)
       iAC2=iAC2+nish(is)*NTOT
       If (NASH(IS).gt.0)
     & Call DGEMM_('N','N',NBAS(IS),NASH(IS),NASH(IS),
     &              1.0d0,Work(ip),NBAS(IS),WORK(IPUNITY+IAC),NAC,
     &              0.0d0,CB(iac2),NTOT)
       iAC=iAC+NASH(IS)*NAC+NASH(IS)
       iAC2=iAC2+(nbas(is)-nish(is))*NTOT+nbas(is)
       ip=ip+nbas(is)*(nbas(is)-nish(is))
      endif
      END DO
*
      ip=1
      i=0
      ii=0
      Do iS=1,nSYM
       call dcopy_(nBAS(is),[0.0d0],0,OCCB(ip),1)
       call dcopy_((nFro(is)+nish(is)),[1.0d0],0,OCCB(ip),1)
       ip=ip+nFro(is)+nish(is)
       Do iA=1,nash(is)
        ii=ii+1
        i=i+ii
        OCCB(ip)=Work(ipDB+i-1)
        ip=ip+1
       End Do
       ip=ip+nbas(is)-nFro(is)-nish(is)-nash(is)
      End Do
      ENDIF
*
*     OK, CLEAN UP
*
      Call GetMem('DS','FREE','REAL',ipDS,NACPAR)
      Call GetMem('DT','FREE','REAL',ipDT,NACPAR)
      Call GetMem('DA','FREE','REAL',ipDA,NACPAR)
      Call GetMem('DB','FREE','REAL',ipDB,NACPAR)
      Call GetMem('CMO','FREE','REAL',ipC,NTOT2)
      Call GetMem('UNITY','FREE','REAL',ipUnity,NAC**2)
      Return
      End
