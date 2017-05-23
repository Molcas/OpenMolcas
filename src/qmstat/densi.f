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
* Here we construct the density matrix given the orbital
* coefficients.
      SUBROUTINE DENSI_MO(DENS,ORBCO,IS,IA,NBAS,IDIM)
      IMPLICIT Real*8 (A-H,O-Z)
      DIMENSION DENS(*),ORBCO(IDIM,*)
      IJ=0
      DO 8 I=1,NBAS
        DO 9 J=1,I
          IJ=IJ+1
          DENS(IJ)=0.0d0
   9    CONTINUE
   8  CONTINUE
      DO 10 I=IS,IS+IA-1
        IJ=0
        DO 11 J=1,NBAS
          DO 12 K=1,J
            IJ=IJ+1
            DENS(IJ)=DENS(IJ)+4.d0*ORBCO(J,I)*ORBCO(K,I)
 12       CONTINUE
          DENS(IJ)=DENS(IJ)-ORBCO(J,I)*ORBCO(J,I)*2.d0
 11     CONTINUE
  10  CONTINUE
      RETURN
      END


* The RASSI-density matrix subroutine.
      Subroutine DensiSt(Dens,StVec,iS,nSt,iDim)
      Implicit Real*8 (a-h,o-z)
      Dimension Dens(*),StVec(iDim,*)

* iS        -        Which state that is occupied.
* Dens        -        The density
* StVec        -        The coefficients for how the new states are expressed
*                with the old.
      kaunt=0
      Do 101, i=1,nSt
        Do 102, j=1,i
          kaunt=kaunt+1
          Dens(kaunt)=0.0d0
102     Continue
101   Continue
      kaunt=0
      Do 112, ii=1,nSt
        Do 113, jj=1,ii
          kaunt=kaunt+1
          If(ii.eq.jj) then
            Dens(kaunt)=1.0d0*StVec(ii,iS)*StVec(jj,iS)
          Else
            Dens(kaunt)=2.0d0*StVec(ii,iS)*StVec(jj,iS)
          Endif
113     Continue
112   Continue
      Return
      End

* MP2 density correction.
      Subroutine DCorrCorr(Dens,DenCorr,Trace_Diff,iOrb,iOcc)
      Implicit Real*8 (a-h,o-z)
      Dimension Dens(*),DenCorr(*)
      Trace_HF=dble(iOcc*2)
      kaunt=0
      T=Trace_HF/(Trace_HF-Trace_Diff)
      Do 183, i=1,iOrb
        Do 184, j=1,i
          kaunt=kaunt+1
          Dens(kaunt)=T*(Dens(kaunt)-DenCorr(kaunt))
184     Continue
183   Continue
*      Trace=0.0d0
*      kaunt=0
*      Do 181, i=1,iOrb
*        Do 182, j=1,i
*          kaunt=kaunt+1
*          If(i.eq.j)Trace=Trace+Dens(kaunt)
*182     Continue
*181   Continue
*      call triprt('KKK',' ',Dens,iorb)
*      write(6,*)'QQQ:',Trace
      Return
      End
