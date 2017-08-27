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
      function ranf(idum)
      integer idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      real*8 ranf,AM,EPS,RNMX
      parameter(IM1=2147483563,IM2=2147483399,AM=1./dble(IM1),
     *IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,
     *IR2=3791,NTAB=32,NDIV=1+int(dble(IMM1)/NTAB),EPS=1.2e-7,
     *RNMX=1.-EPS)
      integer idum2,j,k,iv(NTAB),iy
      save iv,iy,idum2
      data idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
           k=idum/IQ1
           idum=IA1*(idum-k*IQ1)-k*IR1
           if(idum.lt.0) idum=idum+IM1
           if(j.le.NTAB)iv(j)=idum
11      continue
        iy = iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum = idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ranf=min(AM*iy,RNMX)
      return
      end

*      function ranf(iSeed)
*      Implicit Real*4 (a-h,o-z)
*      Parameter (scale=0.5e0**31)
*      Integer fact
*      Data fact,mask/x'00003ED7',x'7FFFFFFF'/
*      iSeed=iAnd(mask,iSeed*fact)
*      Ranf=scale*iSeed
*      Return
*      End
