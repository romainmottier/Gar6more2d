!! Copyright INRIA. Contributors : Julien DIAZ and Abdelaaziz EZZIANI
!! 
!! Julien.Diaz@inria.fr and Abdelaaziz.Ezziani@univ-pau.fr
!! 
!! This software is a computer program whose purpose is to
!! compute the analytical solution of problems of waves propagation in two 
!! layered media such as
!! - acoustic/acoustic
!! - acoustic/elastodynamic
!! - acoustic/porous
!! - porous/porous,
!! based on the Cagniard-de Hoop method.
!! 
!! This software is governed by the CeCILL license under French law and
!! abiding by the rules of distribution of free software.  You can  use, 
!! modify and/ or redistribute the software under the terms of the CeCILL
!! license as circulated by CEA, CNRS and INRIA at the following URL
!! "http://www.cecill.info". 
!! 
!! As a counterpart to the access to the source code and  rights to copy,
!! modify and redistribute granted by the license, users are provided only
!! with a limited warranty  and the software's author,  the holder of the
!! economic rights,  and the successive licensors  have only  limited
!! liability. 
!! 
!! In this respect, the user's attention is drawn to the risks associated
!! with loading,  using,  modifying and/or developing or reproducing the
!! software by the user in light of its specific status of free software,
!! that may mean  that it is complicated to manipulate,  and  that  also
!! therefore means  that it is reserved for developers  and  experienced
!! professionals having in-depth computer knowledge. Users are therefore
!! encouraged to load and test the software's suitability as regards their
!! requirements in conditions enabling the security of their systems and/or 
!! data to be ensured and,  more generally, to use and operate it in the 
!! same conditions as regards security. 
!! 
!! The fact that you are presently reading this means that you have had
!! knowledge of the CeCILL license and that you accept its terms.
!! ========================================================================

SUBROUTINE sub_reflex_ae
  USE m_phys  
  USE m_num
  USE m_source
  USE m_sismo
  USE m_const
  IMPLICIT NONE 
#include "m_interface.h90"
  INTEGER I,J,K,Nmin,Nmax,N1,N2
  REAL*8 t0,r,x,ddt,src,t,tau,F,thead,green,dq,q,q0
  COMPLEX*16::pp,Coeff,k1, kp,ks
  WRITE(6,*) 'Computation of the reflected wave'
  F=1/V1**2
  vmax=MAX(Vp2,Vs2)
  DO I=1,Nx
     x=X1+(I-1)*dx
     r=dsqrt(x**2+(y+h)**2)
     t0=r/V1
     IF (ABS(x/r).GT.V1/Vmax) THEN
        !Computation of the head wave'
        thead=ABS(y+h)*SQRT(1/V1**2-1/Vmax**2)+ABS(x)/Vmax;
        Nmin=CEILING((thead-T1)/dt+1)
        Nmax=FLOOR((t0+2*tdelay-T1)/dt+1)
        N1=FLOOR((t0-T1)/dt+1)
        N2=FLOOR((thead+2*tdelay-T1)/dt+1)
        Nmin=MAX(Nmin,1)
        Nmax=MIN(Nmax,Nt)
        N1=MAX(N1,1)
        N2=MAX(N2,1)
        N1=MIN(N1,Nt)
        N2=MIN(N2,Nt)
        !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,dq,q0,k,tau,src,pp,k1,kp,ks,Coeff,green,q)
        DO J=Nmin,MIN(N2,N1)
           t=T1+(j-1)*dt
           dq=ASIN(t/t0)/Nint
           q0=0
           DO  k=1,Nint
              q=q0+(k-1)*dq+dq/2.
              tau=t0*SIN(q)
              src=source(t-tau)
              pp=-ima*(ABS(y+h)/r**2*t0*COS(q)-ABS(x)*tau/r**2)
              k1=zsqrt(pp**2+1/V1**2)
              IF(AIMAG(pp)**2.GT.1/Vp2**2) THEN
                 kp=ima*zsqrt(-pp**2-1/Vp2**2)
              ELSE
                 kp=zsqrt(pp**2+1/Vp2**2)
              ENDIF
              IF(AIMAG(pp)**2.GT.1/Vs2**2) THEN
                 ks=ima*zsqrt(-pp**2-1/Vs2**2)
              ELSE
                 ks=zsqrt(pp**2+1/Vs2**2)
              ENDIF
              Coeff=(k1*((lambda2*(kp**2-pp**2)+2*mu2*kp**2)*(ks**2+pp**2)-4*mu2*pp**2*ks*kp)-rho1*kp*(ks**2-pp**2))
              Coeff=Coeff/(k1*((lambda2*(kp**2-pp**2)+2*mu2*kp**2)*(ks**2+pp**2)-4*mu2*pp**2*ks*kp)+rho1*kp*(ks**2-pp**2))
              Coeff=Coeff/(2.*pi)*dq
              Ux(i,j)=Ux(i,j)-F*AIMAG(ima*pp*Coeff)*src/rho1
              Uy(i,j)=Uy(i,j)-F*AIMAG(k1*Coeff)*src/rho1
              P(i,j)=P(i,j)-F*AIMAG(Coeff)*src
           END DO
        END DO
        !$OMP END PARALLEL DO 
        !$OMP TASKWAIT
        !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,dq,q0,k,tau,src,pp,k1,kp,ks,Coeff,green,q)
        DO J=MAX(N2+1,Nmin),N1
           t=T1+(j-1)*dt
           dq=(ASIN(t/t0)-ASIN((t-2.*tdelay)/t0))/Nint
           q0=ASIN((t-2.*tdelay)/t0)

           DO  k=1,Nint
              q=q0+(k-1)*dq+dq/2.
              tau=t0*SIN(q)
              src=source(t-tau)
              pp=-ima*(ABS(y+h)/r**2*t0*COS(q)-ABS(x)*tau/r**2)
              k1=zsqrt(pp**2+1/V1**2)
              IF(AIMAG(pp)**2.GT.1/Vp2**2) THEN
                 kp=ima*zsqrt(-pp**2-1/Vp2**2)
              ELSE
                 kp=zsqrt(pp**2+1/Vp2**2)
              ENDIF
              IF(AIMAG(pp)**2.GT.1/Vs2**2) THEN
                 ks=ima*zsqrt(-pp**2-1/Vs2**2)
              ELSE
                 ks=zsqrt(pp**2+1/Vs2**2)
              ENDIF
              Coeff=(k1*((lambda2*(kp**2-pp**2)+2*mu2*kp**2)*(ks**2+pp**2)-4*mu2*pp**2*ks*kp)-rho1*kp*(ks**2-pp**2))
              Coeff=Coeff/(k1*((lambda2*(kp**2-pp**2)+2*mu2*kp**2)*(ks**2+pp**2)-4*mu2*pp**2*ks*kp)+rho1*kp*(ks**2-pp**2))
              Coeff=Coeff/(2.*pi)*dq
              Ux(i,j)=Ux(i,j)-F*AIMAG(ima*pp*Coeff)*src/rho1
              Uy(i,j)=Uy(i,j)-F*AIMAG(k1*Coeff)*src/rho1
              P(i,j)=P(i,j)-F*AIMAG(Coeff)*src
           END DO
        END DO
        !$OMP END PARALLEL DO 
        !$OMP TASKWAIT
        !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,dq,q0,k,tau,src,pp,k1,kp,ks,Coeff,green,q)
        DO J=N1+1,MIN(N2,Nmax)
           t=T1+(j-1)*dt
           dq=pi/2./Nint
           q0=0
           DO  k=1,Nint
              q=q0+(k-1)*dq+dq/2.
              tau=t0*SIN(q)
              src=source(t-tau)
              pp=-ima*(ABS(y+h)/r**2*t0*COS(q)-ABS(x)*tau/r**2)
              k1=zsqrt(pp**2+1/V1**2)
              IF(AIMAG(pp)**2.GT.1/Vp2**2) THEN
                 kp=ima*zsqrt(-pp**2-1/Vp2**2)
              ELSE
                 kp=zsqrt(pp**2+1/Vp2**2)
              ENDIF
              IF(AIMAG(pp)**2.GT.1/Vs2**2) THEN
                 ks=ima*zsqrt(-pp**2-1/Vs2**2)
              ELSE
                 ks=zsqrt(pp**2+1/Vs2**2)
              ENDIF
              Coeff=(k1*((lambda2*(kp**2-pp**2)+2*mu2*kp**2)*(ks**2+pp**2)-4*mu2*pp**2*ks*kp)-rho1*kp*(ks**2-pp**2))
              Coeff=Coeff/(k1*((lambda2*(kp**2-pp**2)+2*mu2*kp**2)*(ks**2+pp**2)-4*mu2*pp**2*ks*kp)+rho1*kp*(ks**2-pp**2))
              Coeff=Coeff/(2.*pi)*dq
              Ux(i,j)=Ux(i,j)-F*AIMAG(ima*pp*Coeff)*src/rho1
              Uy(i,j)=Uy(i,j)-F*AIMAG(k1*Coeff)*src/rho1
              P(i,j)=P(i,j)-F*AIMAG(Coeff)*src
           END DO
        END DO
        !$OMP END PARALLEL DO 
        !$OMP TASKWAIT
        !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,dq,q0,k,tau,src,pp,k1,kp,ks,Coeff,green,q)
        DO J=MAX(N1+1,N2+1),Nmax
           t=T1+(j-1)*dt
           dq=(pi/2.-ASIN((t-2.*tdelay)/t0))/Nint
           q0=ASIN((t-2.*tdelay)/t0)
           DO  k=1,Nint
              q=q0+(k-1)*dq+dq/2.
              tau=t0*SIN(q)
              src=source(t-tau)
              pp=-ima*(ABS(y+h)/r**2*t0*COS(q)-ABS(x)*tau/r**2)
              k1=zsqrt(pp**2+1/V1**2)
              IF(AIMAG(pp)**2.GT.1/Vp2**2) THEN
                 kp=ima*zsqrt(-pp**2-1/Vp2**2)
              ELSE
                 kp=zsqrt(pp**2+1/Vp2**2)
              ENDIF
              IF(AIMAG(pp)**2.GT.1/Vs2**2) THEN
                 ks=ima*zsqrt(-pp**2-1/Vs2**2)
              ELSE
                 ks=zsqrt(pp**2+1/Vs2**2)
              ENDIF
              Coeff=(k1*((lambda2*(kp**2-pp**2)+2*mu2*kp**2)*(ks**2+pp**2)-4*mu2*pp**2*ks*kp)-rho1*kp*(ks**2-pp**2))
              Coeff=Coeff/(k1*((lambda2*(kp**2-pp**2)+2*mu2*kp**2)*(ks**2+pp**2)-4*mu2*pp**2*ks*kp)+rho1*kp*(ks**2-pp**2))
              Coeff=Coeff/(2.*pi)*dq
              Ux(i,j)=Ux(i,j)-F*AIMAG(ima*pp*Coeff)*src/rho1
              Uy(i,j)=Uy(i,j)-F*AIMAG(k1*Coeff)*src/rho1
              P(i,j)=P(i,j)-F*AIMAG(Coeff)*src
           END DO
        END DO
        !$OMP END PARALLEL DO 
        !$OMP TASKWAIT
     END IF
     !! Computation of the volume wave
     Nmin=CEILING((t0-T1)/dt+1)
     Nmax=FLOOR((t0+2*tdelay-T1)/dt+1)
     Nmin=MAX(1,Nmin)
     Nmax=MIN(Nmax,Nt)
     !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,dq,q0,k,tau,src,pp,k1,kp,ks,Coeff,green,q)
     DO J=Nmin,Nmax
        t=T1+(j-1)*dt
        dq=acosh(t/t0)/Nint
        q0=0
        DO k=1,Nint
           q=q0+(k-1)*dq+dq/2.
           tau=t0*COSH(q)
           src=source(t-tau)
           pp=ima*tau*ABS(x)/r**2+ABS(y+h)/r**2*t0*SINH(q)
           k1=zsqrt(pp**2+1/V1**2)
           kp=zsqrt(pp**2+1/Vp2**2)
           ks=zsqrt(pp**2+1/Vs2**2)
           Coeff=(k1*((lambda2*(kp**2-pp**2)+2*mu2*kp**2)*(ks**2+pp**2)-4*mu2*pp**2*ks*kp)-rho1*kp*(ks**2-pp**2))
           Coeff=Coeff/(k1*((lambda2*(kp**2-pp**2)+2*mu2*kp**2)*(ks**2+pp**2)-4*mu2*pp**2*ks*kp)+rho1*kp*(ks**2-pp**2))
           Coeff=Coeff/(2*pi)*dq
           Ux(i,j)=Ux(i,j)+F*REAL(ima*pp*Coeff)*src/rho1
           Uy(i,j)=Uy(i,j)+F*REAL(k1*Coeff)*src/rho1
           P(i,j)=P(i,j)+F*REAL(Coeff)*src
        END DO
     END DO
     !$OMP END PARALLEL DO 
     !$OMP TASKWAIT
     !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,dq,q0,k,q,tau,src,pp,k1,kp,ks,Coeff,green)
     DO J=Nmax+1,Nt
        t=T1+(j-1)*dt
        dq=(acosh(t/t0)-acosh((t-2.*tdelay)/t0))/Nint
        q0=acosh((t-2.*tdelay)/t0)
        DO k=1,Nint
           q=q0+(k-1)*dq+dq/2.
           tau=t0*COSH(q)
           src=source(t-tau)
           pp=ima*tau*ABS(x)/r**2+ABS(y+h)/r**2*t0*SINH(q)
           k1=zsqrt(pp**2+1/V1**2)
           kp=zsqrt(pp**2+1/Vp2**2)
           ks=zsqrt(pp**2+1/Vs2**2)
           Coeff=(k1*((lambda2*(kp**2-pp**2)+2*mu2*kp**2)*(ks**2+pp**2)-4*mu2*pp**2*ks*kp)-rho1*kp*(ks**2-pp**2))
           Coeff=Coeff/(k1*((lambda2*(kp**2-pp**2)+2*mu2*kp**2)*(ks**2+pp**2)-4*mu2*pp**2*ks*kp)+rho1*kp*(ks**2-pp**2))
           Coeff=Coeff/(2*pi)*dq
           Ux(i,j)=Ux(i,j)+F*REAL(ima*pp*Coeff)*src/rho1
           Uy(i,j)=Uy(i,j)+F*REAL(k1*Coeff)*src/rho1
           P(i,j)=P(i,j)+F*REAL(Coeff)*src
        END DO
     END DO
     !$OMP END PARALLEL DO 
     !$OMP TASKWAIT
  END DO
END SUBROUTINE sub_reflex_ae
