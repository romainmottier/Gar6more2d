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

SUBROUTINE sub_transmitpsi_ap
  USE m_phys  
  USE m_num
  USE m_source
  USE m_sismo
  USE m_const
  IMPLICIT NONE 
#include "m_interface.h90"
  INTEGER I,J,K,Nmin,Nmax,N1,N2
  REAL*8 t0,x,ddt,src,t,tau,F,thead,tau0
  COMPLEX*16::pp,Coeff(4),kpsi2,pp0,green
  WRITE(6,*) 'Computation of the transmitted psi wave'
  F=1;
  vmax=MAX(Vs2,Vf2)
  DO I=1,Nx
     x=X1+(I-1)*dx
     t0=arrivaltime(h,x,y,V1,Vpsi2)
     pp0=path(x,y,h,V1,Vpsi2,t0)     
       IF (ABS(AIMAG(pp0)).GT.1/Vmax) THEN
          !Computation of the head wave'
          thead=ABS(y)*SQRT(1/Vpsi2**2-1/vmax**2)+h*SQRT(1/V1**2-1/vmax**2)+ABS(x)/vmax;
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
       !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,ddt,tau0,k,tau,src,pp,kpsi2,Coeff,green)          
          DO J=Nmin,MIN(N2,N1)
             t=T1+(j-1)*dt
             ddt=(t-thead)/Nint
             tau0=thead
             DO  k=1,Nint
                tau=tau0+(k-1)*ddt+ddt/2.
                src=source(t-tau)
                pp=path(-ABS(x),y,h,V1,Vpsi2,tau)  
                pp=ima*AIMAG(pp)   
                kpsi2=zsqrt(pp**2+1/Vpsi2**2)
                Coeff=calccoeff_ap(pp)
                green=derivee(-ABS(x),y,h,V1,Vpsi2,pp)/pi*ddt
                Ux(i,j)=Ux(i,j)+F*REAL(kpsi2*Coeff(4)*green)*src
                Uy(i,j)=Uy(i,j)+F*REAL(ima*pp*Coeff(4)*green)*src
             END DO
          END DO
          !$OMP END PARALLEL DO 
          !$OMP TASKWAIT
       !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,ddt,tau0,k,tau,src,pp,kpsi2,Coeff,green)
          DO J=MAX(N2+1,Nmin),N1
             t=T1+(j-1)*dt
             ddt=2.*tdelay/Nint
             tau0=t-2.*tdelay
             DO  k=1,Nint
                tau=tau0+(k-1)*ddt+ddt/2.
                src=source(t-tau)
                pp=path(-ABS(x),y,h,V1,Vpsi2,tau)  
                pp=ima*AIMAG(pp)   
                kpsi2=zsqrt(pp**2+1/Vpsi2**2)
                Coeff=calccoeff_ap(pp)
                green=derivee(-ABS(x),y,h,V1,Vpsi2,pp)/pi*ddt
                Ux(i,j)=Ux(i,j)+F*REAL(kpsi2*Coeff(4)*green)*src
                Uy(i,j)=Uy(i,j)+F*REAL(ima*pp*Coeff(4)*green)*src
             END DO
          END DO
          !$OMP END PARALLEL DO 
          !$OMP TASKWAIT
        !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,ddt,tau0,k,tau,src,pp,kpsi2,Coeff,green)
          DO J=N1+1,MIN(N2,Nmax)
             t=T1+(j-1)*dt
             ddt=(t0-thead)/Nint
             tau0=thead
             DO  k=1,Nint-1
                tau=tau0+(k-1)*ddt+ddt/2.
                src=source(t-tau)
                pp=path(-ABS(x),y,h,V1,Vpsi2,tau)  
                pp=ima*AIMAG(pp)   
                kpsi2=zsqrt(pp**2+1/Vpsi2**2)
                Coeff=calccoeff_ap(pp)
                green=derivee(-ABS(x),y,h,V1,Vpsi2,pp)/pi*ddt
                Ux(i,j)=Ux(i,j)+F*REAL(kpsi2*Coeff(4)*green)*src
                Uy(i,j)=Uy(i,j)+F*REAL(ima*pp*Coeff(4)*green)*src
             END DO
             tau=t0-ddt/2.
             src=source(t-tau)
             pp=path(-ABS(x),y,h,V1,Vpsi2,tau)   
                pp=ima*AIMAG(pp)  
             kpsi2=zsqrt(pp**2+1/Vpsi2**2)
             Coeff=calccoeff_ap(pp)
             green=SQRT(t0**2-tau**2)*derivee(-ABS(x),y,h,V1,Vpsi2,pp)*ACOS(1-ddt/t0)/pi*ddt
             Ux(i,j)=Ux(i,j)+F*REAL(kpsi2*Coeff(4)*green)*src
             Uy(i,j)=Uy(i,j)+F*REAL(ima*pp*Coeff(4)*green)*src   
          END DO
          !$OMP END PARALLEL DO 
          !$OMP TASKWAIT

        !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,ddt,tau0,k,tau,src,pp,kpsi2,Coeff,green)
          DO J=MAX(N1+1,N2+1),Nmax
            t=T1+(j-1)*dt
            ddt=(t0-t+2.*tdelay)/Nint
            tau0=t-2.*tdelay
            DO  k=1,Nint-1
               tau=tau0+(k-1)*ddt+ddt/2.
               src=source(t-tau)
               pp=path(-ABS(x),y,h,V1,Vpsi2,tau)  
               pp=ima*AIMAG(pp)   
               kpsi2=zsqrt(pp**2+1/Vpsi2**2)
               Coeff=calccoeff_ap(pp)
               green=derivee(-ABS(x),y,h,V1,Vpsi2,pp)/pi*ddt
               Ux(i,j)=Ux(i,j)+F*REAL(kpsi2*Coeff(4)*green)*src
                Uy(i,j)=Uy(i,j)+F*REAL(ima*pp*Coeff(4)*green)*src
             END DO
             tau=t0-ddt/2.
             src=source(t-tau)
             pp=path(-ABS(x),y,h,V1,Vpsi2,tau)   
             pp=ima*AIMAG(pp)  
             kpsi2=zsqrt(pp**2+1/Vpsi2**2)
             Coeff=calccoeff_ap(pp)
             green=SQRT(t0**2-tau**2)*derivee(-ABS(x),y,h,V1,Vpsi2,pp)*ACOS(1-ddt/t0)/pi*ddt
             Ux(i,j)=Ux(i,j)+F*REAL(kpsi2*Coeff(4)*green)*src
             Uy(i,j)=Uy(i,j)+F*REAL(ima*pp*Coeff(4)*green)*src   
          END DO
          !$OMP END PARALLEL DO 
          !$OMP TASKWAIT
 END IF
!! Computation of the volume wave
    Nmin=CEILING((t0-T1)/dt+1)
    Nmax=FLOOR((t0+2*tdelay-T1)/dt+1) 
    Nmin=MAX(1,Nmin)
    Nmax=MIN(Nmax,Nt)
    !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,ddt,tau0,k,tau,src,pp,kpsi2,Coeff,green)
    DO J=Nmin,Nmax
       t=T1+(j-1)*dt
          ddt=(t-t0)/Nint
          tau0=t0
          tau=tau0+ddt/2.
          src=source(t-tau)
          pp=path(-ABS(x),y,h,V1,Vpsi2,tau)   
          kpsi2=zsqrt(pp**2+1/Vpsi2**2)
          Coeff=calccoeff_ap(pp)
          green=SQRT(tau**2-t0**2)*derivee(-ABS(x),y,h,V1,Vpsi2,pp)*dacosh(1+ddt/t0)/pi
          Ux(i,j)=Ux(i,j)+F*REAL(kpsi2*Coeff(4)*green)*src
          Uy(i,j)=Uy(i,j)+F*REAL(ima*pp*Coeff(4)*green)*src
          DO k=2,Nint
             tau=tau0+(k-1)*ddt+ddt/2.
             src=source(t-tau)
             pp=path(-ABS(x),y,h,V1,Vpsi2,tau)   
             kpsi2=zsqrt(pp**2+1/Vpsi2**2)
             Coeff=calccoeff_ap(pp)
             green=derivee(-ABS(x),y,h,V1,Vpsi2,pp)/pi*ddt
             Ux(i,j)=Ux(i,j)+F*REAL(kpsi2*Coeff(4)*green)*src
             Uy(i,j)=Uy(i,j)+F*REAL(ima*pp*Coeff(4)*green)*src
          END DO
       END DO
       !$OMP END PARALLEL DO 
       !$OMP TASKWAIT

       !$OMP PARALLEL DO     DEFAULT(SHARED) PRIVATE(J,t,ddt,tau0,k,tau,src,pp,kpsi2,Coeff,green)
       DO J=Nmax+1,Nt
          t=T1+(j-1)*dt
          ddt=2.*tdelay/Nint
          tau0=t-2.*tdelay
          DO k=1,Nint
             tau=tau0+(k-1)*ddt+ddt/2.
             src=source(t-tau)
             pp=path(-ABS(x),y,h,V1,Vpsi2,tau)   
             kpsi2=zsqrt(pp**2+1/Vpsi2**2)
             Coeff=calccoeff_ap(pp)
             green=derivee(-ABS(x),y,h,V1,Vpsi2,pp)/pi*ddt
             Ux(i,j)=Ux(i,j)+F*REAL(kpsi2*Coeff(4)*green)*src
             Uy(i,j)=Uy(i,j)+F*REAL(ima*pp*Coeff(4)*green)*src
          END DO
       END DO
       !$OMP END PARALLEL DO 
       !$OMP TASKWAIT
    END DO
END SUBROUTINE sub_transmitpsi_ap
