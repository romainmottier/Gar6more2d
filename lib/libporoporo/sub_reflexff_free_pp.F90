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

SUBROUTINE sub_reflexff_free_pp
  Use m_phys  
  Use m_num
  Use m_source
  Use m_sismo
  Use m_const
  implicit none 
#include "m_interface.h90"
  integer I,J,K,Nmin,Nmax,N1,N2,integ(1000)
  real*8 t0,r,x,ddt,src,t,tau,thead,green,F,dq,q,q0
  complex*16::pp,Coeff(3),kf1,test_complex(1000)
  F=P1(1,1)*f1
  write(6,*) 'Computation of the reflected ff wave'
  vmax=max(max(Vf1,Vpsi1),Vs1)
  DO I=1,Nx
     x=X1+(I-1)*dx
     r=dsqrt(x**2+(y+h)**2)
     t0=r/Vf1
       if (abs(x/r).gt.Vf1/Vmax) then
          !Computation of the head wave'
          thead=abs(y+h)*sqrt(1/Vf1**2-1/Vmax**2)+abs(x)/Vmax;
        Nmin=ceiling((thead-T1)/dt+1)
        Nmax=floor((t0+2*tdelay-T1)/dt+1)
        N1=floor((t0-T1)/dt+1)
        N2=floor((thead+2*tdelay-T1)/dt+1)
        Nmin=max(Nmin,1)
        Nmax=min(Nmax,Nt)
        N1=max(N1,1)
        N2=max(N2,1)
        N1=min(N1,Nt)
        N2=min(N2,Nt)
    !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,dq,q0,q,k,tau,src,pp,kf1,Coeff,green)
        DO J=Nmin,min(N2,N1)
           t=T1+(j-1)*dt
           dq=asin(t/t0)/Nint
           q0=0
             DO  k=1,Nint
                q=q0+(k-1)*dq+dq/2.
                tau=t0*sin(q)
                src=source(t-tau)
                pp=-ima*(abs(y+h)/r**2*t0*cos(q)-abs(x)*tau/r**2)
                kf1=zsqrt(pp**2+1/Vf1**2)
                Coeff=F*calccoefff_free_pp(pp)/pi*dq
                Ux(i,j)=Ux(i,j)+aimag(ima*pp*Coeff(1)*kf1)*src
                Uy(i,j)=Uy(i,j)+aimag(kf1**2*Coeff(1))*src
             END DO
          END DO
          !$OMP END PARALLEL DO 
          !$OMP TASKWAIT
          !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,dq,q0,q,k,tau,src,pp,kf1,Coeff,green)
        DO J=max(N2+1,Nmin),N1
            t=T1+(j-1)*dt
           dq=(asin(t/t0)-asin((t-2.*tdelay)/t0))/Nint
           q0=asin((t-2.*tdelay)/t0)
             DO  k=1,Nint
                q=q0+(k-1)*dq+dq/2.
                tau=t0*sin(q)
                src=source(t-tau)
                pp=-ima*(abs(y+h)/r**2*t0*cos(q)-abs(x)*tau/r**2)
                kf1=zsqrt(pp**2+1/Vf1**2)
                Coeff=F*calccoefff_free_pp(pp)/pi*dq
                Ux(i,j)=Ux(i,j)+aimag(ima*pp*Coeff(1)*kf1)*src
                Uy(i,j)=Uy(i,j)+aimag(kf1**2*Coeff(1))*src
             END DO
            END DO
       !$OMP END PARALLEL DO 
       !$OMP TASKWAIT
    !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,dq,q0,q,k,tau,src,pp,kf1,Coeff,green)
        DO J=N1+1,min(N2,Nmax)
           t=T1+(j-1)*dt
           dq=pi/2./Nint
           q0=0
             DO  k=1,Nint
                q=q0+(k-1)*dq+dq/2.
                tau=t0*sin(q)
                src=source(t-tau)
                pp=-ima*(abs(y+h)/r**2*t0*cos(q)-abs(x)*tau/r**2)
                kf1=zsqrt(pp**2+1/Vf1**2)
                Coeff=F*calccoefff_free_pp(pp)/pi*dq
                Ux(i,j)=Ux(i,j)+aimag(ima*pp*Coeff(1)*kf1)*src
                Uy(i,j)=Uy(i,j)+aimag(kf1**2*Coeff(1))*src
             END DO
          END DO
          !$OMP END PARALLEL DO 
          !$OMP TASKWAIT
    !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,dq,q0,q,k,tau,src,pp,kf1,Coeff,green)
            DO J=max(N1+1,N2+1),Nmax
                t=T1+(j-1)*dt
           dq=(pi/2.-asin((t-2.*tdelay)/t0))/Nint
           q0=asin((t-2.*tdelay)/t0)
             DO  k=1,Nint
                q=q0+(k-1)*dq+dq/2.
                tau=t0*sin(q)
                src=source(t-tau)
                pp=-ima*(abs(y+h)/r**2*t0*cos(q)-abs(x)*tau/r**2)
                kf1=zsqrt(pp**2+1/Vf1**2)
                Coeff=F*calccoefff_free_pp(pp)/pi*dq
                Ux(i,j)=Ux(i,j)+aimag(ima*pp*Coeff(1)*kf1)*src
                Uy(i,j)=Uy(i,j)+aimag(kf1**2*Coeff(1))*src
             END DO
          end DO
             !$OMP END PARALLEL DO 
             !$OMP TASKWAIT
       end if
!! Computation of the volume wave
    Nmin=ceiling((t0-T1)/dt+1)
    Nmax=floor((t0+2*tdelay-T1)/dt+1)
    Nmin=max(1,Nmin)
    Nmax=min(Nmax,Nt)
    !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,dq,q0,q,k,tau,src,pp,kf1,Coeff,green)
    DO J=Nmin,Nmax
       t=T1+(j-1)*dt
             dq=acosh(t/t0)/Nint
             q0=0
          do k=1,Nint
             q=q0+(k-1)*dq+dq/2.
             tau=t0*cosh(q)
             src=source(t-tau)
             pp=ima*tau*abs(x)/r**2+abs(y+h)/r**2*t0*sinh(q)
             kf1=zsqrt(pp**2+1/Vf1**2)
             Coeff=F*calccoefff_free_pp(pp)/pi*dq
             Ux(i,j)=Ux(i,j)-real(ima*pp*kf1*Coeff(1))*src
             Uy(i,j)=Uy(i,j)-real(kf1**2*Coeff(1))*src
          end do
    end DO
       !$OMP END PARALLEL DO 
       !$OMP TASKWAIT
   !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,dq,q0,q,k,tau,src,pp,kf1,Coeff,green)
       DO J=Nmax+1,Nt
          t=T1+(j-1)*dt
             dq=(acosh(t/t0)-acosh((t-2.*tdelay)/t0))/Nint
             q0=acosh((t-2.*tdelay)/t0)
          do k=1,Nint
             q=q0+(k-1)*dq+dq/2.
             tau=t0*cosh(q)
             src=source(t-tau)
             pp=ima*tau*abs(x)/r**2+abs(y+h)/r**2*t0*sinh(q)
             kf1=zsqrt(pp**2+1/Vf1**2)
             Coeff=F*calccoefff_free_pp(pp)/pi*dq
             Ux(i,j)=Ux(i,j)-real(ima*pp*kf1*Coeff(1))*src
             Uy(i,j)=Uy(i,j)-real(kf1**2*Coeff(1))*src
          end do
       end DO
       !$OMP END PARALLEL DO 
       !$OMP TASKWAIT
    END DO
END SUBROUTINE sub_reflexff_free_pp
