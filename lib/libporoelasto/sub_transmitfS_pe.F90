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

subroutine sub_transmitfS_pe
  use m_phys  
  use m_num
  use m_source
  use m_sismo
  use m_const
  implicit none 
#include "m_interface.h90"
  integer I,J,K,Nmin,Nmax,N1,N2
  real*8 t0,x,ddt,src,t,tau,F,thead,tau0
  complex*16::pp,Coeff(5),kS2,pp0,green
  write(6,*) 'Computation of the transmitted fS wave'
  F=f1
  vmax=max(max(VP2,Vpsi1),Vs1)
  do I=1,Nx
     x=X1+(I-1)*dx
     t0=arrivaltime(h,x,y,Vf1,VS2)
     pp0=path(x,y,h,Vf1,VS2,t0)     
     if (abs(aimag(pp0)).gt.1/Vmax) then
        !Computation of the head wave'
        thead=abs(y)*sqrt(1/VS2**2-1/vmax**2)+h*sqrt(1/Vf1**2-1/vmax**2)+abs(x)/vmax;
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
        !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,ddt,tau0,k,tau,src,pp,kS2,Coeff,green)
        do J=Nmin,min(N2,N1)
           t=T1+(j-1)*dt
           ddt=(t-thead)/Nint
           tau0=thead
           do  k=1,Nint
              tau=tau0+(k-1)*ddt+ddt/2.
              src=source(t-tau)
              pp=path(-abs(x),y,h,Vf1,VS2,tau)     
              pp=ima*aimag(pp)       
              kS2=zsqrt(pp**2+1/VS2**2)
              Coeff=calccoefff_pe(pp)
              green=derivee(-abs(x),y,h,Vf1,VS2,pp)/pi*ddt
              Ux(i,j)=Ux(i,j)+F*real(kS2*coeff(5)*green)*src
              Uy(i,j)=Uy(i,j)+F*real(ima*pp*coeff(5)*green)*src
           end do
        end do
        !$OMP END PARALLEL DO 
        !$OMP TASKWAIT
        !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,ddt,tau0,k,tau,src,pp,kS2,Coeff,green)
        do J=max(N2+1,Nmin),N1
           t=T1+(j-1)*dt
           ddt=2.*tdelay/Nint
           tau0=t-2.*tdelay
           do  k=1,Nint
              tau=tau0+(k-1)*ddt+ddt/2.
              src=source(t-tau)
              pp=path(-abs(x),y,h,Vf1,VS2,tau)     
              pp=ima*aimag(pp)       
              kS2=zsqrt(pp**2+1/VS2**2)
              Coeff=calccoefff_pe(pp)
              green=derivee(-abs(x),y,h,Vf1,VS2,pp)/pi*ddt
              Ux(i,j)=Ux(i,j)+F*real(kS2*coeff(5)*green)*src
              Uy(i,j)=Uy(i,j)+F*real(ima*pp*coeff(5)*green)*src
           end do
        end do
        !$OMP END PARALLEL DO 
        !$OMP TASKWAIT
        !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,ddt,tau0,k,tau,src,pp,kS2,Coeff,green)
        do J=N1+1,min(N2,Nmax)
           t=T1+(j-1)*dt
           ddt=(t0-thead)/Nint
           tau0=thead
           do  k=1,Nint-1
              tau=tau0+(k-1)*ddt+ddt/2.
              src=source(t-tau)
              pp=path(-abs(x),y,h,Vf1,VS2,tau)     
              pp=ima*aimag(pp)       
              kS2=zsqrt(pp**2+1/VS2**2)
              Coeff=calccoefff_pe(pp)
              green=derivee(-abs(x),y,h,Vf1,VS2,pp)/pi*ddt
              Ux(i,j)=Ux(i,j)+F*real(kS2*coeff(5)*green)*src
              Uy(i,j)=Uy(i,j)+F*real(ima*pp*coeff(5)*green)*src
           end do
           tau=t0-ddt/2.
           src=source(t-tau)
           pp=path(-abs(x),y,h,Vf1,VS2,tau)    
           pp=ima*aimag(pp)        
           kS2=zsqrt(pp**2+1/VS2**2)
           Coeff=calccoefff_pe(pp)
           green=sqrt(t0**2-tau**2)*derivee(-abs(x),y,h,Vf1,VS2,pp)*acos(1-ddt/t0)/pi*ddt
           Ux(i,j)=Ux(i,j)+F*real(kS2*coeff(5)*green)*src
           Uy(i,j)=Uy(i,j)+F*real(ima*pp*coeff(5)*green)*src   
        end do
        !$OMP END PARALLEL DO 
        !$OMP TASKWAIT
        !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,ddt,tau0,k,tau,src,pp,kS2,Coeff,green)
        do J=max(N1+1,N2+1),Nmax
           t=T1+(j-1)*dt
           ddt=(t0-t+2.*tdelay)/Nint
           tau0=t-2.*tdelay
           do  k=1,Nint-1
              tau=tau0+(k-1)*ddt+ddt/2.
              src=source(t-tau)
              pp=path(-abs(x),y,h,Vf1,VS2,tau)     
              pp=ima*aimag(pp)       
              kS2=zsqrt(pp**2+1/VS2**2)
              Coeff=calccoefff_pe(pp)
              green=derivee(-abs(x),y,h,Vf1,VS2,pp)/pi*ddt
              Ux(i,j)=Ux(i,j)+F*real(kS2*coeff(5)*green)*src
              Uy(i,j)=Uy(i,j)+F*real(ima*pp*coeff(5)*green)*src
           end do
           tau=t0-ddt/2.
           src=source(t-tau)
           pp=path(-abs(x),y,h,Vf1,VS2,tau)    
           pp=ima*aimag(pp)        
           kS2=zsqrt(pp**2+1/VS2**2)
           Coeff=calccoefff_pe(pp)
           green=sqrt(t0**2-tau**2)*derivee(-abs(x),y,h,Vf1,VS2,pp)*acos(1-ddt/t0)/pi*ddt
           Ux(i,j)=Ux(i,j)+F*real(kS2*coeff(5)*green)*src
           Uy(i,j)=Uy(i,j)+F*real(ima*pp*coeff(5)*green)*src
        enddo
        !$OMP END PARALLEL DO 
        !$OMP TASKWAIT
     end if
     !! Computation of the volume wave
     Nmin=ceiling((t0-T1)/dt+1)
     Nmax=floor((t0+2*tdelay-T1)/dt+1)
     Nmin=max(1,Nmin)
     Nmax=min(Nmax,Nt)
     !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,ddt,tau0,k,tau,src,pp,kS2,Coeff,green)
     do J=Nmin,Nmax
        t=T1+(j-1)*dt
        ddt=(t-t0)/Nint
        tau0=t0
        tau=tau0+ddt/2.
        src=source(t-tau)
        pp=path(-abs(x),y,h,Vf1,VS2,tau)   
        kS2=zsqrt(pp**2+1/VS2**2)
        Coeff=calccoefff_pe(pp)
        green=sqrt(tau**2-t0**2)*derivee(-abs(x),y,h,Vf1,VS2,pp)*dacosh(1+ddt/t0)/pi
        Ux(i,j)=Ux(i,j)+F*real(kS2*coeff(5)*green)*src
        Uy(i,j)=Uy(i,j)+F*real(ima*pp*coeff(5)*green)*src
        do k=2,Nint
           tau=tau0+(k-1)*ddt+ddt/2.
           src=source(t-tau)
           pp=path(-abs(x),y,h,Vf1,VS2,tau)   
           kS2=zsqrt(pp**2+1/VS2**2)
           Coeff=calccoefff_pe(pp)
           green=derivee(-abs(x),y,h,Vf1,VS2,pp)/pi*ddt
           Ux(i,j)=Ux(i,j)+F*real(kS2*coeff(5)*green)*src
           Uy(i,j)=Uy(i,j)+F*real(ima*pp*coeff(5)*green)*src
        end do
     end do
     !$OMP END PARALLEL DO 
     !$OMP TASKWAIT
     !$OMP PARALLEL DO   DEFAULT(SHARED) PRIVATE(J,t,ddt,tau0,k,tau,src,pp,kS2,Coeff,green)
     do J=Nmax+1,Nt
        t=T1+(j-1)*dt
        ddt=2.*tdelay/Nint
        tau0=t-2.*tdelay
        do k=1,Nint
           tau=tau0+(k-1)*ddt+ddt/2.
           src=source(t-tau)
           pp=path(-abs(x),y,h,Vf1,VS2,tau)   
           kS2=zsqrt(pp**2+1/VS2**2)
           Coeff=calccoefff_pe(pp)
           green=derivee(-abs(x),y,h,Vf1,VS2,pp)/pi*ddt
           Ux(i,j)=Ux(i,j)+F*real(kS2*coeff(5)*green)*src
           Uy(i,j)=Uy(i,j)+F*real(ima*pp*coeff(5)*green)*src
        end do
     end do
     !$OMP END PARALLEL DO 
     !$OMP TASKWAIT
  end do
end subroutine sub_transmitfS_pe
