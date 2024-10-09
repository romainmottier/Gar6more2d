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

program Gar6more2D

  use m_phys  
  use m_num
  use m_source
  use m_sismo
  use m_const
  implicit none 
  integer I,J
  real*8 delta,Fstar(2),Mataux(2,2),P1inv(2,2),A1inv(2,2),A2inv(2,2)
  real*8 kb1,ks1,kf1,kb2,ks2,kf2,x
  pi=2.d0*dasin(1.D0)
  ima=(0,1)
  open(20,FILE="Gar6more2D.dat")
  write(6,*) 'Do you want'
  write(6,*) '- An infinite homogeneous medium (1);'
  write(6,*) '- A semi-infinite homogeneous media with a free surface boundary condition at it bottom (2):'
  write(6,*) '- A semi-infinite homogeneous media with a wall boundary condition at it bottom (3):'
  write(6,*) '- A bilayered infinite media with a plane interface (4):'
  write(6,*) '- A bilayered infinite media with a plane interface  a free surface boundary condition at it top (5):'
  read(20,*)  type_medium
  write(6,*)  type_medium
  if(type_medium.ge.4) then
     write(6,*) 'Is the first medium acoustic (0), elastodynamic (1) or porous (2)?'
     read(20,*) type_medium_1
     write(6,*) type_medium_1
  else
     write(6,*) 'Is the first medium acoustic (0), elastodynamic (1) or porous (2)?'
     read(20,*) type_medium_1
     write(6,*) type_medium_1
  end if
  write(6,*) 'Source time derivative (2 for Ricker)?'
  read(20,*) deriv
  write(6,*) deriv
  write(6,*) 'Frequency of the source?'
  read(20,*) f0
  write(6,*) f0
  if(type_medium_1.eq.0) then
     write(6,*) 'Amplitude of the source?'
     read(20,*) amp
     write(6,*) amp
  elseif(type_medium_1.eq.1) then
     write(6,*) 'Amplitude of the P source?'
     read(20,*) ampP
     write(6,*) ampP
     write(6,*) 'Amplitude of the S source?'
     read(20,*) ampS
     write(6,*) ampS
  elseif(type_medium_1.eq.2) then
     write(6,*) 'Amplitude of the bulk source source (fu)?'
     read(20,*) ampbulk
     write(6,*) ampbulk
     write(6,*) 'Amplitude of the pressure source (fp)?'
     read(20,*) ampP
     write(6,*) ampP
  end if
  write(6,*) 'delay of the source (f(t)=0 if t>2tdelay)?'
  read(20,*) tdelay
  write(6,*) tdelay
  write(6,*) 'Height of the source?'
  read(20,*) h
  write(6,*) h
  write(6,*) 'Height of the line of receivers?'
  read(20,*) Y
  write(6,*) Y
  write(6,*) 'abscissa of the first receiver?'
  read(20,*) X1
  write(6,*) X1
  write(6,*) 'abscissa of the last receiver?'
  read(20,*) X2
  write(6,*) X2
  write(6,*) 'How many receivers?'
  read(20,*) NX
  write(6,*) NX
  if (nx.ge.2) then
     dx=(X2-X1)/(Nx-1)
  else
     dx=0
  end if
  write(6,*) 'dx=',dx
  write(6,*) 'At which time shall the seismo begin?'
  read(20,*) T1
  write(6,*) T1
  write(6,*) 'At which time shall the seismo stop?'
  read(20,*) T2
  write(6,*) T2
  write(6,*) 'What is the time-step?'
  read(20,*) dt
  write(6,*) dt
  Nt=aint((T2-T1)/dt)+1
  write(6,*) Nt
  write(6,*) 'How many intervals are required for the numerical computation of the time convolution?'
  read(20,*) Nint
  write(6,*) Nint
  if (type_medium_1.eq.0) then
     write(6,*) 'mu and rho?'
     read(20,*) mu1, rho1
     write(6,*) mu1, rho1
     V1=dsqrt(mu1/rho1)
     write(6,*) 'V1=', V1
  elseif  (type_medium_1.eq.1) then
     write(6,*) 'mu,lambda,rho'
     read(20,*) mu1,lambda1,rho1
     write(6,*) mu1,lambda1,rho1
     Vp1=dsqrt((lambda1+2*mu1)/rho1)
     Vs1=dsqrt(mu1/rho1)
     write(6,*) 'Vp1, Vs1',Vp1, Vs1
  elseif  (type_medium_1.eq.2) then
     write(6,*) 'kb (Frame bulk modulus), ks (Solid bulk modulus),&
          & kf (Fluid bulk modulus), mu (Frame shear modulus),&
          & rhof (Fluid density), rhos (Solid density), &
          &phi (Porosity), a(Tortuosity)?' 
     read(20,*) kb1,ks1,kf1,mu1,rhof1,rhos1,phi1,a1  
     write(6,*) kb1,ks1,kf1,mu1,rhof1,rhos1,phi1,a1  
     lambda1=kb1-2*mu1/3;
     rho1=phi1*rhof1+(1-phi1)*rhos1;
     rhow1=a1*rhof1/phi1;
     beta1=1-kb1/ks1;
     M1=1/(phi1/kf1+(beta1-phi1)/ks1);
     la1=lambda1+M1*(beta1-phi1)**2;
     R1=M1*phi1**2;
     ga1=M1*phi1*(beta1-phi1);
     S2=la1+2*mu1;
     AA1(1,1)=rho1 
     AA1(1,2)=rhof1
     AA1(2,1)=rhof1
     AA1(2,2)=rhow1;
     A1inv(1,1)=AA1(2,2)
     A1inv(1,2)=-AA1(1,2)
     A1inv(2,1)=-AA1(2,1)
     A1inv(2,2)=AA1(1,1)
     A1inv=A1inv/(AA1(1,1)*AA1(2,2)-AA1(1,2)**2)
     B1(1,1)=lambda1+2*mu1+M1*beta1**2 
     B1(1,2)=M1*beta1
     B1(2,1)=M1*beta1 
     B1(2,2)=M1;
     TT1=matmul(A1inv,B1);
     delta=(TT1(1,1)+TT1(2,2))**2+4*(TT1(1,2)*TT1(2,1)-TT1(1,1)*TT1(2,2))
     D1(1)=((TT1(1,1)+TT1(2,2))+sqrt(Delta))/2.
     D1(2)=((TT1(1,1)+TT1(2,2))-sqrt(Delta))/2.
     if (TT1(2,1).eq.0.d0) then
        if (D1(1).eq.TT1(1,1)) then
           P1(1,1)=1.;P1(2,1)=0.;
           if(TT1(1,2).eq.0.d0) then
              P1(1,2)=0.;P1(2,2)=1.;
           else
              P1(1,2)=-TT1(1,2)/(TT1(1,1)-D1(2));P1(2,2)=1.;
           end if
        else
           P1(1,2)=1.;P1(2,2)=0.
           if(TT1(1,2).eq.0.d0) then
              P1(1,1)=0.;P1(2,1)=1.;
           else
              P1(1,1)=-TT1(1,2)/(TT1(1,1)-D1(1));P1(2,1)=1.;
           end if
        end if
     else
        P1(1,1)=-(TT1(2,2)-D1(1))/TT1(2,1);P1(2,1)=1.
        P1(1,2)=-(TT1(2,2)-D1(2))/TT1(2,1);P1(2,2)=1.
     end if
     P1(:,1)=P1(:,1)/sqrt(sum(P1(:,1)**2))
     P1(:,2)=P1(:,2)/sqrt(sum(P1(:,2)**2))
     Vpsi1=sqrt(mu1*rhow1/(rho1*rhow1-rhof1**2))
     Vf1=sqrt(D1(1))
     Vs1=sqrt(D1(2))
     write(6,*) 'Vpsi1, Vf1, Vs1',Vpsi1,Vf1,Vs1
     write(6,*) 'P1'
     write(6,*) P1(1,1),P1(1,2)
     write(6,*) P1(2,1),P1(2,2)
     Fstar(1)=ampbulk-M1*beta1*ampP
     Fstar(2)=ampbulk-M1*ampP
     P1inv(1,1)=P1(2,2)
     P1inv(1,2)=-P1(1,2)
     P1inv(2,1)=-P1(2,1)
     P1inv(2,2)=P1(1,1)
     P1inv=P1inv/(P1(1,1)*P1(2,2)-P1(1,2)*P1(2,1))
     Mataux=matmul(P1inv,A1inv)
     Fstar=matmul(Mataux,Fstar);
     F1=Fstar(1);
     F2=Fstar(2);
     write(6,*) 'F1, F2', F1, F2
  end if
  if (type_medium.ge.4) then
     write(6,*) 'Is the second medium acoustic (0) elastodynamic (1) or porous (2)?'
     read(20,*) type_medium_2
     write(6,*) type_medium_2
     select case(type_medium_2)
     case(0)
        write(6,*) 'mu and rho?'
        read(20,*) mu2,rho2
        write(6,*) mu2,rho2
        V2=dsqrt(mu2/rho2)
        write(6,*) 'V2=', V2
     case(1)
        write(6,*) 'mu,lambda,rho'
        read(20,*) mu2,lambda2,rho2
        write(6,*) mu2,lambda2,rho2
        Vp2=dsqrt((lambda2+2*mu2)/rho2)
        Vs2=dsqrt(mu2/rho2)
        write(6,*) 'Vp2, Vs2',Vp2, Vs2
     case(2)
        write(6,*) ' kb (Frame bulk modulus), ks (Solid bulk modulus), kf (Fluid bulk modulus), mu (Frame shear modulus),&
             & rhof (Fluid density), rhos (Solid density), phi (Porosity), a(Tortuosity)?' 
        read(20,*) kb2,ks2,kf2,mu2,rhof2,rhos2,phi2,a2  
        write(6,*) kb2,ks2,kf2,mu2,rhof2,rhos2,phi2,a2  
        lambda2=kb2-2*mu2/3;
        rho2=phi2*rhof2+(1-phi2)*rhos2;
        rhow2=a2*rhof2/phi2;
        beta2=1-kb2/ks2;
        M2=1/(phi2/kf2+(beta2-phi2)/ks2);
        la2=lambda2+M2*(beta2-phi2)**2;
        R2=M2*phi2**2;
        ga2=M2*phi2*(beta2-phi2);
        S2=la2+2*mu2;
        AA2(1,1)=rho2 
        AA2(1,2)=rhof2
        AA2(2,1)=rhof2
        AA2(2,2)=rhow2;
        B2(1,1)=lambda2+2*mu2+M2*beta2**2 
        B2(1,2)=M2*beta2
        B2(2,1)=M2*beta2 
        B2(2,2)=M2;
        A2inv(1,1)=AA2(2,2)
        A2inv(1,2)=-AA2(1,2)
        A2inv(2,1)=-AA2(2,1)
        A2inv(2,2)=AA2(1,1)
        A2inv=A2inv/(AA2(1,1)*AA2(2,2)-AA2(1,2)**2)
        TT2=matmul(A2inv,B2);
        delta=(TT2(1,1)+TT2(2,2))**2+4*(TT2(1,2)*TT2(2,1)-TT2(1,1)*TT2(2,2))
        D2(1)=((TT2(1,1)+TT2(2,2))+sqrt(Delta))/2.
        D2(2)=((TT2(1,1)+TT2(2,2))-sqrt(Delta))/2.
        if (TT2(2,1).eq.0.d0) then
           if (D2(1).eq.TT2(1,1)) then
              P2(1,1)=1.;P2(2,1)=0.;
              if(TT2(1,2).eq.0.d0) then
                 P2(1,2)=0.;P2(2,2)=1.;
              else
                 P2(1,2)=-TT2(1,2)/(TT2(1,1)-D2(2));P2(2,2)=1.;
              end if
           else
              P2(1,2)=1.;P2(2,2)=0.
              if(TT2(1,2).eq.0.d0) then
                 P2(1,1)=0.;P2(2,1)=1.;
              else
                 P2(1,1)=-TT2(1,2)/(TT2(1,1)-D2(1));P2(2,1)=1.;
              end if
           end if
        else
           P2(1,1)=-(TT2(2,2)-D2(1))/TT2(2,1);P2(2,1)=1.
           P2(1,2)=-(TT2(2,2)-D2(2))/TT2(2,1);P2(2,2)=1.
        end if
        P2(:,1)=P2(:,1)/sqrt(sum(P2(:,1)**2))
        P2(:,2)=P2(:,2)/sqrt(sum(P2(:,2)**2))
        Vpsi2=sqrt(mu2*rhow2/(rho2*rhow2-rhof2**2))
        Vf2=sqrt(D2(1))
        Vs2=sqrt(D2(2))
        write(6,*) 'Vpsi2, Vf2, Vs2',Vpsi2,Vf2,Vs2
        write(6,*) 'P2'
        write(6,*) P2(1,1),P2(1,2)
        write(6,*) P2(2,1),P2(2,2)
        if (type_medium_1.eq.0) then
           write(6,*) 'Do you want open pore (1) or sealed pore (2) conditions ?'
           read(20,*) open
           write(6,*) open
        end if
     case DEFAULT
        write(6,*) 'Not implemented'
        stop
     end select
  end if
  if(type_medium.ge.5) then
     write(6,*) 'Width of the layer ?'
     read(20,*) width
     write(6,*) width
  endif
  close(20)
  if (type_medium.lt.4) then 
     select case(type_medium_1)
     case(0)
        call acousacous
     case(1)
        call elastoelasto
     case(2)
        call poroporo
     case DEFAULT
        write(6,*) 'Not implemented'
        stop
     end select
  else
     select case(type_medium_1)
     case(0)
        select case(type_medium_2)
        case(0)
           call acousacous
        case(1)
           call acouselasto
        case(2)
           call acousporo
        case  DEFAULT
           write(6,*) 'Not implemented'
           stop
        end select
     case(1)
        select case(type_medium_2)
        case(1)
           call elastoelasto
        case  DEFAULT
           write(6,*) 'Not implemented'
           stop
        end select
     case(2)
        select case(type_medium_2)
     case(1)
        call poroelasto
     case(2)
        call poroporo
     case  DEFAULT
        write(6,*) 'Not implemented'
        stop
     end select
  case DEFAULT
     write(6,*) 'Not implemented'
     stop
  end select
end if
do I=1,Nx
  x=X1+(I-1)*dx
  if (x.gt.0.) Ux(i,:)=-Ux(i,:)
end do
if(type_medium_1.eq.0) then
  open(20,FILE='Ux.dat')
  open(21,FILE='Uy.dat')
  open(22,FILE='P.dat')
  do J=1,Nt
     write(20,*) T1+(j-1)*dt,(Ux(i,j),i=1,Nx)
     write(21,*) T1+(j-1)*dt,(Uy(i,j),i=1,Nx)
     write(22,*) T1+(j-1)*dt,(P(i,j),i=1,Nx)
  end do
  close(20)
  close(21)
  close(22)
else
  open(20,FILE='Ux.dat')
  open(21,FILE='Uy.dat')
  do J=1,Nt
     write(20,*) T1+(j-1)*dt,(Ux(i,j),i=1,Nx)
     write(21,*) T1+(j-1)*dt,(Uy(i,j),i=1,Nx)
  end do
  close(20)
  close(21)
end if
end program
