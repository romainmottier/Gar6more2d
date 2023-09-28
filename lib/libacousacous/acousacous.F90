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

subroutine acousacous
  use m_phys  
  use m_num
  use m_source
  use m_sismo
  implicit none 
!!! Computation of the incident wave
  allocate(Ux(1:nx,1:nt))
  allocate(Uy(1:nx,1:nt))
  allocate(P(1:nx,1:nt))
  Ux=0.d0
  Uy=0.d0
  P=0.d0
  if (type_medium.eq.1) then
!!!Infinite Medium
     write(6,*) 'Computation of the incident wave'
     call sub_incid_aa
  elseif (type_medium.eq.2) then
!!!Free Surface
     if (y.ge.0d0) then
        write(6,*) 'Computation of the incident wave'
        call sub_incid_aa
        write(6,*) 'Computation of the reflected wave'
        call sub_reflex_free_aa
     else
        write(6,*) 'Error, y can not be smaller than 0 in this configuration'
        stop
     end if
  elseif (type_medium.eq.3) then
!!!Wall Boundary
     if (y.ge.0d0) then
        write(6,*) 'Computation of the incident wave'
        call sub_incid_aa
        write(6,*) 'Computation of the reflected wave'
        call sub_reflex_wall_aa
     else
        write(6,*) 'Error, y can not be smaller than 0 in this configuration'
        stop
     end if
  elseif (type_medium.eq.4) then
!!! Bilayered Medium
     if (y.ge.0d0) then
        write(6,*) 'Computation of the incident wave'
        call sub_incid_aa
        write(6,*) 'Computation of the reflected wave'
        call sub_reflex_aa
     else
        write(6,*)  'Computation of the transmitted wave'
        call sub_transmit_aa
     end if
  elseif (type_medium.eq.5) then
!!! Bilayered Medium
     if (y.ge.0d0) then
        write(6,*) 'Computation of the incident wave'
        call sub_incid_aa
        write(6,*) 'Computation of the first reflected wave'
        call sub_reflex_aa
        write(6,*) 'Computation of the second reflected wave'
        h=2*width-h
        P=-P
        Uy=-Uy
        call sub_incid_aa
        write(6,*) 'Computation of the third reflected wave'
       call sub_reflex_aa
             write(6,*) 'Computation of the fourth reflected wave'
        y=2*width-y
        h=2*width-h
        call sub_reflex_aa
        P=-P
        Uy=-Uy
        h=2*width-h
        write(6,*) 'Computation of the fifth reflected wave'
        call sub_reflex_aa
!        P=-P
!        Uy=-Uy
     else
        write(6,*)  'Computation of the first transmitted wave'
        call sub_transmit_aa
        h=-2*width+h
        write(6,*)  'Computation of the first transmitted wave'
        call sub_transmit_aa
     end if
  end if
end subroutine acousacous
