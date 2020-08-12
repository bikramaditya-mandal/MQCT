      SUBROUTINE  INI_ARRAYS !!! remeber to change int( to round(
! This subroutine is written by Alexander Semenov and modified by Bikramaditya Mandal
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA	  
      IMPLICIT NONE
      INTEGER i,st,j_count,j_summ,m_count,j1_count,j2_count
      INTEGER p_count
      INTEGER st1,st2,KRONEKER,p_lim_max_ini,round	 
      EXTERNAL KRONEKER,round	  
      st = 0
      IF(MYID.eq.0) PRINT *, "ARRAY_INI_STARTED"
      ALLOCATE(chann_indx(number_of_channels))
      IF(.not.identical_particles_defined) THEN	  
      p_lim_max = 1
      p_lim_min = 1	
      ELSE
      SELECT CASE(coll_type)		  
      CASE(5)	  

      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
      CASE(6)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))*
     & KRONEKER(v1_ch(chann_ini),v2_ch(chann_ini)) 	  
      CASE(0)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
     & *KRONEKER(ka1_ch(chann_ini),ka2_ch(chann_ini))
     & *KRONEKER(kc1_ch(chann_ini),kc2_ch(chann_ini))	  
      END SELECT		  
      ENDIF	
      IF(identical_particles_defined) THEN
      IF(coll_type.eq.5 .or. coll_type.eq.6 .or. coll_type.eq.0) THEN
      IF(myid.eq.0) THEN
      PRINT*,"EXCHANGE SYMMETRIES FOR INDETICAL PARTICLES:"
      IF(exch_par_w_pl.eq.1d0) THEN
      PRINT*, "POSITIVE EXCHANGE PARITY ONLY"
      ENDIF	  
      IF(exch_par_w_pl.lt.1d0 .and. exch_par_w_pl.gt.0d0) THEN
      PRINT*, "BOTH EXCHANGE PARITIES"
      ENDIF	  
      IF(exch_par_w_pl.eq.0d0) THEN
      PRINT*, "NEGATIVE EXCHANGE PARITY ONLY"
      ENDIF	  	  
      ENDIF	 
      ELSE
      IF(myid.eq.0) 
     & PRINT*,
     & "ERROR: EXCHANGE SYMMETRY DEFINED FOR NON-IDENTICAL PARTICLES"
      STOP	
      ENDIF	  
      ENDIF	  
      chann_indx = 0!!! CREATING CHANNEL -> STATE	  
      DO i=1,number_of_channels
      SELECT CASE(coll_type)
      CASE(1)
      IF(fine_structure_defined .and. LORB_FINE.gt.0 .and.
     & SPIN_FINE.eq.2 ) THEN
      DO j_count = 1, round(j_h_ch(i)*2) + 1	  
      st = st + 1
      IF(j_count*2-2 .eq. round(j_h_ch(i)*2) + 1) chann_indx(i) = st !!!! SLOWLY MODIFIYNG
      ENDDO	
      ELSE
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      IF(j_count.eq.0) chann_indx(i) = st
      ENDDO	  
      ENDIF	  
      CASE(2)
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      IF(j_count.eq.0) chann_indx(i) = st
      ENDDO
      CASE(3)
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      IF(j_count.eq.0) chann_indx(i) = st
      ENDDO
      CASE(4)
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      IF(j_count.eq.0) chann_indx(i) = st
      ENDDO
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i)))
     &	  chann_indx(i) = st
      ENDDO
      ENDDO
      ELSE
      p_lim_min = 1
      p_lim_max = min(2-KRONEKER(j1_ch(i),j2_ch(i)),p_lim_max_ini)	  
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.p_lim_min  )
     &	  chann_indx(i) = st
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i)))
     &	  chann_indx(i) = st
      ENDDO
      ENDDO
      ELSE
      p_lim_min = 1
      p_lim_max = min(2-KRONEKER(j1_ch(i),j2_ch(i))*
     & KRONEKER(v1_ch(i),v2_ch(i)),p_lim_max_ini) 
      DO p_count=p_lim_min,p_lim_max	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and.  p_count .eq.p_lim_min  )
     &	  chann_indx(i) = st
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	  
      CASE(7)
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i)))
     &	  chann_indx(i) = st
      ENDDO
      ENDDO
      CASE(8)
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i)))
     &	  chann_indx(i) = st
      ENDDO
      ENDDO
      CASE(9)
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i)))
     &	  chann_indx(i) = st
      ENDDO
      ENDDO
      CASE(0)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i)))
     &	  chann_indx(i) = st
      ENDDO
      ENDDO
      ELSE
      p_lim_min = 1
      p_lim_max = min(2-KRONEKER(j1_ch(i),j2_ch(i))
     & *KRONEKER(ka1_ch(i),ka2_ch(i))*KRONEKER(kc1_ch(i),kc2_ch(i)),
     & p_lim_max_ini)
      DO p_count=p_lim_min,p_lim_max	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and.  p_count .eq.1  )
     &	  chann_indx(i) = st
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	   	  
      END SELECT	  
      ENDDO
      states_size = st
      IF(MYID.eq.0) PRINT *, "STATES_SIZE_FOUND",states_size
      IF(MYID.eq.0) PRINT *, "number_of_channels",number_of_channels
      IF(MYID.eq.0) PRINT *, "jmax_included",jmax_included
	  
      ALLOCATE(indx_chann(states_size))
      ALLOCATE(j_min_ind(number_of_channels))
      ALLOCATE(j_max_ind(number_of_channels))	  
      ALLOCATE(m12(states_size))
      ALLOCATE(j12(states_size))
      ALLOCATE(Ej(states_size))
      IF(fine_structure_defined .and. LORB_FINE.gt.0 .and.
     & SPIN_FINE.eq.2 ) THEN
      ALLOCATE(j12_h(states_size))
      ALLOCATE(m12_h(states_size))
      ENDIF	 
      ALLOCATE(indx_corr(jmax_included*4+1,
     & jmax_included*2+1,number_of_channels))
      indx_corr = 0	 
      IF(identical_particles_defined) THEN	  
      ALLOCATE(parity_state(states_size))
      ALLOCATE(indx_corr_id(2,4*jmax_included+1,
     & jmax_included*2+1,number_of_channels))
      indx_corr_id = 0	 
      ENDIF
      IF(MYID.eq.0) PRINT *, "BASIS ARRAYS CREATED"
      indx_chann = 0
      m12 = 0
      j12 = 0	  
c      PRINT*,"states_size",states_size
c      PRINT*,"states",chann_indx
      st = 0	  
      DO i=1,number_of_channels
	  
      SELECT CASE(coll_type)
      CASE(1)
      IF(fine_structure_defined .and. LORB_FINE.gt.0 .and.
     & SPIN_FINE.eq.2 ) THEN
      DO j_count = 1, round(j_h_ch(i)*2) + 1	  
      st = st + 1
      indx_chann(st) = i
      j12_h(st) = j_h_ch(i)
      m12_h(st) = j_count-1 - j_h_ch(i)	
      indx_corr(round(m12_h(st) +j_h_ch(i)+1),round(j_h_ch(i)+0.5d0),i)
     &	  = st	!!!!! STOPED HERE 
      j_min_ind(i) = int(j_h_ch(i))
      j_max_ind(i) = int(j_h_ch(i))	
      Ej(st) = E_ch(i)	  
      ENDDO		 
      ELSE	 
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      indx_chann(st) = i
      j12(st) = j_ch(i)
      m12(st) = j_count
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_ch(i)+1,j_ch(i)+1,i) = st
      ENDDO
      j_min_ind(i) = j_ch(i)
      j_max_ind(i) = j_ch(i)
      ENDIF	  
      CASE(2)
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      indx_chann(st) = i
      j12(st) = j_ch(i)
      m12(st) = j_count
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_ch(i)+1,j_ch(i)+1,i) = st
      ENDDO
      j_min_ind(i) = j_ch(i)
      j_max_ind(i) = j_ch(i)	  
      CASE(3)
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      indx_chann(st) = i
      j12(st) = j_ch(i)
      m12(st) = j_count
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_ch(i)+1,j_ch(i)+1,i) = st
      ENDDO
      j_min_ind(i) = j_ch(i)
      j_max_ind(i) = j_ch(i)		  
      CASE(4)
      DO j_count = -j_ch(i),j_ch(i)	  
      st = st + 1
      indx_chann(st) = i
      j12(st) = j_ch(i)
      m12(st) = j_count
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_ch(i)+1,j_ch(i)+1,i) = st
      ENDDO
      j_min_ind(i) = j_ch(i)
      j_max_ind(i) = j_ch(i)	  
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(*,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_summ+1,j_summ+1,i) = st	  
      ENDDO
      ENDDO
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)		  
      ELSE
      p_lim_min = 1
      p_lim_max = min(2-KRONEKER(j1_ch(i),j2_ch(i)),p_lim_max_ini)		  

      DO p_count=p_lim_min,p_lim_max	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
      j12(st) = j_summ
      m12(st) = j_count
	  
      Ej(st) = E_ch(i)
      parity_state(st) = (-1)**(p_count-1)
c      IF((j1_ch(i)-j2_ch(i)).eq.0) THEN
c      parity_state(st) = (-1)**j_summ	  
c      ENDIF
      indx_corr_id(p_count,
     & j_count+j_summ+1,j_summ+1,i) = st
c      PRINT*,st,j12(st),m12(st),      indx_corr_id(p_count,
c     & j_count+j_summ+1,j_summ+1,i),parity_state(st),indx_chann(st),
c     & ej(st)

	 
      ENDDO
      ENDDO
      ENDDO	 
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_summ+1,j_summ+1,i) = st	  
      ENDDO
      ENDDO
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)		  
      ELSE
      p_lim_max = min(2-KRONEKER(j1_ch(i),j2_ch(i))*
     & KRONEKER(v1_ch(i),v2_ch(i)),p_lim_max_ini) 	  

      DO p_count=p_lim_min,p_lim_max	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      parity_state(st) = (-1)**(p_count-1)
c      IF((j1_ch(i)-j2_ch(i)).eq.0) THEN
c      parity_state(st) = (-1)**j_summ	  
c      ENDIF
      indx_corr_id(p_count,
     & j_count+j_summ+1,j_summ+1,i) = st	  
      ENDDO
      ENDDO
      ENDDO	 
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)	  
      ENDIF	 
      CASE(7)
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_summ+1,j_summ+1,i) = st  
      ENDDO
      ENDDO
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)	  
      CASE(8)
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_summ+1,j_summ+1,i) = st  
      ENDDO
      ENDDO
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)	  
      CASE(9)
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_summ+1,j_summ+1,i) = st  
      ENDDO
      ENDDO
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)	  
      CASE(0)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      indx_corr(j_count+j_summ+1,j_summ+1,i) = st	  
      ENDDO
      ENDDO
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)		  
      ELSE
      p_lim_max = min(2-KRONEKER(j1_ch(i),j2_ch(i))
     & *KRONEKER(ka1_ch(i),ka2_ch(i))*KRONEKER(kc1_ch(i),kc2_ch(i)),
     & p_lim_max_ini)	  

      DO p_count=p_lim_min,p_lim_max
      DO j_summ = abs(j1_ch(i)-j2_ch(i)),j1_ch(i)+j2_ch(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      indx_chann(st) = i
c      WRITE(1,*) i,st
      j12(st) = j_summ
      m12(st) = j_count
c      PRINT*,"channel,states=",i,st+m_count+j_summ,j_summ	  
      Ej(st) = E_ch(i)
      parity_state(st) = (-1)**(p_count-1)
c      IF((j1_ch(i)-j2_ch(i)).eq.0) THEN
c      parity_state(st) = (-1)**j_summ	  
c      ENDIF
      indx_corr_id(p_count,
     & j_count+j_summ+1,j_summ+1,i) = st	  
      ENDDO
      ENDDO
      ENDDO	 
      j_min_ind(i) =  abs(j1_ch(i)-j2_ch(i))
      j_max_ind(i) = j1_ch(i)+j2_ch(i)	  
      ENDIF
      END SELECT      
	  
      ENDDO
c      PRINT*, "INITIATION FINISHED"	
	  
c       PRINT*, j_min_ind,j_max_ind   
c       PRINT*,m12,j12
c c      PRINT*,Ej
c       PRINT*,parity_state
c       PRINT*,  indx_corr_id    	  
c      STOP "HERE" 	  
      total_size = 0	  
      DO st1=1,states_size
      DO st2=1,st1
      IF(identical_particles_defined) THEN
      IF(parity_state(st1).ne.parity_state(st2)) CYCLE  
      ENDIF	  
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) THEN
      IF(int(m12_h(st1)).eq.int(m12_h(st2)).and.
     & m12_h(st2)*m12_h(st1).gt.0 ) total_size = total_size + 1	  
      ELSE	  
      IF(m12(st1).eq.m12(st2)) total_size = total_size + 1
      ENDIF	  
      ENDDO
      ENDDO
!      PRINT*, "TOTAL_SIZE",total_size        !IT WAS CHANGED
!      STOP	  
      ALLOCATE(ind_mat(2,total_size))
	  
c      IF(MYID.eq.0) PRINT *, "TOTAL_SIZE_NOW",total_size
c      STOP	  
	  st = 0
      DO st1=1,states_size
      DO st2=1,st1
      IF(identical_particles_defined) THEN
      IF(parity_state(st1).ne.parity_state(st2)) CYCLE  
      ENDIF	  	  
      
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) THEN
      IF(int(m12_h(st1)).eq.int(m12_h(st2)) .and.
     & m12_h(st2)*m12_h(st1).gt.0 )THEN
      st  = st+1	  
      ind_mat(1,st) = st1
      ind_mat(2,st) = st2	  
      ENDIF	  
      ELSE	  
	  IF(m12(st1).eq.m12(st2)) THEN
      st  = st+1	  
      ind_mat(1,st) = st1
      ind_mat(2,st) = st2	  
      ENDIF
      ENDIF	  
      ENDDO
      ENDDO
	  
      IF(st.ne.total_size) STOP"ERROR: INI FAILED"
      IF(mpi_task_defined) THEN
      IF(MYID.EQ.0) THEN	  
      ALLOCATE(Mat_el(n_r_coll,total_size))
      ALLOCATE(Mat_el_der(n_r_coll,total_size))
      ENDIF	  
      ELSE
      ALLOCATE(Mat_el(n_r_coll,total_size))
      ALLOCATE(Mat_el_der(n_r_coll,total_size))	  
      ENDIF	  
      ALLOCATE(R_COM(n_r_coll))
c      IF(identical_particles_defined) THEN	  
c      DO st=1,states_size
c      i = indx_chann(st)   
c      WRITE(*,'(i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
c    & i,st,j1_ch(i),j2_ch(i),j12(st),m12(st),
c     &l_state(st),parity_state(st)	  
c      ENDDO
c      ENDIF	  
c     STOP	
      IF(MYID.eq.0 .and. print_states_defined) THEN  !!! COMMENT FOR A WHILE
      PRINT*, "ARRAYS_INITIALIZED"
      OPEN(4,FILE="STATES.out")

      SELECT CASE(coll_type)
      CASE(1)
      IF(SPIN_FINE.eq.2 .and. LORB_FINE.gt.0) THEN
      WRITE(4,'(a8,1x,a8,1x,a4,1x,a5,1x,1x,a4)')
     & "#STATE","#CHANNEL","J12","M12","J"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,f4.1,1x,f5.1,1x,1x,f4.1)')
     & st,i,j12_h(st),m12_h(st),j_h_ch(i)	  
      ENDDO	  
      ELSE	  
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,1x,i3)')
     & st,i,j12(st),m12(st),j_ch(i)	  
      ENDDO
      ENDIF	  
      CASE(2)
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","V","J"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),v_ch(i),j_ch(i)	  
      ENDDO
      CASE(3)
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J","K","EPS"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j_ch(i),k_ch(i),eps_ch(i)	  
      ENDDO
      CASE(4)
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J","KA","KC"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j_ch(i),ka_ch(i),kc_ch(i)	  
      ENDDO	  
      CASE(5)
      IF(.NOT.identical_particles_defined) THEN 	  
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","J2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i),parity_state(st)
      ENDDO	 
      ENDIF
      CASE(6)
      IF(.NOT.identical_particles_defined) THEN 	  
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","V1","J1","V2","J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),v1_ch(i),j1_ch(i),v2_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      WRITE(4,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","V1","J1","V2","J2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),v1_ch(i),j1_ch(i),v2_ch(i),j2_ch(i),
     & parity_state(st)	 
      ENDDO 
      ENDIF
      CASE(7)
      WRITE(4,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","K1","EPS", "J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),k1_ch(i),eps1_ch(i),j2_ch(i)	  
      ENDDO	  
      CASE(8)
      WRITE(4,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1","J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i),j2_ch(i)	  
      ENDDO	  
      CASE(9)
      WRITE(4,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1",
     & "J2","K2","EPS"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),k2_ch(i),eps2_ch(i)	 
      ENDDO	  
      CASE(0)
      IF(.NOT.identical_particles_defined) THEN  
      WRITE(4,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1",
     & "J2","KA2","KC2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i)	 
      ENDDO
      ELSE
      WRITE(4,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a2)
     & ')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1",
     & "J1","KA2","KC2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(4,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i2)
     & ')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i),parity_state(st)
      ENDDO	 
      ENDIF	 
      END SELECT
      CLOSE(4)
      OPEN(6,FILE="STATE_Mij_INDX.out")
      WRITE(6,"(a9,1x,a4,1x,a4)") "Mij_INDEX", "ST_1", "ST_2"
      DO i=1,total_size
      WRITE(6,"(i9,1x,i4,1x,i4)") i,ind_mat(1,i),ind_mat(2,i)	  
      ENDDO	  
      CLOSE(6)
      ENDIF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      END SUBROUTINE INI_ARRAYS 
      SUBROUTINE MATRIX_MAKE
! This subroutine is written by Alexander Semenov and modified by Bikramaditya Mandal
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      USE MPI_TASK_TRAJECT	  
      USE POT_STORE
      USE CS_MATRIX	 	  
      USE OLD_MIJ	  
      IMPLICIT NONE
      LOGICAL file_exst
      CHARACTER(LEN=22) :: R_GRID_DAT = "USER_DEFINED_RGRID.DAT"	  
      INTEGER ISTAT	  
      INTEGER i,st_1,st_2,k,k_st_mpi,k_fn_mpi,chunk_mpi_size,task_size
      INTEGER L_NJU(8)
      INTEGER mean_size_1,mean_size_2,i_nr_fin,i_nr_ini
      INTEGER status(MPI_STATUS_SIZE)	  
      INTEGER k_check,i_check,n_point_dis_1,n_point_dis_2	  
      REAL*8 intgeral,R_dist,dR_step
      REAL*8 TIME_WRITING_MATRIX
	  real*8 bk_time_bgn, bk_time_fin,bk_time		!Bikram Jan,'19	  
      REAL*8, ALLOCATABLE :: Mat_el_temp(:,:),Mat_el_der_temp(:,:),
     & Mat_el_res(:),Mat_el_der_res(:),Mat_el_resf(:,:)
     &,Mat_el_der_resf(:,:)
      REAL*8,ALLOCATABLE :: Mat_el_r_temp(:,:), Mat_el_der_r_temp(:,:)
      LOGICAL, ALLOCATABLE :: K_SKIPPED_BY_ROUTINE(:)
      LOGICAL :: K_SKIPPED_RES = .FALSE.
      LOGICAL BELONGS
      EXTERNAL BELONGS
      IF(myid.eq.0) WRITE(*,'(a32,1x,f17.3)')
     & "ARRAYS ALLOCATION TOOK TIME, s =",
     & TIME_1_ARRAY	  
      IF(myid.eq.0) PRINT*,"TOTAL_SIZE_OF_M = ",total_size
      IF(total_size.le.0) STOP "ERROR:ZERO DIM Mij MATRIX"
      IF(.not. calc_matrix_defined) THEN
      IF(MYID.eq.0) PRINT*,"PROGRAM STOPED"	  
      CALL MPI_FINALIZE(ierr_mpi)
      STOP      	  
      ENDIF	  
      conv_unit_r = 1d0
      conv_unit_e = 1d0
      mean_size_1 = total_size
      mean_size_2 = 0
      i_nr_ini = max(ir_bgn_exp_pnt,1)
      i_nr_fin = max(min(ir_fin_exp_pnt,n_r_coll),i_nr_ini)
      IF(ir_fin_exp_pnt.lt.0 .or. ir_bgn_exp_pnt.lt.0) THEN
      i_nr_ini = 1
      i_nr_fin = n_r_coll	  
      ENDIF
      IF(i_nr_ini.gt.n_r_coll .or. i_nr_fin.lt.1) THEN
      IF(myid.eq.0) THEN
      PRINT*,
     & "ERROR: DEFINE PROPERLY GRID RANGE"
      PRINT*,"i_nr_ini",i_nr_ini
      PRINT*,"n_r_coll",n_r_coll	  
      ENDIF	
      STOP	  
      ENDIF	 	  
!      CALL READ_USER_TERMS
      IF(angs_unit) conv_unit_r = a_bohr
      IF(cm_unit)  conv_unit_e =1d0/eVtown/autoeV 	  
      IF(klv_unit)	conv_unit_e =1d0/autokelv/autoeV
      IF(kal_unit) conv_unit_e = 1d0/kcal	  
      massAB_M = mass_red*amutoau
      Ej = Ej/eVtown/autoeV
     	  
      IF(expansion_defined) THEN	  
      IF(terms_file_defined ) THEN
      pot_expansion_defined = .FALSE.
      ELSE
      pot_expansion_defined = .TRUE.	  
      ENDIF
      IF(terms_onfly_defined) THEN
      expansion_defined = .FALSE.
      pot_expansion_defined = .TRUE.	  
      ENDIF	  
      ENDIF		  
      IF(.not.expansion_defined .and. pot_expansion_defined .and.
     & .not. mpi_task_defined ) THEN
      DEALLOCATE(Mat_el,Mat_el_der)
      CALL EXCLUDE_STATES	  
      RETURN
      ELSE
      IF(pot_expansion_defined)expansion_defined = .TRUE.	  
      ENDIF	  
      IF(expansion_defined) THEN
      make_grid_file = .FALSE.	  
      ENDIF	  
      ALLOCATE(Ks_have_to_be_skipped(total_size))	  
      IF(n_r_coll.le.1) STOP "ERROR: TOO FEW POINTS FOR COLL GRID"
! EXPANSION NOT DEFINED,	  
      IF(.not. expansion_defined .or. pot_expansion_defined .or.
     & calc_expansion_defined) THEN
      IF(.not.read_r_grid_from_file_defined) THEN	 
      IF(.not.eq_grid_defined) THEN	  
      IF(MYID.eq.0) PRINT*," MAKING AN ""OPTIMIZED"" GRID	"  
      n_point_dis_1 = int(n_r_coll/3)
      n_point_dis_2 = 	int(n_r_coll/2)+1  
      R_COM(1) = R_min_dist
      R_COM(n_r_coll) = R_max_dist	  
      dR_step = (R_max_dist-R_min_dist)/(n_r_coll-1d0)	  
      DO i=1,n_point_dis_1
      R_COM(i) = R_min_dist+ dR_step*(i-1d0)/3d0    
      ENDDO
      DO i=n_point_dis_1+1,n_point_dis_2
      R_COM(i) = R_COM(n_point_dis_1)+ dR_step*(i-n_point_dis_1)    
      ENDDO
      dR_step = (R_max_dist-R_COM(n_point_dis_2))/
     & (n_r_coll-n_point_dis_2)
      DO i=n_point_dis_2,n_r_coll
      R_COM(i) = R_COM(n_point_dis_2)+ dR_step*(i-n_point_dis_2)    
      ENDDO
      ELSE
      IF(MYID.eq.0) PRINT*," MAKING AN EQUIDIST GRID	" 	  
      dR_step = (R_max_dist-R_min_dist)/(n_r_coll-1d0)	  
      DO i=1,n_r_coll
      R_COM(i) = R_min_dist+ dR_step*(i-1d0)   
      ENDDO	  
      ENDIF	  
      ENDIF
	  IF(read_r_grid_from_file_defined) THEN
      INQUIRE( FILE=R_GRID_DAT, EXIST=file_exst)
      IF(file_exst) THEN	  
      IF(MYID.eq.0) THEN
      OPEN(23,FILE=R_GRID_DAT,STATUS="OLD",ACTION="READ")
      READ(23,*,IOSTAT=ISTAT)	  
      DO i=1,n_r_coll
      READ(23,*,IOSTAT=ISTAT) i_old,R_COM(i)
      ENDDO
      R_COM = R_COM/conv_unit_r	  
      CLOSE(23)
      IF(ISTAT.gt.0) STOP "ERROR:USER DEFINED R_GRID FILE NOT FOUND"
      ENDIF	
      CALL MPI_BCAST(R_COM, n_r_coll, MPI_DOUBLE_PRECISION,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      ELSE
      IF(MYID.eq.0) PRINT*," R_GRID FILE NOT FOUND" 	  
      IF(MYID.eq.0) PRINT*," MAKING AN ""OPTIMIZED"" GRID	"  
      n_point_dis_1 = int(n_r_coll/3)
      n_point_dis_2 = 	int(n_r_coll/2)+1  
      R_COM(1) = R_min_dist
      R_COM(n_r_coll) = R_max_dist	  
      dR_step = (R_max_dist-R_min_dist)/(n_r_coll-1d0)	  
      DO i=1,n_point_dis_1
      R_COM(i) = R_min_dist+ dR_step*(i-1d0)/3d0    
      ENDDO
      DO i=n_point_dis_1+1,n_point_dis_2
      R_COM(i) = R_COM(n_point_dis_1)+ dR_step*(i-n_point_dis_1)    
      ENDDO
      dR_step = (R_max_dist-R_COM(n_point_dis_2))/
     & (n_r_coll-n_point_dis_2)
      DO i=n_point_dis_2,n_r_coll
      R_COM(i) = R_COM(n_point_dis_2)+ dR_step*(i-n_point_dis_2)    
      ENDDO	  
      ENDIF	  
      ELSE
      INQUIRE( FILE=R_GRID_DAT, EXIST=file_exst )	  
      IF (MYID.EQ.0 .and. .not.file_exst) THEN
!     IF GRID DOES NOT EXIST CREATE ITS OWN R-GRID FILE	  
      OPEN(23,FILE=R_GRID_DAT,STATUS="NEW",ACTION="WRITE")	  
      WRITE(23,'(a3,1x,a12)') '#iR','R_COM(iR)'	  
      DO i=1,n_r_coll
      WRITE(23,'(i3,1x,e19.12)') i,R_COM(i)
!      IF(i.ne.i_old) STOP "ERROR IN R_GRID_FILE"	  
      ENDDO
      CLOSE(23)	  
      ENDIF	  
      ENDIF

      ENDIF
!     CREATING/READING V-POTENTIAL FILE ON THE TOTAL ANGLE(S)/R-GRID 	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.le.1))
     & THEN	   
      INQUIRE(FILE=potential_file_name,EXIST=grid_file_found)       
      IF(grid_file_found) THEN
      IF(MYID.eq.0)PRINT*,"GRID WILL BE READ FROM FILE" 
      ELSE
      IF(MYID.eq.0)PRINT*,
     & "GRID WILL BE STORED IN FILE FOR FURTHER USE"
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi ) 	   	   
      CALL INI_POT_STORAGE
!!    CREATING A GRID	  
      IF(MYID.eq.0) PRINT*,"GRID HAD BEEN STORED IN FILE"
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi ) 	   
      ENDIF	   
      ENDIF
!      STOP
!     COMPUTING POTENTIAL EXPANSION	  
      IF(calc_expansion_defined) THEN
	  bk_time_bgn=MPI_Wtime()		!Bikram Jan,'19
      IF(MYID.EQ.0) WRITE(*,*)"EXPANSION WILL BE COMPUTED"
      L_NJU  = 0 	  
      L_NJU(1) = L_MAX_EXPAN
      L_NJU(2) = L1_MAX_EXPAN
      L_NJU(3) = L2_MAX_EXPAN
      L_NJU(4) = NJU_MAX_EXPAN
      L_NJU(5) = NJU1_MAX_EXPAN
      L_NJU(6) = NJU2_MAX_EXPAN
      L_NJU(7) = L1_MIN_EXPAN
      L_NJU(8) = L2_MIN_EXPAN	  
      CALL CALC_EXPANSION(coll_type,n_r_coll,L_NJU)
!     CALCULATING EXPANSION AND EXISITING
      bk_time_fin=MPI_Wtime()		!Bikram Jan,'19
	  bk_time=bk_time_fin-bk_time_bgn		!Bikram Jan,'19
	  if(myid.eq.0)print*,"TIME TOOK TO COMPUTE EXPANSION TERMS = ",
     & bk_time,"SEC"	 !Bikram Jan,'19
!	  if (myid.eq.0) then
!	  INQUIRE(FILE="PROGRESS_EXPAN.tmp", EXIST=file_exst )	
!      IF(file_exst) open(1111,FILE="PROGRESS_EXPAN.tmp")
!      close(1111,status='delete')
!	  endif
      IF(MYID.eq.0)PRINT*,"EXPANSION HAS BEEN COMPUTED,PROGRAM HAS BEEN
     & STOPED"
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi ) 	   
      STOP	  
      ENDIF
!!!!   ASSIGNING THE SIZE OF A CHUNK	  
      chunk_mpi_size = total_size/nproc
      k_st_mpi = myid*chunk_mpi_size + 1
      k_fn_mpi = (myid+1)*chunk_mpi_size 	  
!      IF(myid.eq.nproc-1) k_fn_mpi = total_size
      task_size = (k_fn_mpi-k_st_mpi +1)*n_r_coll
      IF(MYID.eq.0 .and. .not.matrix_reading_defined
     & .and. chunk_mpi_size.gt.0 )
     & WRITE(*,'(a47,1x,i7)')
     & "THE SIZE OF THE Mij COMPUTED BY EACH PROCESSOR=",
     & chunk_mpi_size
      IF(chunk_mpi_size.eq.0 .and. myid.eq.0 ) WRITE(*,*)
     & "NUMBER OF CPUS > SIZE OF Mij"	  
       ! MATRIX READING OR CALCULATINGS
      TIME_MAT_START = MPI_Wtime()	   
      IF(.not.matrix_reading_defined) THEN
!!!!  COMPUTING MATRIX	  
      IF(myid.eq.0) PRINT*,"MATRIX_INI_STARTED"
!      PRINT*, "HELLO"!DELETE	  
      IF(chunk_mpi_size.ge.0) THEN
!! ALLOCATING CHUNK_SIZE ARRAYS
      IF(chunk_mpi_size.gt.0) THEN	  
      ALLOCATE(Mat_el_temp(n_r_coll,k_fn_mpi-k_st_mpi +1)
     &, Mat_el_der_temp(n_r_coll,k_fn_mpi-k_st_mpi +1))
      ALLOCATE(K_SKIPPED_BY_ROUTINE(k_fn_mpi-k_st_mpi +1))
      DO k=k_st_mpi,k_fn_mpi
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.FALSE.	  
      ENDDO
      ENDIF	  
	  
!!!!! ALLOCATION OF SKIPPING ARRAYS
!!!! THIS IS THE CASE WHEN FOR LARGE SYSTEMS IT IS NOT POSSIBLE TO
!! OR PRACTICAL TO STORE IN THE MEMORY THE ENTIRE GRID
	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.eq.0))
     & THEN
      IF(coll_type.eq.1) THEN
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_beta,nr_cold
      ALLOCATE(V_2_1(n_beta,nr_cold))
      READ(1)  V_2_1
      CLOSE(1)
      ENDIF	  
      IF(MYID.EQ.0) THEN
      PRINT*, "MATRIX WILL BE COMPUTED USING V-GRID FILE"	  
      SELECT CASE(coll_type)
      CASE(6)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_beta1,n_beta2,n_gamma1,n_r_vib1,
     & n_r_vib2, nr_cold
      ALLOCATE(V_VIB_2_2(n_gamma1,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1,nr_cold)) 
      CASE(7)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_alpha2,n_beta1,n_beta2,n_gamma1,nr_cold 
      ALLOCATE(
     & V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))		  
      CASE(8)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_alpha2,n_beta1,n_beta2,n_gamma1,nr_cold 
      ALLOCATE(
     & V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      CASE(9)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat)
     &	n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,nr_cold 
      ALLOCATE(
     & V_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  	  
      CASE(0)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat)
     & n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,nr_cold 
      ALLOCATE(
     & V_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      END SELECT
      PRINT*,"VIJ_BUFFER_ALLOCATED"	
!!!!! READS AND PASSE V_ij BUFFER FOR FURTHER BROADCASTING,
!!!! ONLY FOR CASE OF ASYMETRIC + ASYMETRIC	
      IF(istat.gt.0) THEN
      CRITICAL_ERROR=.TRUE.
      PRINT*,"ERROR in V_GRID"	  
      ENDIF	  
      ENDIF
      CALL MPI_BCAST(CRITICAL_ERROR, 1, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(CRITICAL_ERROR) STOP	 
      CALL MPI_BCAST(n_alpha2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_beta1, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_beta2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_gamma1, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_gamma2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_r_vib1, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_r_vib2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
      CALL MPI_BCAST(nr_cold, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(nr_cold.lt.n_r_coll) THEN 
      IF(MYID.EQ.0) THEN
      CRITICAL_ERROR=.TRUE.
      PRINT*,"ERROR in V_GRID"
      ENDIF
      STOP	  
      ENDIF	  
	  
      SELECT CASE(coll_type)
      CASE(6)
      ALLOCATE(V_VIB_2_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1)) 
      CASE(7)
      ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	 	  
      CASE(8)
      ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(9)
      ALLOCATE(
     & V_3_3_int_buffer(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(0)
!!!  MINI BUFFERS FOR EACH R-VALUE TO AVOID MEMORY LIMITATIONS	  
      ALLOCATE(
     & V_3_3_int_buffer(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))
!      PRINT*,"BUFFER_ALLOCATED FOR=",MYID
      END SELECT	 

      DO i=1,n_r_coll
      IF(MYID.EQ.0) THEN	  
      SELECT CASE(coll_type)
      CASE(6)

      READ(1) 
     & V_VIB_2_2_int_buffer
      V_VIB_2_2(:,:,:,:,:,i)=
     & V_VIB_2_2_int_buffer	  

      CASE(7)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_2_int_buffer
      V_3_2(:,:,:,:,i)=
     & V_3_2_int_buffer 	  
      CASE(8)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_2_int_buffer
      V_3_2(:,:,:,:,i)=
     & V_3_2_int_buffer 

      CASE(9)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_3_int_buffer
      V_3_3(:,:,:,:,:,i)=
     & V_3_3_int_buffer 
 
      CASE(0)	  

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_3_int_buffer
      V_3_3(:,:,:,:,:,i)=
     & V_3_3_int_buffer	 
      END SELECT
      IF(i.ge.i_nr_ini .and. i.le. i_nr_fin) THEN	  
	  PRINT*,"BUFFER HAS BEEN READ FOR Ri = ",i
      ENDIF	  
	  
      ENDIF
      SELECT CASE(coll_type)
      CASE(6)
      buffer_size_V = n_r_vib1*n_r_vib2*n_beta1*n_beta2*n_gamma1
      CASE(7)
      buffer_size_V = n_gamma1*n_beta2*n_beta1*n_alpha2		  
      CASE(8)
      buffer_size_V = n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(9)
      buffer_size_V = n_gamma2*n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(0)	  
      buffer_size_V = n_gamma2*n_gamma1*n_beta2*n_beta1*n_alpha2
      END SELECT	  
      IF(MYID.EQ.0) PRINT*,"BUFFER OF VIJ SIZE= ",
     & buffer_size_V!,size(V_3_3_int_buffer)	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
!!!!! BROADCASTING MINI _BUFFER	  
      IF(MYID.EQ.0) PRINT*,"BROADCASTING OF VIJ_BUFFER STARTED"
!      IF(MYID.EQ.0) WRITE(3,*)	V_3_3_int_buffer
      SELECT CASE(coll_type)
      CASE(6)
      CALL MPI_BCAST(V_VIB_2_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CASE(7)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(8)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(9)
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(0)  
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      END SELECT	 
      IF(MYID.EQ.0) PRINT*,"THE BUFFER HAS BEEN BROADCASTED"
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i) !!! CALCULATING VALUE OVER INTEGRATION FOR EACH R
!      PRINT*,st_1,st_2,intgeral	  
      Mat_el_temp(i,k-k_st_mpi+1) = intgeral*conv_unit_e
      IF(ABS(Mat_el_temp(1,k-k_st_mpi+1)).LT.MIJ_ZERO_CUT) THEN    !!! EXCLUDE SOME MATRIX ELEMENTS IF THEY ARE BELOW CERTAIN VALUE
      Mat_el_temp(:,k-k_st_mpi+1) = 0d0
      Mat_el_der_temp(:,k-k_st_mpi+1) = 0d0	  
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.TRUE.	  
      ENDIF	  
!      PRINT*,st_1,st_2,	 intgeral 
      ENDDO
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      ENDDO
      ELSE
!!!! IF NOT THE GRID FORM FILE DEFINED COMPUTE MATRIX ELEMENTS
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "COMPUTING MATRIX ELEMENTS STARTED"	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE !!! SKIP SOME ELEMENTS IF NESSEACRY		  
      DO i=1,n_r_coll
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i)
!!! COMPUTING AND STROING IN THE TEMPORARY ARRAY ASSIGNED FOR EACN ID_PROC	  
      Mat_el_temp(i,k-k_st_mpi+1) = intgeral*conv_unit_e
      IF(ABS(Mat_el_temp(1,k-k_st_mpi+1)).LT.MIJ_ZERO_CUT) THEN !!!! THIS IS EXCLUSION CONDITION
      Mat_el_temp(:,k-k_st_mpi+1) = 0d0
      Mat_el_der_temp(:,k-k_st_mpi+1) = 0d0	  
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.TRUE.	  
      ENDIF	  
      ENDDO
      ENDDO
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      ENDIF	
!!!      SPLININING STARTED
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "SPLINING OF Mij STARTED"	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE
      deriv_bgn = (Mat_el_temp(2,k-k_st_mpi+1)
     & - Mat_el_temp(1,k-k_st_mpi+1))/
     & (R_COM(2)-R_COM(1)) 
      deriv_end = (Mat_el_temp(n_r_coll,k-k_st_mpi+1)
     & - Mat_el_temp(n_r_coll-1,k-k_st_mpi+1))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1)) 	  
      CALL  spline(R_COM,Mat_el_temp(:,k-k_st_mpi+1),
     & n_r_coll,deriv_bgn,deriv_end,
     & Mat_el_der_temp(:,k-k_st_mpi+1))	 
      ENDDO	 
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "SPLINING OF Mij FINISHED"	 
!!! GATHERING MATRIX IN Mij.dat IN ONE PROCESSOR	 

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(chunk_mpi_size.gt.0) THEN	  
      CALL MPI_GATHER(Mat_el_temp,task_size,MPI_DOUBLE_PRECISION,Mat_el,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 	  
      CALL MPI_GATHER(Mat_el_der_temp,task_size,MPI_DOUBLE_PRECISION,
     & Mat_el_der,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      ENDIF	 
      ENDIF		 
!!!! COMPUTING THE REST OF THE MATRIX MIJ
      IF(expansion_defined) THEN
      CALL READ_EXPANSION_TERMS	  
      ENDIF
      IF(total_size.gt.nproc*chunk_mpi_size) THEN
      IF(MYID.EQ.0) PRINT*,"COMPUTING THE RESIDUE OF Mij"
      ALLOCATE(Mat_el_res(n_r_coll),Mat_el_der_res(n_r_coll))
      ALLOCATE(Mat_el_resf(n_r_coll,nproc),
     & Mat_el_der_resf(n_r_coll,nproc))
      Mat_el_res = 0d0
      Mat_el_der_res = 0d0
      DO i=1,n_r_coll
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
      IF(K_SKIPPED_RES) CYCLE
!!! AGAIN IF READING FROM GRID DEFINED!	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.eq.0))
     & THEN
      IF(MYID.EQ.0) THEN
      SELECT CASE(coll_type)
      CASE(6)
      V_VIB_2_2_int_buffer(:,:,:,:,:) = 
     & V_VIB_2_2(:,:,:,:,:,i)
      CASE(7)
      V_3_2_int_buffer(:,:,:,:) = V_3_2(:,:,:,:,i)	  
      CASE(8)
      V_3_2_int_buffer(:,:,:,:) = V_3_2(:,:,:,:,i)	  
      CASE(0)	  
      V_3_3_int_buffer(:,:,:,:,:) = V_3_3(:,:,:,:,:,i)
      END SELECT	  
      ENDIF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      SELECT CASE(coll_type)
      CASE(6)
      CALL MPI_BCAST(V_VIB_2_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CASE(7)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(8)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(9)
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(0)	  
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      PRINT*, "CHUNK BUFFER HAS BEEN BROADCASTED FOR IR=",i	 
      END SELECT	 
      ENDIF	  
!      R_dist =  R_COM(i)*conv_unit_r
      IF(myid.lt.total_size-nproc*chunk_mpi_size) THEN
      k=nproc*chunk_mpi_size+1+myid	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i)
	  
!      PRINT*,"chunk test",intgeral,k,i	  
      Mat_el_res(i) = intgeral*conv_unit_e
      IF(ABS(Mat_el_res(1)).LT.MIJ_ZERO_CUT) THEN
      Mat_el_res(:) = 0d0
      Mat_el_der_res(:) = 0d0	  
      K_SKIPPED_RES=.TRUE.	  
      ENDIF
      ENDIF	  
      ENDDO
!!! SPLINING THE REST OF MIJ
      IF(MYID.eq.0) PRINT*,"SPLINING RESIDUE OF Mij"
      IF(.NOT. K_SKIPPED_RES) THEN
      deriv_bgn = (Mat_el_res(2)
     & - Mat_el_res(1))/
     & (R_COM(2)-R_COM(1)) 
      deriv_end = (Mat_el_res(n_r_coll)
     & - Mat_el_res(n_r_coll-1))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1)) 	  	  
      CALL  spline(R_COM,Mat_el_res,n_r_coll,
     & deriv_bgn,deriv_end,
     & Mat_el_der_res)	  
      ENDIF

!!!!!  GATHERING THE REST	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      CALL MPI_GATHER(Mat_el_res,n_r_coll,MPI_DOUBLE_PRECISION,
     & Mat_el_resf,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 	  
      CALL MPI_GATHER(Mat_el_der_res,n_r_coll,MPI_DOUBLE_PRECISION,
     & Mat_el_der_resf,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
!!!!!  GATHERING THE REST      
      IF(MYID.EQ.0) THEN
      DO i=1,total_size-nproc*chunk_mpi_size
      Mat_el(:,i+nproc*chunk_mpi_size) = Mat_el_resf(:,i)
      Mat_el_der(:,i+nproc*chunk_mpi_size) = Mat_el_der_resf(:,i)	  
      ENDDO	  
      ENDIF	  
      ENDIF
	  TIME_1_MAT = MPI_Wtime()
      IF(MYID.EQ.0 .and. print_matrix_defined) THEN
      IF(.not.unformat_defined) THEN	  
      CALL PRINT_MIJ
      WRITE(*,'(a35,1x,a2,a11,a2)')
     & "MATRIX HAS BEEN SAVED INTO THE FILE"
     & ,"""",MATRIX_NAME_MIJ,""""
      ELSE	 
      CALL PRINT_MIJ_USER
      WRITE(*,'(a35,1x,a2,a11,a2)')
     & "MATRIX HAS BEEN SAVED INTO THE FILE"
     & ,"""",MATRIX_NAME_MIJ_UF,""""	  
      ENDIF	  
      ENDIF
	  
	  TIME_2_MAT = MPI_Wtime()
      TIME_WRITING_MATRIX = TIME_2_MAT - TIME_1_MAT		  
!!! WAITING FOR OTHER PROCESSORS
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(.not.run_prog_defined) THEN
      TIME_MAT_FINISH= MPI_Wtime()
      IF(MYID.EQ.0) WRITE(*,'(a57,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij READING/COMPUTING/SAVING ON DISK"
     & ,",s",(TIME_MAT_FINISH-TIME_MAT_START)/nproc
	  IF(print_matrix_defined.and. myid.eq.0) 
     & WRITE(*,'(a47,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij SAVING ON DISK"
     & ,",s",TIME_WRITING_MATRIX	 
      IF(MYID.EQ.0) PRINT*, "ALL WORK ON MATRIX IS DONE"
      STOP	  
      ENDIF

! Bikram Start Dec 2019:	  
	  if(bikram_mij_shift) then
	  DO i=1,n_r_coll
	  Mat_el(i,:) = Mat_el(i,:) - Mat_el(n_r_coll,:)
      ENDDO 
	  endif
! Bikram End.
	  
	  IF(coupled_states_defined .and. myid.eq.0) THEN
      ALLOCATE(mat_buffer(n_r_coll,total_size))
      ALLOCATE(mat_buffer_der(n_r_coll,total_size))
	  mat_buffer = Mat_el
	  mat_buffer_der = Mat_el_der
      ENDIF		  
      IF(.not.mpi_task_defined) THEN
!!!!! MATRIX BROADCASTING	  
      CALL MPI_BCAST(Mat_el, n_r_coll*total_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CALL MPI_BCAST(Mat_el_der, n_r_coll*total_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(R_COM, n_r_coll, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
!!!!!  CHECKING IF Mij.dat HAS BEEN BROADCASTED CORRECTLY	 
      DO  k=k_st_mpi,k_fn_mpi
      DO i=1,n_r_coll	  
      IF(Mat_el_temp(i,k-k_st_mpi+1).ne.Mat_el(i,k)) THEN
      CALL MPI_FINALIZE (ierr_mpi)
      PRINT*,R_COM
!!!!!! IF NOT STOP THE PROGRAM	  
      STOP "ERROR:MATRIX_INI_WRONG"	  
      ENDIF
      IF(Mat_el_der_temp(i,k-k_st_mpi+1).ne.Mat_el_der(i,k)) THEN
      CALL MPI_FINALIZE (ierr_mpi)
!!!!!! IF NOT STOP THE PROGRAM	  
      STOP "ERROR:MATRIX_INI_WRONG_IN_DER"	  
      ENDIF	  
      ENDDO
      ENDDO 
      ENDIF
      IF(mpi_task_defined) THEN
!!!!! IF MPI_TASKS DEFINED THEN IT WILL REALLOCATE Mij.dat	  
      IF(MYID.eq.0) THEN	  
      ALLOCATE(Mat_el_non_zero(n_r_coll,total_size))
      ALLOCATE(Mat_el_non_zero_der(n_r_coll,total_size))
      Mat_el_non_zero = Mat_el
      Mat_el_non_zero_der = Mat_el_der
      DEALLOCATE(Mat_el,Mat_el_der)	  
      ENDIF	  
!      TIME_MAT_START = MPI_Wtime()	  
      IF(myid.eq.0) PRINT*,"MPI TASK PER TRAJECTORY WILL BE USED"
      IF(myid.eq.0)	
     & WRITE(*,'(a53,1x,i4)')
     & "MPI TASKS WHICH ARE ASSOCIATED WITH ONE TRAJECTORY = ",
     & mpi_task_per_traject	 
      traject_roots = nproc/mpi_task_per_traject
      IF(traject_roots*mpi_task_per_traject.ne.nproc) THEN
      STOP "mpi_task_number must be a delimeter of nproc"	  
      ENDIF	  
      ALLOCATE(mpi_traject_roots(traject_roots))
      ALLOCATE(portion_of_MIJ_per_task(2,nproc))
      ALLOCATE(portion_of_state_per_task(2,nproc))
      ALLOCATE(portion_of_work_per_task(2,nproc))	  
      ALLOCATE(mpi_root_belongs(nproc))	  
!!!!!!!!!!!!!   REMAKE	  
      size_mij_chunk_mpi = total_size/mpi_task_per_traject
      residue_mij_mpi = total_size-
     & size_mij_chunk_mpi*mpi_task_per_traject
	 
      size_state_chunk_mpi = states_size/mpi_task_per_traject
      residue_state_mpi = states_size-
     & size_state_chunk_mpi*mpi_task_per_traject
	 
      size_work_chunk_mpi = (2*states_size+8)/mpi_task_per_traject !!!! DO IN FUTURE
      residue_work_mpi = 2*states_size+8-!!!! DO IN FUTURE
     & size_work_chunk_mpi*mpi_task_per_traject	!!!! DO IN FUTURE
	 
      DO k=1,traject_roots
      mpi_traject_roots(k) = (k-1)*mpi_task_per_traject
      DO k_p=1,mpi_task_per_traject	  
      mpi_root_belongs(k_p+mpi_traject_roots(k))=mpi_traject_roots(k)
      ENDDO	  
      ENDDO
      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_mij_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_MIJ_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_mij_chunk_mpi+1)
      portion_of_MIJ_per_task(2,k_mpi_proc) =
     & (k_p)*(size_mij_chunk_mpi+1)
      ENDDO
      total_size_check = total_size_check + size_mij_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_MIJ_per_task(1,k_mpi_proc)
     & = residue_mij_mpi*(size_mij_chunk_mpi+1)+1  
     & + (k_p-1-residue_mij_mpi)*size_mij_chunk_mpi
      portion_of_MIJ_per_task(2,k_mpi_proc)
     & = residue_mij_mpi*(size_mij_chunk_mpi+1)+	  
     & (k_p-residue_mij_mpi)*size_mij_chunk_mpi
      ENDDO	 
      total_size_check = total_size_check + size_mij_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_MIJ_per_task(2,k_mpi_proc).ne.total_size) 
     & STOP "WRONG MATRIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO



      	  
      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_state_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_state_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_state_chunk_mpi+1)
      portion_of_state_per_task(2,k_mpi_proc) =
     & (k_p)*(size_state_chunk_mpi+1)
      ENDDO
      state_size_check = state_size_check + size_state_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_state_per_task(1,k_mpi_proc)
     & = residue_state_mpi*(size_state_chunk_mpi+1)+1  
     & + (k_p-1-residue_state_mpi)*size_state_chunk_mpi
      portion_of_state_per_task(2,k_mpi_proc)
     & = residue_state_mpi*(size_state_chunk_mpi+1)+	  
     & (k_p-residue_state_mpi)*size_state_chunk_mpi
      ENDDO	 
      state_size_check = state_size_check + size_state_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_state_per_task(2,k_mpi_proc).ne.states_size) 
     & STOP "WRONG MATRIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO	


      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_work_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_work_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_work_chunk_mpi+1)
      portion_of_work_per_task(2,k_mpi_proc) =
     & (k_p)*(size_work_chunk_mpi+1)
      ENDDO
      work_size_check = work_size_check + size_work_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_work_per_task(1,k_mpi_proc)
     & = residue_work_mpi*(size_work_chunk_mpi+1)+1  
     & + (k_p-1-residue_work_mpi)*size_work_chunk_mpi
      portion_of_work_per_task(2,k_mpi_proc)
     & = residue_work_mpi*(size_work_chunk_mpi+1)+	  
     & (k_p-residue_work_mpi)*size_work_chunk_mpi
      ENDDO	 
      work_size_check = work_size_check + size_work_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_work_per_task(2,k_mpi_proc).ne.states_size*2+8) 
     & STOP "WRONG MATRIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO

	  
	  
!      PRINT*,"STATES",portion_of_state_per_task	  
      IF(state_size_check.ne.states_size) THEN
      IF(myid.eq.0) WRITE(*,*)state_size_check,states_size	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF	
	  
      IF(total_size_check.ne.total_size) THEN
      IF(myid.eq.0) WRITE(*,*)total_size_check,total_size	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF
		  
      IF(work_size_check.ne.states_size*2+8) THEN
      IF(myid.eq.0) WRITE(*,*)work_size_check,states_size*2+8	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF	  
	  
      total_size_mpi = portion_of_MIJ_per_task(2,myid+1) - 
     & portion_of_MIJ_per_task(1,myid+1) + 1	  
      ALLOCATE(Mat_el(n_r_coll,total_size_mpi))
      ALLOCATE(Mat_el_der(n_r_coll,total_size_mpi))
      CALL MPI_Comm_group(MPI_COMM_WORLD, wrld_group,ierr_mpi)
      ALLOCATE(process_rank_distr(mpi_task_per_traject,traject_roots))
      ALLOCATE(comms_distr(mpi_task_per_traject),
     & groups_distr(mpi_task_per_traject))   	  
      DO i=1,traject_roots	  
      DO k=1,mpi_task_per_traject
      process_rank_distr(k,i) = k - 1 + (i-1)*mpi_task_per_traject     	  
      ENDDO
      ENDDO	

      DO i=1,mpi_task_per_traject      	  
      CALL MPI_Group_incl(wrld_group, traject_roots,
     & process_rank_distr(i,:), groups_distr(i),ierr_mpi)
      CALL MPI_Comm_create(MPI_COMM_WORLD,groups_distr(i),
     & comms_distr(i),ierr_mpi)
      ENDDO	  
!      PRINT*,myid,total_size_mpi	  
      tag1 = 1
      tag2 = 2
      IF(MYID.eq.0) PRINT*,"Mij ROOT PROC DISTRIBUTION STARTED"	  
      IF(MYID.eq.0) THEN
      DO k=2,mpi_task_per_traject
      total_size_mpi= portion_of_MIJ_per_task(2,k) - 
     & portion_of_MIJ_per_task(1,k) + 1	
      task_portion_size = total_size_mpi*n_r_coll	 
      ALLOCATE(buffer_mpi_portion(n_r_coll,total_size_mpi))
      buffer_mpi_portion = Mat_el_non_zero(:,
     & portion_of_MIJ_per_task(1,k):portion_of_MIJ_per_task(2,k)) 	  
      CALL MPI_SEND(buffer_mpi_portion,
     & task_portion_size, MPI_REAL8, k-1, 
     &  tag1, MPI_COMM_WORLD, ierr_mpi)
      buffer_mpi_portion = 	 Mat_el_non_zero_der (:,
     & portion_of_MIJ_per_task(1,k):portion_of_MIJ_per_task(2,k))
      CALL MPI_SEND(buffer_mpi_portion,
     & task_portion_size, MPI_REAL8, k-1, 
     &  tag2, MPI_COMM_WORLD, ierr_mpi)
      DEALLOCATE(buffer_mpi_portion) 
      ENDDO
      total_size_mpi= portion_of_MIJ_per_task(2,1) - 
     & portion_of_MIJ_per_task(1,1) + 1		  
      Mat_el=Mat_el_non_zero(:,
     & portion_of_MIJ_per_task(1,1):portion_of_MIJ_per_task(2,1))
      Mat_el_der=Mat_el_non_zero_der(:,
     & portion_of_MIJ_per_task(1,1):portion_of_MIJ_per_task(2,1))
      DEALLOCATE(Mat_el_non_zero,Mat_el_non_zero_der)	 
      ELSE
      IF(myid.le.mpi_task_per_traject-1) THEN	  
      task_portion_size = total_size_mpi*n_r_coll	  
      CALL MPI_RECV(Mat_el, task_portion_size, MPI_REAL8, 
     & 0, tag1, MPI_COMM_WORLD, status, ierr_mpi)	  
      CALL MPI_RECV(Mat_el_der, task_portion_size, MPI_REAL8, 
     & 0, tag2, MPI_COMM_WORLD, status, ierr_mpi)
      ENDIF	 
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
      IF(myid.eq.0)	 PRINT*, "Mij ROOT PROC DISTRIBUTION DONE"  
      IF(myid.eq.0) PRINT*,"Mij ALL PROC DISTRIBUTION STARTED" 	  
      DO i=1,mpi_task_per_traject
      IF(i-1 .eq. myid - int(myid/mpi_task_per_traject)
     & *mpi_task_per_traject) THEN
      CALL MPI_Comm_rank(comms_distr(i),id_proc_in_group ,ierr_mpi)
      IF(id_proc_in_group.ne.myid/mpi_task_per_traject) PRINT*,
     & "ERROR IN COMUNICATIONS ASSIGNEMNET_1_distr",id_proc_in_group,
     & i-1
      IF(myid.lt.mpi_task_per_traject .and. id_proc_in_group.ne.0 )
     & PRINT*,"ERROR IN COMUNICATIONS ASSIGNEMNET_2_distr"
      task_portion_size = total_size_mpi*n_r_coll	  
      CALL MPI_BCAST(Mat_el,task_portion_size, MPI_REAL8,0,
     &  comms_distr(i),ierr_mpi)
      CALL MPI_BCAST(Mat_el_der,task_portion_size, MPI_REAL8,0,
     &  comms_distr(i),ierr_mpi)
      ENDIF	 
      ENDDO
      IF(myid.eq.0) PRINT*,"Mij ALL PROC DISTRIBUTION DONE"	  
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)	  
      IF(MYID.EQ.0) PRINT *,"MPI_COMMUNICATORS_CREATION STARTED"
      ALLOCATE(comms(traject_roots),groups(traject_roots))     	  
      ALLOCATE(process_rank(traject_roots,mpi_task_per_traject))		  
      DO i=1,traject_roots	  
      DO k=1,mpi_task_per_traject
      process_rank(i,k) = k - 1 + (i-1)*mpi_task_per_traject     	  
      ENDDO
      ENDDO	

      DO i=1,traject_roots      	  
      CALL MPI_Group_incl(wrld_group, mpi_task_per_traject,
     & process_rank(i,:), groups(i),ierr_mpi)
      CALL MPI_Comm_create(MPI_COMM_WORLD,groups(i),comms(i),ierr_mpi)
      ENDDO

      DO i=1,traject_roots
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN
      CALL MPI_Comm_rank(comms(i),id_proc_in_group ,ierr_mpi)	  
      IF(mpi_traject_roots(i)+id_proc_in_group.ne.myid) STOP
     & "WRONG GROUP MPI ASSIGNEMENT"	  
	  
      ENDIF		  
      ENDDO	 
      IF(MYID.EQ.0) PRINT *,"MPI_COMMUNICATORS_CREATION DONE"
      TIME_MAT_FINISH = MPI_Wtime()
      TIME_2_MAT = TIME_MAT_FINISH - TIME_MAT_START
      CALL MPI_BCAST(R_COM, n_r_coll, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(ind_mat, 2*total_size, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      ENDIF		

      ELSE
!!!   READING Mij FROM A FILE	  
!      TIME_MAT_START = MPI_Wtime()	  
      IF (myid.eq.0) THEN
      IF(.not.unformat_defined) THEN
      CALL READ_MIJ !!!! TESTING
      ELSE 
      CALL READ_MIJ_USER	  
      ENDIF	 
      ENDIF
!      STOP	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      TIME_MAT_FINISH = MPI_Wtime()
      TIME_1_MAT = TIME_MAT_FINISH - TIME_MAT_START	  
      CALL MPI_BCAST(CRITICAL_ERROR, 1, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(CRITICAL_ERROR) THEN
      IF(MYID.EQ.0) PRINT*,"CRITICAL ERROR IN MATRIX READING"
      STOP 	  
      ENDIF	  
      CALL MPI_BCAST(total_size_old, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
!      IF(MYID.EQ.0) PRINT*,"TOTAL_OLD_SIZE",total_size_old
!!!   COMPUTING OF NEW MATRIX BEGINS
!!!        EXPANDING ON R_GRID
      IF(total_size.eq.total_size_old) THEN
      IF(ir_bgn_exp_pnt.gt.0 .and. ir_fin_exp_pnt.gt.0 .and.
     &  .not. (ir_bgn_exp_pnt.eq.1 .and. ir_fin_exp_pnt.eq.
     & n_r_coll ) ) THEN
      IF(myid.eq.0) THEN	  
      ALLOCATE(Mat_el_r_temp(n_r_coll,total_size))
      ALLOCATE(Mat_el_der_r_temp(n_r_coll,total_size))
      ENDIF	  
!!!!  COMPUTING MATRIX	  
      IF(myid.eq.0) PRINT*,"MATRIX_INI__R_ADD_STARTED"
!      PRINT*, "HELLO"!DELETE	  
      IF(chunk_mpi_size.ge.0) THEN
!! ALLOCATING CHUNK_SIZE ARRAYS
      IF(chunk_mpi_size.gt.0) THEN	  
      ALLOCATE(Mat_el_temp(n_r_coll,k_fn_mpi-k_st_mpi +1)
     &, Mat_el_der_temp(n_r_coll,k_fn_mpi-k_st_mpi +1))
      ALLOCATE(K_SKIPPED_BY_ROUTINE(k_fn_mpi-k_st_mpi +1))
      DO k=k_st_mpi,k_fn_mpi
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.FALSE.	  
      ENDDO
      ENDIF	  
	  
!!!!! ALLOCATION OF SKIPPING ARRAYS
!!!! THIS IS THE CASE WHEN FOR LARGE SYSTEMS IT IS NOT POSSIBLE TO
!! OR PRACTICAL TO STORE IN THE MEMORY THE ENTIRE GRID
	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.eq.0))
     & THEN
      IF(MYID.EQ.0) THEN
      PRINT*, "MATRIX WILL BE COMPUTED FROM GRID"	  
      SELECT CASE(coll_type)
      CASE(6)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_beta1,n_beta2,n_gamma1,n_r_vib1,
     & n_r_vib2, nr_cold
      ALLOCATE(V_VIB_2_2(n_gamma1,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1,nr_cold)) 
      CASE(7)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_alpha2,n_beta1,n_beta2,n_gamma1,nr_cold 
      ALLOCATE(
     & V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))		  
      CASE(8)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat) n_alpha2,n_beta1,n_beta2,n_gamma1,nr_cold 
      ALLOCATE(
     & V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      CASE(9)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat)
     &	n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,nr_cold 
      ALLOCATE(
     & V_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  	  
      CASE(0)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1,IOSTAT=istat)
      READ(1,IOSTAT=istat)
     & n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,nr_cold 
      ALLOCATE(
     & V_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      END SELECT
      PRINT*,"VIJ_BUFFER_ALLOCATED"	
!!!!! READS AND PASSE V_ij BUFFER FOR FURTHER BROADCASTING,
!!!! ONLY FOR CASE OF ASYMETRIC + ASYMETRIC	
      IF(istat.gt.0) THEN
      CRITICAL_ERROR=.TRUE.
      PRINT*,"ERROR in V_GRID"	  
      ENDIF	  
      ENDIF
      CALL MPI_BCAST(CRITICAL_ERROR, 1, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(CRITICAL_ERROR) STOP	 
      CALL MPI_BCAST(n_alpha2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_beta1, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_beta2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_gamma1, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_gamma2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_r_vib1, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_r_vib2, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
      CALL MPI_BCAST(nr_cold, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(nr_cold.lt.n_r_coll) THEN 
      IF(MYID.EQ.0) THEN
      CRITICAL_ERROR=.TRUE.
      PRINT*,"ERROR in V_GRID"
      ENDIF
      STOP	  
      ENDIF	  
	  
      SELECT CASE(coll_type)
      CASE(6)
      ALLOCATE(V_VIB_2_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1)) 
      CASE(7)
      ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	 	  
      CASE(8)
      ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(9)
      ALLOCATE(
     & V_3_3_int_buffer(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(0)
!!!  MINI BUFFERS FOR EACH R-VALUE TO AVOID MEMORY LIMITATIONS	  
      ALLOCATE(
     & V_3_3_int_buffer(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))
!      PRINT*,"BUFFER_ALLOCATED FOR=",MYID
      END SELECT	 

      DO i=1,n_r_coll
      IF(MYID.EQ.0) THEN	  
      SELECT CASE(coll_type)
      CASE(6)

      READ(1) 
     & V_VIB_2_2_int_buffer
      V_VIB_2_2(:,:,:,:,:,i)=
     & V_VIB_2_2_int_buffer	  

      CASE(7)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_2_int_buffer
      V_3_2(:,:,:,:,i)=
     & V_3_2_int_buffer 	  
      CASE(8)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_2_int_buffer
      V_3_2(:,:,:,:,i)=
     & V_3_2_int_buffer 

      CASE(9)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_3_int_buffer
      V_3_3(:,:,:,:,:,i)=
     & V_3_3_int_buffer 
 
      CASE(0)	  

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_3_int_buffer
      V_3_3(:,:,:,:,:,i)=
     & V_3_3_int_buffer	 
      END SELECT	  
      PRINT*,"BUFFER HAS BEEN READ FOR Ri = ",i
	  
      ENDIF
      SELECT CASE(coll_type)
      CASE(6)
      buffer_size_V = n_r_vib1*n_r_vib2*n_beta1*n_beta2*n_gamma1
      CASE(7)
      buffer_size_V = n_gamma1*n_beta2*n_beta1*n_alpha2		  
      CASE(8)
      buffer_size_V = n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(9)
      buffer_size_V = n_gamma2*n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(0)	  
      buffer_size_V = n_gamma2*n_gamma1*n_beta2*n_beta1*n_alpha2
      END SELECT	  
      IF(MYID.EQ.0) PRINT*,"BUFFER OF VIJ SIZE= ",
     & buffer_size_V!,size(V_3_3_int_buffer)	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
!!!!! BROADCASTING MINI _BUFFER	  
      IF(MYID.EQ.0) PRINT*,"BROADCASTING OF VIJ_BUFFER STARTED"
!      IF(MYID.EQ.0) WRITE(3,*)	V_3_3_int_buffer
      SELECT CASE(coll_type)
      CASE(6)
      CALL MPI_BCAST(V_VIB_2_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CASE(7)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(8)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(9)
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(0)  
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      END SELECT	 
      IF(MYID.EQ.0) PRINT*,"THE BUFFER HAS BEEN BROADCASTED"
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i) !!! CALCULATING VALUE OVER INTEGRATION FOR EACH R
!      PRINT*,st_1,st_2,intgeral	  
      Mat_el_temp(i,k-k_st_mpi+1) = intgeral*conv_unit_e
      IF(ABS(Mat_el_temp(1,k-k_st_mpi+1)).LT.MIJ_ZERO_CUT) THEN    !!! EXCLUDE SOME MATRIX ELEMENTS IF THEY ARE BELOW CERTAIN VALUE
      Mat_el_temp(:,k-k_st_mpi+1) = 0d0
      Mat_el_der_temp(:,k-k_st_mpi+1) = 0d0	  
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.TRUE.	  
      ENDIF	  
!      PRINT*,st_1,st_2,	 intgeral 
      ENDDO
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      ENDDO
      ELSE
!!!! IF NOT THE GRID FORM FILE DEFINED COMPUTE MATRIX ELEMENTS
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "COMPUTING MATRIX ELEMENTS STARTED"	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE !!! SKIP SOME ELEMENTS IF NESSEACRY		  
      DO i=1,n_r_coll 
      IF(i.gt.i_nr_fin .or. i.lt.i_nr_ini) CYCLE	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i)
!!! COMPUTING AND STROING IN THE TEMPORARY ARRAY ASSIGNED FOR EACN ID_PROC	  
      Mat_el_temp(i,k-k_st_mpi+1) = intgeral*conv_unit_e
      IF(ABS(Mat_el_temp(1,k-k_st_mpi+1)).LT.MIJ_ZERO_CUT) THEN !!!! THIS IS EXCLUSION CONDITION
      Mat_el_temp(:,k-k_st_mpi+1) = 0d0
      Mat_el_der_temp(:,k-k_st_mpi+1) = 0d0	  
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.TRUE.	  
      ENDIF	  
      ENDDO
      ENDDO
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      ENDIF	
!!!      SPLININING STARTED
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "SPLINING OF Mij STARTED"	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE
      deriv_bgn = (Mat_el_temp(2,k-k_st_mpi+1)
     & - Mat_el_temp(1,k-k_st_mpi+1))/
     & (R_COM(2)-R_COM(1)) 
      deriv_end = (Mat_el_temp(n_r_coll,k-k_st_mpi+1)
     & - Mat_el_temp(n_r_coll-1,k-k_st_mpi+1))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1)) 	  
      CALL  spline(R_COM,Mat_el_temp(:,k-k_st_mpi+1),
     & n_r_coll,deriv_bgn,deriv_end,
     & Mat_el_der_temp(:,k-k_st_mpi+1))	 
      ENDDO	 
      IF(MYID.EQ.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "SPLINING OF Mij FINISHED"	 
!!! GATHERING MATRIX IN Mij.dat IN ONE PROCESSOR	 

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(chunk_mpi_size.gt.0) THEN	  
      CALL MPI_GATHER(Mat_el_temp,task_size,MPI_DOUBLE_PRECISION,
     & Mat_el_r_temp,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 	  
      CALL MPI_GATHER(Mat_el_der_temp,task_size,MPI_DOUBLE_PRECISION,
     & Mat_el_der_r_temp,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      ENDIF	 
      ENDIF
      IF(MYID.eq.0) THEN	  
      Mat_el(i_nr_ini:i_nr_fin,:) = Mat_el_r_temp(i_nr_ini:i_nr_fin,:)
      Mat_el_der(i_nr_ini:i_nr_fin,:) = 
     & Mat_el_der_r_temp(i_nr_ini:i_nr_fin,:) 
      DEALLOCATE(Mat_el_r_temp,Mat_el_der_r_temp)
      ENDIF	  
!!!! COMPUTING THE REST OF THE MATRIX MIJ
      IF(expansion_defined) THEN
      CALL READ_EXPANSION_TERMS	  
      ENDIF
      IF(total_size.gt.nproc*chunk_mpi_size) THEN
      IF(MYID.EQ.0) PRINT*,"COMPUTING THE RESIDUE OF Mij"
      ALLOCATE(Mat_el_res(n_r_coll),Mat_el_der_res(n_r_coll))
      ALLOCATE(Mat_el_resf(n_r_coll,nproc),
     & Mat_el_der_resf(n_r_coll,nproc))
      Mat_el_res = 0d0
      Mat_el_der_res = 0d0
      DO i=1,n_r_coll
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
      IF(K_SKIPPED_RES) CYCLE
!!! AGAIN IF READING FROM GRID DEFINED!	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.eq.0))
     & THEN
      IF(MYID.EQ.0) THEN
      SELECT CASE(coll_type)
      CASE(6)
      V_VIB_2_2_int_buffer(:,:,:,:,:) = 
     & V_VIB_2_2(:,:,:,:,:,i)
      CASE(7)
      V_3_2_int_buffer(:,:,:,:) = V_3_2(:,:,:,:,i)	  
      CASE(8)
      V_3_2_int_buffer(:,:,:,:) = V_3_2(:,:,:,:,i)	  
      CASE(0)	  
      V_3_3_int_buffer(:,:,:,:,:) = V_3_3(:,:,:,:,:,i)
      END SELECT	  
      ENDIF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      SELECT CASE(coll_type)
      CASE(6)
      CALL MPI_BCAST(V_VIB_2_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CASE(7)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(8)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(9)
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(0)	  
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      PRINT*, "CHUNK BUFFER HAS BEEN BROADCASTED FOR IR=",i	 
      END SELECT	 
      ENDIF	  
!      R_dist =  R_COM(i)*conv_unit_r
      IF(myid.lt.total_size-nproc*chunk_mpi_size) THEN
      k=nproc*chunk_mpi_size+1+myid	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i)
	  
!      PRINT*,"chunk test",intgeral,k,i	  
      Mat_el_res(i) = intgeral*conv_unit_e
      IF(ABS(Mat_el_res(1)).LT.MIJ_ZERO_CUT) THEN
      Mat_el_res(:) = 0d0
      Mat_el_der_res(:) = 0d0	  
      K_SKIPPED_RES=.TRUE.	  
      ENDIF
      ENDIF	  
      ENDDO
!!! SPLINING THE REST OF MIJ
      IF(MYID.eq.0) PRINT*,"SPLINING RESIDUE OF Mij"
      IF(.NOT. K_SKIPPED_RES) THEN
      deriv_bgn = (Mat_el_res(2)
     & - Mat_el_res(1))/
     & (R_COM(2)-R_COM(1)) 
      deriv_end = (Mat_el_res(n_r_coll)
     & - Mat_el_res(n_r_coll-1))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1)) 	  	  
      CALL  spline(R_COM,Mat_el_res,n_r_coll,
     & deriv_bgn,deriv_end,
     & Mat_el_der_res)	  
      ENDIF

!!!!!  GATHERING THE REST	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      CALL MPI_GATHER(Mat_el_res,n_r_coll,MPI_DOUBLE_PRECISION,
     & Mat_el_resf,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 	  
      CALL MPI_GATHER(Mat_el_der_res,n_r_coll,MPI_DOUBLE_PRECISION,
     & Mat_el_der_resf,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
!!!!!  GATHERING THE REST      
      IF(MYID.EQ.0) THEN
      DO i=1,total_size-nproc*chunk_mpi_size
      Mat_el(i_nr_ini:i_nr_fin,i+nproc*chunk_mpi_size)
     & = Mat_el_resf(i_nr_ini:i_nr_fin,i)
      Mat_el_der(i_nr_ini:i_nr_fin,i+nproc*chunk_mpi_size) =
     & Mat_el_der_resf(i_nr_ini:i_nr_fin,i)	  
      ENDDO	  
      ENDIF	  
      ENDIF	 
	  TIME_1_MAT = MPI_Wtime() 
      IF(MYID.EQ.0 .and. print_matrix_defined) THEN
      IF(.not.unformat_defined) THEN	  
      CALL PRINT_MIJ
      PRINT*,
     &  "MATRIX HAS BEEN SAVED INTO THE FILE ""MTRX.DAT"" "	  
      ELSE	  
      CALL PRINT_MIJ_USER
      PRINT*,"MATRIX HAS BEEN SAVED INTO THE FILE ""MTRX_UF.DAT"" "	  
      ENDIF		  
      ENDIF		 
	  TIME_2_MAT = MPI_Wtime()
      TIME_WRITING_MATRIX = TIME_2_MAT - TIME_1_MAT	  
	 
! Bikram Start Dec 2019:	  
	  if(bikram_mij_shift) then
	  DO i=1,n_r_coll
	  Mat_el(i,:) = Mat_el(i,:) - Mat_el(n_r_coll,:)
      ENDDO 
	  endif
! Bikram End.

      ENDIF
      ENDIF  
      IF(total_size.gt.total_size_old) THEN
      ALLOCATE(Mat_rest(n_r_coll,total_size-total_size_old))
      ALLOCATE(Mat_rest_der(n_r_coll,total_size-total_size_old))	  
      chunk_mpi_size = (total_size-total_size_old)/nproc
!      IF(MYID.EQ.0) PRINT*,"CHUNK SIZE",chunk_mpi_size	  
      k_st_mpi = myid*chunk_mpi_size + 1+total_size_old
      k_fn_mpi = (myid+1)*chunk_mpi_size+total_size_old
!      PRINT*,"INTEGRATION RANGE",myid,k_st_mpi,k_fn_mpi 	  
!      IF(myid.eq.nproc-1) k_fn_mpi = total_size
      task_size = (k_fn_mpi-k_st_mpi +1)*n_r_coll	  
      IF(myid.eq.0) PRINT*,"MATRIX_ADD_STARTED"
      IF(chunk_mpi_size.ge.0) THEN
!      PRINT*, "HELLO"!DELETE
      IF(chunk_mpi_size.gt.0)	THEN  
      ALLOCATE(Mat_el_temp(n_r_coll,k_fn_mpi-k_st_mpi +1)
     &, Mat_el_der_temp(n_r_coll,k_fn_mpi-k_st_mpi +1))
      ALLOCATE(K_SKIPPED_BY_ROUTINE(k_fn_mpi-k_st_mpi +1))
      DO k=k_st_mpi,k_fn_mpi
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.FALSE.	  
      ENDDO
      ENDIF	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.eq.0))
     & THEN
      IF(MYID.EQ.0) THEN
      PRINT*, "MATRIX WILL BE COMPUTED FROM GRID"	  
      SELECT CASE(coll_type)
      CASE(6)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1)
      READ(1) n_beta1,n_beta2,n_gamma1,n_r_vib1,
     & n_r_vib2, nr_cold
      ALLOCATE(V_VIB_2_2(n_gamma1,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1,nr_cold)) 
      CASE(7)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1)
      READ(1) n_alpha2,n_beta1,n_beta2,n_gamma1,nr_cold 
      ALLOCATE(
     & V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      CASE(8)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1)
      READ(1) n_alpha2,n_beta1,n_beta2,n_gamma1,nr_cold 
      ALLOCATE(
     & V_3_2(n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      CASE(9)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1)
      READ(1) n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,nr_cold 
      ALLOCATE(
     & V_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  	  
      CASE(0)
      OPEN(1,FILE=potential_file_name,
     & FORM="UNFORMATTED",STATUS="OLD",ACTION="READ")
      READ(1)
      READ(1) n_alpha2,n_beta1,n_beta2,n_gamma1,n_gamma2,nr_cold 
      ALLOCATE(
     & V_3_3(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2,nr_cold))	  
      END SELECT
      PRINT*,"VIJ_BUFFER_ALLOCATED"	
!!!!! READS AND PASSE V_ij BUFFER FOR FURTHER BROADCASTING,
!!!! ONLY FOR CASE OF ASYMETRIC + ASYMETRIC	  
      ENDIF
      CALL MPI_BCAST(n_alpha2, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_beta1, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_beta2, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_gamma1, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_gamma2, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_r_vib1, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(n_r_vib2, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	 
      CALL MPI_BCAST(nr_cold, 1, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(nr_cold.lt.n_r_coll) THEN 
      IF(MYID.EQ.0)
     & PRINT*,"ERROR: WRONG FILE V_IJ GRID"
!!!! BROADCASTING ENDED NOW CHECKING	  
      STOP
      ENDIF
      SELECT CASE(coll_type)
      CASE(6)
      ALLOCATE(V_VIB_2_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_r_vib2,
     & n_r_vib1)) 
      CASE(7)
      ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(8)
      ALLOCATE(V_3_2_int_buffer(n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(9)
      ALLOCATE(
     & V_3_3_int_buffer(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))	  
      CASE(0)
!!!  MINI BUFFERS FOR EACH R-VALUE TO AVOID MEMORY LIMITATIONS	  
      ALLOCATE(
     & V_3_3_int_buffer(n_gamma2,n_gamma1,n_beta2,n_beta1,n_alpha2))
!      PRINT*,"BUFFER_ALLOCATED FOR=",MYID
      END SELECT	 

      DO i=1,n_r_coll
      IF(MYID.EQ.0) THEN	  
      SELECT CASE(coll_type)
      CASE(6)

      READ(1) 
     & V_VIB_2_2_int_buffer
      V_VIB_2_2(:,:,:,:,:,i)=
     & V_VIB_2_2_int_buffer	  

      CASE(7)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_2_int_buffer
      V_3_2(:,:,:,:,i)=
     & V_3_2_int_buffer	 	 
  
      CASE(8)
 
!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_2_int_buffer
      V_3_2(:,:,:,:,i)=
     & V_3_2_int_buffer	 

      CASE(9)

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_3_int_buffer
      V_3_3(:,:,:,:,:,i)=
     & V_3_3_int_buffer		 
	  
      CASE(0)	  

!!!!!!    READING MINI BUFFER	  
      READ(1) 
     & V_3_3_int_buffer
      V_3_3(:,:,:,:,:,i)=
     & V_3_3_int_buffer	 
      END SELECT	  
      PRINT*,"BUFFER HAS BEEN READ FOR Ri = ",i
	  
      ENDIF
      SELECT CASE(coll_type)
      CASE(6)
      buffer_size_V = n_r_vib1*n_r_vib2*n_beta1*n_beta2*n_gamma1
      CASE(7)
      buffer_size_V = n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(8)
      buffer_size_V = n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(9)
      buffer_size_V = n_gamma2*n_gamma1*n_beta2*n_beta1*n_alpha2	  
      CASE(0)	  
      buffer_size_V = n_gamma2*n_gamma1*n_beta2*n_beta1*n_alpha2
      END SELECT	  
      IF(MYID.EQ.0) PRINT*,"BUFFER OF VIJ SIZE= ",
     & buffer_size_V!,size(V_3_3_int_buffer)	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
!!!!! BROADCASTING MINI _BUFFER	  
      IF(MYID.EQ.0) PRINT*,"BROADCASTING OF VIJ_BUFFER STARTED"
!      IF(MYID.EQ.0) WRITE(3,*)	V_3_3_int_buffer
      SELECT CASE(coll_type)
      CASE(6)
      CALL MPI_BCAST(V_VIB_2_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CASE(7)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(8)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(9)
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(0)  
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      END SELECT	 
      IF(MYID.EQ.0) PRINT*,"THE BUFFER HAS BEEN BROADCASTED"
      IF(MYID.eq.0 .and. chunk_mpi_size.gt.0 .and. i.eq.1 )
     & PRINT*, "ADD_MATRIX IS BEING COMPUTED ON V_GRID"	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i) !!! CALCULATING VALUE OVER INTEGRATION FOR EACH R
!      PRINT*,st_1,st_2,intgeral	  
      Mat_el_temp(i,k-k_st_mpi+1) = intgeral*conv_unit_e
      IF(ABS(Mat_el_temp(1,k-k_st_mpi+1)).LT.MIJ_ZERO_CUT) THEN    !!! EXCLUDE SOME MATRIX ELEMENTS IF THEY ARE BELOW CERTAIN VALUE
      Mat_el_temp(:,k-k_st_mpi+1) = 0d0
      Mat_el_der_temp(:,k-k_st_mpi+1) = 0d0	  
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.TRUE.	  
      ENDIF	  
!      PRINT*,st_1,st_2,	 intgeral 
      ENDDO
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	  
      ENDDO	  
      ELSE	  
!!!   REGULAR COMPUTING
      IF(MYID.eq.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "ADD_MATRIX IS BEING COMPUTED"	  
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE 	  
      DO i=1,n_r_coll
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
	  
c      R_dist =  R_COM(i)*conv_unit_r
c      PRINT*,st_1,st_2,R_dist
c      STOP	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i)
c      PRINT*,st_1,st_2,intgeral	  
      Mat_el_temp(i,k-k_st_mpi+1) = intgeral*conv_unit_e
      IF(ABS(Mat_el_temp(1,k-k_st_mpi+1)).LT.MIJ_ZERO_CUT) THEN
      Mat_el_temp(:,k-k_st_mpi+1) = 0d0
      Mat_el_der_temp(:,k-k_st_mpi+1) = 0d0	  
      K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)=.TRUE.	  
      ENDIF	  
c      PRINT*,st_1,st_2,	 intgeral 
      ENDDO
  
      ENDDO
      ENDIF	  
!!!   END REGULAR COMPUTING
!! SPLINING	 
      IF(MYID.eq.0 .and. chunk_mpi_size.gt.0)
     & PRINT*, "ADD_MATRIX IS BEING SPLINED" 
      DO  k=k_st_mpi,k_fn_mpi
      IF(K_SKIPPED_BY_ROUTINE(k-k_st_mpi +1)) CYCLE
      deriv_bgn = (Mat_el_temp(2,k-k_st_mpi+1)
     & - Mat_el_temp(1,k-k_st_mpi+1))/
     & (R_COM(2)-R_COM(1)) 
      deriv_end = (Mat_el_temp(n_r_coll,k-k_st_mpi+1)
     & - Mat_el_temp(n_r_coll-1,k-k_st_mpi+1))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1))	  
      CALL  spline(R_COM,Mat_el_temp(:,k-k_st_mpi+1),
     & n_r_coll,deriv_bgn,deriv_end,
     & Mat_el_der_temp(:,k-k_st_mpi+1))		  
      ENDDO	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      IF(chunk_mpi_size.gt.0) THEN	  
      CALL MPI_GATHER(Mat_el_temp,task_size,MPI_DOUBLE_PRECISION,
     & Mat_rest,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 	  
      CALL MPI_GATHER(Mat_el_der_temp,task_size,MPI_DOUBLE_PRECISION,
     & Mat_rest_der,
     &	  task_size,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      ENDIF	 
      ENDIF		 
!!!   COMPUTING CHUNKS RESIDUE
      IF(expansion_defined) THEN
      CALL READ_EXPANSION_TERMS	  
      ENDIF
      IF(total_size-total_size_old
     & .gt.nproc*chunk_mpi_size) THEN
      IF(myid.eq.0) PRINT*,"REST OF THE CHUNK COMPUTING STARTED"
      IF(myid.eq.0) PRINT*,"THE SIZE OF CHUNK PER EACH PROCESSOR= ",
     & chunk_mpi_size	  
      ALLOCATE(Mat_el_res(n_r_coll),Mat_el_der_res(n_r_coll))
      ALLOCATE(Mat_el_resf(n_r_coll,nproc),
     & Mat_el_der_resf(n_r_coll,nproc))
      Mat_el_res = 0d0
      Mat_el_der_res = 0d0	  
      DO i=1,n_r_coll
      IF(i.lt.i_nr_ini .or. i.gt. i_nr_fin) CYCLE	  
      IF(make_grid_file .and. (coll_type.gt.5 .or. coll_type.eq.0))
     & THEN
      IF(MYID.EQ.0) THEN
      SELECT CASE(coll_type)
      CASE(6)
      V_VIB_2_2_int_buffer(:,:,:,:,:) = 
     & V_VIB_2_2(:,:,:,:,:,i)
      CASE(7)
      V_3_2_int_buffer(:,:,:,:) = V_3_2(:,:,:,:,i)	  
      CASE(8)
      V_3_2_int_buffer(:,:,:,:) = V_3_2(:,:,:,:,i)	  
      CASE(0)	  
      V_3_3_int_buffer(:,:,:,:,:) = V_3_3(:,:,:,:,:,i)
      END SELECT	  
      ENDIF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      SELECT CASE(coll_type)
      CASE(6)
      CALL MPI_BCAST(V_VIB_2_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CASE(7)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(8)
      CALL MPI_BCAST(V_3_2_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(9)
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CASE(0)	  
      CALL MPI_BCAST(V_3_3_int_buffer, buffer_size_V, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      PRINT*, "CHUNK BUFFER HAS BEEN BROADCASTED FOR IR=",i	 
      END SELECT	 
      ENDIF	  	  
      IF(K_SKIPPED_RES) CYCLE	  
c      R_dist =  R_COM(i)*conv_unit_r
      IF(myid.lt.total_size-total_size_old
     & -nproc*chunk_mpi_size) THEN
      k=nproc*chunk_mpi_size+1+myid+total_size_old 	  
      CALL	INTEGRATOR_TEMP(intgeral,k,i)
      Mat_el_res(i) = intgeral*conv_unit_e
      IF(ABS(Mat_el_res(1)).LT.MIJ_ZERO_CUT) THEN
      Mat_el_res(:) = 0d0
      Mat_el_der_res(:) = 0d0	  
      K_SKIPPED_RES=.TRUE.	  
      ENDIF
      ENDIF	  
      ENDDO
c      PRINT*,st_1,st_2,intgeral	  
      IF(.NOT.K_SKIPPED_RES) THEN
      deriv_bgn = (Mat_el_res(2)
     & - Mat_el_res(1))/
     & (R_COM(2)-R_COM(1)) 
      deriv_end = (Mat_el_res(n_r_coll)
     & - Mat_el_res(n_r_coll-1))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1))	 	  
      CALL  spline(R_COM,Mat_el_res,n_r_coll,
     & deriv_bgn,deriv_end,
     & Mat_el_der_res)	  
      ENDIF

      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
!      STOP		  
!!!!!!!!!!!!! SEND RECIEVE
      CALL MPI_GATHER(Mat_el_res,n_r_coll,MPI_DOUBLE_PRECISION,
     & Mat_el_resf,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi) 	  
      CALL MPI_GATHER(Mat_el_der_res,n_r_coll,MPI_DOUBLE_PRECISION,
     & Mat_el_der_resf,
     &	  n_r_coll,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      IF(MYID.EQ.0) THEN
      DO i=1,total_size-total_size_old-nproc*chunk_mpi_size
      Mat_rest(:,i+nproc*chunk_mpi_size) = Mat_el_resf(:,i)
      Mat_rest_der(:,i+nproc*chunk_mpi_size) = Mat_el_der_resf(:,i)	  
      ENDDO	  
      ENDIF	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )	 
	 
!!!!!!!!! SEND RECIEVE	 
	  
      ENDIF
      IF(MYID.EQ.0) THEN	  
      DO i=total_size_old+1,total_size
      Mat_el(:,i) = Mat_rest(:,i-total_size_old)
      Mat_el_der(:,i) = Mat_rest_der(:,i-total_size_old)	  
      ENDDO
      ENDIF	  
	  
	  
      IF(MYID.EQ.0 .and. print_matrix_defined) THEN
      IF(.not.unformat_defined) THEN	  
      CALL PRINT_MIJ
      PRINT*,
     &  "MATRIX HAS BEEN SAVED INTO THE FILE ""MTRX.DAT"" "	  
      ELSE	  
      CALL PRINT_MIJ_USER
      PRINT*,"MATRIX HAS BEEN SAVED INTO THE FILE ""MTRX_UF.DAT"" "	  
      ENDIF		  
      ENDIF	  
!      IF(MYID.EQ.0) CALL PRINT_MREST
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)	  
	  
! Bikram Start Dec 2019:	  	  
	  if(bikram_mij_shift) then
	  DO i=1,n_r_coll
	  Mat_el(i,:) = Mat_el(i,:) - Mat_el(n_r_coll,:)
      ENDDO 
	  endif
! Bikram End.
	  
      ENDIF
	  IF(run_prog_defined .and. 
     & coupled_states_defined .and. myid.eq.0) THEN
      ALLOCATE(mat_buffer(n_r_coll,total_size))
      ALLOCATE(mat_buffer_der(n_r_coll,total_size))
	  mat_buffer = Mat_el
	  mat_buffer_der = Mat_el_der
      ENDIF	
!      TIME_MAT_START = MPI_Wtime()	 
      IF(MYID.EQ.0) THEN
      DO k=1,total_size 
      IF(ABS(Mat_el(1,k)).lt.MIJ_ZERO_CUT) THEN
      Ks_have_to_be_skipped(k) = .TRUE.
      mean_size_1 = mean_size_1 - 1	  
!   	  IF(ind_mat(1,k).eq.ind_mat(2,k)) THEN
 !     PRINT*, "ERROR GLOBAL IN MIJ",k
!      run_prog_defined = .FALSE.	  
!      CALL MPI_BCAST(run_prog_defined,1, MPI_LOGICAL,0,
!     &  MPI_COMM_WORLD,ierr_mpi)	  
!      ENDIF	  
	  
      ELSE
      Ks_have_to_be_skipped(k) = .FALSE.     
      ENDIF	 
      ENDDO	 
      ENDIF	
      IF(MYID.EQ.-1) THEN
      OPEN(345,FILE="MAT_CHECK.out")	  
      DO i=1,n_r_coll	  
      WRITE(345,'(e19.12,1x,e19.12)')R_COM(i),Mat_el(i,1)
      ENDDO
      CLOSE(345)	  
      ENDIF
      IF(MYID.EQ.0)WRITE(*,'(a35,1x,i9)')
     & "THE SIZE OF NON_ZERO PART OF Mij = ",mean_size_1	  
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
      IF(.not. run_prog_defined) THEN
      TIME_MAT_FINISH= MPI_Wtime()
      IF(MYID.EQ.0) WRITE(*,'(a57,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij READING/COMPUTING/SAVING ON DISK"
     & ,",s",(TIME_MAT_FINISH-TIME_MAT_START)/nproc
	  IF(print_matrix_defined.and. myid.eq.0) 
     & WRITE(*,'(a47,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij SAVING ON DISK"
     & ,",s",TIME_WRITING_MATRIX	 
      IF(MYID.eq.0)PRINT*, "ALL WORK ON MATRIX IS DONE"	  
      CALL MPI_FINALIZE (ierr_mpi)	  
      STOP
      ELSE
      IF(MYID.EQ.0) THEN	  
      ALLOCATE(Mat_el_non_zero(n_r_coll,mean_size_1))
      ALLOCATE(Mat_el_non_zero_der(n_r_coll,mean_size_1))	  
      ALLOCATE(ind_mat_non_zero(2,mean_size_1))
      k_non_zero = 0	  
      DO k=1,total_size
      IF(Ks_have_to_be_skipped(k)) CYCLE
      k_non_zero = k_non_zero + 1
      ind_mat_non_zero(:,k_non_zero) = ind_mat(:,k)
      Mat_el_non_zero(:,k_non_zero) = Mat_el(:,k)
      Mat_el_non_zero_der(:,k_non_zero) = Mat_el_der(:,k)	  
      ENDDO 
      IF(k_non_zero.ne.mean_size_1) CRITICAL_ERROR = .TRUE.
      total_size = mean_size_1	  
      DEALLOCATE(Mat_el)
      DEALLOCATE(Mat_el_der)	  
      DEALLOCATE(ind_mat)
      IF(.NOT.mpi_task_defined) THEN	  
      ALLOCATE(Mat_el(n_r_coll,total_size))	  
      ALLOCATE(Mat_el_der(n_r_coll,total_size))
      ENDIF	  
      ALLOCATE(ind_mat(2,total_size))
      IF(.NOT.mpi_task_defined) THEN  
      Mat_el = Mat_el_non_zero
      Mat_el_der = Mat_el_non_zero_der
      ENDIF	  
      ind_mat = ind_mat_non_zero
      	  
      IF(.NOT.mpi_task_defined) THEN
      DEALLOCATE(Mat_el_non_zero,Mat_el_non_zero_der,ind_mat_non_zero)
      ENDIF 	  
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)	  
	  
      CALL MPI_BCAST(total_size, 1, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  
      CALL MPI_BCAST(CRITICAL_ERROR, 1, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      IF(CRITICAL_ERROR) STOP "CRITICAL ERROR"
      IF(MYID.NE.0) THEN
      DEALLOCATE(ind_mat)
      ALLOCATE(ind_mat(2,total_size))
      IF(.NOT.mpi_task_defined) THEN	  
      DEALLOCATE(Mat_el)
      DEALLOCATE(Mat_el_der)
      ALLOCATE(Mat_el(n_r_coll,total_size))	  
      ALLOCATE(Mat_el_der(n_r_coll,total_size))
      ENDIF
      ENDIF	  
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)	  
      TIME_MAT_FINISH = MPI_Wtime()
      TIME_3_MAT = TIME_MAT_FINISH - TIME_MAT_START 	  
      IF(mpi_task_defined) THEN
!      TIME_MAT_START = MPI_Wtime()	  
      IF(myid.eq.0) PRINT*,"MPI TASK PER TRAJECTORY WILL BE USED"
      IF(myid.eq.0)	
     & WRITE(*,'(a53,1x,i4)')
     & "MPI TASKS WHICH ARE ASSOCIATED WITH ONE TRAJECTORY = ",
     & mpi_task_per_traject	 
      traject_roots = nproc/mpi_task_per_traject
      IF(traject_roots*mpi_task_per_traject.ne.nproc) THEN
      STOP "mpi_task_number must be a delimeter of nproc"	  
      ENDIF	  
      ALLOCATE(mpi_traject_roots(traject_roots))
      ALLOCATE(portion_of_MIJ_per_task(2,nproc))
      ALLOCATE(portion_of_state_per_task(2,nproc))
      ALLOCATE(portion_of_work_per_task(2,nproc))	  
      ALLOCATE(mpi_root_belongs(nproc))	  
!!!!!!!!!!!!!   REMAKE	  
      size_mij_chunk_mpi = total_size/mpi_task_per_traject
      residue_mij_mpi = total_size-
     & size_mij_chunk_mpi*mpi_task_per_traject
      size_state_chunk_mpi = states_size/mpi_task_per_traject
      residue_state_mpi = states_size-
     & size_state_chunk_mpi*mpi_task_per_traject
      size_work_chunk_mpi = (2*states_size+8)/mpi_task_per_traject !!!! DO IN FUTURE
      residue_work_mpi = 2*states_size+8-!!!! DO IN FUTURE
     & size_work_chunk_mpi*mpi_task_per_traject	!!!! DO IN FUTURE
	 
      DO k=1,traject_roots
      mpi_traject_roots(k) = (k-1)*mpi_task_per_traject
      DO k_p=1,mpi_task_per_traject	  
      mpi_root_belongs(k_p+mpi_traject_roots(k))=mpi_traject_roots(k)
      ENDDO	  
      ENDDO
      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_mij_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_MIJ_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_mij_chunk_mpi+1)
      portion_of_MIJ_per_task(2,k_mpi_proc) =
     & (k_p)*(size_mij_chunk_mpi+1)
      ENDDO
      total_size_check = total_size_check + size_mij_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_MIJ_per_task(1,k_mpi_proc)
     & = residue_mij_mpi*(size_mij_chunk_mpi+1)+1  
     & + (k_p-1-residue_mij_mpi)*size_mij_chunk_mpi
      portion_of_MIJ_per_task(2,k_mpi_proc)
     & = residue_mij_mpi*(size_mij_chunk_mpi+1)+	  
     & (k_p-residue_mij_mpi)*size_mij_chunk_mpi
      ENDDO	 
      total_size_check = total_size_check + size_mij_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_MIJ_per_task(2,k_mpi_proc).ne.total_size) 
     & STOP "WRONG MATIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO	
      IF(total_size_check.ne.total_size) THEN
      IF(myid.eq.0) WRITE(*,*)total_size_check,total_size	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF	  


      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_state_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_state_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_state_chunk_mpi+1)
      portion_of_state_per_task(2,k_mpi_proc) =
     & (k_p)*(size_state_chunk_mpi+1)
      ENDDO
      state_size_check = state_size_check + size_state_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_state_per_task(1,k_mpi_proc)
     & = residue_state_mpi*(size_state_chunk_mpi+1)+1  
     & + (k_p-1-residue_state_mpi)*size_state_chunk_mpi
      portion_of_state_per_task(2,k_mpi_proc)
     & = residue_state_mpi*(size_state_chunk_mpi+1)+	  
     & (k_p-residue_state_mpi)*size_state_chunk_mpi
      ENDDO	 
      state_size_check = state_size_check + size_state_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_state_per_task(2,k_mpi_proc).ne.states_size) 
     & STOP "WRONG MATRIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO	  
	  
!      PRINT*,"STATES",portion_of_state_per_task	  
      IF(state_size_check.ne.states_size) THEN
      IF(myid.eq.0) WRITE(*,*)state_size_check,states_size	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF	

      DO k_p=1,mpi_task_per_traject
      IF(k_p.le.residue_work_mpi) THEN
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k) 	  
      portion_of_work_per_task(1,k_mpi_proc)
     & = 1 + (k_p-1)*(size_work_chunk_mpi+1)
      portion_of_work_per_task(2,k_mpi_proc) =
     & (k_p)*(size_work_chunk_mpi+1)
      ENDDO
      work_size_check = work_size_check + size_work_chunk_mpi+1	  
      ELSE
      DO k=1,traject_roots
      k_mpi_proc = k_p + mpi_traject_roots(k)  	  
      portion_of_work_per_task(1,k_mpi_proc)
     & = residue_work_mpi*(size_work_chunk_mpi+1)+1  
     & + (k_p-1-residue_work_mpi)*size_work_chunk_mpi
      portion_of_work_per_task(2,k_mpi_proc)
     & = residue_work_mpi*(size_work_chunk_mpi+1)+	  
     & (k_p-residue_work_mpi)*size_work_chunk_mpi
      ENDDO	 
      work_size_check = work_size_check + size_work_chunk_mpi
      IF(k_p.eq.mpi_task_per_traject) THEN
      IF(portion_of_work_per_task(2,k_mpi_proc).ne.states_size*2+8) 
     & STOP "WRONG MATRIX ASSIGNEMENT"  
      ENDIF

	  
      ENDIF
	  
      ENDDO	  
	  
      IF(work_size_check.ne.states_size*2+8) THEN
      IF(myid.eq.0) WRITE(*,*)work_size_check,states_size*2+8	  
      STOP! "WRONG ASSIGNMENT"
      ENDIF
	  
      total_size_mpi = portion_of_MIJ_per_task(2,myid+1) - 
     & portion_of_MIJ_per_task(1,myid+1) + 1	  
      ALLOCATE(Mat_el(n_r_coll,total_size_mpi))
      ALLOCATE(Mat_el_der(n_r_coll,total_size_mpi))
      CALL MPI_Comm_group(MPI_COMM_WORLD, wrld_group,ierr_mpi)
      ALLOCATE(process_rank_distr(mpi_task_per_traject,traject_roots))
      ALLOCATE(comms_distr(mpi_task_per_traject),
     & groups_distr(mpi_task_per_traject))   	  
      DO i=1,traject_roots	  
      DO k=1,mpi_task_per_traject
      process_rank_distr(k,i) = k - 1 + (i-1)*mpi_task_per_traject     	  
      ENDDO
      ENDDO	

      DO i=1,mpi_task_per_traject      	  
      CALL MPI_Group_incl(wrld_group, traject_roots,
     & process_rank_distr(i,:), groups_distr(i),ierr_mpi)
      CALL MPI_Comm_create(MPI_COMM_WORLD,groups_distr(i),
     & comms_distr(i),ierr_mpi)
      ENDDO	  
!      PRINT*,myid,total_size_mpi	  
      tag1 = 1
      tag2 = 2
      IF(MYID.eq.0) PRINT*,"Mij ROOT PROC DISTRIBUTION STARTED"	  
      IF(MYID.eq.0) THEN
      DO k=2,mpi_task_per_traject
      total_size_mpi= portion_of_MIJ_per_task(2,k) - 
     & portion_of_MIJ_per_task(1,k) + 1	
      task_portion_size = total_size_mpi*n_r_coll	 
      ALLOCATE(buffer_mpi_portion(n_r_coll,total_size_mpi))
      buffer_mpi_portion = Mat_el_non_zero(:,
     & portion_of_MIJ_per_task(1,k):portion_of_MIJ_per_task(2,k)) 	  
      CALL MPI_SEND(buffer_mpi_portion,
     & task_portion_size, MPI_REAL8, k-1, 
     &  tag1, MPI_COMM_WORLD, ierr_mpi)
      buffer_mpi_portion = 	 Mat_el_non_zero_der (:,
     & portion_of_MIJ_per_task(1,k):portion_of_MIJ_per_task(2,k))
      CALL MPI_SEND(buffer_mpi_portion,
     & task_portion_size, MPI_REAL8, k-1, 
     &  tag2, MPI_COMM_WORLD, ierr_mpi)
      DEALLOCATE(buffer_mpi_portion) 
      ENDDO
      total_size_mpi= portion_of_MIJ_per_task(2,1) - 
     & portion_of_MIJ_per_task(1,1) + 1		  
      Mat_el=Mat_el_non_zero(:,
     & portion_of_MIJ_per_task(1,1):portion_of_MIJ_per_task(2,1))
      Mat_el_der=Mat_el_non_zero_der(:,
     & portion_of_MIJ_per_task(1,1):portion_of_MIJ_per_task(2,1))
      DEALLOCATE(Mat_el_non_zero,Mat_el_non_zero_der)	 
      ELSE
      IF(myid.le.mpi_task_per_traject-1) THEN	  
      task_portion_size = total_size_mpi*n_r_coll	  
      CALL MPI_RECV(Mat_el, task_portion_size, MPI_REAL8, 
     & 0, tag1, MPI_COMM_WORLD, status, ierr_mpi)	  
      CALL MPI_RECV(Mat_el_der, task_portion_size, MPI_REAL8, 
     & 0, tag2, MPI_COMM_WORLD, status, ierr_mpi)
      ENDIF	 
      ENDIF
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
      IF(myid.eq.0)	 PRINT*, "Mij ROOT PROC DISTRIBUTION DONE"  
      IF(myid.eq.0) PRINT*,"Mij ALL PROC DISTRIBUTION STARTED" 	  
      DO i=1,mpi_task_per_traject
      IF(i-1 .eq. myid - int(myid/mpi_task_per_traject)
     & *mpi_task_per_traject) THEN
      CALL MPI_Comm_rank(comms_distr(i),id_proc_in_group ,ierr_mpi)
      IF(id_proc_in_group.ne.myid/mpi_task_per_traject) PRINT*,
     & "ERROR IN COMUNICATIONS ASSIGNEMNET_1_distr",id_proc_in_group,
     & i-1
      IF(myid.lt.mpi_task_per_traject .and. id_proc_in_group.ne.0 )
     & PRINT*,"ERROR IN COMUNICATIONS ASSIGNEMNET_2_distr"
      task_portion_size = total_size_mpi*n_r_coll	  
      CALL MPI_BCAST(Mat_el,task_portion_size, MPI_REAL8,0,
     &  comms_distr(i),ierr_mpi)
      CALL MPI_BCAST(Mat_el_der,task_portion_size, MPI_REAL8,0,
     &  comms_distr(i),ierr_mpi)
      ENDIF	 
      ENDDO
      IF(myid.eq.0) PRINT*,"Mij ALL PROC DISTRIBUTION DONE"	  
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)	  
      IF(MYID.EQ.0) PRINT *,"MPI_COMMUNICATORS_CREATION STARTED"
      ALLOCATE(comms(traject_roots),groups(traject_roots))     	  
      ALLOCATE(process_rank(traject_roots,mpi_task_per_traject))		  
      DO i=1,traject_roots	  
      DO k=1,mpi_task_per_traject
      process_rank(i,k) = k - 1 + (i-1)*mpi_task_per_traject     	  
      ENDDO
      ENDDO	

      DO i=1,traject_roots      	  
      CALL MPI_Group_incl(wrld_group, mpi_task_per_traject,
     & process_rank(i,:), groups(i),ierr_mpi)
      CALL MPI_Comm_create(MPI_COMM_WORLD,groups(i),comms(i),ierr_mpi)
      ENDDO

      DO i=1,traject_roots
      IF(BELONGS(myid,process_rank(i,:),mpi_task_per_traject)) THEN
      CALL MPI_Comm_rank(comms(i),id_proc_in_group ,ierr_mpi)	  
      IF(mpi_traject_roots(i)+id_proc_in_group.ne.myid) STOP
     & "WRONG GROUP MPI ASSIGNEMENT"	  
	  
      ENDIF		  
      ENDDO	 
      IF(MYID.EQ.0) PRINT *,"MPI_COMMUNICATORS_CREATION DONE"
      TIME_MAT_FINISH = MPI_Wtime()
      TIME_2_MAT = TIME_MAT_FINISH - TIME_MAT_START	  
      ENDIF	  
! mpi_task_per_traject
! mpi_task_defined
      IF(.NOT.mpi_task_defined) THEN 	  
      CALL MPI_BCAST(Mat_el, n_r_coll*total_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(Mat_el_der, n_r_coll*total_size, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(Ks_have_to_be_skipped,total_size, MPI_LOGICAL,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      ENDIF
      CALL MPI_BCAST(R_COM, n_r_coll, MPI_REAL8,0,
     &  MPI_COMM_WORLD,ierr_mpi)
      CALL MPI_BCAST(ind_mat, 2*total_size, MPI_INTEGER,0,
     &  MPI_COMM_WORLD,ierr_mpi)	  

      ENDIF
      IF(states_to_exclude_defined) THEN
      ALLOCATE(stts_to_excl(states_size))	  
      IF(myid.eq.0) THEN
      CALL EXCLUDE_STATES	  
      ENDIF	  
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      	  
      CALL MPI_BCAST(stts_to_excl,states_size,MPI_LOGICAL,0,
     & MPI_COMM_WORLD,ierr_mpi)
      ENDIF	 	  
!      ALLOCATE(dq_dt_mpi(size_mij_chunk_mpi,states_size*2+8))
      TIME_MAT_FINISH= MPI_Wtime()
      IF(myid.eq.0) THEN	  
      WRITE(*,'(a57,1x,a2,f10.1)') 
     & "TIME SPENT ON MATRIX Mij READING/COMPUTING/SAVING ON DISK"
     & ,",s",TIME_MAT_FINISH-TIME_MAT_START
	  IF(print_matrix_defined.and. myid.eq.0) 
     & WRITE(*,'(a47,1x,a2,f10.1)')  
     & "TIME SPENT ON MATRIX Mij SAVING ON DISK"
     & ,",s",TIME_WRITING_MATRIX
      ENDIF	 
      IF(.not.run_prog_defined) STOP	  
      CALL MPI_BARRIER(MPI_COMM_WORLD,ierr_mpi)
!      STOP	  
      END SUBROUTINE MATRIX_MAKE
      SUBROUTINE PRINT_MIJ
! This subroutine is written by Alexander Semenov and modified by Bikramaditya Mandal
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      IMPLICIT NONE
      INTEGER st,i,k	  
      PRINT*,"SYSTEM_SETUP_DONE"
      OPEN(1,FILE=MATRIX_NAME_MIJ,ACTION="WRITE")
      WRITE(1,'(a15,1x,i1)') "COLLISION_TYPE=",coll_type
      WRITE(1,'(a17,x,i4)') "NUMBER_OF_CHANLS=",number_of_channels
      WRITE(1,'(a17,1x,i6)') "NUMBER_OF_STATES=",states_size
      WRITE(1,'(a12,1x,i9)') "MATRIX_SIZE=",total_size
      WRITE(1,'(a12,1x,i4)') "R_GRID_SIZE=",n_r_coll	  
      SELECT CASE(coll_type)
      CASE(1)
      IF(.not.fine_structure_defined) THEN	  
      WRITE(1,'(a8,1x,a8,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j_ch(i)	  
      ENDDO
      ELSE
      IF(SPIN_FINE.ne.2) THEN       	  
      WRITE(1,'(a8,1x,a8,1x,a3,1x,a3,1x,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J","F","N"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i3,1x,i3,1x,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j_ch(i),2+f_ch(i),j_ch(i)+f_ch(i)	  
      ENDDO
      ELSE
      WRITE(1,'(a8,1x,a8,1x,a4,1x,a5,1x,1x,a4,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J","F","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,f4.1,1x,f5.1,1x,1x,f4.1,1x,i3,1x,i3)')
     & st,i,j12_h(st),m12_h(st),j_h_ch(i),f_ch(i),par_lorb_ch(i)	  !!!! TO MODIFY
      ENDDO	  
      ENDIF	  
	  
      ENDIF	  
      CASE(2)
      WRITE(1,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","V","J"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),v_ch(i),j_ch(i)	  
      ENDDO
      CASE(3)
      WRITE(1,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J","K", "EPS"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j_ch(i),k_ch(i),eps_ch(i)	  
      ENDDO
      CASE(4)
      WRITE(1,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J","KA", "KC"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j_ch(i),ka_ch(i),kc_ch(i)	  
      ENDDO      	  
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      WRITE(1,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      WRITE(1,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","J2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i),parity_state(i)	  
      ENDDO	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      WRITE(1,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","V1","V2","J1","J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),v1_ch(i),v2_ch(i),j1_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      WRITE(1,'(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","V1","V2","J1","J2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),v1_ch(i),v2_ch(i),j1_ch(i),j2_ch(i),
     & parity_state(i)	 
      ENDDO	  
      ENDIF	  
      CASE(7)
      WRITE(1,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","K","EPS",
     & "J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),k1_ch(i),eps1_ch(i)
     & ,j2_ch(i) 
      ENDDO	  
      CASE(8)
      WRITE(1,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","K1A","K1C",
     & "J2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i) 
      ENDDO		  
      CASE(9)
      WRITE(1,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","K1A","K1C",
     & "J2","K2","EPS"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),k2_ch(i),eps2_ch(i) 
      ENDDO		  
      CASE(0)
      IF(.NOT.identical_particles_defined) THEN  
      WRITE(1,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3)')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1",
     & "J2","KA2","KC2"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i)	 
      ENDDO
      ELSE
      WRITE(1,
     & '(a8,1x,a8,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a3,1x,a2)
     & ')
     & "#STATE","#CHANNEL","J12","M12","J1","KA1","KC1",
     & "J2","KA2","KC2","P"
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(1,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i2)
     & ')
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i),parity_state(st)
      ENDDO	 
      ENDIF	 
      END SELECT
      WRITE(1,"(a9,1x,a4,1x,a4)") "MIJ_INDEX", "ST_1", "ST_2"
      DO i=1,total_size
      WRITE(1,"(i9,1x,i4,1x,i4)") i,ind_mat(1,i),ind_mat(2,i)	  
      ENDDO		  
      WRITE(1,'(a5,1x,a19)') "#iR","R_COM(iR)"	  
      DO i=1,n_r_coll
      WRITE(1,'(i5,1x,e19.12)')i,R_COM(i)	  
      ENDDO
      WRITE(1,'(a8,1x,a3,1x,a19)') "#k_matel","#iR",
     & "Mij(iR,k_matel)"!,
!     & "Mij_der(iR,k_matel)"	  
      DO k=1,total_size
      DO i=1,n_r_coll	  
      WRITE(1,'(i8,1x,i3,1x,e19.12)')k,i,Mat_el(i,k)!,
!     & Mat_el_der(i,k)	  
      ENDDO
      ENDDO
      CLOSE(1)
!      OPEN(333,file="LAST_MIJ.DAT")	   !!!!! DELETE AFTER ALL
!	  DO k=1,total_size !!!!! DELETE AFTER ALL
!	  IF(ind_mat(1,k).eq.ind_mat(2,k)) THEN
!	  WRITE(333,*) k,Mat_el(n_r_coll,k)
!      ENDIF
!      ENDDO	
!      CLOSE(333)	  
      IF(test_expansion_defined) THEN 
      CALL PRINT_ELASTIC_MIJ
      PRINT*, "ELASTIC_TERMS_FILE_CREATED"  
      ENDIF	  
      END SUBROUTINE PRINT_MIJ

      SUBROUTINE READ_MIJ !!! TO DO
! This subroutine is written by Alexander Semenov and modified by Bikramaditya Mandal
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      USE OLD_MIJ
      IMPLICIT NONE
      REAL*8 Mij_buffer	  
      INTEGER i,k,st,tot_siz_read,st_siz_read,j_count,j_summ,p_count
      LOGICAL file_exst
      INTEGER KRONEKER,p_lim_max_ini	  
      EXTERNAL KRONEKER	  
      INQUIRE( FILE=MATRIX_NAME_MIJ, EXIST=file_exst )
      IF(.not.file_exst) THEN 
      PRINT*, "ERROR: MTRX.DAT NOT FOUND"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF
      IF(.not.identical_particles_defined) THEN	  
      p_lim_max = 1
      p_lim_min = 1	
      ELSE
      SELECT CASE(coll_type)		  
      CASE(5)	  

      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
      CASE(6)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))*
     & KRONEKER(v1_ch(chann_ini),v2_ch(chann_ini)) 	  
      CASE(0)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
     & *KRONEKER(ka1_ch(chann_ini),ka2_ch(chann_ini))
     & *KRONEKER(kc1_ch(chann_ini),kc2_ch(chann_ini))	  
      END SELECT		  
      ENDIF	  
      OPEN(1,FILE=MATRIX_NAME_MIJ,STATUS="OLD",ACTION="READ")
      PRINT*,"MTRX.DAT READING STARTED"	  
      READ(1,'(a15,1x,i1)') buffer_word_2,coll_type_old
      IF(coll_type_old.ne.coll_type) THEN
      PRINT*, "ERROR: SYSTEM IS DIFFERENT"
      CRITICAL_ERROR = .TRUE.	  
      ENDIF	 
      READ(1,'(a17,x,i4)') buffer_word_3,number_of_channels_old	 
      READ(1,'(a17,1x,i6)') buffer_word_3,states_size_old
      READ(1,'(a12,1x,i9)') buffer_word_1,total_size_old	  
      READ(1,'(a12,1x,i4)') buffer_word_1,n_r_coll_old
      IF(n_r_coll_old.ne.n_r_coll) THEN
      PRINT*, "ERROR:WRONG GRID"
      CRITICAL_ERROR = .TRUE.	  
      RETURN	  
      ENDIF	  
      ALLOCATE(j12_old(states_size_old),m12_old(states_size_old))
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) THEN
      ALLOCATE(j12_h_old(states_size_old),m12_h_old(states_size_old))
      ENDIF	  
      ALLOCATE(indx_chann_old(states_size_old))
      IF(identical_particles_defined) 
     & ALLOCATE(parity_states_old(states_size_old))  
      SELECT CASE(coll_type)
      CASE(1)
      READ(1,*)
      ALLOCATE(j_ch_old(number_of_channels_old))
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined)
     & ALLOCATE(j_h_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      IF(fine_structure_defined) THEN	!!!! MODIFY
      IF(SPIN_FINE.ne.2) THEN	  
	  READ(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,f_old_b,n_old_b
      ELSE
	  READ(1,'(i8,1x,i8,1x,f4.1,1x,f5.1,1x,f4.1,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_h_old(st),m12_h_old(st),j_h_old_b,f_old_b,par_old_b	  
!!!! TO ADD	  
      ENDIF
      ELSE
	  READ(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) THEN
      j_h_ch_old(indx_chann_old(st_old)) = j_h_old_b	
      ELSE	  
      j_ch_old(indx_chann_old(st_old)) = j_old_b  
      ENDIF
      ENDDO		  
      CASE(2)
      READ(1,*)
      ALLOCATE(v_ch_old(number_of_channels_old),
     & j_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),v_old_b,j_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      v_ch_old(indx_chann_old(st_old)) = v_old_b	  
      ENDDO
      CASE(3)
      READ(1,*)
      ALLOCATE(j_ch_old(number_of_channels_old),
     & k_ch_old(number_of_channels_old),
     & eps_ch_old(number_of_channels_old) )	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,k_old_b,eps_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      k_ch_old(indx_chann_old(st_old)) = k_old_b
      eps_ch_old(indx_chann_old(st_old)) = eps_old_b	  
      ENDDO	  
      CASE(4)
      READ(1,*)
      ALLOCATE(j_ch_old(number_of_channels_old),
     & ka_ch_old(number_of_channels_old),
     & kc_ch_old(number_of_channels_old) )	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,ka_old_b,kc_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      ka_ch_old(indx_chann_old(st_old)) = ka_old_b
      kc_ch_old(indx_chann_old(st_old)) = kc_old_b	  
      ENDDO	  	  
      CASE(5)
      READ(1,*)
      IF(.not.identical_particles_defined) THEN	  
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & j1_old_b,j2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Mij FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b		 
      ENDDO
      ELSE
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1xi3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & j1_old_b,j2_old_b,par_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG Muj FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      parity_states_old(st_old) = par_old_b	  
      ENDDO	  
      ENDIF	  
      CASE(6)
      READ(1,*)
      IF(.not.identical_particles_defined) THEN	  
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))
      ALLOCATE(v1_ch_old(number_of_channels_old),
     & v2_ch_old(number_of_channels_old))	 
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & v1_old_b,v2_old_b,j1_old_b,j2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      v1_ch_old(indx_chann_old(st_old)) = v1_old_b
      v2_ch_old(indx_chann_old(st_old)) = v2_old_b		  
      ENDDO
      ELSE
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & v1_old_b,v2_old_b,j1_old_b,j2_old_b,par_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      v1_ch_old(indx_chann_old(st_old)) = v1_old_b
      v2_ch_old(indx_chann_old(st_old)) = v2_old_b		  
      parity_states_old(st_old) = par_old_b	  
      ENDDO	  
      ENDIF	  
      CASE(7)
      READ(1,*)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & k1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),	 
     & eps1_ch_old(number_of_channels_old) )	  
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j1_old_b,k1_old_b,eps1_old_b,j2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b	  
      k1_ch_old(indx_chann_old(st_old)) = k1_old_b
      eps1_ch_old(indx_chann_old(st_old)) = eps1_old_b
      ENDDO	  
      CASE(8)
      READ(1,*)	  
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old))
      DO st=1,states_size_old
      READ(1,'(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j1_old_b,ka1_old_b,kc1_old_b,
     & j2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
	  
      ENDDO	 	 
      CASE(9)
      READ(1,*)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old)
     & ,k2_ch_old(number_of_channels_old),
     & eps2_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old
      READ(1,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,k2_old_b,eps2_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      k2_ch_old(indx_chann_old(st_old)) = k2_old_b
      eps2_ch_old(indx_chann_old(st_old)) = eps2_old_b
      ENDDO	  
      CASE(0)
      READ(1,*)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old)
     & ,ka2_ch_old(number_of_channels_old),
     & kc2_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old
      IF(.NOT.identical_particles_defined) THEN	  
      READ(1,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3)')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b
      ELSE
      READ(1,
     & '(i8,1x,i8,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i3,1x,i2)
     & ')
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b,par_old_b	  
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      ka2_ch_old(indx_chann_old(st_old)) = ka2_old_b
      kc2_ch_old(indx_chann_old(st_old)) = kc2_old_b
      IF(identical_particles_defined)
     & parity_states_old(st_old) = par_old_b	  
      ENDDO	 
      END SELECT
      st = 0    
!!!! TO BE CONTINUED	  
      DO i=1,min(number_of_channels_old,number_of_channels)
      SELECT CASE(coll_type)
      CASE(1)
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) THEN
      DO j_count = 1, int(j_h_ch_old(i)*2) + 1	  
      st = st + 1
      IF(j_count*2 .eq. int(j_h_ch_old(i)*2)+1
     & .and. chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF	 
      ENDDO	  
      ELSE	  
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	  
      ENDDO
      ENDIF	  
	  
      CASE(2)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(3)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(4)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
      ENDDO
      ENDDO
      ELSE
      p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch(i),j2_ch(i)))	  
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	
      ENDDO
      ENDDO
      ELSE
      p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch(i),j2_ch(i))*
     & KRONEKER(v1_ch(i),v2_ch(i))) 	  
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO 
  

	  
      ENDIF	  
      CASE(7)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(8)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(9)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(0)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      ELSE
      p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i)))	  
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
!     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
!     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     &	  abs(j1_ch_old(i)-j2_ch_old(i))
     & .and.  p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	   	  
      END SELECT	  
      ENDDO
! HERE WE CHECK THE ORDER !! FOR FINE SPIN=2 WE JUST OVERJUM
      IF(SPIN_FINE.eq.2 .and. fine_structure_defined) GOTO 1994
      DO st=1,min(states_size_old,states_size)
      IF(j12(st).ne.j12_old(st)) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: J12s ARE DIFFERENT"
      RETURN	  
      ENDIF	  
      IF(m12(st).ne.m12_old(st)) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR IN INI: CHECK M12"
      RETURN	  
      ENDIF	  
      SELECT CASE(coll_type)
      CASE(1)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*,"ERROR IN INI : CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      CASE(2)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      IF(v_ch_old(indx_chann_old(st)).ne.v_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI : CHECK V"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	    
      CASE(3)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      IF(k_ch_old(indx_chann_old(st)).ne.k_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK K"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      IF(eps_ch_old(indx_chann_old(st)).ne.eps_ch(indx_chann(st)))
     & THEN
      PRINT*, "ERROR IN INI: CHECK EPS"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      CASE(4)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka_ch_old(indx_chann_old(st)).ne.ka_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc_ch_old(indx_chann_old(st)).ne.kc_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(5)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(6)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(v1_ch_old(indx_chann_old(st)).ne.v1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK V1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(v2_ch_old(indx_chann_old(st)).ne.v2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK V2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(7)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(k1_ch_old(indx_chann_old(st)).ne.k1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK K1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(eps1_ch_old(indx_chann_old(st)).ne.eps1_ch(indx_chann(st)))
     & THEN
      PRINT*, "ERROR IN INI: CHECK EPS1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(8)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(9)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(k2_ch_old(indx_chann_old(st)).ne.k2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK K2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(eps2_ch_old(indx_chann_old(st)).ne.eps2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK EPS2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(0)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka2_ch_old(indx_chann_old(st)).ne.ka2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc2_ch_old(indx_chann_old(st)).ne.kc2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(identical_particles_defined) THEN	 
      IF(parity_states_old(st)
     & .ne.parity_state(st))
     & STOP "ERROR IN INI: CHECK PARITY"
      ENDIF	 
      END SELECT	  
	  
	  
      ENDDO	 	  

! HERE WE CHECK THE ORDER	  
!      PRINT*, "WE HAVE CKECED"	  
1994      READ(1,*)
!      PRINT*, buffer_word_3
      DO i=1,total_size_old
      READ(1,"(i9,1x,i4,1x,i4)") i_old,ind_mat_old_1,ind_mat_old_2
      IF(i.ne.i_old) CRITICAL_ERROR = .TRUE.
      IF(i.le.total_size) THEN	  
      IF(ind_mat_old_1.ne.ind_mat(1,i)) CRITICAL_ERROR = .TRUE.
      IF(ind_mat_old_2.ne.ind_mat(2,i))	CRITICAL_ERROR = .TRUE.
      ENDIF	  
      ENDDO		  
      READ(1,*)	  
      DO i=1,n_r_coll
      READ(1,'(i5,1x,e19.12)')i_old,R_COM(i)	  
      ENDDO
      READ(1,*)	  
      DO  k=1,min(total_size,total_size_old)
      DO i=1,n_r_coll
      READ(1,'(i8,1x,i3,1x,e19.12)')k_old,i_old,
     & Mat_el(i,k) 	  
      IF(k_old.ne.k) PRINT*,"ERROR IN MIJ: k is wrong"
      IF(i.ne.i_old) PRINT*,"ERROR IN MIJ : ir is wrong"
      ENDDO

! Bikram Start Dec 2019:
	  if(bikram_mij_shift) then
	  do i=1,n_r_coll
	  Mat_el(i,k) = Mat_el(i,k) - Mat_el(n_r_coll,k)		!Bikram
	  enddo
	  endif
! Bikram End.

      deriv_bgn = (Mat_el(2,k)
     & - Mat_el(1,k))/
     & (R_COM(2)-R_COM(1)) 
      deriv_end = (Mat_el(n_r_coll,k)
     & - Mat_el(n_r_coll-1,k))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1)) 	 
      CALL  spline(R_COM,Mat_el(:,k),n_r_coll,
     & deriv_bgn,deriv_end,
     & Mat_el_der(:,k))
 
      ENDDO
      CLOSE(1)
	  
      PRINT*,"FORMATTED READING DONE"
      IF(test_expansion_defined)CALL PRINT_ELASTIC_MIJ
      IF(CRITICAL_ERROR) PRINT*,"ERROR IN ASSIGNMENT"	  
      RETURN	  
      END SUBROUTINE READ_MIJ	  

      SUBROUTINE EXCLUDE_STATES
! This subroutine is written by Alexander Semenov
      USE VARIABLES
      USE OLD_MIJ	  
      IMPLICIT NONE
      INTEGER i_state,st_user_exclude,file_stat
      LOGICAL exst 	  
      CHARACTER(LEN=21)::STATES_TO_EXCLUDE="STATES_TO_EXCLUDE.DAT"
      INQUIRE( FILE=STATES_TO_EXCLUDE, EXIST=exst )
      DO i_state = 1,states_size
      stts_to_excl(i_state) = .FALSE.	  
      ENDDO	  
      IF(.not.exst) RETURN	  
      OPEN(UNIT=234,FILE=STATES_TO_EXCLUDE,STATUS="OLD",
     & ACTION = "READ")	  
      READ(234,*,IOSTAT=file_stat) !HEADER
      IF(file_stat.ne.0) THEN
      CLOSE(234)
      RETURN
      ENDIF	  
      DO WHILE(file_stat.eq.0)
      READ(234,*,IOSTAT=file_stat) st_user_exclude
      IF(st_user_exclude.gt.states_size) THEN
      CLOSE(234)	  
      RETURN
      ENDIF	   
      stts_to_excl(st_user_exclude) = .TRUE.	  
      ENDDO	  
      RETURN	  
      END SUBROUTINE EXCLUDE_STATES	  
      SUBROUTINE READ_MIJ_USER
! This subroutine is written by Alexander Semenov and modified by Bikramaditya Mandal
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      USE OLD_MIJ
      IMPLICIT NONE
      REAL*8 Mij_buffer	  
      INTEGER i,k,st,tot_siz_read,st_siz_read,j_count,j_summ,p_count
      INTEGER istat,p_lim_max_ini	  
      LOGICAL file_exst
      INTEGER KRONEKER
      CHARACTER (LEN=12) ::
     & MIJ_FILE_NAME_N ="Mij_UF_N.dat"
      CHARACTER (LEN=10) ::
     & MIJ_FILE_NAME ="Mij_UF.dat"       	  
      EXTERNAL KRONEKER	
      INQUIRE( FILE=MATRIX_NAME_MIJ_UF, EXIST=file_exst )
      IF(.not.file_exst) THEN
      PRINT*,"ERROR: MTRX_UF.DAT NOT FOUND"
      CRITICAL_ERROR = .TRUE.	  
      RETURN	  
      ENDIF
      IF(.not.identical_particles_defined) THEN	  
      p_lim_max = 1
      p_lim_min = 1	
      ELSE
      SELECT CASE(coll_type)		  
      CASE(5)	  

      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
      CASE(6)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))*
     & KRONEKER(v1_ch(chann_ini),v2_ch(chann_ini)) 	  
      CASE(0)
      p_lim_min = 1
      p_lim_max_ini = 2-KRONEKER(j1_ch(chann_ini),j2_ch(chann_ini))
     & *KRONEKER(ka1_ch(chann_ini),ka2_ch(chann_ini))
     & *KRONEKER(kc1_ch(chann_ini),kc2_ch(chann_ini))	  
      END SELECT		  
      ENDIF		  
      OPEN(1,FORM="UNFORMATTED",
     % FILE=MATRIX_NAME_MIJ_UF,STATUS="OLD",ACTION="READ")
      PRINT*,"MTRX_UF.DAT READING STARTED"	  
      READ(1,IOSTAT=istat) coll_type_old
      IF(coll_type_old.ne.coll_type) THEN 
      PRINT*, "ERROR: SYSTEM IS DIFFERENT"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      READ(1,IOSTAT=istat) number_of_channels_old	 
      READ(1,IOSTAT=istat) states_size_old
      READ(1,IOSTAT=istat) total_size_old	  
      READ(1,IOSTAT=istat)n_r_coll_old
      IF(n_r_coll_old.ne.n_r_coll) THEN
      CRITICAL_ERROR = .TRUE.
      PRINT*,"WRONG GRID"
      RETURN	  
      ENDIF	  
      IF(istat.ne.0) THEN
      CRITICAL_ERROR = .TRUE.
      PRINT*,"ERROR in ",MATRIX_NAME_MIJ_UF
      RETURN	  
      ENDIF	  
!      IF(total_size_old.lt.total_size) THEN
!      CALL SYSTEM('cp '//
!     & MIJ_FILE_NAME//' '//MIJ_FILE_NAME_N)	  
!      ENDIF
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) GOTO 1992	  
      ALLOCATE(j12_old(states_size_old),m12_old(states_size_old))
      ALLOCATE(indx_chann_old(states_size_old))
      IF(identical_particles_defined) 
     & ALLOCATE(parity_states_old(states_size_old))  
      SELECT CASE(coll_type)
      CASE(1)
      ALLOCATE(j_ch_old(number_of_channels_old))     
      DO st=1,states_size_old
      READ(1)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	  
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      ENDDO		  
      CASE(2)
      ALLOCATE(v_ch_old(number_of_channels_old),
     & j_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),v_old_b,j_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      v_ch_old(indx_chann_old(st_old)) = v_old_b	  
      ENDDO
      CASE(3)
      ALLOCATE(k_ch_old(number_of_channels_old),
     & j_ch_old(number_of_channels_old),
     & eps_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,k_old_b,eps_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      k_ch_old(indx_chann_old(st_old)) = k_old_b
      eps_ch_old(indx_chann_old(st_old)) = eps_old_b
      ENDDO	  
      CASE(4)	  
      ALLOCATE(ka_ch_old(number_of_channels_old),
     & j_ch_old(number_of_channels_old),
     & kc_ch_old(number_of_channels_old))	  
      DO st=1,states_size_old
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),j_old_b,ka_old_b,kc_old_b
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF		 
      j_ch_old(indx_chann_old(st_old)) = j_old_b
      ka_ch_old(indx_chann_old(st_old)) = ka_old_b
      kc_ch_old(indx_chann_old(st_old)) = kc_old_b
      ENDDO	  
      CASE(5)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))
      DO st=1,states_size_old
      IF(.not.identical_particles_defined) THEN	 
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & j1_old_b,j2_old_b
      ELSE
      READ(1,IOSTAT=istat)	  
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & j1_old_b,j2_old_b,par_old_b
      parity_states_old(st_old) = par_old_b	 
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b		 
      ENDDO
      CASE(6)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old))
      ALLOCATE(v1_ch_old(number_of_channels_old),
     & v2_ch_old(number_of_channels_old))	 
      DO st=1,states_size_old
      IF(.not.identical_particles_defined) THEN	 
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & v1_old_b,j1_old_b,v2_old_b,j2_old_b
      ELSE
      READ(1,IOSTAT=istat)	  
     & st_old,indx_chann_old(st),j12_old(st),m12_old(st),
     & v1_old_b,j1_old_b,v2_old_b,j2_old_b,par_old_b
      parity_states_old(st_old) = par_old_b	 
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF			 
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      v1_ch_old(indx_chann_old(st_old)) = v1_old_b
      v2_ch_old(indx_chann_old(st_old)) = v2_old_b	  
      ENDDO
      CASE(7)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & k1_ch_old(number_of_channels_old),
     & eps1_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old

      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,k1_old_b,eps1_old_b
     & ,j2_old_b
	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      k1_ch_old(indx_chann_old(st_old)) = k1_old_b
      eps1_ch_old(indx_chann_old(st_old)) = eps1_old_b
      ENDDO		  
      CASE(8)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old

      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b
	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      ENDDO	 	  	  
      CASE(9)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old)
     & ,k2_ch_old(number_of_channels_old),
     & eps2_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old

      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,k2_old_b,eps2_old_b
	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      k2_ch_old(indx_chann_old(st_old)) = k2_old_b
      eps2_ch_old(indx_chann_old(st_old)) = eps2_old_b
      ENDDO	 	  
      CASE(0)
      ALLOCATE(j1_ch_old(number_of_channels_old),
     & j2_ch_old(number_of_channels_old),
     & ka1_ch_old(number_of_channels_old),
     & kc1_ch_old(number_of_channels_old)
     & ,ka2_ch_old(number_of_channels_old),
     & kc2_ch_old(number_of_channels_old))		  
      DO st=1,states_size_old
      IF(.NOT.identical_particles_defined) THEN	  
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b
      ELSE
      READ(1,IOSTAT=istat)
     & st_old,indx_chann_old(st),
     & j12_old(st),m12_old(st),
     & j1_old_b,ka1_old_b,kc1_old_b
     & ,j2_old_b,ka2_old_b,kc2_old_b,par_old_b	  
      ENDIF	  
      IF(st_old.ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: WRONG MIJ FILE"
      RETURN	  
      ENDIF	
      j1_ch_old(indx_chann_old(st_old)) = j1_old_b
      j2_ch_old(indx_chann_old(st_old)) = j2_old_b
      ka1_ch_old(indx_chann_old(st_old)) = ka1_old_b
      kc1_ch_old(indx_chann_old(st_old)) = kc1_old_b
      ka2_ch_old(indx_chann_old(st_old)) = ka2_old_b
      kc2_ch_old(indx_chann_old(st_old)) = kc2_old_b
      IF(identical_particles_defined)
     & parity_states_old(st_old) = par_old_b	  
      ENDDO	 
      END SELECT
      st = 0       
      DO i=1,min(number_of_channels_old,number_of_channels)
      SELECT CASE(coll_type)
      CASE(1)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	  
      ENDDO
      CASE(2)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(3)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(4)
      DO j_count = -j_ch_old(i),j_ch_old(i)	  
      st = st + 1
      IF(j_count.eq.0) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF
      ENDDO
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
      ENDDO
      ENDDO
      ELSE
      p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch(i),j2_ch(i)))	  
	  
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.1  ) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	
      ENDDO
      ENDDO
      ELSE
      p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch(i),j2_ch(i))*
     & KRONEKER(v1_ch(i),v2_ch(i))) 	  
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch(i),j2_ch(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq. abs(j1_ch(i)-j2_ch(i))
     & .and. p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*,
     & "ERROR IN ASSINGNMENT CHANNEL -  STATE"
      RETURN	 
      ENDIF
      ENDIF	 
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO 
  

	  
      ENDIF	  
      CASE(7)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(8)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(9)
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      CASE(0)
      IF(.not.identical_particles_defined) THEN	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ	  
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     & abs(j1_ch_old(i)-j2_ch_old(i))) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
      ENDDO
      ENDDO
      ELSE
      p_lim_max=min(p_lim_max_ini,2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i)))	  
      DO p_count=p_lim_min,p_lim_max!2-KRONEKER(j1_ch_old(i),j2_ch_old(i))
!     & *KRONEKER(ka1_ch_old(i),ka2_ch_old(i))
!     & *KRONEKER(kc1_ch_old(i),kc2_ch_old(i))	  
      DO j_summ = abs(j1_ch_old(i)-j2_ch_old(i)),
     & j1_ch_old(i)+j2_ch_old(i)
      DO j_count = -j_summ,j_summ
      st = st + 1
      IF(j_count.eq.0 .and. j_summ .eq.
     &	  abs(j1_ch_old(i)-j2_ch_old(i))
     & .and.  p_count .eq.p_lim_min  ) THEN
      IF(chann_indx(i).ne.st) PRINT*,
     & "ERROR IN ASSINGNMENT CHANN INDEX",
     & "chann_indx(i)","st",chann_indx(i),st,"i",i
      ENDIF	
c      WRITE(*,'(a2,a3,1x,i3,1x,i3,a3,i2,a3,i2,a3,i3,a3,i3)')
c     & "st"," ch",st,i,
c     & "j1",j1_ch(i),"j2",j2_ch(i),"lc",l_count,"pc",p_count
      ENDDO
      ENDDO
	  
      ENDDO	  
  

	  
      ENDIF	   	  
      END SELECT	  
      ENDDO
! HERE WE CHECK THE ORDER 
      DO st=1,min(states_size_old,states_size)
      IF(j12(st).ne.j12_old(st)) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR: J12s ARE DIFFERENT"
      RETURN	  
      ENDIF	  
      IF(m12(st).ne.m12_old(st)) THEN
      CRITICAL_ERROR = .TRUE.	  
      PRINT*, "ERROR IN INI: CHECK M12"
      RETURN	  
      ENDIF	  
      SELECT CASE(coll_type)
      CASE(1)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*,"ERROR IN INI : CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      CASE(2)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      IF(v_ch_old(indx_chann_old(st)).ne.v_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI : CHECK V"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	    
      CASE(3)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      IF(k_ch_old(indx_chann_old(st)).ne.k_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK K"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      IF(eps_ch_old(indx_chann_old(st)).ne.eps_ch(indx_chann(st)))
     & THEN
      PRINT*, "ERROR IN INI: CHECK EPS"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	 
      CASE(4)
      IF(j_ch_old(indx_chann_old(st)).ne.j_ch(indx_chann(st))) THEN
      PRINT*, "ERROR IN INI: CHECK J"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka_ch_old(indx_chann_old(st)).ne.ka_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc_ch_old(indx_chann_old(st)).ne.kc_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(5)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(6)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(v1_ch_old(indx_chann_old(st)).ne.v1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK V1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(v2_ch_old(indx_chann_old(st)).ne.v2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK V2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(7)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(k1_ch_old(indx_chann_old(st)).ne.k1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK K1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(eps1_ch_old(indx_chann_old(st)).ne.eps1_ch(indx_chann(st)))
     & THEN
      PRINT*, "ERROR IN INI: CHECK EPS1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(8)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(9)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(k2_ch_old(indx_chann_old(st)).ne.k2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK K2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(eps2_ch_old(indx_chann_old(st)).ne.eps2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK EPS2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      CASE(0)
      IF(j1_ch_old(indx_chann_old(st)).ne.j1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(j2_ch_old(indx_chann_old(st)).ne.j2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK J2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka1_ch_old(indx_chann_old(st)).ne.ka1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(ka2_ch_old(indx_chann_old(st)).ne.ka2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KA2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc1_ch_old(indx_chann_old(st)).ne.kc1_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(kc2_ch_old(indx_chann_old(st)).ne.kc2_ch(indx_chann(st)))THEN
      PRINT*, "ERROR IN INI: CHECK KC2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF		  
      IF(identical_particles_defined) THEN	 
      IF(parity_states_old(st)
     & .ne.parity_state(st))
     & STOP "ERROR IN INI: CHECK PARITY"
      ENDIF	 
      END SELECT	  
	  
	  
      ENDDO	  

! HERE WE CHECK THE ORDER	  
1992  DO i=1,total_size_old
      READ(1) i_old,ind_mat_old_1,ind_mat_old_2
      IF(i.ne.i_old) STOP "ERROR IN INI: CHECK FILE"
      IF(i.le.total_size) THEN	  
      IF(ind_mat_old_1.ne.ind_mat(1,i)) THEN
      PRINT*, "ERROR IN INI:IND_MAT_1"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      IF(ind_mat_old_2.ne.ind_mat(2,i))	THEN
      PRINT*,"ERROR IN INI:IND_MAT_2"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
      ENDIF	  
      ENDDO		  
      DO i=1,n_r_coll
      READ(1)i_old,R_COM(i)	  
      ENDDO
      DO  k=1,min(total_size,total_size_old)
      READ(1) k_old
!      DO i=1,n_r_coll
      READ(1) Mat_el(:,k) 	  
      IF(k_old.ne.k) THEN
      PRINT*,"ERROR IN MIJ: k is wrong"
      CRITICAL_ERROR = .TRUE.
      RETURN	  
      ENDIF	  
!      IF(i.ne.i_old) PRINT*,"ERROR IN MIJ : ir is wrong"
!      ENDDO
! SPLINING	  


! Bikram Start Dec 2019:
	  if(bikram_mij_shift) then
	  do i=1,n_r_coll
	  Mat_el(i,k) = Mat_el(i,k) - Mat_el(n_r_coll,k)		!Bikram
	  enddo
	  endif
! Bikram End.

      deriv_bgn = (Mat_el(2,k)
     & - Mat_el(1,k))/
     & (R_COM(2)-R_COM(1)) 
      deriv_end = (Mat_el(n_r_coll,k)
     & - Mat_el(n_r_coll-1,k))/
     & (R_COM(n_r_coll)-R_COM(n_r_coll-1))
!      deriv_end = 0d0
!      deriv_bgn = 0d0	  
      CALL  spline(R_COM,Mat_el(:,k),n_r_coll,
     & deriv_bgn,deriv_end,
     & Mat_el_der(:,k))
      ENDDO
      CLOSE(1)
      IF(test_expansion_defined)CALL PRINT_ELASTIC_MIJ	  
      PRINT*,"UNFORMATTED READING DONE"
      !!! TESTING
!      CALL SYSTEM('rm '// MIJ_FILE_NAME_N)
!      OPEN(333,file="LAST_MIJ.DAT")	            !!!!! DELETE AFTER ALL
!	  DO k=1,total_size
!	  IF(ind_mat(1,k).eq.ind_mat(2,k)) THEN
!	  WRITE(333,*) k,Mat_el(n_r_coll,k)
!      ENDIF
!      ENDDO	
!      CLOSE(333)		  
      RETURN	  
      !!! TESTING	  
      END SUBROUTINE READ_MIJ_USER 
!!!!!!!!!!!!! PREPARE IT FOR MPI TASK PER TRAJECTORY
     
      SUBROUTINE PRINT_MIJ_USER
! This subroutine is written by Alexander Semenov and modified by Bikramaditya Mandal
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      USE OLD_MIJ
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      IMPLICIT NONE
      INTEGER st,i,k
      CHARACTER (LEN=12) ::
     & MIJ_FILE_NAME_N ="Mij_UF_N.dat"
      CHARACTER (LEN=10) ::
     & MIJ_FILE_NAME ="Mij_UF.dat"    	  
      IF(MYID.ne.0) RETURN	  
      PRINT*,"SYSTEM_SETUP_DONE"
      OPEN(111,FORM="UNFORMATTED",FILE=MATRIX_NAME_MIJ_UF,
     & ACTION="WRITE")
      WRITE(111)coll_type
      WRITE(111)number_of_channels
      WRITE(111)states_size
      WRITE(111)total_size
      WRITE(111)n_r_coll
      IF(fine_structure_defined .and. SPIN_FINE.eq.2) GOTO 1993	  
      SELECT CASE(coll_type)
      CASE(1)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j_ch(i)	  
      ENDDO		  
      CASE(2)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),v_ch(i),j_ch(i)
      ENDDO	
      CASE(3)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j_ch(i),k_ch(i),eps_ch(i)
      ENDDO
      CASE(4)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j_ch(i),ka_ch(i),kc_ch(i)
      ENDDO		  
      CASE(5)
      IF(.not.identical_particles_defined) THEN	  
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),j2_ch(i),parity_state(st)	  
      ENDDO	  
      ENDIF	  
      CASE(6)
      IF(.not.identical_particles_defined)	THEN  
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),v1_ch(i),j1_ch(i),v2_ch(i),j2_ch(i)	  
      ENDDO
      ELSE
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),v1_ch(i),j1_ch(i),v2_ch(i),j2_ch(i),
     & parity_state(st)	  
      ENDDO	  
      ENDIF
      CASE(7)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),k1_ch(i),eps1_ch(i)
     & ,j2_ch(i)	 
      ENDDO
      CASE(8)	  
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i)	 
      ENDDO
      CASE(9)
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),kc2_ch(i),eps2_ch(i)	 
      ENDDO	  
      CASE(0)
      IF(.NOT.identical_particles_defined) THEN  
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i)	 
      ENDDO
      ELSE
      DO st=1,states_size
      i = indx_chann(st)   	  
      WRITE(111)
     & st,i,j12(st),m12(st),j1_ch(i),ka1_ch(i),kc1_ch(i)
     & ,j2_ch(i),ka2_ch(i),kc2_ch(i),parity_state(st)
      ENDDO	 
      ENDIF	 
      END SELECT
1993  DO i=1,total_size
      WRITE(111) i,ind_mat(1,i),ind_mat(2,i)	  
      ENDDO		  
      DO i=1,n_r_coll
      WRITE(111)i,R_COM(i)	  
      ENDDO
      DO k=1,total_size
      WRITE(111) k	  
!      DO i=1,n_r_coll	  
      WRITE(111) Mat_el(:,k)!,
!      ENDDO
      ENDDO
      CLOSE(111)
!      CALL SYSTEM('rm '// MIJ_FILE_NAME_N)
!      OPEN(333,file="LAST_MIJ.DAT")	            !!!!! DELETE AFTER ALL
!	  DO k=1,total_size
!	  IF(ind_mat(1,k).eq.ind_mat(2,k)) THEN
!	  WRITE(333,*) k,Mat_el(n_r_coll,k)
!      ENDIF
!      ENDDO	
!      CLOSE(333)	
      IF(test_expansion_defined)CALL PRINT_ELASTIC_MIJ  
      END SUBROUTINE PRINT_MIJ_USER
      SUBROUTINE PRINT_ELASTIC_MIJ
! This subroutine is written by Alexander Semenov
      USE CONSTANTS	  
      USE VARIABLES
      USE MPI	  
      USE MPI_DATA
      USE OLD_MIJ
      IMPLICIT NONE
      INTEGER i,j,k
      INTEGER non_zero_size	  
      REAL*8, ALLOCATABLE :: Mij_elast(:,:),
     & Mij_elast_cs(:,:)
      INTEGER, ALLOCATABLE :: index_elastic_corr(:)
      CHARACTER(LEN=22) :: elast_mij_out = "ELASTIC_ELEMENTS  .out"	  
      ALLOCATE(Mij_elast(n_r_coll,states_size))
      non_zero_size = 0	  
      DO k=1,total_size
      IF(ind_mat(1,k).eq.ind_mat(2,k)) THEN
      IF(m12(ind_mat(1,k)).eq.m_elastic_proj_print)
     & non_zero_size = non_zero_size +1	  
      Mij_elast(:,ind_mat(1,k)) = Mat_el(:,k)	  
      ENDIF	  
      ENDDO
      ALLOCATE(Mij_elast_cs(n_r_coll,non_zero_size))
	  ALLOCATE(index_elastic_corr(non_zero_size))
      non_zero_size = 0	  
      DO k=1,total_size
      IF(ind_mat(1,k).eq.ind_mat(2,k)) THEN
      IF(m12(ind_mat(1,k)).eq.m_elastic_proj_print) THEN
      non_zero_size = non_zero_size +1	  
      Mij_elast_cs(:,non_zero_size) = Mat_el(:,k)
      index_elastic_corr(non_zero_size) = k	  
      ENDIF	  
      ENDIF	  
      ENDDO	  
      WRITE(elast_mij_out(17:18),'(i2.2)')m_elastic_proj_print	  
      OPEN(234,FILE=elast_mij_out)
       DO j=1,non_zero_size
      IF(j.eq.1) WRITE(234,"(a6)",ADVANCE="NO") "R_COM" 
      k = 	index_elastic_corr(j)  
      WRITE(234,"(2x,i4,2x,i3,1x,i3,2x)",ADVANCE="NO") k,
     & ind_mat(1,k),ind_mat(2,k)
      ENDDO
      WRITE(234,*)	  
	  
      DO i=1,n_r_coll
      DO j=1,non_zero_size
      IF(j.eq.1) WRITE(234,"(f6.2)",ADVANCE="NO") R_COM(i) 	  
      WRITE(234,"(e17.5)",ADVANCE="NO") Mij_elast_cs(i,j)*autoeV*eVtown
      ENDDO
      WRITE(234,*)	  
      ENDDO	  
      CLOSE(234)	  
	  
      END SUBROUTINE PRINT_ELASTIC_MIJ	  