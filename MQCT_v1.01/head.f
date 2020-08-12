      MODULE MPI_DATA
! This module is written by Alexander Semenov
      IMPLICIT NONE	  
      INTEGER myid,nproc,ierr_mpi,id_cnt
      REAL*8 TIME_START,TIME_FINISH,TIME_SPEND
      REAL*8,ALLOCATABLE ::	  TOTAL_TIME(:), TOTAL_TIME_INI(:)
      REAL*8 TIME_BEGIN,TIME_CURRENT,TIME_1_MAT,TIME_2_MAT,TIME_3_MAT
      REAL*8 TIME_1_ARRAY
      REAL*8 TIME_INI,TIME_PROPAGATE,TIME_DOTRAJECT	  
      END MODULE MPI_DATA
      PROGRAM MQCT_HEAD
! This program is written by Alexander Semenov and modified by Bikramaditya Mandal
      USE MPI_DATA
      USE MPI
      CALL MPI_INIT( ierr_mpi ) !/MPI INITIALIZATION
      IF(ierr_mpi.ne.0) THEN 
      WRITE(*,'(a22,1x,i5)')'ERROR: MPI_INI_FAIL = ',ierr_mpi 
      STOP
      ENDIF
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, nproc, ierr_mpi ) !/MPI NUMBER OF PROCESSORS INI
      IF(ierr_mpi.ne.0) THEN 
      WRITE(*,'(a22,1x,i5)')'ERROR: MPI_INI_FAIL = ',ierr_mpi 
      STOP
      ENDIF
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr_mpi ) !/ MPI PROCESSORS ID INI
      IF(ierr_mpi.ne.0) THEN 
      WRITE(*,'(a22,1x,i5)')'ERROR: MPI_INI_FAIL = ',ierr_mpi 
      STOP
      ENDIF
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )     !/ MPI CHECKING IF ALL PROCESSORS READY
      IF(ierr_mpi.ne.0) THEN 
      WRITE(*,'(a22,1x,i5)')'ERROR: MPI_INI_FAIL = ',ierr_mpi 
      STOP
      ENDIF
      IF(MYID.EQ.0) PRINT*, "MPI INITIALIZED"
c/ READING INPUT AND SETUP MATRIX Mij(R)
      IF(MYID.EQ.0) PRINT*, 
     & "PARSERING OF USER INPUT"
      TIME_BEGIN = MPI_Wtime()
      CALL INITIALIZATION
!      CALL CPU_TIME(TIME_START)	  
      TIME_START = 	MPI_Wtime()  ! TIME START
      TIME_INI = TIME_START-TIME_BEGIN ! TIME SPENT ON INITIALIZATION
      ALLOCATE(TOTAL_TIME_INI(nproc))
      CALL MPI_GATHER(TIME_INI,1,MPI_DOUBLE_PRECISION,
     & TOTAL_TIME_INI,
     &	  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)	  
c/ SAMPLING AND PROPAGATION OF MQCT EQUATIONS
      TIME_PROPAGATE = MPI_Wtime() 	  
      CALL PROPAGATE
!      CALL CPU_TIME(TIME_SPEND)	  
      TIME_FINISH = MPI_Wtime()  ! TIME FINISH
      TIME_SPEND = TIME_FINISH-TIME_START	  
      ALLOCATE(TOTAL_TIME(nproc)) ! GATHERING TIMES FROM ALL TASKS
      CALL MPI_BARRIER( MPI_COMM_WORLD, ierr_mpi )
      CALL MPI_GATHER(TIME_SPEND,1,MPI_DOUBLE_PRECISION,
     & TOTAL_TIME,
     &	  1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr_mpi)
      CALL PRINT_OUTPUT ! PRINTING OUTPUT
      DEALLOCATE(TOTAL_TIME)	 !! MEMORY FREE
      CALL MPI_BARRIER(MPI_COMM_WORLD, ierr_mpi )	  
      CALL MPI_FINALIZE (ierr_mpi)		  
      END
      INCLUDE "iotest.f"
      INCLUDE "propagator.f"
      INCLUDE "print_output.f"	  
      INCLUDE "asym_top.f"
      INCLUDE "vib_diatomic.f"
      INCLUDE "pot_comp.f"	  
      INCLUDE "matrix_ini.f"
      INCLUDE "mqct_eqtns.f"	  
      INCLUDE "pes_sys_type.f"	  
      INCLUDE "suppl_routines.f"
	  

