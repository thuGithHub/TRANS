SUBROUTINE TURBTWO_OUTPUT_Residual
    USE Global
    IMPLICIT NONE

    !INTEGER,PARAMETER:: myid=0 !Not used until MPI


    OPEN(UNIT=19,FILE=ThisBlock%FLHISk, MODE='WRITE', STATUS='UNKNOWN',ACCESS='APPEND')  !,SHARE='DENYNONE')
    OPEN(UNIT=20,FILE=ThisBlock%FLHISo, MODE='WRITE', STATUS='UNKNOWN',ACCESS='APPEND')  !,SHARE='DENYNONE')
  
!
	   write(*,*) " Residuals of k-w-TM in blk: ",ThisBlock%id_present_blk

         WRITE(19,41)myid,KI,Itr_Chld,ThisBlock%DRSk,ThisBlock%IMk,ThisBlock%JMk,ThisBlock%KMk,ThisBlock%DUk
         WRITE(20,41)myid,KI,Itr_Chld,ThisBlock%DRSo,ThisBlock%IMo,ThisBlock%JMo,ThisBlock%KMo,ThisBlock%DUo
         WRITE(*,41)myid,KI,Itr_Chld,ThisBlock%DRSk,ThisBlock%IMk,ThisBlock%JMk,ThisBlock%KMk,ThisBlock%DUk
         WRITE(*,41)myid,KI,Itr_Chld,ThisBlock%DRSo,ThisBlock%IMo,ThisBlock%JMo,ThisBlock%KMo,ThisBlock%DUo
	   write(*,*) 

41      FORMAT(2X,I5,2X,I8,2X,I5,2x,E12.6,2X,3I5,2X,E12.6)
!
	close(19)
	close(20)

END SUBROUTINE TURBTWO_OUTPUT_Residual