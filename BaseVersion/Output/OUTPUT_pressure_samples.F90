SUBROUTINE OUTPUT_pressure_samples
    USE Global
    IMPLICIT NONE
    
    INTEGER:: iBlock
    
    DO iBlock = 1, Max_Block
        CALL GLOBAL_SetPointers(iBlock)
        CALL OUTPUT_pressure_sample
    ENDDO

END SUBROUTINE OUTPUT_pressure_samples


SUBROUTINE OUTPUT_pressure_sample
    USE Global
    IMPLICIT NONE
    
    INTEGER:: Ip,Jp,Kp,L,M,N,KZ,ii
    REAL::    timee
    character (LEN=100):: sample_name
    
!	  if (Ki>= Ki_aver) then
        KZ=MOD(KI,KQ)
        do N=1,ThisBlock%Num_sample
            IP          =   ThisBlock%Isam(N)
            JP          =   ThisBlock%Jsam(N)
            KP          =   ThisBlock%Ksam(N)
            sample_name =   ThisBlock%FL_sam(N)

          ThisBlock%SAMPLE(1,N,KZ)=V(1,Ip,Jp,Kp)  !den
          ThisBlock%SAMPLE(2,N,KZ)=V(2,Ip,Jp,Kp)  !ru
          ThisBlock%SAMPLE(3,N,KZ)=V(3,Ip,Jp,Kp)  !rv
          ThisBlock%SAMPLE(4,N,KZ)=V(4,Ip,Jp,Kp)   ! rw
        ! ThisBlock%SAMPLE(2,N,KZ)=V(6,Ip,Jp,Kp) / V(1,Ip,Jp,Kp)  !vK
          ThisBlock%SAMPLE(5,N,KZ)= PP(Ip,Jp,Kp)  !pre
          ThisBlock%SAMPLE(6,N,KZ)=  T(Ip,Jp,Kp)   ! temp
!
        if (KZ == 0) then
            OPEN(30,FILE=trim(FilePathPrefix)//sample_name,STATUS='UNKNOWN',ACCESS='APPEND')
            do ii=1,KQ-1
                timee=DT_LUSGS_Main*float(KI-KQ+ii)
                write(30,101)KI-KQ+ii,timee,ThisBlock%SAMPLE(1,N,ii),ThisBlock%SAMPLE(2,N,ii),ThisBlock%SAMPLE(3,N,ii),ThisBlock%SAMPLE(4,N,ii),ThisBlock%SAMPLE(5,N,ii),ThisBlock%SAMPLE(6,N,ii)
            enddo
            timee=DT_LUSGS_Main*float(KI)
            write(30,101)KI,timee,ThisBlock%SAMPLE(1,N,0),ThisBlock%SAMPLE(2,N,0),ThisBlock%SAMPLE(3,N,0),ThisBlock%SAMPLE(4,N,0),ThisBlock%SAMPLE(5,N,0),ThisBlock%SAMPLE(6,N,0)
            close(30)
        endif
    
    ! write(*,*)KI,timee,ThisBlock%SAMPLE(1,N,0),ThisBlock%SAMPLE(2,N,0),ThisBlock%SAMPLE(3,N,0),ThisBlock%SAMPLE(4,N,0),ThisBlock%SAMPLE(5,N,0),ThisBlock%SAMPLE(6,N,0)
    !pause       
     
     enddo
       
101  format(1x,I6,2x,7(E16.8,2x))
END SUBROUTINE OUTPUT_pressure_sample
