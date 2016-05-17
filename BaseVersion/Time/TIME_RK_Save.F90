SUBROUTINE TIME_RK_Save
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L
    INTEGER:: NI1, NJ1, NK1
    INTEGER:: MLtemp

    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1
    
    
    MLtemp=ML
    
    IF(Itr_Chld == 1)then
        
        DO L=1,MLtemp
           DO J=1,NJ1
           DO K=1,NK1
           DO I=1,NI1
	        Wt2(L,I,J,K)=  V(L,I,J,K)
	     ENDDO
	     ENDDO
	     ENDDO
        ENDDO
    ENDIF
    
!    IF(Itr_Chld == 2)then
        
!        DO L=1,MLtemp
!           DO J=1,NJ1
!           DO K=1,NK1
!           DO I=1,NI1!
!	        Wt1(L,I,J,K)=  V(L,I,J,K)
!	     ENDDO
!	     ENDDO
!	     ENDDO
!        ENDDO
!    ENDIF
    
ENDSUBROUTINE 