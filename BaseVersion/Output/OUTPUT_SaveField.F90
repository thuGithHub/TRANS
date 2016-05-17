SUBROUTINE OUTPUT_SaveFields
    USE Global
    IMPLICIT NONE
    INTEGER:: iBlock
    INTEGER:: NI1,NJ1,NK1
    INTEGER:: I,J,K,L,LL
    
    DO iBlock = 1, Max_Block
        CALL GLOBAL_SetPointers(iBlock)
        NI1=ThisBlock%NI1
        NJ1=ThisBlock%NJ1
        NK1=ThisBlock%NK1

        DO L=1,2
          IF(L == 1)THEN
            OPEN(UNIT=15, FILE=ThisBlock%FLDAT, MODE='WRITE', STATUS='UNKNOWN', FORM="UNFORMATTED")  !, SHARE='DENYWR')
          ELSE !Datb
            OPEN(UNIT=15, FILE=ThisBlock%FLDATb, MODE='WRITE', STATUS='UNKNOWN', FORM="UNFORMATTED")  !, SHARE='DENYWR')
          ENDIF
      
            WRITE(15) KI
            write(15) ((((V(LL,I,J,K),K=1,NK1),J=1,NJ1),I=1,NI1),LL=1,ML)
            write(15) (((Rmiu(i,j,k),K=1,NK1),j=1,NJ1),I=1,NI1)
            if(Kind_dual==1)then
                write(15) ((((Wt0(LL,i,j,k),k=1,Nk1),j=1,NJ1),i=1,NI1),LL=1,ML)
            elseif(Kind_dual==2)then
                write(15) ((((Wt0(LL,i,j,k),k=1,Nk1),j=1,NJ1),i=1,NI1),LL=1,ML)
                write(15) ((((Wt1(LL,i,j,k),k=1,Nk1),j=1,NJ1),i=1,NI1),LL=1,ML)
            endif
            close(15)
        ENDDO
    enddo
END SUBROUTINE




