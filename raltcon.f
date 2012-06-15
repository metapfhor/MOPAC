      SUBROUTINE RALTCON(IREAD,LOPT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      COMMON /CNSTR / ICONXN, APPLIED, INFINT
      LOGICAL LEADSP, APPLIED
      CHARACTER (LEN=MAXCHAR) :: LINE
      CHARACTER (LEN=MAXCHAR) :: CHUNK
      DOUBLE PRECISION TMPC(8)
      INTEGER ICONXN(6,NUMATM),CONX(4),II,JJ,EMPTY(1),TSIZE
      INTEGER, ALLOCATABLE :: TATMS(:)
      DIMENSION ISTART(MAXCHAR/2), LOPT(3,NUMATM)
      CHARACTER SPACE*1,COMMA*1,LSQB*1,RSQB*1,LCRB*1,RCRB*1
      DATA COMMA,SPACE,LSQB,RSQB,LCRB,RCRB,II,JJ
     1 /',',' ','[',']','{','}',0,0/


   10 READ(IREAD,'(A)',END=130)LINE
      IF(LINE.EQ.' ') GOTO 130


      LEADSP=.TRUE.
      NVALUE=0
      DO 20 I=1,MAXCHAR
         IF (LEADSP.AND.LINE(I:I).NE.SPACE) THEN
            NVALUE=NVALUE+1
            ISTART(NVALUE)=I
         END IF
         LEADSP=(LINE(I:I).EQ.SPACE)
   20 CONTINUE

      II=0
      DO WHILE(II.LE.MAXCHAR)
        II=II+1
        PRINT *, II
        IF(LINE(II:II).EQ.LSQB)THEN
            JJ=II
            DO WHILE(JJ.LE.MAXCHAR)
            JJ=JJ+1
                IF(LINE(JJ:JJ).EQ.COMMA)THEN
                    LINE(JJ:JJ)=';'
                ELSEIF(LINE(JJ:JJ).EQ.RSQB)THEN
                    II=JJ+1
                    JJ=MAXCHAR
                ENDIF
           END DO
        ENDIF
      END DO

      TMPC(7)=0
      DO 30 I=1,6
       TMPC(I)=0
       TMPC(I)=READN(LINE,ISTART(I))
       IF(TMPC(7).EQ.0.AND.TMPC(I).EQ.0)TMPC(7)=I

   30 CONTINUE


      IF(TMPC(7).EQ.3) THEN
       GOTO 40
      ELSEIF(TMPC(7).EQ.4) THEN
       GOTO 50
      ELSEIF(TMPC(7).EQ.5) THEN
       GOTO 60
      ENDIF

C     FIRST INFO ON THIS LINE IS A TRANSLATION
   40   CONX(1)=INT(TMPC(1))
        LINE=LINE(ISTART(2)-1:)
      DO I=1,MAXCHAR/2
           IF(LINE.EQ.' ')EXIT
          CALL SPLIT(LINE,',',CHUNK)
          LEADSP=.TRUE.
          NVALUE=0
          DO II=1,MAXCHAR/2
            ISTART(II)=0
          END DO
          DO 70 II=1,MAXCHAR
             IF (LEADSP.AND.CHUNK(II:II).NE.SPACE) THEN
                NVALUE=NVALUE+1
                ISTART(NVALUE)=II
             END IF
             LEADSP=(CHUNK(II:II).EQ.SPACE)
   70     CONTINUE

          TMPC(7)=0
          DO 80 II=1,6
           TMPC(II)=0
           TMPC(II)=READN(CHUNK,ISTART(II))
           IF(TMPC(7).EQ.0.AND.TMPC(II).EQ.0)TMPC(7)=II
   80     CONTINUE
          IF(ALLOCATED(TATMS).EQ.FALSE)THEN
            ALLOCATE(TATMS(1))
          ELSE
            DEALLOCATE(TATMS)
            ALLOCATE(TATMS(1))
          ENDIF
          TATMS=EMPTY(:)
          IF(TMPC(7).EQ.2) THEN
            CHUNK=CHUNK(ISTART(4):)
            CONX(2)=INT(TMPC(2))
            CALL GETRANGE(CHUNK,LEN(CHUNK),CONX(2),TATMS,TSIZE)
           CALL ADDTRANS(CONX(1),CONX(2),TATMS,TSIZE)
          ELSEIF(TMPC(7).EQ.3) THEN
           GOTO 50
          ELSEIF(TMPC(7).EQ.4) THEN
           GOTO 60
          ENDIF



      END DO

      IF(ICONXN(1,CONX(2)).EQ.0.OR.ICONXN(1,CONX(2)).EQ.CONX(1))THEN
        ICONXN(1,CONX(2))=CONX(1)
        ICONXN(4,CONX(2))=1
        IF(INDEX(LINE,'F').NE.0)THEN
            LOPT(1,CONX(2))=0
C            IF(TMPC(4).NE.0) IVAL(1,CONX(2))=TMPC(4)
        ELSEIF(INDEX(LINE,'S').NE.0.AND.TMPC(4).NE.0)THEN
C            IVAL(1,CONX(2))=TMPC(4)
        ENDIF
      ELSE
        GOTO 330
      ENDIF
      GOTO 10

C       Laurent: Reordering still needs to be implemented
   50 CONX(1)=INT(TMPC(1))
      CONX(2)=INT(TMPC(2))
      CONX(3)=INT(TMPC(3))
      IF((ICONXN(2,CONX(3)).EQ.0.OR.ICONXN(2,CONX(3)).EQ.CONX(1))
     1.AND.(ICONXN(1,CONX(3)).EQ.0.OR.ICONXN(1,CONX(3)).EQ.CONX(2)))THEN
        ICONXN(2,CONX(3))=CONX(1)
        ICONXN(1,CONX(3))=CONX(2)
        ICONXN(5,CONX(3))=1
        IF(INDEX(LINE,'F').NE.0)THEN
            LOPT(2,CONX(3))=0
C            IF(TMPC(5).NE.0)IVAL(2,CONX(3))=TMPC(5)
        ELSEIF(INDEX(LINE,'S').NE.0.AND.TMPC(5).NE.0)THEN
C             IVAL(2,CONX(3))=TMPC(5)
        ENDIF
      ELSE
        GOTO 330
      ENDIF

      GOTO 10
C       Laurent: Reordering still needs to be implemented
   60 CONX(1)=INT(TMPC(1))
      CONX(2)=INT(TMPC(2))
      CONX(3)=INT(TMPC(3))
      CONX(4)=INT(TMPC(4))
      IF((ICONXN(3,CONX(4)).EQ.0.OR.ICONXN(3,CONX(4)).EQ.CONX(1))
     1.AND.(ICONXN(2,CONX(4)).EQ.0.OR.ICONXN(2,CONX(4)).EQ.CONX(2))
     2.AND.(ICONXN(1,CONX(4)).EQ.0.OR.ICONXN(1,CONX(4)).EQ.CONX(3)))THEN
        ICONXN(3,CONX(4))=CONX(1)
        ICONXN(2,CONX(4))=CONX(2)
        ICONXN(1,CONX(4))=CONX(3)
        ICONXN(6,CONX(4))=1
        IF(INDEX(LINE,'F').NE.0)THEN
            LOPT(3,CONX(4))=0
C            IF(TMPC(6).NE.0)IVAL(3,CONX(4))=TMPC(6)
        ELSEIF(INDEX(LINE,'S').NE.0.AND.TMPC(6).NE.0)THEN
C            IVAL(3,CONX(4))=TMPC(6)
        ENDIF
      ELSE
        GOTO 330
      ENDIF
      GOTO 10

      GOTO 10

  330 WRITE(6,'(A)')'Impossible Gemoetric Constraints:', LINE
      STOP
C       End Laurent
  130 RETURN
      END

      SUBROUTINE SPLIT(STR,DELIMS,BEFORE)

        ! ROUTINE FINDS THE FIRST INSTANCE OF A CHARACTER FROM 'DELIMS' IN THE
        ! THE STRING 'STR'. THE CHARACTERS BEFORE THE FOUND DELIMITER ARE
        ! OUTPUT IN 'BEFORE'. THE CHARACTERS AFTER THE FOUND DELIMITER ARE
        ! OUTPUT IN 'STR'. THE OPTIONAL OUTPUT CHARACTER 'SEP' CONTAINS THE
        ! FOUND DELIMITER. A DELIMITER IN 'STR' IS TREATED LIKE AN ORDINARY
        ! CHARACTER IF IT IS PRECEDED BY A BACKSLASH (\). IF THE BACKSLASH
        ! CHARACTER IS DESIRED IN 'STR', THEN PRECEDE IT WITH ANOTHER BACKSLASH.

        !Modified version of algorithm written by Benthien, George
        !http://gbenthien.net/strings/index.html

      CHARACTER(LEN=*) :: STR,DELIMS,BEFORE
      LOGICAL :: PRES
      CHARACTER :: CH,CHA

      STR=ADJUSTL(STR)
      CALL COMPACT(STR)
      LENSTR=LEN_TRIM(STR)
      IF(LENSTR.EQ.0) RETURN        ! STRING STR IS EMPTY
      K=0
      IBSL=0                        ! BACKSLASH INITIALLY INACTIVE
      BEFORE=' '
      DO I=1,LENSTR
         CH=STR(I:I)
         IF(IBSL.EQ.1) THEN          ! BACKSLASH ACTIVE
            K=K+1
            BEFORE(K:K)=CH
            IBSL=0
            CYCLE
         END IF
         IF(CH.EQ.'\') THEN          ! BACKSLASH WITH BACKSLASH INACTIVE
            K=K+1
            BEFORE(K:K)=CH
            IBSL=1
            CYCLE
         END IF
         IPOS=INDEX(DELIMS,CH)
         IF(IPOS.EQ. 0) THEN          ! CHARACTER IS NOT A DELIMITER
            K=K+1
            BEFORE(K:K)=CH
            CYCLE
         END IF
         IF(CH.NE.' ') THEN          ! CHARACTER IS A DELIMITER THAT IS NOT A SPACE
            STR=STR(I+1:)
            EXIT
         END IF
         CHA=STR(I+1:I+1)            ! CHARACTER IS A SPACE DELIMITER
         IPOSA=INDEX(DELIMS,CHA)
         IF(IPOSA.GT.0) THEN          ! NEXT CHARACTER IS A DELIMITER
            STR=STR(I+2:)
            EXIT
         ELSE
            STR=STR(I+1:)
            EXIT
         END IF
      END DO
      IF(I.GE.LENSTR) STR=''
      STR=ADJUSTL(STR)              ! REMOVE INITIAL SPACES


      END

      SUBROUTINE COMPACT(STR)

      ! CONVERTS MULTIPLE SPACES AND TABS TO SINGLE SPACES; DELETES CONTROL CHARACTERS;
      ! REMOVES INITIAL SPACES.

      CHARACTER(LEN=*):: STR
      CHARACTER(LEN=1):: CH
      CHARACTER(LEN=LEN_TRIM(STR)):: OUTSTR

      STR=ADJUSTL(STR)
      LENSTR=LEN_TRIM(STR)
      OUTSTR=' '
      ISP=0
      K=0

      DO I=1,LENSTR
        CH=STR(I:I)
        ICH=IACHAR(CH)

        SELECT CASE(ICH)

          CASE(9,32)     ! SPACE OR TAB CHARACTER
            IF(ISP.EQ.0) THEN
              K=K+1
              OUTSTR(K:K)=' '
            END IF
            ISP=1

          CASE(33:)      ! NOT A SPACE, QUOTE, OR CONTROL CHARACTER
            K=K+1
            OUTSTR(K:K)=CH
            ISP=0

        END SELECT

      END DO

      STR=ADJUSTL(OUTSTR)

      END

      SUBROUTINE FOO(I)
      INTEGER I
      I=2+3
      RETURN
      END

      SUBROUTINE ADDVAL(A,V,N)
      INTEGER N
      INTEGER, ALLOCATABLE :: B(:),A(N)
      ALLOCATE(B(N+1))
      B(1:N)=A(1:N)
      B(N+1)=V
      DEALLOCATE(A)
      ALLOCATE(A(N+1))
      A(1:N+1)=B(1:N+1)
      END

      SUBROUTINE GETRANGE(STRING,LENSTR,I,ATMS,RANGESIZE)
      CHARACTER STRING*(*),TMP*(LEN(STRING)),CHUNK*(LEN(STRING))
      INTEGER I,LENSTR,ATMS(*),RANGESIZE,J,K,L
      TMP=STRING(2:LEN(STRING)-1)
      ATMS(1)=I
      RANGESIZE=1
      DO I=0,LENSTR
        IF(TMP.EQ.' ')EXIT
        CALL SPLIT(TMP,";",CHUNK)
        IF(INDEX(CHUNK,'-').NE.0)THEN
             J=INT(READA(CHUNK,1))
             K=INT(READA(CHUNK,INDEX(CHUNK,'-')))
             L=MAX(J,K)
             K=MIN(J,K)
             DO J=K,L
                IF(NINARR(ATMS,J,RANGESIZE).EQ.TRUE)THEN
                  CALL ADDVAL(ATMS,J,RANGESIZE)
                ENDIF
             END DO
        ELSE
            J=INT(READA(CHUNK,1))
            IF(NINARR(ATMS,J,RANGESIZE).EQ.TRUE)THEN
               CALL ADDVAL(ATMS,J,RANGESIZE)
            ENDIF
        ENDIF
      END DO
      END

      LOGICAL FUNCTION NINARR(A,I,N)
      INTEGER A(*),I,N
      NINARR=.TRUE.
      DO J=1,N
        IF(A(J).EQ.I)THEN
            NINARR=.FALSE.
            RETURN
        ENDIF
      END DO
      RETURN
      END


      SUBROUTINE ADDTRANS(I,J,ATMS,L)
      INTEGER I,J,ATMS(L)
      END
