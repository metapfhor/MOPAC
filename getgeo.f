      SUBROUTINE GETGEO(IREAD,LABELS,GEO,LOPT,NA,NB,NC,AMS,NATOMS,INTR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INCLUDE 'SIZES'
      DIMENSION GEO(3,*),NA(*),NB(*),NC(*),AMS(*), LOPT(3,*)
     1,LABELS(*)
      LOGICAL INTR
************************************************************************
*
*   GETGEO READS IN THE GEOMETRY. THE ELEMENT IS SPECIFIED BY IT'S
*          CHEMICAL SYMBOL, OR, OPTIONALLY, BY IT'S ATOMIC NUMBER.
*
*  ON INPUT   IREAD  = CHANNEL NUMBER FOR READ, NORMALLY 5
*             AMS    = DEFAULT ATOMIC MASSES.
*
* ON OUTPUT LABELS = ATOMIC NUMBERS OF ALL ATOMS, INCLUDING DUMMIES.
*           GEO    = INTERNAL COORDINATES, IN ANGSTROMS, AND DEGREES.
*           LOPT   = INTEGER ARRAY, A '1' MEANS OPTIMIZE THIS PARAMETER,
*                    '0' MEANS DO NOT OPTIMIZE, AND A '-1' LABELS THE
*                    REACTION COORDINATE.
*           NA     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
*           NB     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
*           NC     = INTEGER ARRAY OF ATOMS (SEE DATA INPUT)
*           ATMASS = ATOMIC MASSES OF ATOMS.
************************************************************************
      COMMON /PATH  / IDUM(2),REACT(3,66), DUMM1,DUMM2
      COMMON /SIMBOL/ SIMBOL(MAXPAR)
      COMMON /ATMASS/ ATMASS(NUMATM)
      COMMON /ATOMTX/ LTXT, TXTATM(NUMATM)
      COMMON /KEYWRD/ KEYWRD
C       Laurent Modification
      COMMON /CNSTR / ICONXN, IVAL, APPLIED, INFINT
C       End Laurent
      DIMENSION ISTART(40), XYZ(3,NUMATM), VALUE(4)
C       Laurent Modification
      INTEGER ICONXN(6,NUMATM),CONX(4)
      DOUBLE PRECISION IVAL(3,NUMATM),INFINT
C       End Laurent
      LOGICAL LEADSP, IRCDRC, APPLIED
      CHARACTER KEYWRD*241, TXTATM*8, SIMBOL*10, LTXT*1
      CHARACTER ELEMNT(107)*2, LINE*80, SPACE*1, NINE*1,ZERO*1,
     1TAB*1, COMMA*1, STRING*80, ELE*2, TURN*1
      DOUBLE PRECISION TMPC(8)
      SAVE ELEMNT, COMMA, SPACE, NINE, ZERO
      DATA (ELEMNT(I),I=1,107)/'H','HE',
     1 'LI','BE','B','C','N','O','F','NE',
     2 'NA','MG','AL','SI','P','S','CL','AR',
     3 'K','CA','SC','TI','V','CR','MN','FE','CO','NI','CU',
     4 'ZN','GA','GE','AS','SE','BR','KR',
     5 'RB','SR','Y','ZR','NB','MO','TC','RU','RH','PD','AG',
     6 'CD','IN','SN','SB','TE','I','XE',
     7 'CS','BA','LA','CE','PR','ND','PM','SM','EU','GD','TB','DY',
     8 'HO','ER','TM','YB','LU','HF','TA','W','RE','OS','IR','PT',
     9 'AU','HG','TL','PB','BI','PO','AT','RN',
     1 'FR','RA','AC','TH','PA','U','NP','PU','AM','CM','BK','CF','XX',
     2 'FM','MD','CB','++','+','--','-','TV'/
      DATA COMMA,SPACE,NINE,ZERO,INFINT/',',' ','9','0',
     1    z'7FF0000000000000'/
      TAB=CHAR(9)
      IRCDRC=(INDEX(KEYWRD,'IRC')+INDEX(KEYWRD,'DRC') .NE.0)
      ILOWA = ICHAR('a')
      ILOWZ = ICHAR('z')
      ICAPA = ICHAR('A')
      ICAPZ = ICHAR('Z')
      MAXTXT=0
      NATOMS=0
      NUMAT=0
      ISERR=0
C             Laurent Modification
      APPLIED=.FALSE.
      DO 400 I=1,NUMATM
        DO 410 J=1,3
            IVAL(J,I)=INFINT
  410   CONTINUE
  400 CONTINUE
C           END LAURENT
      DO 10 I=1,MAXPAR
   10 SIMBOL(I)= '---'
   20 READ(IREAD,'(A)',END=130,ERR=230)LINE
C       Laurent Modification: To allow for cartesian coordinates to be read in with internal constraints
C      IF(LINE.EQ.' ') GOTO 130
C       Added
      IF(LINE.EQ.' ')THEN
       IF(INDEX(KEYWRD,'ALTCON').EQ.0)THEN
       GOTO 130
       ELSE
       GOTO 270
       ENDIF
      ENDIF
C       End Laurent
      NATOMS=NATOMS+1
C
C   SEE IF TEXT IS ASSOCIATED WITH THIS ELEMENT
C
      I=INDEX(LINE,'(')
      IF(I.NE.0)THEN
C
C  YES, ELEMENT IS LABELLED.
C
         K=INDEX(LINE,')')
         TXTATM(NATOMS)=LINE(I:K)
         MAXTXT=MAX(MAXTXT,K-I+1)
         STRING=LINE(1:I-1)//LINE(K+1:)
         LINE=STRING
      ELSE
         TXTATM(NATOMS)=' '
      ENDIF
*   CLEAN THE INPUT DATA
************************************************************************
      DO 30 I=1,80
         ILINE=ICHAR(LINE(I:I))
         IF(ILINE.GE.ILOWA.AND.ILINE.LE.ILOWZ) THEN
            LINE(I:I)=CHAR(ILINE+ICAPA-ILOWA)
         ENDIF
   30 CONTINUE
************************************************************************
      ICOMMA=ICHAR(COMMA)
      ITAB=ICHAR(TAB)
      DO 40 I=1,80
         KHAR=ICHAR(LINE(I:I))
         IF(KHAR.EQ.ICOMMA.OR.KHAR.EQ.ITAB)LINE(I:I)=SPACE
   40 CONTINUE
*
*   INITIALIZE ISTART TO INTERPRET BLANKS AS ZERO'S
      DO 50 I=1,10
   50 ISTART(I)=80
*
* FIND INITIAL DIGIT OF ALL NUMBERS, CHECK FOR LEADING SPACES FOLLOWED
*     BY A CHARACTER AND STORE IN ISTART
      LEADSP=.TRUE.
      NVALUE=0
      DO 60 I=1,80
         IF (LEADSP.AND.LINE(I:I).NE.SPACE) THEN
            NVALUE=NVALUE+1
            ISTART(NVALUE)=I
         END IF
         LEADSP=(LINE(I:I).EQ.SPACE)
C       Laurent Modification: THis makes it so the first line is completely ignored
C       Removed:
C      IF(LEADSP.AND.NATOMS.EQ.1.AND.NVALUE.EQ.1)LINE(I:)=' '
C       End Laurent
   60 CONTINUE
*
* ESTABLISH THE ELEMENT'S NAME AND ISOTOPE, CHECK FOR ERRORS OR E.O.DATA
*
      WEIGHT=0.D0
      STRING=LINE(ISTART(1):ISTART(2)-1)
      IF( STRING(1:1) .GE. ZERO .AND. STRING(1:1) .LE. NINE) THEN
*  ATOMIC NUMBER USED: NO ISOTOPE ALLOWED
         LABEL=READA(STRING,1)
         IF (LABEL.EQ.0) GO TO 120
         IF (LABEL.LT.0.OR.LABEL.GT.107) THEN
            WRITE(6,'(''  ILLEGAL ATOMIC NUMBER'')')
            GO TO 240
         END IF
         GO TO 80
      END IF
*  ATOMIC SYMBOL USED
      REAL=ABS(READA(STRING,1))
      IF (REAL.LT.1.D-15) THEN
*   NO ISOTOPE
         ELE=STRING(1:2)
      ELSE
         WEIGHT=REAL
         IF( STRING(2:2) .GE. ZERO .AND. STRING(2:2) .LE. NINE) THEN
            ELE=STRING(1:1)
         ELSE
            ELE=STRING(1:2)
         END IF
      END IF
*   CHECK FOR ERROR IN ATOMIC SYMBOL
      IF(ELE(1:1).EQ.'-'.AND.ELE(2:2).NE.'-')ELE(2:2)=' '
      DO 70 I=1,107
         IF(ELE.EQ.ELEMNT(I)) THEN
            LABEL=I
            GO TO 80
         END IF
   70 CONTINUE
      IF(ELE(1:1).EQ.'X')THEN
         LABEL=99
         GOTO 80
      ENDIF
      WRITE(6,'(''  UNRECOGNIZED ELEMENT NAME: ('',A,'')'')')ELE
      GOTO 240
*
* ALL O.K.
*
   80 IF (LABEL.NE.99) NUMAT=NUMAT+1
      IF(WEIGHT.NE.0.D0)THEN
         WRITE(6,'('' FOR ATOM'',I4,''  ISOTOPIC MASS:''
     1    ,F15.5)')NATOMS, WEIGHT
         ATMASS(NUMAT)=WEIGHT
      ELSE
         IF(LABEL .NE. 99)  ATMASS(NUMAT)=AMS(LABEL)
      ENDIF
      IF(NATOMS.GT.NUMATM)THEN
         WRITE(6,'(//10X,''****  MAX. NUMBER OF ATOMS ALLOWED:'',I4)')
     1NUMATM
         STOP
      ENDIF
             LABELS(NATOMS)   =LABEL
      IF(INDEX(KEYWRD,'ALTCON').EQ.0)THEN
       GEO(1,NATOMS)    =READA(LINE,ISTART(2))
       GEO(2,NATOMS)    =READA(LINE,ISTART(4))
       GEO(3,NATOMS)    =READA(LINE,ISTART(6))
      ELSE
       GEO(1,NATOMS)    =READA(LINE,ISTART(2))
       GEO(2,NATOMS)    =READA(LINE,ISTART(3))
       GEO(3,NATOMS)    =READA(LINE,ISTART(4))
      ENDIF
      IF(IRCDRC)THEN
         TURN=LINE(ISTART(3):ISTART(3))
         IF(TURN.EQ.'T')THEN
            LOPT(1,NATOMS)=1
            IF(NATOMS.EQ.1)WRITE(6,'(A)')' IN DRC MONITOR POTENTIAL ENER
     1GY'//' TURNING POINTS'
         ELSE
            LOPT(1,NATOMS)=0
         ENDIF
         TURN=LINE(ISTART(5):ISTART(5))
         IF(TURN.EQ.'T')THEN
            LOPT(2,NATOMS)=1
         ELSE
            LOPT(2,NATOMS)=0
         ENDIF
         TURN=LINE(ISTART(7):ISTART(7))
         IF(TURN.EQ.'T')THEN
            LOPT(3,NATOMS)=1
         ELSE
            LOPT(3,NATOMS)=0
         ENDIF
      ELSE
         IF(INDEX(KEYWRD,'ALTCON').EQ.0)THEN
          LOPT(1,NATOMS)   =READA(LINE,ISTART(3))
          LOPT(2,NATOMS)   =READA(LINE,ISTART(5))
          LOPT(3,NATOMS)   =READA(LINE,ISTART(7))
         ELSE
          LOPT(1,NATOMS)=1
          LOPT(2,NATOMS)=1
          LOPT(3,NATOMS)=1
         ENDIF
         DO 90 I=3,7,2
            IF(ICHAR(LINE(ISTART(I):ISTART(I))).GE.ICAPA.AND.
     1ICHAR(LINE(ISTART(I):ISTART(I))).LE.ICAPZ)ISERR=1
   90    CONTINUE
      ENDIF
      NA(NATOMS)       =READA(LINE,ISTART(8))
      NB(NATOMS)       =READA(LINE,ISTART(9))
      NC(NATOMS)       =READA(LINE,ISTART(10))
C
C  SPECIAL CASE OF USERS FORGETTING TO ADD DIHEDRAL DATA FOR ATOM 3
C
      IF(NATOMS.EQ.3)THEN
         IF(LOPT(3,3).EQ.2)THEN
            NA(3)=1
            NB(3)=2
            GEO(3,3)=0.D0
            LOPT(3,3)=0
         ELSEIF(LOPT(3,3).EQ.1.AND.ABS(GEO(3,3)-2.D0).LT.1.D-4)THEN
            NA(3)=2
            NB(3)=1
            GEO(3,3)=0.D0
            LOPT(3,3)=0
         ENDIF
      ENDIF
      IF(LOPT(1,NATOMS).GT.1.OR.LOPT(2,NATOMS).GT.1.OR.
     1LOPT(3,NATOMS).GT.1)ISERR=1
      IF(ISERR.EQ.1) THEN
C
C  MUST BE GAUSSIAN GEOMETRY INPUT
C
         DO 110 I=2,NATOMS
            DO 110 K=1,3
               J=GEO(K,I)+0.4D0
               IF(ABS(GEO(K,I)-J).GT.1.D-5)THEN
C
C   GEOMETRY CANNOT BE GAUSSIAN
C
                  WRITE(6,'(A)')' GEOMETRY IS FAULTY.  GEOMETRY READ IN
     1IS'
                  CONST=3.141592653598D0/180.D0
                  DO 100 L=1,NATOMS
                     GEO(2,L)=GEO(2,L)*CONST
  100             GEO(3,L)=GEO(3,L)*CONST
                  CALL GEOUT(6)
                  STOP
               ENDIF
  110    CONTINUE
         NATOMS=-1
         RETURN
      ENDIF
      GOTO 20
C       Laurent Modification: Reading-in of alternate constraints
  270 READ(IREAD,'(A)',END=130,ERR=230)LINE
      IF(LINE.EQ.' ') GOTO 130


      LEADSP=.TRUE.
      NVALUE=0
      DO 280 I=1,80
         IF (LEADSP.AND.LINE(I:I).NE.SPACE) THEN
            NVALUE=NVALUE+1
            ISTART(NVALUE)=I
         END IF
         LEADSP=(LINE(I:I).EQ.SPACE)
  280 CONTINUE

      TMPC(7)=0
      DO 290 I=1,6
       TMPC(I)=0
       TMPC(I)=READN(LINE,ISTART(I))
       IF(TMPC(7).EQ.0.AND.TMPC(I).EQ.0)TMPC(7)=I

  290 CONTINUE


      IF(TMPC(7).EQ.3) THEN
       GOTO 300
      ELSEIF(TMPC(7).EQ.4) THEN
       GOTO 310
      ELSEIF(TMPC(7).EQ.5) THEN
       GOTO 320
      ENDIF


  300 IF(TMPC(1).GT.TMPC(2))THEN
      TMPC(3)=TMPC(2)
      TMPC(2)=TMPC(1)
      TMPC(1)=TMPC(3)
      ENDIF
      CONX(1)=INT(TMPC(1))
      CONX(2)=INT(TMPC(2))
      IF(ICONXN(1,CONX(2)).EQ.0.OR.ICONXN(1,CONX(2)).EQ.CONX(1))THEN
        ICONXN(1,CONX(2))=CONX(1)
        ICONXN(4,CONX(2))=1
        IF(INDEX(LINE,'F').NE.0)THEN
            LOPT(1,CONX(2))=0
            IF(TMPC(4).NE.0) IVAL(1,CONX(2))=TMPC(4)
        ELSEIF(INDEX(LINE,'S').NE.0.AND.TMPC(4).NE.0)THEN
            IVAL(1,CONX(2))=TMPC(4)
        ENDIF
      ELSE
        GOTO 330
      ENDIF
      GOTO 270

C       Laurent: Reordering still needs to be implemented
  310 CONX(1)=INT(TMPC(1))
      CONX(2)=INT(TMPC(2))
      CONX(3)=INT(TMPC(3))
      IF((ICONXN(2,CONX(3)).EQ.0.OR.ICONXN(2,CONX(3)).EQ.CONX(1))
     1.AND.(ICONXN(1,CONX(3)).EQ.0.OR.ICONXN(1,CONX(3)).EQ.CONX(2)))THEN
        ICONXN(2,CONX(3))=CONX(1)
        ICONXN(1,CONX(3))=CONX(2)
        ICONXN(5,CONX(3))=1
        IF(INDEX(LINE,'F').NE.0)THEN
            LOPT(2,CONX(3))=0
            IF(TMPC(5).NE.0)IVAL(2,CONX(3))=TMPC(5)
        ELSEIF(INDEX(LINE,'S').NE.0.AND.TMPC(5).NE.0)THEN
             IVAL(2,CONX(3))=TMPC(5)
        ENDIF
      ELSE
        GOTO 330
      ENDIF

      GOTO 270
C       Laurent: Reordering still needs to be implemented
  320 CONX(1)=INT(TMPC(1))
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
            IF(TMPC(6).NE.0)IVAL(3,CONX(4))=TMPC(6)
        ELSEIF(INDEX(LINE,'S').NE.0.AND.TMPC(6).NE.0)THEN
            IVAL(3,CONX(4))=TMPC(6)
        ENDIF
      ELSE
        GOTO 330
      ENDIF
      GOTO 270

      GOTO 270

  330 WRITE(6,'(A)')'Impossible Gemoetric Constraints:', LINE
      STOP
C       End Laurent
*
* ALL DATA READ IN, CLEAN UP AND RETURN
*
  120 NATOMS=NATOMS-1
  130 NA(2)=1
      LTXT=CHAR(MAXTXT)
      IF(NATOMS.GT.3)THEN
         INTR=(NA(4).NE.0)
      ELSE
         IF(GEO(2,3).LT.10.AND.NATOMS.EQ.3)
     1WRITE(6,'(//10X,'' WARNING: INTERNAL COORDINATES ARE ASSUMED -'',/
     210X,'' FOR THREE-ATOM SYSTEMS '',//)')
         INTR=.TRUE.
      ENDIF
      IF(INTR)GEO(2,2)=0
C
C     READ IN VELOCITY VECTOR, IF PRESENT
C
      IF(INDEX(KEYWRD,'VELO').GT.0)THEN
         IF(INTR)THEN
            WRITE(6,'(A)')' COORDINATES MUST BE CARTESIAN WHEN VELOCITY'
     1//' VECTOR IS USED.'
            STOP
         ENDIF
C#      WRITE(6,'(/10X,A)')'INITIAL VELOCITY VECTOR FOR DRC'
         DO 150 I=1,NATOMS
            READ(5,'(A)') LINE
            CALL NUCHAR(LINE,VALUE,NDMY)
            IF(NDMY.NE.3)THEN
               WRITE(6,'(/10X,A)')
     1'  THERE MUST BE EXACTLY THREE VELOCITY DATA PER LINE'
               STOP
            ENDIF
            DO 140 J=1,3
  140       REACT(J,I+2)=VALUE(J)
C#      WRITE(6,'(2X,A2,2X,3F13.5)')ELEMNT(LABELS(I)),(VALUE(J),J=1,3)
  150    CONTINUE
         DO 160 I=1,3
            DO 160 J=1,2
  160    REACT(I,J)=GEO(I,J+1)-GEO(I,1)
C
C  NOW TO ROTATE VELOCITY VECTOR TO SUIT INTERNAL COORDINATE DEFINITION
C
C
C   ROTATE AROUND THE 1-2 X-AXIS TO AS TO ELIMINATE REACT(3,2)
C   (PUT ATOM 2 IN X-Y PLANE)
         SA=REACT(3,1)/SQRT(REACT(2,1)**2+REACT(3,1)**2+1.D-20)
         CA=SIGN(SQRT(1.D0-SA**2),REACT(2,1))
C#      LABELS(NATOMS+1)=1
C#      LABELS(NATOMS+2)=1
C#      WRITE(6,*)' FIRST ROTATION, ABOUT 1-2 X-AXIS'
         DO 170 I=1,NATOMS+2
            TEMP1= REACT(2,I)*CA+REACT(3,I)*SA
            TEMP2=-REACT(2,I)*SA+REACT(3,I)*CA
            REACT(2,I)=TEMP1
            REACT(3,I)=TEMP2
C#      WRITE(6,'(2X,A2,2X,3F13.5)')ELEMNT(LABELS(I)),(REACT(J,I),J=1,3)
  170    CONTINUE
C   ROTATE AROUND THE 1-2 Z-AXIS TO AS TO ELIMINATE REACT(2,2)
C   (PUT ATOM 2 ON X AXIS)
         CA=REACT(1,1)/SQRT(REACT(2,1)**2+REACT(1,1)**2+1.D-20)
         SA=SIGN(SQRT(1.D0-CA**2),REACT(2,1))
C#      WRITE(6,*)' SECOND ROTATION, ABOUT 1-2 Z-AXIS'
         DO 180 I=1,NATOMS+2
            TEMP1= REACT(1,I)*CA+REACT(2,I)*SA
            TEMP2=-REACT(1,I)*SA+REACT(2,I)*CA
            REACT(1,I)=TEMP1
            REACT(2,I)=TEMP2
C#      WRITE(6,'(2X,A2,2X,3F13.5)')ELEMNT(LABELS(I)),(REACT(J,I),J=1,3)
  180    CONTINUE
C   ROTATE AROUND THE 2-3 X-AXIS TO AS TO ELIMINATE REACT(3,3)
C   (PUT ATOM 3 ON X-Y PLANE)
         SA=REACT(3,2)/SQRT(REACT(2,2)**2+REACT(3,2)**2+1.D-20)
         CA=SIGN(SQRT(1.D0-SA**2),REACT(2,2))
C#      WRITE(6,*)' THIRD ROTATION, ABOUT 2-3 X-AXIS'
         DO 190 I=1,NATOMS+2
            TEMP1= REACT(2,I)*CA+REACT(3,I)*SA
            TEMP2=-REACT(2,I)*SA+REACT(3,I)*CA
            REACT(2,I)=TEMP1
            REACT(3,I)=TEMP2
C#      WRITE(6,'(2X,A2,2X,3F13.5)')ELEMNT(LABELS(I)),(REACT(J,I),J=1,3)
  190    CONTINUE
C
C  STRIP OFF FIRST TWO COORDINATES; THESE WERE THE COORDINATE AXIS
C  DEFINITIONS
C
         DO 200 I=1,NATOMS
            DO 200 J=1,3
  200    REACT(J,I)=REACT(J,I+2)
      ENDIF
      IF(  .NOT. INTR) THEN
         DO 210 I=1,NATOMS
            DO 210 J=1,3
  210    XYZ(J,I)=GEO(J,I)
         DEGREE=90.D0/ASIN(1.D0)
         CALL XYZINT(XYZ,NATOMS,NA,NB,NC,DEGREE,GEO)
         IF(INDEX(KEYWRD,' XYZ').EQ.0)THEN
C
C  UNCONDITIONALLY SET FLAGS FOR INTERNAL COORDINATES
C
            DO 220 I=1,3
               DO 220 J=I,3
  220       LOPT(J,I)=0
         ENDIF
         IF(ABS(GEO(2,3)-180.D0).LT.1.D-4.OR.ABS(GEO(2,3)).LT.1.D-4)
     1THEN
            WRITE(6,'(A)')' DUE TO PROGRAM BUG, THE FIRST THREE ATOMS MU
     1ST NOT LIE IN A STRAIGHT LINE.'
            STOP
         ENDIF
      ELSEIF (.NOT.IRCDRC) THEN
         LOPT(2,2)=0
         IF(LOPT(1,1)+LOPT(2,1)+LOPT(3,1)+LOPT(3,2)+
     1        LOPT(3,3) .GT. 0)THEN
            LOPT(1,1)=0
            LOPT(2,1)=0
            LOPT(3,1)=0
            LOPT(3,2)=0
            LOPT(3,3)=0
            WRITE(6,'(//10X,'' AN UNOPTIMIZABLE GEOMETRIC PARAMETER HAS'
     1',/10X,'' BEEN MARKED FOR OPTIMIZATION. THIS IS A NON-FATAL ''
     2,''ERROR'')')
         ENDIF
      ENDIF
      IF(NA(3).EQ.0) THEN
         NB(3)=1
         NA(3)=2
      ENDIF
      RETURN
* ERROR CONDITIONS
  230 IF(IREAD.EQ.5) THEN
         WRITE(6,'( '' ERROR DURING READ AT ATOM NUMBER '', I3 )')NATOMS
      ELSE
         NATOMS=0
         RETURN
      ENDIF
  240 J=NATOMS-1
      WRITE(6,'('' DATA CURRENTLY READ IN ARE '')')
      DO 250 K=1,J
  250 WRITE(6,260)LABELS(K),(GEO(JJ,K),LOPT(JJ,K),JJ=1,3),
     1NA(K),NB(K),NC(K)
  260 FORMAT(I4,2X,3(F10.5,2X,I2,2X),3(I2,1X))
      STOP
      END