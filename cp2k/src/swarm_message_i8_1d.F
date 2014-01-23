! *****************************************************************************
!> \brief Addes an entry from a swarm-message.
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_add_i8_1d(msg, key, value)
   TYPE(swarm_message_type), INTENT(INOUT)   :: msg
   CHARACTER(LEN=*), INTENT(IN)              :: key
   INTEGER(KIND=int_8), DIMENSION(:) , INTENT(IN)                        :: value

   TYPE(message_entry_type), POINTER :: new_entry

   IF(swarm_message_haskey(msg, key)) &
      CALL mp_abort("swarm_message_add_i8_1d: key already exists: "//TRIM(key))

   ALLOCATE(new_entry)
   new_entry%key = key

   ALLOCATE(new_entry%value_i8_1d(SIZE(value)))

   new_entry%value_i8_1d = value

   !WRITE (*,*) "swarm_message_add_i8_1d: key=",key, " value=",new_entry%value_i8_1d

   IF(.NOT. ASSOCIATED(msg%root)) THEN
      msg%root => new_entry
   ELSE
      new_entry%next => msg%root
      msg%root => new_entry
   ENDIF

 END SUBROUTINE swarm_message_add_i8_1d


! *****************************************************************************
!> \brief Returns an entry from a swarm-message.
!> \author Ole Schuett
! *****************************************************************************
 SUBROUTINE swarm_message_get_i8_1d(msg, key, value)
   TYPE(swarm_message_type), INTENT(IN)  :: msg
   CHARACTER(LEN=*), INTENT(IN)          :: key
   INTEGER(KIND=int_8), DIMENSION(:), POINTER                             :: value

   TYPE(message_entry_type), POINTER :: curr_entry
   !WRITE (*,*) "swarm_message_get_i8_1d: key=",key

   IF(ASSOCIATED(value)) CALL mp_abort("swarm_message_get_i8_1d: value already associated")

   curr_entry => msg%root
   DO WHILE(ASSOCIATED(curr_entry))
      IF(TRIM(curr_entry%key) == TRIM(key)) THEN
         IF(.NOT. ASSOCIATED(curr_entry%value_i8_1d)) &
            CALL mp_abort("swarm_message_get_i8_1d: value not associated key: "//TRIM(key))
         ALLOCATE(value(SIZE(curr_entry%value_i8_1d)))
         value = curr_entry%value_i8_1d
         !WRITE (*,*) "swarm_message_get_i8_1d: value=",value
         RETURN
      ENDIF
      curr_entry => curr_entry%next
   END DO
   CALL mp_abort("swarm_message_get: key not found: "//TRIM(key))
 END SUBROUTINE swarm_message_get_i8_1d

!EOF
