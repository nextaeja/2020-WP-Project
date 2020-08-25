#include "helper.h"

long long pointer_to_int(double *pointer) {
	// This union is used to assign the pointer value to a retun array to be freed later
    // Type (double *) takes 8 bytes, so long long is used in the conversion
	union {long long integer; double *pointer;} int_pointer_conversion;

	int_pointer_conversion.integer = 0;
	int_pointer_conversion.pointer = pointer;

	return int_pointer_conversion.integer;
}

double *int_to_pointer(long long integer) {
	// This union is used to assign the pointer value to a retun array to be freed later
    // Type (double *) takes 8 bytes, so long long is used in the conversion
	union {long long integer; double *pointer;} int_pointer_conversion;

	int_pointer_conversion.integer = integer;
	return int_pointer_conversion.pointer;
}
