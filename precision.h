/*===============================;
 *
 * File: precision.h
 * Content: Header file for arbitrary
 * precision arithmetic.
 * Date: 14/11/2025
 * Author: Andy Oldham
 *
 **********************************/

#ifndef __PRECISION_H__
#define __PRECISION_H__

struct num 
{
	unsigned char sign; // 0 for positive, 1 for negative
	unsigned char *binary_digits; // binary digit string to contain the
				      // binary information
	long int len; // number of bytes in binary_digits
        long int radix_point_index; // 0 starting index 
						// into binary_digits
	// radix point is between the least significant bit of 
	// binary_digits[radix_point_index] and most significant bit of
	// binary_digits[radix_point_index + 1]
};

// right shifts n by bit_shift, positive bit_shift is multiply,
//                        negative bit_dhift id divide,
//                        zero bit_shift returns a copy of n
struct num *rshift_num(struct num *n, long int bit_shift);

// make a copy of n1 and return the address of the numb allocated on the heap
struct num *copy_num(struct num *n1);

// compares n1 and n2 returns the number with larger absolute value
// or null pointer if they are equal
struct num *grt_num(struct num *n1, struct num *n2);

// constructs and returns the number 0.
// len is the len field of the returned number and rp_index is the 
// radix_pointz_index of the returned number. sign is positive
struct num *new_num(long int len, long int rp_index);

// add n1 to n2 and return the result. n1 and n2 remain the same
struct num *add_num(struct num *n1, struct num *n2);

// return the result of n1 - n2. n1 and n2 remain the same
struct num *sub_num(struct num *n1, struct num *n2);

// returns the product of n1 and n2, leaving n1 and n2 the same
struct num *mul_num(struct num *n1, struct num *n2);

// returns n1 / n2, leaving n1 and n2 unchanged
struct num *div_num(struct num* n1, struct num *n2);

// finds the ijnverse of n, within precision given by prec
struct num *inverse_num(struct num *n, struct num *prec);

// sets n in place to the integer given by val of sign s where s==0 for positive
// , 1 for negative
void set_num_int(struct num **n, unsigned long int val, unsigned char s);

// returns a null terminated decimal digit string representation of n
char *num_to_char(struct num *n);
// free the number n and all associated memory
void free_num(struct num *n);
#endif
