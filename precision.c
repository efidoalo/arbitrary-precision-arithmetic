#include "precision.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

void free_num(struct num *n)
{
	free(n->binary_digits);
	free(n);
}

struct num *new_num(long int len, long int rp_index)
{
	if (len == 0) {
		printf("Error creating num with len requested as zero.\n");
		exit(EXIT_FAILURE);
	}
	struct num *new_num = (struct num *)malloc(sizeof(struct num));
	if (new_num == 0) {
		printf("Error allocating num structure on heap.%s.\n",
			strerror(errno));
		exit(EXIT_FAILURE);
	}

	new_num->sign = 0;
	new_num->len = len;
	new_num->binary_digits = (unsigned char *)malloc(len);
	if ((new_num->binary_digits) == 0) {
		printf("Error allocating memory on heap for new number binary_digits string.%s.\n", strerror(errno));
		exit(EXIT_FAILURE);
	}
	new_num->radix_point_index = rp_index;
	memset(new_num->binary_digits, 0, new_num->len); // set all bits to zero
	return new_num;
}

// returns the number of bytes constituting the fractional part of n
long int fractional_byte_length(struct num *n)
{
	return (n->len) - 1 - (n->radix_point_index);
}

// returns the number of bytes that caonstitutes the integer component of n
long int integer_byte_length(struct num *n)
{
	return (n->radix_point_index) + 1;
}


unsigned char *obtain_next_integer_temp_array(unsigned char **temp_integer_array,
		                              long int len,
					      long int *returned_len)
{
	if ( (*temp_integer_array == 0) && (len == 0) ) {
		*temp_integer_array = (unsigned char *)malloc(1);
		if (*temp_integer_array == 0) {
			printf("Error allocating temp_integer_array on heap via malloc.%s.\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
		(*temp_integer_array)[0] = 1;
		*returned_len = 1;
		return *temp_integer_array;
	}
	unsigned char carry = 0, temp_carry = 0;
        int temp_val = 0;
	for (int i=(len-1); i>=0; --i) {
		temp_val = (2*(*temp_integer_array)[i]) + carry;
		if (temp_val > 9) {
			carry = 1;
			temp_val -= 10;
		}
		else {
			carry = 0;
		}
		(*temp_integer_array)[i] = temp_val;		
	}
	if (carry > 0) {
		*temp_integer_array = (unsigned char *)realloc(
						*temp_integer_array,
						len+1);
		if (*temp_integer_array == 0) {
			printf("Error allocating memory on heap fro temporary integer array. %s.\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
		for (int i=len; i>0;--i) {
			(*temp_integer_array)[i] = (*temp_integer_array)[i-1];
		}
		(*temp_integer_array)[0] = carry;	
		*returned_len = (len+1);
	}
	else {
		*returned_len = len;
	}
	return *temp_integer_array;
}

// temp_fractional_array has length len. We return a new temp_fractional_array
// half the value of length len+1
unsigned char *obtain_next_temp_fractional_array(unsigned char **temp_fractional_array, long int len)
{
	
	unsigned char *new_temp = (unsigned char *)realloc(*temp_fractional_array, len+1);
	if (new_temp == 0) {
		printf("Error allocating memory on heap for temp fractional array.%s.\n", strerror(errno));
		exit(EXIT_FAILURE);
	}
	new_temp[len] = 0;
	unsigned char carry = 0;
	if (len == 0) {
		carry = 5;
	}
	unsigned char temp_carry = 0;
	for (int i=0; i<=len; ++i) {
		if (new_temp[i] % 2) {
			temp_carry = 5;
		}
		else {
			temp_carry = 0;
		}
		new_temp[i] /= 2;
		new_temp[i] += carry;
		carry = temp_carry;
	}
	return new_temp;
}



void add_fractional_array(unsigned char **fractional_array, 
			  long int *fractional_array_len,
			  unsigned char **temp_fractional_array,
			  long int temp_fractional_array_len)
{
	if ((*fractional_array == 0) && (fractional_array_len == 0)) {
		*fractional_array = (unsigned char *)malloc(temp_fractional_array_len);
		if (*fractional_array == 0) {
			printf("Error allocating memory for fractional array. %s.\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
		for (int i=0; i<temp_fractional_array_len; ++i) {
			(*fractional_array)[i] = (*temp_fractional_array)[i];
		}
		*fractional_array_len = temp_fractional_array_len;
	}
	else {
		*fractional_array = (unsigned char *)realloc(*fractional_array, temp_fractional_array_len);
		if (*fractional_array == 0) {
			printf("Error allocating memory for fractional array on heap.%s.\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
		for (int i=(*fractional_array_len); i<temp_fractional_array_len; ++i) {
			(*fractional_array)[i] = 0;
		}
		unsigned char carry = 0;
		for (int i=(temp_fractional_array_len-1); i>=1; --i) {
			unsigned char temp_val = (*fractional_array)[i] + (*temp_fractional_array)[i] + carry;
			if (temp_val > 9) {
				carry = 1;
				temp_val -= 10;
			}
			else {
				carry = 0;
			}
			(*fractional_array)[i] = temp_val;
		}	
		(*fractional_array)[0] += (carry + (*temp_fractional_array)[0]);
		*fractional_array_len = temp_fractional_array_len;
	}
}

void add_integer_array(unsigned char **integer_array,
                       long int integer_array_len,
                       unsigned char **temp_integer_array,
                       long int temp_integer_array_len,
		       long int *returned_integer_array_len)
{
	unsigned char temp_carry = 0, carry = 0;
	int temp_val = 0;
	int no_of_extra_digits = temp_integer_array_len - integer_array_len;
	for (int i=integer_array_len-1; i>=0; --i) {
		temp_val = (*integer_array)[i] + (*temp_integer_array)[i+no_of_extra_digits] + carry;
		if (temp_val > 9) {
			carry = 1;
			temp_val -= 10;
		}
		else {
			carry = 0;
		}
		(*integer_array)[i] = temp_val;
	}
	for (int i=(no_of_extra_digits-1);
	     i>=0; --i) {
		(*integer_array) = (unsigned char *)realloc(*integer_array,
							 integer_array_len+1);
		if (*integer_array == 0) {
			printf("Error allocating memory on heap for integer array.%s.\n", strerror(errno));
			exit(EXIT_FAILURE);
		}
		for (int j=integer_array_len; j>0; --j) {
			(*integer_array)[j] = (*integer_array)[j-1];
		}
		temp_val = carry + (*temp_integer_array)[i];
		if (temp_val > 9) {
			carry = 1;
			temp_val -= 10;
		}
		else {
			carry = 0;
		}
		(*integer_array)[0] = temp_val;
		++integer_array_len;
	}
	if (carry == 1) {
		(*integer_array) = (unsigned char *)realloc(*integer_array,
                                                         integer_array_len+1);
                if (*integer_array == 0) {
                        printf("Error allocating memory on heap for integer array.%s.\n", strerror(errno));
                        exit(EXIT_FAILURE);
                }
		for (int j=integer_array_len; j>0; --j) {
                        (*integer_array)[j] = (*integer_array)[j-1];
                }
		(*integer_array)[0] = carry;
		++integer_array_len;	
	}
	*returned_integer_array_len = integer_array_len;
}

char *num_to_char(struct num *n)
{
	unsigned char *fractional_array = 0, *integer_array = 0, 
		      *temp_fractional_array = 0, *temp_integer_array=0,
		      fractional_decimal_digit_index = 0, carry_val = 0, 
		      temp_val = 0, first_fractional_digit_encountered = 0;
	long int temp_fractional_array_len = 0, 
		          fractional_array_len = 0, temp_integer_array_len = 0,
			  integer_array_len=0;
	for (int i=((n->radix_point_index)+1); i<(n->len); ++i) {
		for (int j=0; j<8; ++j) {
			temp_fractional_array = 
				obtain_next_temp_fractional_array(
						&temp_fractional_array, 
						temp_fractional_array_len);
			++temp_fractional_array_len;
			if ((n->binary_digits)[i] & (1 << (7-j))) {
				add_fractional_array(&fractional_array, 
						&fractional_array_len, 
						&temp_fractional_array, 
						temp_fractional_array_len);
				fractional_array_len = temp_fractional_array_len;
			}	
		}	
	}
	for (int i=(n->radix_point_index); i>=0; --i) {
		for (int j=0; j<8; ++j) {
			temp_integer_array = obtain_next_integer_temp_array(
					      &temp_integer_array,
					      temp_integer_array_len,
					      &temp_integer_array_len);
			if ( (n->binary_digits)[i] & (1<<j) ) {
				add_integer_array(&integer_array, 
						  integer_array_len,
						  &temp_integer_array,
						  temp_integer_array_len,
						  &integer_array_len);
			}
		}
	}
	char *num_string = (char *)malloc(integer_array_len + 2 + 
					  fractional_array_len);
	if (num_string == 0) {
		printf("Error allocating memory on heap for string representation of num.%s.\n", strerror(errno));
		exit(EXIT_FAILURE);
	}	
	for (int i=0; i<integer_array_len; ++i) {
		num_string[i] = '0' + integer_array[i];
	}
	num_string[integer_array_len] = '.';
	for (int i=0; i<fractional_array_len; ++i) {
		num_string[integer_array_len+1+i] = '0' + fractional_array[i];
	}
	num_string[integer_array_len + fractional_array_len + 1] = 0;
	if (n->sign == 1) {
		num_string = (char *)realloc(num_string, integer_array_len + 3 + 
				                         fractional_array_len);
		if (num_string == 0) {
			printf("Error reallocating num_sring to prepend minus sign.%s.\n",
				strerror(errno));
			exit(EXIT_FAILURE);
		}
		for (int i=integer_array_len+2+fractional_array_len; i>0; --i) {
			num_string[i] = num_string[i-1];
		}
		num_string[0] = '-';
	}
	return num_string;
}

struct num *grt_num(struct num *n1, struct num *n2)
{
	if (n1->radix_point_index > n2->radix_point_index) {
		long int radix_difference = n1->radix_point_index - n2->radix_point_index;
		for (int i=0; i<radix_difference; ++i) {
			if ((n1->binary_digits)[i]) {
				return n1;
			}
		}
		long int n1_len = (n1->len) - radix_difference;
		long int common_len = 0;
	        if (n1_len > (n2->len)) { 
			common_len = n2->len;
		}
		else {
			common_len = n1_len;
		}
		for (int i=0; i<common_len; ++i) {
			if ( ((n1->binary_digits)[radix_difference + i]) >
			     (n2->binary_digits)[i] ) {
				return n1;
			}
			else if ( (n2->binary_digits)[i] > 
				  ((n1->binary_digits)[radix_difference + i]) ) {
				return n2;
			}
		}
		if (n1_len > (n2->len)) {
			for (int i=radix_difference+common_len; i<(n1->len); ++i) {
				if ((n1->binary_digits)[i]) {
					return n1;
				}
			}
			return 0; // they are identical
		}
		else {
			for (int i=0; i<(n2->len - common_len); ++i) {
				if ((n2->binary_digits)[common_len+i]) {
					return n2;
				}
			}
			return 0; // they are identical
		}
	}		
	else {
		long int radix_difference = n2->radix_point_index - n1->radix_point_index;
                for (int i=0; i<radix_difference; ++i) {
                        if ((n2->binary_digits)[i]) {
                                return n2;
                        }
                }
                long int n2_len = (n2->len) - radix_difference;
                long int common_len = 0; 
                if (n2_len > (n1->len)) {
                        common_len = n1->len;
                }
                else {
                        common_len = n2_len;
                }
                for (int i=0; i<common_len; ++i) {
                        if ( ((n2->binary_digits)[radix_difference + i]) >
                             (n1->binary_digits)[i] ) {
                                return n2;
                        }
                        else if ( (n1->binary_digits)[i] > 
                                  ((n2->binary_digits)[radix_difference + i]) ) {
                                return n1;
                        }
                }
                if (n2_len > (n1->len)) {
                        for (int i=radix_difference+common_len; i<(n2->len); ++i) {
                                if ((n2->binary_digits)[i]) {
                                        return n2;
                                }
                        }
                        return 0; // they are identical
                }
                else {
                        for (int i=0; i<(n1->len - common_len); ++i) {
                                if ((n1->binary_digits)[common_len+i]) {
                                        return n1;
                                }
                        }
                        return 0; // they are identical
                }
	}
}

struct num *copy_num(struct num *n)
{
	struct num *ncopy = new_num(n->len, n->radix_point_index);
	ncopy->sign = n->sign;
	memcpy(ncopy->binary_digits, n->binary_digits, n->len);
	return ncopy;
}

struct num *add_num(struct num *n1, struct num *n2)
{
	unsigned char n2_is_zero = 1;
	for (int i=0; i<n2->len; ++i) {
		if (n2->binary_digits[i] != 0) {
			n2_is_zero = 0; 
		}
	}
	if (n2_is_zero) {
		struct num *res = copy_num(n1);
		return res;
	}
	unsigned char n1_is_zero = 1;
	for (int i=0; i<n1->len; ++i) {
                if (n1->binary_digits[i] != 0) {
                        n1_is_zero = 0;
                }
        }
        if (n1_is_zero) {
                struct num *res = copy_num(n2);
                return res;
        }

	if (n1->sign != n2->sign) {
		//
		if (grt_num(n1,n2)==n1) {
			struct num *sum = copy_num(n1);
			long int sum_bytes_before_index = sum->radix_point_index + 1;
			long int n2_bytes_before_index = n2->radix_point_index + 1;
		        long int sum_byte_index = 0;
			for (int i=0; i<(n2->len); ++i) {
				sum_byte_index = ((long int)(sum_bytes_before_index) - 
						 (long int)(n2_bytes_before_index)) + i;
				if ((sum_byte_index+1) > ((long int)(sum->len))) {
					sum->binary_digits = (unsigned char *)realloc(sum->binary_digits, sum_byte_index+1);
					if ((sum->binary_digits) == 0) {
						printf("Error reallocating memory on heap for"
						       " binary digits during add function.%s."
						       "\n", strerror(errno));
						exit(EXIT_FAILURE);
					}
					for (int j=sum->len; j<sum_byte_index+1; ++j) {
						sum->binary_digits[j] = 0;
					}
					sum->len = sum_byte_index+1;
				}
				if (sum_byte_index >= 0) { 
					for (int j=0; j<8; ++j) {
						if ( ((n2->binary_digits)[i]) & (1 << (7-j)) ) {
							if ( ((sum->binary_digits)[sum_byte_index]) & (1 << (7-j)) ) {
								unsigned char mask = 0;
								for (int k=0; k<8; ++k) {
									if (k != j) {
										mask += (1 << (7-k));
									}
								}	
								(sum->binary_digits)[sum_byte_index] &= mask;
							}
							else {
								unsigned char sub_handled = 0;
								unsigned char mask = 0;
								for (int k=j-1; k>=0; --k) {
									if ((sum->binary_digits)[sum_byte_index] & (1<<(7-k))) {
										sub_handled = 1;
										for (int l=j; l>k; --l) {
											(sum->binary_digits)[sum_byte_index] += (1 << (7-l));
										}
										mask = 0;
										for (int l=7; l>=0; --l) {
											if (l!=k) {
												mask += (1 << (7-l));
											}
										}
										(sum->binary_digits)[sum_byte_index] &= mask;
									}
									if (sub_handled) {
										break;
									}
								}
								if (sub_handled == 0) {
									unsigned char mask = 0;
									for (int k=0; k<=j; ++k) {
										(sum->binary_digits)[sum_byte_index] += (1 << 7-k);
									}	
									long int curr_sum_byte_index = sum_byte_index - 1;
									while ((sum->binary_digits)[curr_sum_byte_index] == 0) {
										(sum->binary_digits)[curr_sum_byte_index] = 255;
										--curr_sum_byte_index;
									}
									mask = 0;
									for (int k=0; k<8; ++k) {
										if (( (~((sum->binary_digits)[curr_sum_byte_index])) & (1 << k) )) {
											(sum->binary_digits)[curr_sum_byte_index] += (1 << k);		
										}
										else {
											for (int l=0; l<8; ++l) {
												if (l!=k) {
													mask += (1 << l);
												}
											}
											(sum->binary_digits)[curr_sum_byte_index] &= mask;
											break;
										}
									}
								}
							}
						}
					}
				}
			}
			return sum;		
		}	
		else if (grt_num(n1,n2)==n2) {
			struct num *sum = add_num(n2,n1);
			return sum;
		}	
	}	
	else {
		//
		long int n1_fractional_byte_len = fractional_byte_length(n1);
		long int n2_fractional_byte_len = fractional_byte_length(n2);
		long int sum_fractional_byte_len = 0;
		if (n1_fractional_byte_len > n2_fractional_byte_len) {
			sum_fractional_byte_len = n1_fractional_byte_len;
		}
		else {
			sum_fractional_byte_len = n2_fractional_byte_len;
		}
		long int n1_integer_byte_len = integer_byte_length(n1);
		long int n2_integer_byte_len = integer_byte_length(n2);
		long int sum_integer_byte_len = 1;
		if (n1_integer_byte_len > n2_integer_byte_len) {
			sum_integer_byte_len += n1_integer_byte_len;
		}
		else {
			sum_integer_byte_len += n2_integer_byte_len;
		}
		long int sum_len = sum_integer_byte_len + sum_fractional_byte_len;
		long int sum_radix_point_index = sum_integer_byte_len - 1;
		// zero with the right memory allocated with positive sign byte
		struct num *sum = new_num(sum_len, sum_radix_point_index);

		if (sum_fractional_byte_len == n1_fractional_byte_len) {
			unsigned char curr_byte = 0;
			unsigned char carry_bit = 0;
			long int curr_byte_pos = n1_fractional_byte_len; 
			for (int i=0; i<sum_len; ++i) {
				unsigned int curr_byte_val = 0;
				curr_byte_val += carry_bit;
				if (curr_byte_pos > 0) {
					curr_byte_val += (n1->binary_digits)[n1_integer_byte_len+curr_byte_pos-1];
					if (curr_byte_pos <= n2_fractional_byte_len) {
						curr_byte_val += (n2->binary_digits)[n2_integer_byte_len + curr_byte_pos-1];
					}
				}
				else {
					long int curr_integer_pos = -curr_byte_pos + 1;
					if (curr_integer_pos <= n1_integer_byte_len) {
						curr_byte_val += (n1->binary_digits)[n1_integer_byte_len-curr_integer_pos];
					}
					if (curr_integer_pos <= n2_integer_byte_len) {
						curr_byte_val += (n2->binary_digits)[n2_integer_byte_len-curr_integer_pos];
					}
				}
				if (curr_byte_val > 255) {
					curr_byte_val -= 256;
					carry_bit = 1;
				}
				else {
					carry_bit = 0;
				}
				--curr_byte_pos;
				(sum->binary_digits)[sum_len-1-i] = curr_byte_val;
			}
		}
		else {
			unsigned char curr_byte = 0;
                        unsigned char carry_bit = 0;
                        long int curr_byte_pos = n2_fractional_byte_len;
                        for (int i=0; i<sum_len; ++i) {
                                unsigned int curr_byte_val = 0;
                                curr_byte_val += carry_bit;
                                if (curr_byte_pos > 0) {
                                        curr_byte_val += (n2->binary_digits)[n2_integer_byte_len+curr_byte_pos-1];
                                        if (curr_byte_pos <= n1_fractional_byte_len) {
                                                curr_byte_val += (n1->binary_digits)[n1_integer_byte_len + curr_byte_pos-1];
                                        }
                                }
                                else {
                                        long int curr_integer_pos = (-curr_byte_pos) + 1;
                                        if (curr_integer_pos <= n1_integer_byte_len) {
                                                curr_byte_val += (n1->binary_digits)[n1_integer_byte_len-curr_integer_pos];
                                        }
                                        if (curr_integer_pos <= n2_integer_byte_len) {
                                                curr_byte_val += (n2->binary_digits)[n2_integer_byte_len-curr_integer_pos];
                                        }
                                }
                                if (curr_byte_val > 255) {
                                        curr_byte_val -= 256;
                                        carry_bit = 1;
                                }
				else {
					carry_bit = 0;
				}
                                --curr_byte_pos;
                        	(sum->binary_digits)[sum_len - 1 -i] = curr_byte_val;
			}
		}
		sum->sign = n1->sign;
		return sum;
	}	
}

struct num *sub_num(struct num *n1, struct num *n2)
{
	struct num *minus_n2 = copy_num(n2);
	if (minus_n2->sign == 0) {
		minus_n2->sign = 1;
	}
	else {
		minus_n2->sign = 0;
	}
	struct num *result = add_num(n1, minus_n2);
	free_num(minus_n2);
	return result;
}

struct num *rshift_num(struct num *n, long int bit_shift)
{
	if (bit_shift == 0) {
		struct num *result = copy_num(n);
		return result;
	}
	else if (bit_shift > 0) {
		// move decimal point to the right bit_shift number of times
		if ((bit_shift % 8) ==0) {
			long int byte_shift = bit_shift/8;
			if (byte_shift <= (n->len - n->radix_point_index - 1)) {
				struct num *result = copy_num(n);
				result->radix_point_index += byte_shift;
				return result;
			}
			else {
				unsigned long extra_bytes = byte_shift - ((n->len) - (n->radix_point_index) - 1);
				struct num *result = new_num(n->len + extra_bytes, n->len + extra_bytes - 1);
				memcpy(result->binary_digits, 
				       n->binary_digits,
				       n->len);
				result->sign = n->sign;
				return result;
				
			}
		}
		else {
			if ( ((n->radix_point_index) + (bit_shift/8)) < (n->len - 1) ) {
				long int extra_bytes = 1;
				long int new_radix_point = n->radix_point_index + (bit_shift/8) + 1;
				long int new_len = n->len + extra_bytes;
				struct num *res = new_num(new_len,
							  new_radix_point);
				res->sign = n->sign;
				int upper_bit_count = (8 - (bit_shift % 8)),
				    lower_bit_count = bit_shift % 8;
				unsigned char upper_mask = 0,
					      lower_mask = 0;
				for (int i=0; i<upper_bit_count; ++i) {
					upper_mask += (1 << i);
				}
				for (int i=0; i<lower_bit_count; ++i) {
					lower_mask += (1 << (7-i));
				}
				unsigned char curr_byte = 0;
				for (int i=(n->len-1); i>=0; --i) {
					curr_byte = 0;
					if (i == ((n->len)-1) ) {
						curr_byte += ((n->binary_digits)[i] & upper_mask) << lower_bit_count;
					}
					else {
						curr_byte += (((n->binary_digits)[i] & upper_mask) << lower_bit_count) + (((n->binary_digits)[i+1] & lower_mask) >> (8-lower_bit_count));
					}
					res->binary_digits[i+1] = curr_byte;
				}
				res->binary_digits[0] = (((n->binary_digits)[0] & lower_mask) >> (8-lower_bit_count)); 
				return res;
			}
			else {
				long int extra_bytes = (((long int)(n->radix_point_index)) + (bit_shift/8)) - (((long int)(n->len)) - 1) + 1;
				long int new_radix_point_index =
				       		n->radix_point_index + 
						((bit_shift)/8) + 1;
						
				struct num *res = new_num(n->len + extra_bytes,
							new_radix_point_index);
				res->sign = n->sign;
				int upper_bit_count = (8 - (bit_shift % 8)),
                                    lower_bit_count = bit_shift % 8;
                                unsigned char upper_mask = 0,
                                              lower_mask = 0;
                                for (int i=0; i<upper_bit_count; ++i) {
                                        upper_mask += (1 << i);
                                }
                                for (int i=0; i<lower_bit_count; ++i) {
                                        lower_mask += (1 << (7-i));
                                }
                                unsigned char curr_byte = 0;

				for (int i=(n->len); i>=0; --i) {
					if (i== (n->len)) {
						curr_byte = ((n->binary_digits)[i-1] & upper_mask) << lower_bit_count;	
					}
					else if (i>0) {
						curr_byte = (((n->binary_digits)[i] & lower_mask) >> (8-lower_bit_count)) + (((n->binary_digits)[i-1] & upper_mask) << lower_bit_count);
					}
					else {
						curr_byte = ((n->binary_digits)[i] & lower_mask) >> (8-lower_bit_count);
					}
					(res->binary_digits)[i] = curr_byte;
				}
				return res;			
			}	
		}
	}
	else {
		// bit_shift < 0
		bit_shift = -bit_shift;
		if ((bit_shift % 8)==0) {
			if ((bit_shift/8) <= (n->radix_point_index)) {
				struct num *res = copy_num(n);
				res->radix_point_index -= ((bit_shift)/8);
				return res;
			}
			else {
				long int extra_bytes = ((bit_shift/8) - (n->radix_point_index));
				long int new_len = n->len + extra_bytes;
				struct num *res = new_num(new_len,
							  0);
				for (int i=extra_bytes; i<new_len; ++i) {
					(res->binary_digits)[i] =
					(n->binary_digits)[i-extra_bytes];
				}
				res->sign = n->sign;
				return res;
			}
		}
		else {
			long int extra_bytes = 1;
			if ((bit_shift/8) > (n->radix_point_index)) {
				extra_bytes += ( (bit_shift/8) - (n->radix_point_index) );
			}
			long int res_len = n->len + extra_bytes;
			int long temp_radix_point = n->radix_point_index - (bit_shift/8);
			long int new_radix_point_index = 0;
			if (temp_radix_point > 0) {
				new_radix_point_index = temp_radix_point;
			}
			struct num *result = new_num(res_len, new_radix_point_index);
			long int curr_index = 0;
			long int len = (bit_shift/8) - (long int)(n->radix_point_index);
			for (int i=0; i<len; ++i) {
				(result->binary_digits)[i] = 0;
				++curr_index;
			}
			int lower_bit_count = 8-(bit_shift % 8),
			    upper_bit_count = (bit_shift % 8);
			unsigned char lower_mask = 0, upper_mask = 0;
			for (int i=0; i<lower_bit_count; ++i) {
				lower_mask += (1 << (7-i));
			}
			for (int i=0; i<upper_bit_count; ++i) {
				upper_mask += (1 << i);
			}
			unsigned char curr_byte = 0;
			for (int i=curr_index; i<res_len; ++i) {
				if (i==curr_index) {
					curr_byte = ((n->binary_digits)[i-curr_index] & lower_mask) >> upper_bit_count;
				}
				else if (i < (res_len-1))  {
					curr_byte = (((n->binary_digits)[i-curr_index-1] & upper_mask) << lower_bit_count) + 
						    (((n->binary_digits)[i-curr_index] & lower_mask) >> upper_bit_count);
				}
				else {
					curr_byte = (((n->binary_digits)[i-curr_index-1] & upper_mask) << lower_bit_count);
				}
				(result->binary_digits)[i] = curr_byte;
			}
			return result;	
		}
	}
}

struct num *mul_num(struct num *n1, struct num *n2)
{	
	struct num *result = new_num(1,0);
	long int shift = ((n2->radix_point_index + 1)*8) - 1;
	for (int i=0; i<n2->len; ++i) {
		for (int j=0; j<8; ++j) {
			if ((n2->binary_digits)[i] & (1 << (7-j))) {
				struct num *temp = rshift_num(n1, shift);
				struct num *temp_result = result;
				result = add_num(temp_result, temp);
				free_num(temp);
				free_num(temp_result);
			}
			--shift;
		}
	}
	return result;
}

// precision must be positive
struct num *inverse_num(struct num *n, struct num *precision)
{
	long int precision_bytes_after_radix = (precision->len - 1) - (precision->radix_point_index);
	
	char *n_str = num_to_char(n);
	long int n_str_len = strlen(n_str);
	if (n_str_len > 10) {
		n_str[10] = 0;
		n_str_len = 10;
	}
	struct num *one = new_num(1,0);
	one->binary_digits[0] = 1;
	struct num *initial_guess = 0;
	struct num *cmp_num = grt_num(one, n);
	if (cmp_num == 0) {
		return one;
	}
	else {
		long int largest_power = 0;
		for (int i=0; i<n->len; ++i) {
			for (int j=0; j<8; ++j) {
				if (n->binary_digits[i] & (1 << (7-j))) {
					if (i <= n->radix_point_index) {
						largest_power = (n->radix_point_index - i)*8 + (7-j);
						goto power_found;
					}
					else {
						largest_power = - (((n->radix_point_index - i)-1)*8 +(j+1));
						goto power_found;	
					}
				}
			}
		}
		power_found:
		largest_power *= -1;
		struct num *temp = new_num(1,0);
		temp->binary_digits[0] = 1;
		initial_guess = rshift_num(temp, largest_power);
	}
	struct num *two = new_num(1,0); 
	(two->binary_digits)[0] = 2;
	struct num *curr_approx = initial_guess;
	int iter_no = 0;
	int display_digit_precision = 10;
	while (1) {
		printf("iter no: %d\n", iter_no);
		struct num *temp1 = mul_num(curr_approx, n);
		struct num *temp2 = sub_num(two, temp1);
		free_num(temp1);
		struct num *next_approx = mul_num(curr_approx, temp2);
		char *next_approx_str = num_to_char(next_approx);
		printf("next approximation is: %s\n", next_approx_str);
	        free_num(temp2);
		free_num(curr_approx);
		curr_approx = next_approx;
		free(next_approx_str);
		struct num *prod = mul_num(curr_approx, n);
	        struct num *error = sub_num(prod, one);
		error->sign = 0;
		free_num(prod);
		cmp_num = grt_num(precision, error);
		if (cmp_num == precision) {
			free_num(error);
			free_num(one);
			free_num(two);
			return curr_approx;
		}
		char *curr_error_str = num_to_char(error);
		int curr_error_str_len = strlen(curr_error_str);
		printf("current  error:%s\n", curr_error_str);
		free(curr_error_str);
		free_num(error);
		if (((curr_approx->len-1) - curr_approx->radix_point_index) > (precision_bytes_after_radix+1)) {
			curr_approx->binary_digits = realloc(curr_approx->binary_digits, curr_approx->radix_point_index + precision_bytes_after_radix + 2);
			curr_approx->len = curr_approx->radix_point_index + precision_bytes_after_radix + 2;

		}
		++iter_no;	
	}
}

struct num *div_num(struct num *n1, struct num *n2)
{
	struct num *precision = new_num(n2->len+1,0);
	precision->binary_digits[n2->len] = 1;
	struct num *n2_reciprical = inverse_num(n2, precision);
	free_num(precision);
	struct num *result = mul_num(n1, n2_reciprical);
	free_num(n2_reciprical);
	return result;
}

void set_num_int(struct num **n, unsigned long int val, unsigned char s) 
{
	int val_size = sizeof(unsigned long int);
	if ((*n)==0) {
		*n = (struct num *)malloc(sizeof(struct num));
		if ((*n)==0) {
			printf("Error allocating memory for number n. %s.\n",
				strerror(errno));
			exit(EXIT_FAILURE);
		}
	}
	else {
		free((*n)->binary_digits);
	}
	(*n)->binary_digits = (unsigned char *)malloc(val_size);	
	if (((*n)->binary_digits) == 0) {
		printf("Error allocating memory for binary_digits array.%s.\n",
			strerror(errno));
		exit(EXIT_FAILURE);
	}
	(*n)->len = val_size;
	unsigned char *int_bytes = (unsigned char *)&val;
	for (int i=0; i<val_size; ++i) {
		((*n)->binary_digits)[val_size-1-i] = *int_bytes;
	       	++int_bytes;	
	}	
	(*n)->sign = s;
	(*n)->radix_point_index = val_size-1;
}
