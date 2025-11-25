#include "precision.h"
#include <stdlib.h>
#include <stdio.h>

int main()
{
	struct num *n1 = new_num(2,1);
	n1->binary_digits[0] = 0;
	n1->binary_digits[1] = 2; 
	n1->sign = 0;
	char *n1_string = num_to_char(n1);
	printf("n1: %s\n", n1_string);

	struct num *n2 = new_num(4,0);
	n2->binary_digits[0] = 1;
	n2->binary_digits[1] = 128;
	n2->binary_digits[2] = 0;
	n2->binary_digits[3] = 0;
	n2->sign = 0;
	char *n2_string = num_to_char(n2);
	printf("n2: %s = \n", n2_string);
	
	printf("%s * %s =\n", n1_string, n2_string);
	struct num *n3 = mul_num(n1, n2);
	printf("%s\n", num_to_char(n3));
	int rshift = -18;
	free_num(n3);
	n3 = rshift_num(n1, rshift);
	
	char *n3_string = num_to_char(n3);
	printf("n1 rshifted %d times: %s\n", rshift ,n3_string);
	
	/*struct num *n4 = div_num(n1, n2);
	char *n4_string = num_to_char(n4);
	printf("%s\n", n4_string);*/

	struct num *n5 = new_num(4,0);
	n5->binary_digits[0] = 1;
	n5->binary_digits[1] = 128;
	n5->binary_digits[2] = 0;
	n5->binary_digits[3] = 0;
	struct num *precision = new_num(6,0);
	precision->binary_digits[5] = 1;
	struct num *n6 = inverse_num(n5, precision);
	printf("inverse of %s is %s\n", num_to_char(n5), num_to_char(n6));
}
